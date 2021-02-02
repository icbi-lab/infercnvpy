from typing import Tuple, Union, Sequence
import numpy as np
from scanpy import logging
from anndata import AnnData
import scipy.ndimage
import scipy.sparse
import re


def _natural_sort(l: Sequence):
    """Natural sort without third party libraries.

    Adapted from: https://stackoverflow.com/a/4836734/2340703
    """

    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


def _running_mean(
    x: Union[np.ndarray, scipy.sparse.spmatrix], n: int = 50, step: int = 10
) -> np.ndarray:
    """
    Compute the running mean along rows of an array.

    Densifies the matrix. Use `step` and `chunksize` to save memory.

    Parameters
    ----------
    n
        Length of the running window
    step
        only compute running windows ever `step` columns, e.g. if step is 10
        0:100, 10:110, 20:120 etc. Saves memory.
    """
    return scipy.ndimage.uniform_filter1d(x, size=n, axis=1)[
        :,
        # remove edges. Should be n/2, but for some reason this still leads to an artifact at the left side.
        np.arange(n, x.shape[1] - n, step),
    ]


def _running_mean_by_chromosome(
    expr, var, window_size, step
) -> Tuple[dict, np.ndarray]:
    """Compute the running mean for each chromosome independently. Stack
    the resulting arrays ordered by chromosome.

    Parameters
    ----------
    expr
        A gene expression matrix, appropriately preprocessed
    var
        The var data frame of the associated AnnData object
    window_size
    step

    Returns
    -------
    chr_start_pos
        A Dictionary mapping each chromosome to the index of running_mean where
        this chromosome begins.
    running_mean
        A numpy array with the smoothed gene expression, ordered by chromosome
        and genomic position
    """
    chromosomes = _natural_sort(var["chromosome"].unique())

    def _running_mean_for_chromosome(chr):
        genes = var.loc[var["chromosome"] == chr].sort_values("start").index.values
        tmp_x = expr[:, var.index.get_indexer(genes)]
        return _running_mean(tmp_x, n=window_size, step=step)

    running_means = [_running_mean_for_chromosome(chr) for chr in chromosomes]

    chr_start_pos = {
        chr: i
        for chr, i in zip(chromosomes, np.cumsum([x.shape[1] for x in running_means]))
    }

    return chr_start_pos, np.hstack(running_means)


def _get_reference(
    adata: AnnData,
    reference_key: Union[str, None],
    reference_cat: Union[None, str, Sequence[str]],
    reference: Union[np.ndarray, None],
) -> np.ndarray:
    """Parameter validation extraction of reference gene expression"""
    if reference is None:
        if reference_key is None or reference_cat is None:
            logging.warning(
                "Using mean of all cells as reference. For better results, "
                "provide either `reference`, or both `reference_key` and `reference_cat`. "
            )  # type: ignore
            reference_expr = adata.X
        else:
            if isinstance(reference_cat, str):
                reference_cat = [reference_cat]
            reference_expr = adata.X[adata.obs[reference_key].isin(reference_cat), :]

        reference = np.mean(reference_expr, axis=0)

    if reference.size != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference


def infercnv(
    adata: AnnData,
    *,
    reference_key: Union[str, None] = None,
    reference_cat: Union[None, str, Sequence[str]] = None,
    reference: Union[np.ndarray, None] = None,
    lfc_cap: float = 3,
    window_size: int = 100,
    step: int = 10,
    dynamic_threshold: Union[float, None] = 1.5,
    median_filter: Union[int, None] = 5,
    chunksize: int = 5000,
    inplace: bool = True,
    key_added: str = "cnv",
) -> Union[None, Tuple[dict, scipy.sparse.csr_matrix]]:
    """
    Infer Copy Number Variation (CNV) by averaging gene expression over genomic
    regions.

    adata should already be filtered for low-quality cells. adata.X needs to be
    normalized and log-transformed. The method should be
    fairly robust to different normalization methods (`normalize_total`, scran, etc.).

    The method requires a "reference" value to which the expression of genomic
    regions is compared. If your dataset contains different cell types and includes
    both tumor and normal cells, the average of all cells can be used as reference.
    This is the default.

    If you already know which cells are "normal", you can provide a column
    from adata.obs to `reference_key` that contains the annotation. `reference_cat`
    specifies one or multiple values in `reference_key` that refer to normal cells.

    Alternatively, you can specify a numpy array with average gene expression values
    (for instance, derived from a different dataset), to `reference` which will be used
    as reference instead.


    Parameters
    ----------
    adata
        annotated data matrix
    reference

    reference_key
        Column name in adata.obs that contains tumor/normal annotations.
        If this is set to None, the average of all cells is used as reference.
    reference_cat
        One or multiple values in `adata.obs[reference_key]` that annotate
        normal cells.


    """
    if not adata.var_names.is_unique:
        raise ValueError("Ensure your var_names are unique!")
    if {"chromosome", "start", "end"} - set(adata.var.columns) != set():
        raise ValueError(
            "Genomic positions not found. There need to be `chromosome`, `start`, and "
            "`end` columns in `adata.var`. "
        )

    reference = _get_reference(adata, reference_key, reference_cat, reference)

    # Step 1 - compute log fold change
    x_centered = adata.X - reference
    # Step 2 - clip log fold changes
    x_clipped = np.clip(x_centered, -lfc_cap, lfc_cap)
    # Step 3 - smooth by genomic position
    chr_pos, x_smoothed = _running_mean_by_chromosome(
        x_clipped, adata.var, window_size=window_size, step=step
    )
    # Step 4 - center by cell
    x_cell_centered = x_smoothed - np.median(x_smoothed, axis=1)[:, np.newaxis]

    # Noise filtering
    x_res = x_cell_centered
    # step 5 - standard deviation based noise filtering
    if dynamic_threshold is not None:
        noise_thres = dynamic_threshold * np.std(x_res)
        x_res[np.abs(x_res) < noise_thres] = 0
    # step 6 - median filter
    if median_filter is not None:
        x_res = np.apply_along_axis(
            lambda row: scipy.ndimage.median_filter(row, size=median_filter),
            axis=1,
            arr=x_res,
        )

    x_res = scipy.sparse.csr_matrix(x_res)

    if inplace:
        adata.obsm[f"X_{key_added}"] = x_res
        adata.uns[key_added] = {"chr_pos": chr_pos}

    else:
        return chr_pos, x_res

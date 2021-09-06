from typing import Tuple, Union, Sequence
import numpy as np
from scanpy import logging
from anndata import AnnData
import scipy.ndimage
import scipy.sparse
import re
from .._util import _ensure_array, tqdm
from tqdm.contrib.concurrent import process_map
import itertools
from multiprocessing import cpu_count


def cnv_score(
    adata: AnnData,
    *,
    obs_key: str = "cnv_leiden",
    use_rep: str = "cnv",
    key_added: str = "cnv_score",
    inplace: bool = True,
):
    """Assign each cnv cluster a CNV score.

    Clusters with a high score are likely affected by copy number abberations.
    Based on this score, cells can be divided into tumor/normal cells.

    Ths score is currently simply defined as the mean of result of
    :func:`infercnvpy.tl.infercnv` for each cluster.

    Parameters
    ----------
    adata
        annotated data matrix
    obs_key
        Key under which the clustering is stored in adata.obs. Usually
        the result of :func:`infercnvpy.tl.leiden`, but could also be
        other clusters, e.g. obtained from transcriptomics data.
    use_rep
        Key under which the result of :func:`infercnvpy.tl.infercnv` is stored
        in adata.
    key_added
        Key under which the score will be stored in `adata.obs`.
    inplace
        If True, store the result in adata, otherwise return it.

    Returns
    -------
    Depending on the value of `inplace`, either returns `None` or a vector
    with scores.
    """
    if obs_key not in adata.obs.columns and obs_key == "cnv_leiden":
        raise ValueError(
            "`cnv_leiden` not found in `adata.obs`. Did you run `tl.leiden`?"
        )
    cluster_score = {
        cluster: np.mean(
            np.abs(adata.obsm[f"X_{use_rep}"][adata.obs[obs_key] == cluster, :])
        )
        for cluster in adata.obs[obs_key].unique()
    }
    score_array = np.array([cluster_score[c] for c in adata.obs[obs_key]])

    if inplace:
        adata.obs[key_added] = score_array
    else:
        return score_array


def infercnv(
    adata: AnnData,
    *,
    reference_key: Union[str, None] = None,
    reference_cat: Union[None, str, Sequence[str]] = None,
    reference: Union[np.ndarray, None] = None,
    lfc_clip: float = 3,
    window_size: int = 100,
    step: int = 10,
    dynamic_threshold: Union[float, None] = 1.5,
    exclude_chromosomes: Union[Sequence[str], None] = ("chrX", "chrY"),
    chunksize: int = 5000,
    n_jobs: Union[int, None] = None,
    inplace: bool = True,
    layer: Union[str, None] = None,
    key_added: str = "cnv",
) -> Union[None, Tuple[dict, scipy.sparse.csr_matrix]]:
    """
    Infer Copy Number Variation (CNV) by averaging gene expression over genomic
    regions.

    This method is heavily inspired by `infercnv <https://github.com/broadinstitute/inferCNV/>`_
    but more computationally efficient. The method is described in more detail
    in on the :ref:`infercnv-method` page.

    There, you can also find instructions on how to :ref:`prepare input data <input-data>`.

    Parameters
    ----------
    adata
        annotated data matrix
    reference_key
        Column name in adata.obs that contains tumor/normal annotations.
        If this is set to None, the average of all cells is used as reference.
    reference_cat
        One or multiple values in `adata.obs[reference_key]` that annotate
        normal cells.
    reference
        Directly supply an array of average normal gene expression. Overrides
        `reference_key` and `reference_cat`.
    lfc_clip
        Clip log fold changes at this value
    window_size
        size of the running window (number of genes in to include in the window)
    step
        only compute every nth running window where n = `step`. Set to 1 to compute
        all windows.
    dynamic_threshold
        Values `< dynamic threshold * STDDEV` will be set to 0, where STDDEV is
        the stadard deviation of the smoothed gene expression. Set to `None` to disable
        this step.
    exclude_chromosomes
        List of chromosomes to exclude. The default is to exclude genosomes.
    chunksize
        Process dataset in chunks of cells. This allows to run infercnv on
        datasets with many cells, where the dense matrix would not fit into memory.
    n_jobs
        Number of jobs for parallel processing. Default: use all cores.
        Data will be submitted to workers in chunks, see `chunksize`.
    inplace
        If True, save the results in adata.obsm, otherwise return the CNV matrix.
    layer
        Layer from adata to use. If `None`, use `X`.
    key_added
        Key under which the cnv matrix will be stored in adata if `inplace=True`.
        Will store the matrix in `adata.obsm["X_{key_added}"] and additional information
        in `adata.uns[key_added]`.

    Returns
    -------
    Depending on inplace, either return the smoothed and denoised gene expression
    matrix sorted by genomic position, or add it to adata.
    """
    if not adata.var_names.is_unique:
        raise ValueError("Ensure your var_names are unique!")
    if {"chromosome", "start", "end"} - set(adata.var.columns) != set():
        raise ValueError(
            "Genomic positions not found. There need to be `chromosome`, `start`, and "
            "`end` columns in `adata.var`. "
        )

    var_mask = adata.var["chromosome"].isnull()
    if np.sum(var_mask):
        logging.warning(
            f"Skipped {np.sum(var_mask)} genes because they don't have a genomic position annotated. "
        )  # type: ignore
    if exclude_chromosomes is not None:
        var_mask = var_mask | adata.var["chromosome"].isin(exclude_chromosomes)
    tmp_adata = adata[:, ~var_mask]

    expr = tmp_adata.X if layer is None else tmp_adata.layers[layer]

    if scipy.sparse.issparse(expr):
        expr = expr.tocsr()

    reference = _get_reference(tmp_adata, reference_key, reference_cat, reference)

    var = tmp_adata.var.loc[:, ["chromosome", "start", "end"]]  # type: ignore

    chr_pos, chunks = zip(
        *process_map(
            _infercnv_chunk,
            [expr[i : i + chunksize, :] for i in range(0, adata.shape[0], chunksize)],
            itertools.repeat(var),
            itertools.repeat(reference),
            itertools.repeat(lfc_clip),
            itertools.repeat(window_size),
            itertools.repeat(step),
            itertools.repeat(dynamic_threshold),
            tqdm_class=tqdm,
            max_workers=cpu_count() if n_jobs is None else n_jobs,
        )
    )
    res = scipy.sparse.vstack(chunks)

    if inplace:
        adata.obsm[f"X_{key_added}"] = res
        adata.uns[key_added] = {"chr_pos": chr_pos[0]}

    else:
        return chr_pos[0], res


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
    Compute a pyramidially weighted running mean.

    Densifies the matrix. Use `step` and `chunksize` to save memory.

    Parameters
    ----------
    n
        Length of the running window
    step
        only compute running windows ever `step` columns, e.g. if step is 10
        0:100, 10:110, 20:120 etc. Saves memory.
    """
    r = np.arange(1, n)
    pyramid = np.minimum(r, r[::-1])
    smoothed_x = np.apply_along_axis(
        lambda row: np.convolve(row, pyramid, mode="same"), axis=1, arr=x
    ) / np.sum(pyramid)
    return smoothed_x[:, np.arange(0, smoothed_x.shape[1], step)]


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
    chromosomes = _natural_sort(
        [x for x in var["chromosome"].unique() if x.startswith("chr") and x != "chrM"]
    )

    def _running_mean_for_chromosome(chr):
        genes = var.loc[var["chromosome"] == chr].sort_values("start").index.values
        tmp_x = expr[:, var.index.get_indexer(genes)]
        return _running_mean(tmp_x, n=window_size, step=step)

    running_means = [_running_mean_for_chromosome(chr) for chr in chromosomes]

    chr_start_pos = {
        chr: i
        for chr, i in zip(
            chromosomes, np.cumsum([0] + [x.shape[1] for x in running_means])
        )
    }

    return chr_start_pos, np.hstack(running_means)


def _get_reference(
    adata: AnnData,
    reference_key: Union[str, None],
    reference_cat: Union[None, str, Sequence[str]],
    reference: Union[np.ndarray, None],
) -> np.ndarray:
    """Parameter validation extraction of reference gene expression.

    If multiple reference categories are given, compute the mean per
    category.

    Returns a 2D array with reference categories in rows, cells in columns.
    If there's just one category, it's still a 2D array.
    """
    if reference is None:
        if reference_key is None or reference_cat is None:
            logging.warning(
                "Using mean of all cells as reference. For better results, "
                "provide either `reference`, or both `reference_key` and `reference_cat`. "
            )  # type: ignore
            reference = np.mean(adata.X, axis=0)

        else:
            obs_col = adata.obs[reference_key]
            if isinstance(reference_cat, str):
                reference_cat = [reference_cat]
            reference_cat = np.array(reference_cat)
            reference_cat_in_obs = np.isin(reference_cat, obs_col)
            if not np.all(reference_cat_in_obs):
                raise ValueError(
                    "The following reference categories were not found in "
                    "adata.obs[reference_key]: "
                    f"{reference_cat[~reference_cat_in_obs]}"
                )

            reference = np.vstack(
                [np.mean(adata.X[obs_col == cat, :], axis=0) for cat in reference_cat]
            )

    if reference.ndim == 1:
        reference = reference[np.newaxis, :]

    if reference.shape[1] != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference


def _infercnv_chunk(
    tmp_x, var, reference, lfc_cap, window_size, step, dynamic_threshold
):
    """The actual infercnv work is happening here.

    Process chunks of serveral thousand genes independently since this
    leads to (temporary) densification of the matrix.

    Parameters see `infercnv`.
    """
    # Step 1 - compute log fold change. This densifies the matrix.
    # Per default, use "bounded" difference calculation, if multiple references
    # are available. This mitigates cell-type specific biases (HLA, IGHG, ...)
    if reference.shape[0] == 1:
        x_centered = tmp_x - reference[0, :]
    else:
        ref_min = np.min(reference, axis=0)
        ref_max = np.max(reference, axis=0)
        # entries that are between the two "bounds" are considered having a logFC of 0.
        x_centered = np.zeros(tmp_x.shape, dtype=tmp_x.dtype)
        above_max = tmp_x > ref_max
        below_min = tmp_x < ref_min
        x_centered[above_max] = _ensure_array(tmp_x - ref_max)[above_max]
        x_centered[below_min] = _ensure_array(tmp_x - ref_min)[below_min]

    x_centered = _ensure_array(x_centered)
    # Step 2 - clip log fold changes
    x_clipped = np.clip(x_centered, -lfc_cap, lfc_cap)
    # Step 3 - smooth by genomic position
    chr_pos, x_smoothed = _running_mean_by_chromosome(
        x_clipped, var, window_size=window_size, step=step
    )
    # Step 4 - center by cell
    x_cell_centered = x_smoothed - np.median(x_smoothed, axis=1)[:, np.newaxis]

    # Noise filtering
    x_res = x_cell_centered
    # step 5 - standard deviation based noise filtering
    if dynamic_threshold is not None:
        noise_thres = dynamic_threshold * np.std(x_res)
        x_res[np.abs(x_res) < noise_thres] = 0

    x_res = scipy.sparse.csr_matrix(x_res)

    return chr_pos, x_res

"""Read in result files from scevan."""

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from infercnvpy._logging import log, log_result
from infercnvpy._util import _copy_for_inplace_false, _deprecated_positional


def _read_rdata_frame(path: Path, name: str) -> pd.DataFrame:
    """Read object `name` from an `.RData` file as a DataFrame (keeping row/column names)."""
    import rdata  # imported here, so that `import infercnvpy` does not load it

    obj = rdata.read_rda(path, default_encoding="utf-8")[name]
    if isinstance(obj, pd.DataFrame):
        return obj
    # matrices with dimnames can be returned as xarray DataArray
    if hasattr(obj, "to_pandas"):
        return obj.to_pandas()
    return pd.DataFrame(obj)


def _get_chr_pos_from_array(chr_pos_array):
    """Get starting position for each chromosome.

    Assumes that chr_pos_array is an array with chromosomes ascendingly (like
    in the scevan `count_mtx_annot` table).
    """
    chr_pos = {}
    for i, sn in enumerate(chr_pos_array):
        chr_name = f"chr{int(sn)}"
        if chr_name not in chr_pos:
            chr_pos[chr_name] = i
    return chr_pos


@_deprecated_positional("scevan_res_table")
def read_scevan(
    adata: AnnData,
    scevan_res_dir: str | Path,
    *,
    scevan_res_table: str | Path | None = None,
    subclones: bool = True,
    inplace: bool = True,
    subset: bool = True,
    key_added: str = "scevan",
    copy: bool = False,
) -> AnnData | None:
    """Load results from SCEVAN :cite:`DeFalco2021` for downstream analysis with infercnvpy.

    Requires that the cell barcodes used for SCEVAN and `adata.obs_names` match,
    but the order is irrelevant.

    .. versionchanged:: 0.7.0
        `scevan_res_table` is keyword-only; passing it positionally still works with a
        :class:`FutureWarning` until 0.8.0. New argument `copy` to get a modified copy of adata.

    .. deprecated:: 0.7.0
        `inplace=False` returns a modified copy of adata, with a :class:`FutureWarning`. From 0.8.0, it
        will return the results without changing adata, as the other functions of infercnvpy do. Use
        `copy=True` for a modified copy.

    Parameters
    ----------
    adata
        adata object to which the SCEVAN results shall be added
    scevan_res_table
        The results of `SCEVAN::pipelineCNA` saved as CSV file. Will add
        the columns `{key_added}_class`, `{key_added}_confident_normal` and, if SCEVAN was
        ran with subclone calling, `{key_added}_subclone` to `adata.obs`. This parameter can
        be omitted if you only want to load the CNV matrix.
    scevan_res_dir
        The output directory created by SCEVAN. Must only contain results for a single
        sample. Will read the files `*_CNAmtx.RData`, `*_CNAmtxSubclones.RData` (if available),
        and `*_count_mtx_annot.RData`.
    subclones
        Load the separate subclone segmentation, if available.
    subset
        If `True` (the default), subset anndata to the cells that have not been filtered
        out by the SCEVAN analysis. Otherwise the CNV matrix may contain `nan`-values and
        downstream analysis such as PCA will not work.
    key_added
        Used as prefix for the columns added to `adata.obs` and will add the CNV matrix
        as `X_{key_added}` to `adata.obsm`, and chromosome indices to `adata.uns[key_added]`.
    inplace
        If `True`, modify adata in place. `False` is deprecated, see above.
    copy
        Return a modified copy of adata and leave adata unchanged.

    Returns
    -------
    A modified copy of adata if `copy=True` (or `inplace=False`), otherwise `None`.
    Sets the following fields:

    `adata.obs["{key_added}_class"]`, `adata.obs["{key_added}_confident_normal"]`
        Only with `scevan_res_table`: tumour/normal class of each cell and whether it is a
        confident normal cell.
    `adata.obs["{key_added}_subclone"]`
        Only with `scevan_res_table` from a run with subclone calling: the subclone of each cell.
    `adata.obsm["X_{key_added}"]` : :class:`numpy.ndarray`
        CNV matrix, cells x segments. `NaN` for cells SCEVAN filtered out (with `subset=False`).
    `adata.uns[key_added]["chr_pos"]` : :class:`dict`
        Maps each chromosome to the column of `X_{key_added}` where it starts.
    """
    copy = _copy_for_inplace_false("read_scevan", inplace, copy, future="the CNV matrix, the obs columns and chr_pos")
    scevan_res_dir = Path(scevan_res_dir)
    scevan_res_file = list(scevan_res_dir.glob("*_CNAmtx.RData"))
    scevan_subclones_file = list(scevan_res_dir.glob("*_CNAmtxSubclones.RData"))
    scevan_anno_file = list(scevan_res_dir.glob("*_count_mtx_annot.RData"))

    if len(scevan_res_file) != 1 or len(scevan_subclones_file) > 1 or len(scevan_anno_file) != 1:
        raise ValueError(
            "There must be exactly one CNA_mtx and count_mtx_annot file and at most one "
            "CNAmtxSubclones file in the result directory!"
        )

    # read data
    if scevan_res_table is not None:
        tumor_normal_call = pd.read_csv(scevan_res_table, index_col=0)
    else:
        tumor_normal_call = None
        log.debug("read_scevan: no `scevan_res_table` given, so no tumor/normal classification is added")
    scevan_res = _read_rdata_frame(scevan_res_file[0], "CNA_mtx_relat").T
    scevan_anno = _read_rdata_frame(scevan_anno_file[0], "count_mtx_annot")
    scevan_subclone_res = None
    if subclones and len(scevan_subclones_file):
        scevan_subclone_res = _read_rdata_frame(scevan_subclones_file[0], "results.com").T

    if copy:
        adata = adata.copy()

    # add obs annotation
    if tumor_normal_call is not None:
        adata.obs[f"{key_added}_class"] = tumor_normal_call["class"]
        adata.obs[f"{key_added}_confident_normal"] = tumor_normal_call["confidentNormal"]
        if "subclone" in tumor_normal_call.columns:
            adata.obs[f"{key_added}_subclone"] = tumor_normal_call["subclone"].apply(
                lambda x: f"{int(x)}" if not pd.isnull(x) else np.nan
            )

    if subset:
        n_before = adata.n_obs
        # anndata has no public in-place subset; `inplace=True` has to modify the caller's object
        adata._inplace_subset_obs(scevan_res.index.values)
        n_dropped = n_before - adata.n_obs
        log_result(
            "read_scevan",
            "kept %d of %d cells, dropped %d that SCEVAN filtered out",
            *(adata.n_obs, n_before, n_dropped),
            dropped=n_dropped,
            total=n_before,
            remedy="pass subset=False to keep all cells (their CNV values are NaN)",
        )

    cna = scevan_res.reindex(adata.obs_names)
    if scevan_subclone_res is not None:
        # cells SCEVAN dropped are not in adata unless subset=True; columns are the same segments
        idx = scevan_subclone_res.index.intersection(cna.index)
        cna.loc[idx, :] = scevan_subclone_res.loc[idx].to_numpy()
    adata.obsm[f"X_{key_added}"] = cna.values
    adata.uns[key_added] = {"chr_pos": _get_chr_pos_from_array(scevan_anno["seqnames"])}

    return adata if copy else None

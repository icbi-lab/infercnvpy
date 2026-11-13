"""Read in result files from scevan."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rdata
from anndata import AnnData
from anndata.acc import A
from scverse_misc import deprecated_arg

from infercnvpy._acc import DEPR_INPLACE_COPY, DEPR_KEY_ADDED, _key, _legacy, _obsm_key
from infercnvpy._logging import log


def _read_rdata_frame(path: Path, name: str) -> pd.DataFrame:
    """Read object `name` from an `.RData` file as a DataFrame (keeping row/column names)."""
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


@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE_COPY)
def read_scevan(
    adata: AnnData,
    scevan_res_dir: str | Path,
    scevan_res_table: str | Path | None = None,
    *,
    out: Any = False,
    subclones: bool = True,
    subset: bool = True,
    inplace: bool | None = None,
    key_added: str | None = None,
) -> AnnData | None:
    """Load results from SCEVAN :cite:`DeFalco2021` for downstream analysis with infercnvpy.

    Requires that the cell barcodes used for SCEVAN and `adata.obs_names` match,
    but the order is irrelevant.

    Parameters
    ----------
    adata
        adata object to which the SCEVAN results shall be added
    scevan_res_table
        The results of `SCEVAN::pipelineCNA` saved as CSV file. Will add
        the columns `scevan_class`, `scevan_confident_normal` and, if SCEVAN was
        ran with subclone calling, `scevan_subclone` to `adata.obs` (named after `out`). This parameter can
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
    out
        Where to write the CNV matrix, e.g. `A.obsm["X_scevan"]` (the default).
        Its name without the `X_` prefix is used for the chromosome positions in
        `adata.uns` and as prefix of the columns added to `adata.obs`.
    inplace
        Deprecated, copy `adata` before calling this function instead.
    key_added
        Deprecated, use `out`.

    Returns
    -------
    `None`, or a modified copy of adata if the deprecated `inplace=False` is passed.
    """
    out = _legacy(out, key_added, A.obsm[f"X_{key_added}"], name="key_added", default=False)
    key = _obsm_key("scevan") if out is False else _key(out, "obsm")
    name = key.removeprefix("X_")
    scevan_res_dir = Path(scevan_res_dir)
    scevan_res_file = list(scevan_res_dir.glob("*_CNAmtx.RData"))
    scevan_subclones_file = list(scevan_res_dir.glob("*_CNAmtxSubclones.RData"))
    scevan_anno_file = list(scevan_res_dir.glob("*_count_mtx_annot.RData"))

    if len(scevan_res_file) != 1 or len(scevan_subclones_file) > 1 or len(scevan_anno_file) != 1:
        raise ValueError(
            "There must be exactely one CNA_mtx and count_mtx_annot file and at most one "
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

    if inplace is False:
        adata = adata.copy()

    # add obs annotation
    if tumor_normal_call is not None:
        adata.obs[f"{name}_class"] = tumor_normal_call["class"]
        adata.obs[f"{name}_confident_normal"] = tumor_normal_call["confidentNormal"]
        if "subclone" in tumor_normal_call.columns:
            adata.obs[f"{name}_subclone"] = tumor_normal_call["subclone"].apply(
                lambda x: f"{int(x)}" if not pd.isnull(x) else np.nan
            )

    if subset:
        adata._inplace_subset_obs(scevan_res.index.values)

    cna = scevan_res.reindex(adata.obs_names)
    if scevan_subclone_res is not None:
        # cells SCEVAN dropped are not in adata unless subset=True; columns are the same segments
        idx = scevan_subclone_res.index.intersection(cna.index)
        cna.loc[idx, :] = scevan_subclone_res.loc[idx].to_numpy()
    adata.obsm[key] = cna.values
    adata.uns[name] = {"chr_pos": _get_chr_pos_from_array(scevan_anno["seqnames"])}

    return adata if inplace is False else None

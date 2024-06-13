"""Read in result files from scevan."""

from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from anndata import AnnData
from scanpy import logging


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


def read_scevan(
    adata: AnnData,
    scevan_res_dir: str | Path,
    scevan_res_table: str | Path | None = None,
    *,
    subclones: bool = True,
    inplace: bool = True,
    subset: bool = True,
    key_added: str = "scevan",
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
        If `True`, modify anndata inplace. Otherwise, return a copy.

    Returns
    -------
    Depending on the value of `inplace` returns a modified AnnData object or returns
    `None` and modifies adata inplace.
    """
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
        logging.warning("No `scevan_res_table` specified. Will not add tumor/normal classification.")
    scevan_res = pyreadr.read_r(scevan_res_file[0])["CNA_mtx_relat"].T
    scevan_anno = pyreadr.read_r(scevan_anno_file[0])["count_mtx_annot"]
    scevan_subclone_res = None
    if subclones and len(scevan_subclones_file):
        scevan_subclone_res = pyreadr.read_r(scevan_subclones_file[0])["results.com"].T

    if not inplace:
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
        adata._inplace_subset_obs(scevan_res.index.values)

    adata.obsm[f"X_{key_added}"] = scevan_res.reindex(adata.obs_names).values
    if scevan_subclone_res is not None:
        adata.obsm[f"X_{key_added}"].loc[scevan_subclone_res.index, :] = scevan_subclone_res
    adata.uns[key_added] = {"chr_pos": _get_chr_pos_from_array(scevan_anno["seqnames"])}

    if not inplace:
        return adata

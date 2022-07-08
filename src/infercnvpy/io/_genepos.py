from os import PathLike
from typing import Literal, Union
from pathlib import Path
from anndata import AnnData
import pandas as pd
import gtfparse
from scanpy import logging
import numpy as np


def genomic_position_from_gtf(
    gtf_file: Union[Path, str],
    adata: Union[AnnData, None] = None,
    *,
    gtf_gene_id: Literal["gene_id", "gene_name"] = "gene_name",
    adata_gene_id: Union[str, None] = None,
    inplace: bool = True,
) -> Union[pd.DataFrame, None]:
    """Get genomic gene positions from a GTF file.

    The GTF file needs to match the genome annotation used for your single cell dataset.

    .. warning::
        Currently only tested with GENCODE GTFs.

    Parameters
    ----------
    gtf_file
        Path to the GTF file
    adata
        Adds the genomic positions to `adata.var`. If adata is None, returns
        a data frame with the genomic positions instead.
    gtf_gene_id
        Use this GTF column to match it to anndata
    adata_gene_id
        Match this column to the gene ids from the GTF file. Default: use
        adata.var_names.
    inplace
        If True, add the annotations directly to adata, otherwise return a dataframe.
    """
    gtf = gtfparse.read_gtf(
        gtf_file, usecols=["seqname", "feature", "start", "end", "gene_id", "gene_name"]
    )
    gtf = (
        gtf.loc[
            gtf["feature"] == "gene",
            ["seqname", "start", "end", "gene_id", "gene_name"],
        ]
        .drop_duplicates()
        .rename(columns={"seqname": "chromosome"})
    )

    gene_ids_adata = (
        adata.var_names if adata_gene_id is None else adata.var[adata_gene_id]
    ).values
    gtf = gtf.loc[gtf[gtf_gene_id].isin(gene_ids_adata), :]

    missing_from_gtf = len(set(gene_ids_adata) - set(gtf[gtf_gene_id].values))
    if missing_from_gtf:
        logging.warning(
            f"GTF file misses annotation for {missing_from_gtf} genes in adata."
        )

    duplicated_symbols = np.sum(gtf["gene_name"].duplicated())
    if duplicated_symbols:
        logging.warning(
            f"Skipped {duplicated_symbols} genes because of duplicate identifiers in GTF file."
        )
        gtf = gtf.loc[~gtf[gtf_gene_id].duplicated(keep=False), :]

    tmp_var = adata.var.copy()
    orig_index_name = tmp_var.index.name
    TMP_INDEX_NAME = "adata_var_index"
    tmp_var.index.name = TMP_INDEX_NAME
    tmp_var.reset_index(inplace=True)
    var_annotated = tmp_var.merge(
        gtf,
        how="left",
        left_on=TMP_INDEX_NAME if adata_gene_id is None else adata_gene_id,
        right_on=gtf_gene_id,
        validate="one_to_one",
    )
    var_annotated.set_index(TMP_INDEX_NAME, inplace=True)
    var_annotated.index.name = orig_index_name

    if inplace:
        adata.var = var_annotated
    else:
        return var_annotated

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import scanpy.queries
from anndata import AnnData
from scanpy import logging


def genomic_position_from_biomart(
    adata: AnnData | None = None,
    *,
    adata_gene_id: str | None = None,
    biomart_gene_id="ensembl_gene_id",
    species: str = "hsapiens",
    inplace: bool = True,
    **kwargs,
):
    """
    Get genomic gene positions from ENSEMBL Biomart.

    Parameters
    ----------
    adata
        Adds the genomic positions to `adata.var`. If adata is None, returns
        a data frame with the genomic positions instead.
    adata_gene_id
        Column in `adata.var` that contains (ENSMBL) gene IDs. If not specified,
        use `adata.var_names`.
    biomart_gene_id
        The biomart column to use as gene identifier. Typically this would be `ensembl_gene_id` or `hgnc_symbol`,
        but could be different for other species.
    inplace
        If True, add the annotations directly to adata, otherwise return a dataframe.
    **kwargs
        passed on to :func:`scanpy.queries.biomart_annotations`
    """
    biomart_annot = (
        scanpy.queries.biomart_annotations(
            species,
            [
                biomart_gene_id,
                "start_position",
                "end_position",
                "chromosome_name",
            ],
            **kwargs,
        )
        .rename(
            columns={
                "start_position": "start",
                "end_position": "end",
                "chromosome_name": "chromosome",
            }
        )
        # use chr prefix for chromosome
        .assign(chromosome=lambda x: "chr" + x["chromosome"])
    )

    gene_ids_adata = (adata.var_names if adata_gene_id is None else adata.var[adata_gene_id]).values
    missing_from_biomart = len(set(gene_ids_adata) - set(biomart_annot[biomart_gene_id].values))
    if missing_from_biomart:
        logging.warning(
            f"Biomart misses annotation for {missing_from_biomart} genes in adata. Did you use ENSEMBL ids?"
        )

    duplicated_symbols = np.sum(biomart_annot[biomart_gene_id].duplicated())
    if duplicated_symbols:
        logging.warning(f"Skipped {duplicated_symbols} genes because of duplicate identifiers in GTF file.")
        biomart_annot = biomart_annot.loc[~biomart_annot[biomart_gene_id].duplicated(keep=False), :]

    tmp_var = adata.var.copy()
    orig_index_name = tmp_var.index.name
    TMP_INDEX_NAME = "adata_var_index"
    tmp_var.index.name = TMP_INDEX_NAME
    tmp_var.reset_index(inplace=True)
    var_annotated = tmp_var.merge(
        biomart_annot,
        how="left",
        left_on=TMP_INDEX_NAME if adata_gene_id is None else adata_gene_id,
        right_on=biomart_gene_id,
        validate="one_to_one",
    )
    var_annotated.set_index(TMP_INDEX_NAME, inplace=True)
    var_annotated.index.name = orig_index_name

    if inplace:
        adata.var = var_annotated
    else:
        return var_annotated


def genomic_position_from_gtf(
    gtf_file: Path | str,
    adata: AnnData | None = None,
    *,
    gtf_gene_id: Literal["gene_id", "gene_name"] = "gene_name",
    adata_gene_id: str | None = None,
    inplace: bool = True,
) -> pd.DataFrame | None:
    """
    Get genomic gene positions from a GTF file.

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
    try:
        import gtfparse
    except ImportError:
        raise ImportError(
            "genomic_position_from_gtf requires gtfparse as optional dependency. Please install it using `pip install gtfparse`."
        ) from None
    gtf = gtfparse.read_gtf(
        gtf_file, usecols=["seqname", "feature", "start", "end", "gene_id", "gene_name"]
    ).to_pandas()
    gtf = (
        gtf.loc[
            gtf["feature"] == "gene",
            ["seqname", "start", "end", "gene_id", "gene_name"],
        ]
        .drop_duplicates()
        .rename(columns={"seqname": "chromosome"})
    )
    # remove ensembl versions
    gtf["gene_id"] = gtf["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    gene_ids_adata = (adata.var_names if adata_gene_id is None else adata.var[adata_gene_id]).values
    gtf = gtf.loc[gtf[gtf_gene_id].isin(gene_ids_adata), :]

    missing_from_gtf = len(set(gene_ids_adata) - set(gtf[gtf_gene_id].values))
    if missing_from_gtf:
        logging.warning(f"GTF file misses annotation for {missing_from_gtf} genes in adata.")

    duplicated_symbols = np.sum(gtf["gene_name"].duplicated())
    if duplicated_symbols:
        logging.warning(f"Skipped {duplicated_symbols} genes because of duplicate identifiers in GTF file.")
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

    # if not a gencode GTF, let's add 'chr' prefix:
    if np.all(~var_annotated["chromosome"].dropna().str.startswith("chr")):
        var_annotated["chromosome"] = "chr" + var_annotated["chromosome"]

    if inplace:
        adata.var = var_annotated
    else:
        return var_annotated

import os
from multiprocessing import cpu_count

import pandas as pd
from anndata import AnnData
from scanpy import logging
from scipy.sparse import issparse


def copykat(
    adata: AnnData,
    gene_ids: str = "S",
    organism: str = "human",
    segmentation_cut: float = 0.1,
    distance: str = "euclidean",
    s_name: str = "copykat_result",
    min_genes_chr: int = 5,
    key_added: str = "cnv",
    inplace: bool = True,
    layer: str | None = None,
    n_jobs: int | None = None,
    norm_cell_names: str = "",
    cell_line="no",
    window_size=25,
) -> (pd.DataFrame, pd.Series):
    """Inference of genomic copy number and subclonal structure.

    Runs CopyKAT (Copynumber Karyotyping of Tumors) :cite:`Gao2021` based on integrative
    Bayesian approaches to identify genome-wide aneuploidy at 5MB resolution
    in single cells to separate tumor cells from normal cells, and tumor
    subclones using high-throughput sc-RNAseq data.

    Note on input data from the original authors:

        The matrix values are often the count of unique molecular identifier (UMI)
        from nowadays high througput single cell RNAseq data. The early generation of
        scRNAseq data may be summarized as TPM values or total read counts,
        which should also work.

    This means that unlike for :func:`infercnvpy.tl.infercnv` the input data
    should not be log-transformed.

    CopyKAT also does NOT require running :func:`infercnvpy.io.genomic_position_from_gtf`,
    it infers the genomic position from the gene symbols in `adata.var_names`.

    You can find more info on GitHub: https://github.com/navinlabcode/copykat

    Parameters
    ----------
    adata
        annotated data matrix
    key_added
        Key under which the copyKAT scores will be stored in `adata.obsm` and `adata.uns`.
    inplace
        If True, store the result in adata, otherwise return it.
    layer
        AnnData layer to use for running copykat
    gene_ids
        gene id type: Symbol ("S") or Ensemble ("E").
    segmentation_cut
        segmentation parameters, input 0 to 1; larger looser criteria.
    distance
        distance methods include "euclidean", and correlation coverted distance include "pearson" and "spearman".
    s_name
        sample (output file) name.
    min_genes_chr
        minimal number of genes per chromosome for cell filtering.
    norm_cell_names:
        cell barcodes (`adata.obs.index`) indicate normal cells
    n_jobs
        Number of cores to use for copyKAT analysis. Per default, uses all cores
        available on the system. Multithreading does not work on Windows and this
        value will be ignored.
    organism
        Runs methods for calculating copy numbers from: "human" or "mouse" scRNAseq data (default: "human")
    cell_line
        if the data are from pure cell line (ie. not a mixture of tumor and normal), put "yes" to use a synthetic baseline (default: "no")
    window_size
        Sets a minimal window size for segmentation

    Returns
    -------
    Depending on the value of `inplace`, either returns `None` or a tuple (`CNV Matrix`,`CopyKat prediction`)
    """
    if n_jobs is None:
        n_jobs = cpu_count()
    if os.name != "posix":
        n_jobs = 1

    try:
        from rpy2 import robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr
    except ImportError:
        raise ImportError("copyKAT requires rpy2 to be installed. ") from None

    try:
        importr("copykat")
        importr("stringr")
    except ImportError:
        raise ImportError(
            "copyKAT requires a valid R installation with the following packages: copykat, stringr"
        ) from None

    logging.info("Preparing R objects")
    with localconverter(ro.default_converter + numpy2ri.converter):
        expr = adata.X if layer is None else adata.layers[layer]
        if issparse(expr):
            expr = expr.T.toarray()
        else:
            expr = expr.T
        ro.globalenv["expr_r"] = ro.conversion.py2rpy(expr)
    ro.globalenv["gene_names"] = ro.conversion.py2rpy(list(adata.var.index))
    ro.globalenv["cell_IDs"] = ro.conversion.py2rpy(list(adata.obs.index))
    ro.globalenv["n_jobs"] = ro.conversion.py2rpy(n_jobs)
    ro.globalenv["gene_ids"] = ro.conversion.py2rpy(gene_ids)
    ro.globalenv["segmentation_cut"] = ro.conversion.py2rpy(segmentation_cut)
    ro.globalenv["distance"] = ro.conversion.py2rpy(distance)
    ro.globalenv["s_name"] = ro.conversion.py2rpy(s_name)
    ro.globalenv["min_gene_chr"] = ro.conversion.py2rpy(min_genes_chr)
    ro.globalenv["norm_cell_names"] = ro.conversion.py2rpy(norm_cell_names)
    ro.globalenv["organism"] = ro.conversion.py2rpy(organism)
    ro.globalenv["cell_line"] = ro.conversion.py2rpy(cell_line)
    ro.globalenv["window_size"] = ro.conversion.py2rpy(window_size)

    logging.info("Running copyKAT")
    ro.r(
        """
        rownames(expr_r) <- gene_names
        colnames(expr_r) <- cell_IDs
        if (organism == "mouse"){
            copyKAT_run <- copykat(rawmat = expr_r, id.type = gene_ids, ngene.chr = min_gene_chr, win.size = window_size,
                                KS.cut = segmentation_cut, sam.name = s_name, distance = distance, norm.cell.names = norm_cell_names, cell.line = cell_line,
                                n.cores = n_jobs, output.seg = FALSE, genome = 'mm10')
        } else {
            copyKAT_run <- copykat(rawmat = expr_r, id.type = gene_ids, ngene.chr = min_gene_chr, win.size = window_size,
                                KS.cut = segmentation_cut, sam.name = s_name, distance = distance, norm.cell.names = norm_cell_names, cell.line = cell_line,
                                n.cores = n_jobs, output.seg = FALSE)
        }
        copyKAT_result <- data.frame(copyKAT_run$CNAmat)
        colnames(copyKAT_result) <- str_replace_all(colnames(copyKAT_result), "\\\\.", "-")
        copyKAT_pred <- data.frame(copyKAT_run$prediction)
        if(dim(copyKAT_result)[2] != length(cell_IDs)){
            missing_cells <- setdiff(cell_IDs,colnames(copyKAT_result))
            na_mtrx <- data.frame(matrix(ncol=length(missing_cells),nrow=nrow(copyKAT_result)))
            new_colnames <- c(colnames(copyKAT_result),missing_cells)
            copyKAT_result <- cbind(copyKAT_result,na_mtrx)
            colnames(copyKAT_result) <- new_colnames
        }
        """
    )

    with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
        copyKAT_result = ro.conversion.rpy2py(ro.globalenv["copyKAT_result"])
        copyKAT_pred = ro.conversion.rpy2py(ro.globalenv["copyKAT_pred"])

    chrom_pos = {
        "chr_pos": {
            f"chr{chrom}": int(pos) for pos, chrom in copyKAT_result.loc[:, ["chrom"]].drop_duplicates().itertuples()
        }
    }

    # Drop cols
    new_cpkat = copyKAT_result.drop(["chrom", "chrompos", "abspos"], axis=1)
    # align cells
    new_cpkat = new_cpkat.loc[:, adata.obs.index]
    copyKAT_pred = adata.obs.merge(copyKAT_pred, left_index=True, right_index=True, how="left")["copykat.pred"]
    # transpose
    new_cpkat_trans = new_cpkat.T.values

    if inplace:
        adata.uns[key_added] = chrom_pos
        adata.obsm[f"X_{key_added}"] = new_cpkat_trans
        adata.obs[key_added] = copyKAT_pred
    else:
        return new_cpkat_trans, copyKAT_pred

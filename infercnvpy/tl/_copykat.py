import infercnvpy as cnv
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import issparse
from anndata import AnnData
from scanpy import logging, AnnData

def copykat(
    adata: AnnData,
    key_added: str = "copyKAT",
    inplace: bool = True,
    layer: str = None
)-> pd.DataFrame:
    """Inference of genomic copy number and subclonal structure.

    Runs CopyKAT (Copynumber Karyotyping of Tumors) based on sing integrative 
    Bayesian approaches to identify genome-wide aneuploidy at 5MB resolution 
    in single cells to separate tumor cells from normal cells, and tumor 
    subclones using high-throughput sc-RNAseq data.
    
    Note: The matrix values are often the count of unique molecular identifier (UMI) 
    from nowadays high througput single cell RNAseq data. The early generation of 
    scRNAseq data may be summarized as TPM values or total read counts, 
    which should also work.

    You can find more info on GitHub: https://github.com/navinlabcode/copykat

    :cite:`Gao2021`

    Parameters
    ----------
    adata
        annotated data matrix
    key_added
        Key under which the copyKAT scores will be stored in `adata.obs`.
    inplace
        If True, store the result in adata, otherwise return it.

    Returns
    -------
    Depending on the value of `inplace`, either returns `None` or a vector
    with scores.
    """
    try:
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2 import robjects as ro
    except ImportError:
        raise ImportError("copyKAT requires rpy2 to be installed. ")

    try:
        copyKAT = importr("copykat")
    except ImportError:
        raise ImportError(
            "copyKAT requires a valid R installation with the following packages: "
            "copykat"
        )
    
    logging.info("Preparing R objects")
    with localconverter(ro.default_converter + numpy2ri.converter):
        expr = adata.X if layer is None else tmp_adata.layers[layer]
        if issparse(expr):
            expr = expr.T.toarray()
        else:
            expr = expr.T
        ro.globalenv["expr_r"] = ro.conversion.py2rpy(expr)
    ro.globalenv["gene_names"] = ro.conversion.py2rpy(list(adata.var.index))
    ro.globalenv["cell_IDs"] = ro.conversion.py2rpy(list(adata.obs.index))
    
    logging.info("Running copyKAT")
    ro.r(
        """
        rownames(expr_r) <- gene_names
        colnames(expr_r) <- cell_IDs
        copyKAT_run <- copykat(rawmat = expr_r, id.type = "S", ngene.chr = 5, win.size = 25, 
                                KS.cut = 0.1, sam.name = "test", distance = "euclidean", norm.cell.names = "", 
                                n.cores = 12, output.seg = FALSE)
        copyKAT_result <- copyKAT_run$CNAmat
        colnames(copyKAT_result) <- c('chrom', 'chrompos', 'abspos', cell_IDs)
        """
    )

    with localconverter(
        ro.default_converter + numpy2ri.converter + pandas2ri.converter
    ):
        copyKAT_result = ro.conversion.rpy2py(ro.globalenv["copyKAT_result"])

    adata.uns[key_added] = {
        "chr_pos": {
            f"chr{chrom}": int(pos)
            for pos, chrom in copyKAT_result.loc[:, ["chrom"]].drop_duplicates().itertuples()
        }
    }

    # Drop cols
    new_cpkat = copyKAT_result.drop(["chrom", "chrompos", "abspos"], axis=1).values

    # transpose
    new_cpkat_trans = new_cpkat.T

    if inplace:
        adata.obs["X_%s" % key_added] = new_cpkat_trans
    else:
        return new_cpkat_trans
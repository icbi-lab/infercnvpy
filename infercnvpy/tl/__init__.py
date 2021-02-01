from typing import Union, Sequence
import numpy as np
import scipy.stats


def leiden(adata, *args, **kwargs):
    pass


def pca(adata):
    pass


def umap(adata):
    pass


# def infercnv(
#     adata,
#     *,
#     reference_key: Union[str, None] = None,
#     reference_cat: Union[str, Sequence[str]] = None,
#     reference: np.ndarray = None,
# ):
#     """
#     Infer Copy Number Variation (CNV) by averaging gene expression over genomic
#     regions.

#     adata should already be filtered for low-quality cells. adata.X needs to be
#     normalized and log-transformed. The method should be
#     fairly robust to different normalization methods (`normalize_total`, scran, etc.).

#     The method requires a "reference" value to which the expression of genomic
#     regions is compared. If your dataset contains different cell types and includes
#     both tumor and normal cells, the average of all cells can be used as reference.
#     This is the default.

#     If you already know which cells are "normal", you can provide a column
#     from adata.obs to `reference_key` that contains the annotation. `reference_cat`
#     specifies one or multiple values in `reference_key` that refer to normal cells.

#     Alternatively, you can specify a numpy array with average gene expression values
#     (for instance, derived from a different dataset), to `reference` which will be used
#     as reference instead.


#     Parameters
#     ----------
#     adata
#         annotated data matrix
#     reference_key
#         Column name in adata.obs that contains tumor/normal annotations.
#         If this is set to None, the average of all cells is used as reference.
#     reference_cat
#         One or multiple values in `adata.obs[reference_key]` that annotate
#         normal cells.


#     """
#     pass

#     # Follows approximately the inferCNV website:
#     # https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV#infercnv-step-by-step-exploratory-execution

#     # Step 1 - compute running mean per chromosome

#     # Step 2 -

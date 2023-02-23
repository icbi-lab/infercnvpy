import numpy as np
import infercnvpy as cnv

def subclonal_clustering(
    adata,
    reference_key,
    reference_cat,
    obs_key_used,
    obsm_key_added,
    inplace : bool = True,
):
    
    adata_malignant_subset_indicator = adata.obs[reference_key].isin(reference_cat) == False
    malignant_indices = np.argwhere(adata_malignant_subset_indicator.to_numpy() == True).reshape(-1)
    adata_malignant_subset = adata[adata.obs[reference_key].isin(reference_cat) == False]
    cnv.tl.pca(adata_malignant_subset)
    cnv.pp.neighbors(adata_malignant_subset)
    cnv.tl.leiden(adata_malignant_subset, key_added = "subclone_cnv_leiden")

    subclone_labels = -1 * np.ones(len(adata.obs), dtype=np.int32)
    subclone_labels[malignant_indices] = adata_malignant_subset.obs['subclone_cnv_leiden'].to_numpy()
    if inplace:
        adata.obs[obsm_key_added] = subclone_labels
    else:
        return subclone_labels

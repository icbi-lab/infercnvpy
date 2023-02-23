from typing import Tuple, Union, Sequence
import numpy as np
from scanpy import logging
from anndata import AnnData
import scipy.ndimage
import scipy.sparse
import matplotlib.pyplot as plt
from scipy.stats import norm
from hmmlearn import hmm
from scipy.sparse import issparse
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import anndata

def get_hmm_model(means : np.array, 
                  variances : np.array,
                  HMM_transition_prob : float = 1e-6,
                  learned : bool = True, 
                  n_iter : int = 100) -> hmm.GaussianHMM:
    assert means.size == 3
    assert variances.size == 3
    assert HMM_transition_prob >= 0
    assert HMM_transition_prob <= 1
    if learned:
        model = hmm.GaussianHMM(n_components=3, 
                                covariance_type="full", 
                                init_params='t', 
                                params='t',
                                n_iter=n_iter)
        model.means_ = means.reshape(3,1)
        return model
    else:
        t = HMM_transition_prob
        model = hmm.GaussianHMM(n_components=3,
                                covariance_type='full',
                                n_iter=100) 
        model.transmat_ = np.array([[1-2*t,  t,     t],
                                    [t,     1-2*t,  t],
                                    [t,      t,     1-2*t]])
        model.startprob_ = np.array([t,     1-2*t,  t])
        model.means_ = means.reshape(3,1)
        model.covars_ = variances.reshape((3,1,1))
        return model

def std_helper(A : Union[np.ndarray, scipy.sparse.csr.csr_matrix]) -> float:
    # numpy
    if type(A) == np.ndarray:
        return A.std()
    elif type(A) == scipy.sparse.csr.csr_matrix:
        return np.sqrt((A.power(2)).mean() - A.mean()**2)
    else:
        raise Exception('Unknown type', type(A))

def get_state_parameters_infercnv(
        nonmalignant_cnv : Union[np.ndarray, scipy.sparse.csr.csr_matrix],
        p_val : float = 0.01, 
        halve_mean_delta : bool = True) -> Tuple[np.array, np.array]:
    ref_mean = nonmalignant_cnv.mean()
    ref_sigma = std_helper(nonmalignant_cnv)
    mean_delta = np.abs(norm.ppf(p_val, loc=0, scale=ref_sigma))
    # Divide mean_delta in half to get similar mean_delta to InferCNV - should ideally not be there
    if halve_mean_delta:
        mean_delta /= 2
    sorted_means = np.array([ref_mean - mean_delta, ref_mean, ref_mean + mean_delta]).reshape(-1,1)
    variance = ref_sigma ** 2
    sorted_variances = np.array([variance, variance, variance]).reshape(3,1,1)
    return sorted_means, sorted_variances

def run_hmm_transition_not_learned(adata, subclone_key, means, variances, key_used, log1p_base='e', exponentiate=False): # p_val=0.01, log1p_base='e', trans_prob=1e-6):
    # Make matrix of CNV values
    X = get_X_from_adata(adata, key_used=key_used)
    if exponentiate:
        X = np.exp(X) if log1p_base == 'e' else np.power(log1p_base, np.array(np.asarray(adata.obsm['X_cnv'].todense())))

    subclone_names = np.unique(adata.obs[subclone_key].to_numpy())
    # Initialize matrix that will have 0,1,2 for del, neutral, gain hidden states
    Z_matrix = -1 * np.ones_like(X)
    chr_pos_dict = dict(sorted(adata.uns[key_used]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())
    # List of tuples (chr_start, chr_end)
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [X.shape[1]]))
    model = get_hmm_model(means, variances, learned=False) # mu,sigma,mean_delta,HMM_transition_prob=trans_prob)

    for subclone_name in subclone_names:
        # Which cells are in this subclone?
        cluster_idx = np.argwhere(adata.obs[subclone_key].to_numpy() == subclone_name).reshape(-1)
        cluster_expr = X[cluster_idx]
        # Observed values for hmm is the mean CNV values across the cells in the subclone
        hmm_X = cluster_expr.mean(axis=0).reshape(-1,1)
        # Iterate over each chromosome
        for i, (chr_start_pos, chr_end_pos) in enumerate(var_group_positions):
            # The gene indices for this chromosome
            chr_idx = np.arange(chr_start_pos, chr_end_pos)
            # Get the most likely state for each CNV value
            Z_d = model.predict(hmm_X[chr_idx])
            # Assign the most likely states to all the cells in the subclone
            for cluster_id in cluster_idx:
                assert Z_matrix[cluster_id, chr_idx].mean() == -1., Z_matrix[cluster_id, chr_idx]
                Z_matrix[cluster_id, chr_idx] = Z_d

    return Z_matrix

def get_X_from_adata(adata: AnnData, key_used : str = 'cnv') -> Union[np.ndarray, scipy.sparse.csr.csr_matrix]:
    if type(adata.obsm[f'X_{key_used}']) == anndata._core.views.ArrayView or type(adata.obsm[f'X_{key_used}']) == np.ndarray:
        X = np.array(adata.obsm[f'X_{key_used}'])
    elif type(adata.obsm[f'X_{key_used}']) == scipy.sparse.csr.csr_matrix or type(adata.obsm[f'X_{key_used}']) == anndata._core.views.SparseCSRView:
        X = np.array(adata.obsm[f'X_{key_used}'].todense())
    else:
        raise ValueError(f'X_{key_used} is not a numpy array or a scipy sparse matrix')
    return X

def hmm_denoising(
    adata: AnnData,
    *,
    reference_key: str,
    reference_cat: Union[None, str, Sequence[str]] = None,
    subclone_key='cell_type',
    key_used: str = 'cnv',
    key_added: str = 'hmm',
    iterations: int = 1,
    exponentiate: bool = False,
    inplace : bool = True
) -> Union[None, np.ndarray]:
    """Applies HMM for denoising of CNV value matrix.

    Parameters:
    adata
        annotated data matrix
    reference_key
        Column name in adata.obs that contains tumor/normal annotations.
        used for grouping cells into malignant and non-malignant.
    reference_cat
        One or multiple values in `adata.obs[reference_key]` that annotate
        normal cells.
    subclone_key
        Column in adata.obs by which malgnant cells are grouped into subclones. HMM is applied
        each subclone, where the observed sequence is the sequence of mean CNV values across cells in the subclone.
    key_used
        Key under which the output from tl.infercnv (non-denoised CNV values) is stored in adata.
    key_added
        Key under which the denoised matrix (hidden states inferred by HMM) will be stored in adata if `inplace=True`.
        Will store the matrix in `adata.obsm["X_{key_added}"] and additional information
        in `adata.uns[key_added]`.
    iterations
        If iterations is 1, the HMM transmission probabilies (means and variances of the three Gaussians) are estimated only once
        and HMM is applied once.
        If more than one iteration, the means and variances for the Gaussians are reestimated at the end of each iteration based on the
        hidden states inferred by the HMM, stopping at convergence or number of iterations specified.
    inplace
        If True, save the results in adata.obsm, otherwise return the CNV matrix.

    Returns:
    int:Returning value

    """
    
    assert iterations > 0, 'iterations must be greater than 0'

    # Get the malignant cells as matrix
    adata_nonmalignant_subset = adata[adata.obs[reference_key].isin(reference_cat) == True]
    # Get the malignant cells as matrix
    nonmalignant_cnv = get_X_from_adata(adata_nonmalignant_subset, key_used=key_used)
    # Find the means for three CNV states using K-means
    means, variances = get_state_parameters_infercnv(nonmalignant_cnv)
    # Z-matrix, one row per cell, one column per gene, values are 0,1,2 for del, neutral, gain
    Z_matrix = run_hmm_transition_not_learned(adata = adata,
                                                        subclone_key = subclone_key,
                                                        means = means,
                                                        variances = variances,
                                                        key_used=key_used,
                                                        exponentiate=exponentiate)
    converged = False
    for i in range(1,iterations):
        new_means = np.array([adata.obsm[f'X_{key_used}'][Z_matrix == 0.].mean(), adata.obsm[f'X_{key_used}'][Z_matrix == 1.].mean(), adata.obsm[f'X_{key_used}'][Z_matrix == 2.].mean()]).reshape(-1,1)
        new_variances = np.array([adata.obsm[f'X_{key_used}'][Z_matrix == 0.].std() ** 2, adata.obsm[f'X_{key_used}'][Z_matrix == 1.].std() ** 2, adata.obsm[f'X_{key_used}'][Z_matrix == 2.].std() ** 2]).reshape(-1,1,1)
        Z_matrix_new = run_hmm_transition_not_learned(adata = adata,
                                                                    subclone_key = subclone_key,
                                                                    means = new_means,
                                                                    variances= new_variances,
                                                                    key_used=key_used,
                                                                    exponentiate=exponentiate)
        if(np.array_equal(Z_matrix, Z_matrix_new)):
            converged = True
            break
        else:
            Z_matrix = Z_matrix_new

    # chromosome heatmap can only be displayed if both postive and negative CNV values
    # are present - or equivalently both positive and negative CNV states inferred
    if Z_matrix.min() != 0.:
        print('WARNING, found no deletions, adding a gain on gene 0 of cell 0 for visualization to work')
        Z_matrix[0][0] = 0.
    if Z_matrix.max() != 2.:
        print('WARNING, found no gains, adding a gain on gene 1 of cell 0 for visualization to work')
        Z_matrix[0][1] = 2.

    if inplace:
        # Subtract 1 from Z_matrix to get -1,0,1 for del, neutral, gain
        adata.obsm[f'X_{key_added}'] = Z_matrix - 1
        # For making pl.chromosome_heatmap work, we copy the chromosome positions from the cnv key
        adata.uns[key_added] = adata.uns[key_used]
    else:
        return Z_matrix
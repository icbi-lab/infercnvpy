from typing import Tuple, Union, Sequence
import numpy as np
from anndata import AnnData
import scipy.ndimage
import scipy.sparse
from scipy.stats import norm
from hmmlearn import hmm
from scipy.sparse import issparse
from scipy import sparse
import warnings


def get_X_from_adata(adata: AnnData, key_used: str = 'cnv') -> np.ndarray:
    return adata.obsm[f'X_{key_used}']


def std_helper(A: Union[np.ndarray, scipy.sparse.spmatrix]) -> float:
    """
    Compute the standard deviation of a matrix or sparse matrix

    Args:
        A: matrix or sparse matrix

    Returns:
        float: standard deviation
    """
    if issparse(A):
        return np.sqrt((A.power(2)).mean() - A.mean()**2)
    else:
        return A.std()


def get_hmm_model(means: np.array,
                  variances: np.array,
                  HMM_transition_prob: float = 1e-6,
                  n_iter: int = 100) -> hmm.GaussianHMM:
    """
    Initialize a Gaussian HMM model with 3 hidden states

    Args:
        means: means of the 3 hidden states
        variances: variances of the 3 hidden states
        HMM_transition_prob: Probability of transitioning to a different state for each different state
        n_iter: Number of iterations to run the inference algorithm

    Returns:
        hmm.GaussianHMM: Initialized Gaussian HMM model
    """
    assert means.size == 3
    assert variances.size == 3
    assert HMM_transition_prob >= 0
    assert HMM_transition_prob <= 1

    t = HMM_transition_prob
    model = hmm.GaussianHMM(n_components=3,
                            covariance_type='full',
                            n_iter=n_iter)
    model.transmat_ = np.array([[1 - 2 * t, t, t],
                                [t, 1 - 2 * t, t],
                                [t, t, 1 - 2 * t]])
    model.startprob_ = np.array([t, 1 - 2 * t, t])
    model.means_ = means.reshape(3, 1)
    model.covars_ = variances.reshape((3, 1, 1))
    return model


def get_state_parameters_qnorm(
        nonmalignant_cnv: Union[np.ndarray, scipy.sparse.spmatrix],
        p_val: float = 0.01) -> Tuple[np.array, np.array]:
    """
    Compute the mean and variance of the 3 hidden states based on statistics of the non-malignant cells,
    using the inverse normal CDF.

    Args:
        nonmalignant_cnv: matrix of CNV values for non-malignant cells
        p_val: p-value to use for the inverse normal CDF

    Returns:
        np.ndarray: means of the 3 hidden states
        np.ndarray: variances of the 3 hidden states
    """
    mean = nonmalignant_cnv.mean()
    sigma = std_helper(nonmalignant_cnv)
    variance = sigma ** 2
    mean_delta = np.abs(norm.ppf(p_val, loc=0, scale=sigma))
    # Divide mean_delta in half to get similar mean_delta to InferCNV - should ideally not be there
    mean_delta /= 2
    means = np.array([mean - mean_delta, mean, mean + mean_delta]).reshape(-1, 1)
    variances = np.array([variance, variance, variance]).reshape(3, 1, 1)
    return means, variances


def Z_dict_to_dense(adata, subclone_key, key_used, Z_dict) -> np.ndarray:
    """
    Compute a full matrix of hidden states from a dictionary of hidden states for each subclone

    Args:
        adata: AnnData object
        subclone_key: key in adata.obs that contains the subclone names
        key_used: Key under which the output from tl.infercnv (non-denoised CNV values) is stored in adata.
        Z_dict: Dictionary of hidden states for each subclone
    Returns:
        float: standard deviation
    """
    X = get_X_from_adata(adata, key_used=key_used)
    subclone_names = np.unique(adata.obs[subclone_key].to_numpy())
    Z_matrix = -2 * np.ones(X.shape, dtype=np.int32)
    for subclone_name in subclone_names:
        cluster_idx = np.argwhere(adata.obs[subclone_key].to_numpy() == subclone_name).reshape(-1)
        Z_t_dense = Z_dict[subclone_name].toarray()
        for cluster_id in cluster_idx:
            Z_matrix[cluster_id] = Z_t_dense
    assert np.all(Z_matrix != -2), 'Z_matrix not filled'
    return Z_matrix


def apply_hmm(adata: AnnData,
              subclone_key: str,
              model: hmm.GaussianHMM,
              key_used: str) -> dict:
    """
    Run the HMM model on each subclone, where the sequence of observations for a subclone is attained by averaging across cells

    Args:
        adata: AnnData object
        subclone_key: key in adata.obs that contains the subclone names
        model: HMM model to apply
        key_used: Key under which the output from tl.infercnv (non-denoised CNV values) is stored in adata.
    Returns:
        dict: a dictionary of subclone names to the hidden states for each cell in the subclone stored as a sparse array
    """

    # Make matrix of CNV values
    X = get_X_from_adata(adata, key_used=key_used)
    # Get the names of all distinct subclones
    subclone_names = np.unique(adata.obs[subclone_key].to_numpy())
    chr_pos_dict = dict(sorted(adata.uns[key_used]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())
    # List of tuples (chr_start, chr_end)
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [X.shape[1]]))

    subclones_Z_dict = {}
    for subclone_name in subclone_names:
        # Which cells are in this subclone?
        cluster_idx = np.argwhere(adata.obs[subclone_key].to_numpy() == subclone_name).reshape(-1)
        # Matrix of CNV values for only the cells in this subclone
        cluster_expr = X[cluster_idx]
        # Observed values for hmm is the mean CNV values across the cells in the subclone
        hmm_X = cluster_expr.mean(axis=0).reshape(-1, 1)
        # Iterate over each chromosome
        hmm_Z = -2 * np.ones(hmm_X.shape, dtype=np.int32).reshape(-1)
        for _, (chr_start_pos, chr_end_pos) in enumerate(var_group_positions):
            # The gene indices for this chromosome
            # this is to avoid cases where they are not consecutive, but probably not needed
            chr_idx = np.arange(chr_start_pos, chr_end_pos)
            # Observed values for this chromosome
            X_t = hmm_X[chr_idx]
            # Get the most likely state for each CNV value
            Z_t = model.predict(X_t).astype(np.int32) - 1
            # Assign the most likely states to all the cells in the subclone
            hmm_Z[chr_idx] = Z_t
        assert hmm_Z.min() != -2, 'Did not predict hidden state for all genes'
        hmm_Z_sparse = sparse.csr_matrix(hmm_Z)
        subclones_Z_dict[subclone_name] = hmm_Z_sparse

    return subclones_Z_dict


def hmm_denoising(
    adata: AnnData,
    reference_key: str,
    reference_cat: Union[None, str, Sequence[str]] = None,
    subclone_key='cell_type',
    key_used: str = 'cnv',
    key_added: str = 'hmm',
    p_val: float = 0.01,
    iterations: int = 1,
    inplace: bool = True,
    transition_prob: float = 1e-6
) -> Tuple[bool, np.ndarray, dict]:
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
    p_val
        p-value used for determining the means for the three Gaussians used in the HMM.
    HMM_transition_prob
        Parameter for transition matrix of HMM. HMM_transition_prob: Probability of transitioning to a different state (for each different state), so half the probability of moving to another state
    iterations
        If iterations is 1, the HMM transmission probabilies (means and variances of the three Gaussians) are estimated only once
        and HMM is applied once.
        If more than one iteration, the means and variances for the Gaussians are reestimated at the end of each iteration based on the
        hidden states inferred by the HMM, stopping at convergence or number of iterations specified.
    inplace
        If True, save the results in adata.obsm, otherwise return the CNV matrix.

    Returns:
        converged (bool): convergence status, only relevant when multiple iterations are used
        Z_matrix (np.ndarray): hidden states inferred by HMM as a dense matrix
        Z_dict (dict): hidden states inferred by HMM as a dictionary of sparse matrices, one for each subclone
    """
    warnings.filterwarnings('ignore', category=FutureWarning)

    assert iterations > 0, 'iterations must be greater than 0'

    # Get the malignant cells
    adata_nonmalignant_subset = adata[adata.obs[reference_key].isin(reference_cat) == True]
    # Get the malignant cells as matrix
    nonmalignant_cnv = get_X_from_adata(adata_nonmalignant_subset, key_used=key_used)
    # Find the means for three CNV states using K-means
    means, variances = get_state_parameters_qnorm(nonmalignant_cnv, p_val=p_val)
    # Construct HMM model
    model = get_hmm_model(means, variances, HMM_transition_prob=transition_prob)

    # Z-matrix, one row per cell, one column per gene, values are 0,1,2 for del, neutral, gain
    Z_dict = apply_hmm(adata=adata,
                       subclone_key=subclone_key,
                       model=model,
                       key_used=key_used)
    Z_matrix = Z_dict_to_dense(adata, subclone_key, key_used, Z_dict)

    converged = False
    for _ in range(1, iterations):
        # Get CNV values associated with each inferred hidden state in previous step
        XZ0 = np.asarray(adata.obsm[f'X_{key_used}'][Z_matrix == -1])
        XZ1 = np.asarray(adata.obsm[f'X_{key_used}'][Z_matrix == 0])
        XZ2 = np.asarray(adata.obsm[f'X_{key_used}'][Z_matrix == 1])
        # Recompute means and variances for each Gaussian
        new_means = np.array([XZ0.mean(), XZ1.mean(), XZ2.mean()]).reshape(-1, 1)
        new_variances = np.array([XZ0.std() ** 2, XZ1.std() ** 2, XZ2.std() ** 2]).reshape(-1, 1, 1)
        new_model = get_hmm_model(new_means, new_variances, HMM_transition_prob=transition_prob)
        # Infer hidden states using the new HMM
        Z_dict_new = apply_hmm(adata=adata,
                               subclone_key=subclone_key,
                               model=new_model,
                               key_used=key_used)
        Z_matrix_new = Z_dict_to_dense(adata, subclone_key, key_used, Z_dict_new)
        if(np.array_equal(Z_matrix, Z_matrix_new)):
            converged = True
            break
        else:
            Z_matrix = Z_matrix_new
            Z_dict = Z_dict_new

    if inplace:
        adata.obsm[f'X_{key_added}'] = Z_matrix
        # For making pl.chromosome_heatmap work, we copy the chromosome positions from the cnv key
        adata.uns[key_added] = adata.uns[key_used]

    return converged, Z_matrix, Z_dict


def fix_vmin_vmax_for_plotting(adata: AnnData, key_used: str = 'hmm', inplace: bool = True):
    """
    Makes sure results of HMM can be plotted by pl.chromosome_heatmap, ensuring that the Z_matrix has both positive and negative values.

    Parameters:
        adata: AnnData
        key_used: the key under which the HMM results are stored in adata.obsm
        inplace: if True, the modified Z_matrix is stored in adata.obsm
    """

    Z_matrix = get_X_from_adata(adata, key_used=key_used)
    # chromosome heatmap can only be displayed if both postive and negative CNV values
    # are present - or equivalently both positive and negative CNV states inferred
    if Z_matrix.min() >= 0:
        print('WARNING, found no deletions, adding a gain on gene 0 of cell 0 for visualization to work')
        Z_matrix[0][0] = -0.00001
    if Z_matrix.max() <= 0:
        print('WARNING, found no gains, adding a gain on gene 1 of cell 0 for visualization to work')
        Z_matrix[0][1] = 0.00001

    if inplace:
        adata.obsm[f'X_{key_used}'] = Z_matrix

    return Z_matrix

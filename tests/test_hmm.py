from infercnvpy.tl._hmm_denoising import *
import numpy as np
from scipy.stats import norm
from scipy import sparse
import infercnvpy as cnv
import warnings


def test_get_model():
    means = np.array([1, 2, 3])
    variances = np.array([4, 5, 6])
    t = 0.1
    transmat = np.array([[1 - 2 * t, t, t],
                         [t, 1 - 2 * t, t],
                         [t, t, 1 - 2 * t]])
    model = get_hmm_model(means,
                          variances,
                          HMM_transition_prob=t)
    assert np.allclose(model.means_, means.reshape(3, 1))
    assert np.allclose(model._covars_, variances.reshape(3, 1, 1))
    assert np.allclose(model.transmat_, transmat)


def test_get_state_parameters_qnorm():
    intended_mean = 10.
    nonmalignant_cnv = np.array([9, 10., 11])
    intended_std = nonmalignant_cnv.std()
    intended_variance = intended_std**2
    mean_delta = np.abs(norm.ppf(0.01, loc=0, scale=intended_std))
    sorted_means, sorted_variances = get_state_parameters_qnorm(nonmalignant_cnv,
                                                                p_val=0.01)
    assert sorted_means[0] == 10. - mean_delta / 2
    assert sorted_means[1] == 10.
    assert sorted_means[2] == 10. + mean_delta / 2
    assert sorted_variances[0] == intended_variance
    assert sorted_variances[1] == intended_variance
    assert sorted_variances[2] == intended_variance


def test_get_state_parameters_infercnv_csr():
    intended_mean = 10.
    nonmalignant_cnv = np.array([9, 10., 11])
    intended_std = nonmalignant_cnv.std()
    intended_variance = intended_std**2
    mean_delta = np.abs(norm.ppf(0.01, loc=0, scale=intended_std))
    sorted_means, sorted_variances = get_state_parameters_qnorm(sparse.csr_matrix(nonmalignant_cnv),
                                                                p_val=0.01)
    assert np.isclose(sorted_means[0], 10. - mean_delta / 2)
    assert np.isclose(sorted_means[1], 10.)
    assert np.isclose(sorted_means[2], 10. + mean_delta / 2)
    assert np.isclose(sorted_variances[0], intended_variance)
    assert np.isclose(sorted_variances[1], intended_variance)
    assert np.isclose(sorted_variances[2], intended_variance)


def test_std_helper_numpy():
    A = np.array([[1, 2, 0],
                  [0, 0, 3],
                  [1, 0, 4]])
    assert std_helper(A) == A.std()


def test_std_helper_csr():
    A = np.array([[1, 2, 0],
                  [0, 0, 3],
                  [1, 0, 4]])
    A_csr = sparse.csr_matrix(A)
    assert std_helper(A_csr) == A.std()


def test_workflow_hmm(adata_oligodendroma):
    # This first call results in 36 warnings, which are also present in test_tools.py/test_workflow()
    cnv.tl.infercnv(adata_oligodendroma,
                    reference_key="cell_type",
                    reference_cat=["Oligodendrocytes (non-malignant)", "Microglia/Macrophage"],
                    window_size=101,
                    step=1,
                    n_jobs=1,)
    converged, Z_matrix, Z_dict = cnv.tl.hmm_denoising(adata_oligodendroma,
                                                       reference_key='cell_type',
                                                       reference_cat=[
                                                           "Oligodendrocytes (non-malignant)", "Microglia/Macrophage"],
                                                       subclone_key='cell_type',
                                                       key_used='cnv',
                                                       key_added='hmm',
                                                       iterations=2)
    adata_oligodendroma.obsm['X_cnv'] = adata_oligodendroma.obsm['X_cnv'].toarray()
    converged, Z_matrix2, Z_dict2 = cnv.tl.hmm_denoising(adata_oligodendroma,
                                                         reference_key='cell_type',
                                                         reference_cat=[
                                                             "Oligodendrocytes (non-malignant)", "Microglia/Macrophage"],
                                                         subclone_key='cell_type',
                                                         key_used='cnv',
                                                         key_added='hmm',
                                                         iterations=2)
    # Testing that the results are the same for both sparse and dense CNV matrices
    assert np.allclose(Z_matrix, Z_matrix2)

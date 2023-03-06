import pandas as pd
from hmmlearn import hmm
import matplotlib.pyplot as plt
import numpy as np
from os.path import join, exists
from os import listdir
from rpy2 import robjects as ro
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import infercnvpy as cnv

from rpy2.robjects.conversion import localconverter
from itertools import chain
import anndata as anndata


def get_df_obs():
    idx_name_tuples = []
    for key, value in chain(ro.conversion.rpy2py(ro.globalenv['r_ref_idx']).items(), ro.conversion.rpy2py(ro.globalenv['r_obs_idx']).items()):
        for idx in value:
            idx_name_tuple = {'cell_index': idx - 1,
                              'cell_type': key}
            idx_name_tuples.append(idx_name_tuple)
    sorted_tuples = sorted(idx_name_tuples, key=lambda x: x['cell_index'])
    df_obs = pd.DataFrame(sorted_tuples)
    colnames = np.array(list(ro.globalenv['r_colnames']))
    df_obs['cell_id'] = df_obs.apply(lambda row: colnames[row['cell_index']], axis=1)
    df_obs = df_obs.set_index('cell_id')
    return df_obs


def get_df_var():
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_var = ro.conversion.rpy2py(ro.globalenv['r_gene_order'])
    df_var = df_var.rename(columns={"chr": "chromosome", 'stop': 'end'})
    rownames = np.array(list(ro.globalenv['r_rownames']))
    for a, b in zip(rownames, df_var.index.to_numpy()):
        assert a == b, 'Mismatch between rownames from expr.data and gene_order'
    return df_var


def get_X(df_var, df_obs):
    n_genes = len(df_var)
    n_cells = len(df_obs)
    X = np.array(list(ro.globalenv['r_expr_data']), dtype=np.float32).reshape(n_cells, n_genes)
    return X


def add_X_to_adata_according_to_stage(adata, X, stage):
    if stage == 3:
        adata.X = X
    elif stage == 15:
        adata.obsm['X_cnv'] = np.log(X)
    elif stage == 17:
        # Hidden states will be 1,2,3 transformed to -1,0,1
        adata.obsm['X_cnv'] = X - 2
    else:
        raise Exception(f'Adding X to adata at stage {stage} is not implemented')


def make_anndata(X, df_obs, df_var, add_cnv_field, stage):
    adata = anndata.AnnData(X, df_obs, df_var)
    if add_cnv_field:
        add_X_to_adata_according_to_stage(adata, X, stage)
        chr_array = adata.var['chromosome'].to_numpy()
        chr_pos = {}
        for chromosome_name in adata.var['chromosome'].unique():
            start = np.argwhere(chr_array == chromosome_name).min()
            stop = np.argwhere(chr_array == chromosome_name).max()
            chr_pos[chromosome_name] = start
        adata.uns['cnv'] = {'chr_pos': chr_pos}
        #adata.obsm['X_cnv_nolog'] = X
        #adata.uns['cnv_nolog'] = {'chr_pos': chr_pos}
    return adata


def infercnvobj_to_anndata(infercnvobj_filename, stage):
    ro.r(
        """
    library(infercnv)
    input_filename <- "{infercnvobj_filename}"
    infercnv_obj <- readRDS(input_filename)
    r_expr_data <- infercnv_obj@expr.data
    r_gene_order <- infercnv_obj@gene_order
    r_ref_idx <- infercnv_obj@reference_grouped_cell_indices
    r_obs_idx <- infercnv_obj@observation_grouped_cell_indices
    r_rownames <- rownames(infercnv_obj@expr.data)
    r_colnames <- colnames(infercnv_obj@expr.data)
    """.format(infercnvobj_filename=infercnvobj_filename)
    )
    df_obs = get_df_obs()
    df_var = get_df_var()
    X = get_X(df_var, df_obs)
    adata = make_anndata(X, df_obs, df_var, add_cnv_field=True, stage=stage)
    return adata


def evaluate_example(example_output_dir, reference_cat, step_size=10, full_pipeline_n_iters=10):
    quickstart_filename = join(example_output_dir, '15_no_subclusteringHMMi3.infercnv_obj')
    quick15_adata = infercnvobj_to_anndata(quickstart_filename, stage=15)
    print('A) Plot of InferCNV-R at stage 15 before it applies HMM:')
    cnv.pl.chromosome_heatmap(quick15_adata, groupby="cell_type", dendrogram=False)
    cnv.tl.hmm_denoising(
        adata=quick15_adata,
        reference_key="cell_type",  # key to know the healthy cells
        reference_cat=reference_cat,  # actual fields for healthy cells,
        subclone_key='cell_type',
        key_used='cnv',  # X_cnv
        key_added='hmm',  # X_hmm,
        #state_means_setting = 'infercnv',
        iterations=1,
        exponentiate=False
    )
    quick17_filename = join(example_output_dir, '17_HMM_predHMMi3.hmm_mode-samples.infercnv_obj')
    quick17_adata = infercnvobj_to_anndata(quick17_filename, stage=17)
    print('B) Plot of InferCNV-R at stage 17 after it has applied HMM:')

    cnv.pl.chromosome_heatmap(quick17_adata, use_rep='cnv', groupby="cell_type", dendrogram=False)
    print('C) Plot of InferCNV-Py HMM applied to (A), should look similar to (B):')
    cnv.pl.chromosome_heatmap(quick15_adata, use_rep='hmm', groupby="cell_type")

    quick3_filename = join(example_output_dir, '03_normalized_by_depthHMMi3.infercnv_obj')
    quick3_adata = infercnvobj_to_anndata(quick3_filename, stage=3)
    cnv.tl.infercnv(
        quick3_adata,
        reference_key="cell_type",
        reference_cat=reference_cat,
        window_size=101,
        step=step_size,
        n_jobs=1,
    )
    print(
        f'D) Taking InferCNV-R output at early stage 3 (before log-normalization), and using InferCNV-Py for the rest of the pipeline with windowing step-size {step_size} before HMM:')
    cnv.pl.chromosome_heatmap(quick3_adata, use_rep='cnv', groupby="cell_type")
    cnv.tl.hmm_denoising(
        adata=quick3_adata,
        reference_key='cell_type',
        reference_cat=reference_cat,
        subclone_key='cell_type',
        key_used='cnv',
        key_added='hmm',
        iterations=full_pipeline_n_iters,
    )
    print(f'E) Applying InferCNV-Py HMM to (D)')
    cnv.pl.chromosome_heatmap(quick3_adata, use_rep='hmm', groupby="cell_type")

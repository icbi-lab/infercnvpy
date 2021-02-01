import scanpy as sc


def umap(adata, *args, **kwargs):
    return sc.pl.umap(adata, "umap_cnv", **kwargs)


def chromosome_heatmap(adata):
    pass

(large-data)=

# Large, lazy and spatial datasets

`infercnv` processes the cells in chunks (see the `chunksize` parameter) and never needs the whole expression
matrix in memory. It therefore works with

- in-memory `numpy` arrays and `scipy` sparse matrices,
- [dask](https://www.dask.org/) arrays (also with sparse `scipy` blocks), e.g. the lazily loaded tables of
  [spatialdata-io](https://spatialdata.scverse.org/projects/io/),
- backed AnnData objects (`anndata.read_h5ad(..., backed="r")`).

Install the optional dependency with `pip install "infercnvpy[dask]"` to use dask arrays.

## How the data is read

1. The reference expression is computed in one streaming pass over the data (mean of all cells, or of every
   reference category).
2. The noise threshold (`dynamic_threshold`) is determined once for the whole dataset and the CNV matrix is
   computed chunk by chunk. For datasets with up to `noise_sample_size` cells (default: 50,000) the standard
   deviation is exact and the chunks computed for it are reused, so the data is read exactly twice in total
   (three times with `calculate_gene_values=True`, where the chunks are not kept).
   For larger datasets, it is estimated from at least ten blocks of `chunksize` cells spread evenly over the
   data, and only these blocks are read an additional time.

Therefore, the result does not depend on `chunksize`. Use `noise_sd="chunk"` to get the behaviour of infercnvpy
0.6 and earlier, where the standard deviation was computed separately for every chunk.

Sparse chunks with less than 10% non-zero entries are not densified: the log fold change of all zeros is identical
for every cell of a gene, so it is smoothed once and only the non-zero entries are processed. Denser chunks are
densified one at a time. Peak memory is therefore determined by `chunksize` and the size of the result, not by the
size of the input.

:::{note}
The result, `adata.obsm["X_cnv"]`, is an in-memory sparse matrix. In our benchmarks (synthetic data) about 10% of
its entries survive the noise filtering, which amounts to roughly 1.5 GB per million cells.
:::

:::{note}
`infercnv` expects log-normalised expression (see {ref}`input-data`). For lazy data, `scanpy.pp.normalize_total`
and `scanpy.pp.log1p` work on dask arrays and stay lazy. Backed AnnData objects cannot be normalised in place by
scanpy; normalise them in chunks beforehand (e.g. when writing the file).
:::

## Scores and plots for many cells

- {func}`infercnvpy.tl.cnv_score` and {func}`infercnvpy.pl.chromosome_heatmap_summary` make a single pass over the
  result and scale to millions of cells.
- {func}`infercnvpy.tl.ithgex` and {func}`infercnvpy.tl.ithcna` are based on the cell x cell correlation matrix
  of each group, which is quadratic in the group size. Groups with more than `max_cells` (default: 5,000) cells are
  randomly subsampled (`random_state` makes this reproducible), so the score is an estimate for such groups.
  Use `max_cells=None` for the exact value. Groups with a single cell get `NaN` in `adata.obs`.
- {func}`infercnvpy.pl.chromosome_heatmap` plots at most `max_cells_per_group` (default: 2,000) randomly selected
  cells per group; the color scale is always computed from all cells. Use `max_cells_per_group=None` to plot all.

## Example: spatial data (SpatialData)

Spatial technologies often produce tables with a lazy dask `X` (CSR blocks), indexed by Ensembl ID, that are part of
a [SpatialData](https://spatialdata.scverse.org/) object. Next to the genes, such a table can contain control
probes or codewords, which don't have a genomic position and have to be removed before normalisation.

```python
import scanpy as sc
import infercnvpy as cnv

sdata = ...  # a SpatialData object (e.g. from `spatialdata-io`); X of its table is a dask array
table = sdata.tables["table"]

# 1. genomic positions: the table is indexed by Ensembl ID, so match the GTF's `gene_id`
cnv.io.genomic_position_from_gtf("gencode.v44.annotation.gtf", table, gtf_gene_id="gene_id")
# (alternatively: cnv.io.genomic_position_from_biomart(table))

# 2. keep gene features only, normalise lazily
genes = table[:, (table.var["feature_type"] == "Gene Expression").to_numpy()].copy()
sc.pp.normalize_total(genes, target_sum=1e4)
sc.pp.log1p(genes)

# 3. infer CNVs; `cell_type` / the categories are placeholders for an annotation of normal cells in `obs`
cnv.tl.infercnv(genes, reference_key="cell_type", reference_cat=["Fibroblast", "Immune"], chunksize=20_000)

# 4. downstream analysis, and map the scores back to the table of the SpatialData object
cnv.tl.pca(genes)
cnv.pp.neighbors(genes)
cnv.tl.leiden(genes)
cnv.tl.cnv_score(genes)
table.obs["cnv_score"] = genes.obs["cnv_score"].reindex(table.obs_names)
```

Notes:

- Which genes are used depends on your panel: genes on chromosomes without a position (and `chrX`/`chrY` by
  default) are skipped, and chromosomes with fewer genes than `window_size` are summarised by a single window.
  Targeted panels with few genes per chromosome carry much less CNV signal than whole-transcriptome data;
  consider reducing `window_size`.
- Use a `chunksize` such that `n_jobs * chunksize * n_genes` values fit into memory when the data is dense; for
  sparse data much larger chunks are fine.

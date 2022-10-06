.. _infercnv-method:

The inferCNV method
===================

Essentially, this package is a Python reimplementation of
`infercnv <https://github.com/broadinstitute/inferCNV/>`_. It mostly follows the computation steps
outlined `here <https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV>`_,
with minor modifications. The computation steps are outlined below.
By making use of `numpy`, `scipy`, and sparse matrices,
it is a lot more computationally efficient.

Computation steps
-----------------

The function parameters are documented at :func:`infercnvpy.tl.infercnv`.

1. Subtract the reference gene expression from all cells. Since the data is in log
   space, this effectively computes the log fold change. If references for
   multiple categories are available (i.e. multiple values are specified to
   `reference_cat`), the log fold change is "bounded":

      * Compute the mean gene expression for each category separately.
      * Values that are within the minimum and the maximum of the mean of all
        references, receive a log fold change of 0, since they are not considered
        different from the background.
      * From values smaller than the minimum of the mean of all references, subtract that minimum.
      * From values larger than the maximum of the mean of all references, subtract that maximum.

   This procedure avoids calling false positive CNV regions due to cell-type specific
   expression of clustered gene regions (e.g. Immunoglobulin- or HLA genes in different
   immune cell types).
2. Clip the fold changes at `-lfc_cap` and `+lfc_cap`.
3. Smooth the gene expression by genomic position. Computes the average over a
   running window of length `window_size`. Compute only every nth window
   to save time & space, where n = `step`.
4. Center the smoothed gene expression by cell, by subtracting the median of each cell
   from each cell.
5. Perform noise filtering. Values `< dynamic_theshold * STDDEV` are set to 0,
   where `STDDEV` is the standard deviation of the smoothed gene expression
6. Smooth the final result using a median filter.

.. _input-data:

Preparing input data
--------------------

The :class:`anndata.AnnData` object should already be filtered for low-quality cells.
`adata.X` needs to be normalized and log-transformed. The method should be
fairly robust to different normalization methods (:func:`scanpy.pp.normalize_total`, scran, etc.).

The method requires a "reference" value to which the expression of genomic
regions is compared. If your dataset contains different cell types and includes
both tumor and normal cells, the average of all cells can be used as reference.
This is the default.

If you already know which cells are "normal", you can provide a column
from `adata.obs` to `reference_key` that contains the annotation. `reference_cat`
specifies one or multiple values in `reference_key` that refer to normal cells.

Alternatively, you can specify a numpy array with average gene expression values
(for instance, derived from a different dataset), to `reference` which will be used
as reference instead. This array needs to align to the `var` axis of adata.

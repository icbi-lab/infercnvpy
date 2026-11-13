# Settings

infercnvpy has no settings object of its own.
It logs through the `scverse.infercnvpy` logger, which shares its level with the other scverse packages.
The default level is `"warning"`, so infercnvpy prints only warnings about the data, e.g. when half or more genes have no genomic position.

To set the level for all scverse packages, use

```python
from scverse_misc.logging import config

config.verbosity = "info"
```

At `"info"`, each infercnvpy function logs one line with its result, e.g. how many windows {func}`infercnvpy.tl.infercnv` computed and how many genes it dropped, and {func}`~infercnvpy.tl.infercnv` shows a progress bar.
A line about dropped input also names a fix, which shows when the level is more detailed than the line's own: at `"info"` for a warning, at `"debug"` for an info line.
To set the level for infercnvpy alone, use `logging.getLogger("scverse.infercnvpy").setLevel("INFO")`.

The thin wrappers around scanpy ({func}`infercnvpy.tl.pca`, {func}`~infercnvpy.tl.umap`, {func}`~infercnvpy.tl.tsne`, {func}`~infercnvpy.tl.leiden`, {func}`infercnvpy.pp.neighbors`, {func}`infercnvpy.pl.umap` and {func}`~infercnvpy.pl.tsne`) call scanpy, which logs through its own logger.
Their messages still follow {attr}`scanpy.settings.verbosity <scanpy.settings.verbosity>`, so set both until scanpy moves to the shared logger:

```python
import scanpy as sc

sc.settings.verbosity = "info"
```

Problems with a call, such as a deprecated argument or a missing reference in {func}`~infercnvpy.tl.infercnv`, are Python warnings and are controlled with the {mod}`warnings` module, not the log level.

"""What infercnvpy logs and warns, and at which level."""

import logging

import numpy as np
import pytest
import scanpy as sc
from scverse_misc.logging import config as log_config

import infercnvpy as cnv

LOGGER = "scverse.infercnvpy"


def _messages(caplog, level=logging.INFO):
    return [(r.levelno, r.getMessage()) for r in caplog.records if r.name == LOGGER and r.levelno >= level]


def _with_unplaced(adata, n):
    adata.var.loc[adata.var_names[:n], "chromosome"] = None
    return adata


@pytest.mark.parametrize(
    ("n_unplaced", "level"),
    [(1, logging.INFO), (5, logging.WARNING)],
)
def test_infercnv_reports_genes_without_position(adata_full_mock, caplog, n_unplaced, level):
    adata = _with_unplaced(adata_full_mock, n_unplaced)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        cnv.tl.infercnv(adata, reference=np.ones(adata.n_vars), window_size=2, step=1, n_jobs=1)
    [(levelno, msg)] = _messages(caplog)
    assert levelno == level
    assert msg.startswith("infercnv: 4 cells x ")
    assert msg.endswith(f"dropped {n_unplaced} of 10 genes without a genomic position")


def test_infercnv_silent_at_default_verbosity(adata_full_mock, caplog, capsys):
    adata = _with_unplaced(adata_full_mock, 1)
    cnv.tl.infercnv(adata, reference=np.ones(adata.n_vars), window_size=2, step=1, n_jobs=1)
    assert _messages(caplog, logging.DEBUG) == []
    assert "it/s" not in capsys.readouterr().err  # no progress bar


def test_infercnv_progress_bar_at_info(adata_full_mock, capsys):
    log_config.verbosity = "info"
    cnv.tl.infercnv(adata_full_mock, reference=np.ones(adata_full_mock.n_vars), window_size=2, step=1, n_jobs=1)
    assert "it/s" in capsys.readouterr().err


def test_infercnv_mean_reference_warns(adata_full_mock):
    with pytest.warns(UserWarning, match="mean of all cells as reference") as record:
        cnv.tl.infercnv(adata_full_mock, window_size=2, step=1, n_jobs=1)
    assert record[0].filename == __file__


@pytest.mark.parametrize("func", [cnv.pp.neighbors, cnv.tl.tsne])
def test_auto_pca_is_info(adata_oligodendroma, caplog, func):
    cnv.tl.infercnv(adata_oligodendroma, reference_key="cell_type", reference_cat="Microglia/Macrophage")
    with caplog.at_level(logging.INFO, logger=LOGGER):
        func(adata_oligodendroma)
    assert _messages(caplog) == [
        (logging.INFO, f"{func.__name__}: computed X_cnv_pca with default parameters, as it was not in adata.obsm")
    ]
    assert "X_cnv_pca" in adata_oligodendroma.obsm


def test_cnv_score_obs_key_is_deprecated(adata_ithgex):
    with pytest.warns(FutureWarning, match="renamed to `groupby`"):
        cnv.tl.cnv_score(adata_ithgex, obs_key="group", inplace=False)


@pytest.mark.parametrize(
    ("n_unknown", "level"),
    [(2, logging.INFO), (6, logging.WARNING)],
)
def test_gtf_reports_unannotated_genes(testdata, caplog, n_unknown, level):
    pytest.importorskip("gtfparse")
    gtf = testdata / "chr21_gencode.gtf"
    known = cnv.io.genomic_position_from_gtf(gtf)["gene_name"].drop_duplicates(keep=False)
    # "Y_RNA" has several entries in the GTF file, so it cannot be placed either
    genes = ["Y_RNA", *known[: 9 - n_unknown], *(f"unknown{i}" for i in range(n_unknown))]
    adata = sc.AnnData(np.zeros((1, len(genes))), var={"gene": genes}).copy()
    adata.var_names = genes

    with caplog.at_level(logging.INFO, logger=LOGGER):
        cnv.io.genomic_position_from_gtf(gtf, adata)

    assert _messages(caplog) == [
        (
            level,
            f"genomic_position_from_gtf: annotated {9 - n_unknown} of 10 genes, "
            f"{n_unknown} with no entry in the GTF file, 1 with more than one entry for their gene_name",
        )
    ]
    assert adata.var["chromosome"].notna().sum() == 9 - n_unknown


@pytest.mark.parametrize(
    ("dropped", "total", "level"),
    [(0, 3, logging.INFO), (1, 3, logging.INFO), (2, 3, logging.WARNING), (1, 2, logging.WARNING)],
)
def test_log_result_warns_from_half_dropped(caplog, dropped, total, level):
    from infercnvpy._logging import log_result

    with caplog.at_level(logging.INFO, logger=LOGGER):
        log_result("f", "dropped %d of %d", dropped, total, dropped=dropped, total=total)
    [(levelno, msg)] = _messages(caplog)
    assert levelno == level
    assert msg == f"f: dropped {dropped} of {total}"

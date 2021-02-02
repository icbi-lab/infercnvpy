from .fixtures import adata_oligodendroma
import infercnvpy as cnv
import pytest


@pytest.mark.parametrize(
    "reference_key,reference_cat",
    [
        (None, None),
        ("cell_type", ["Microglia/Macrophage", "Oligodendrocytes (non-malignant)"]),
    ],
)
def test_infercnv(adata_oligodendroma, reference_key, reference_cat):
    cnv.tl.infercnv(
        adata_oligodendroma, reference_key=reference_key, reference_cat=reference_cat
    )

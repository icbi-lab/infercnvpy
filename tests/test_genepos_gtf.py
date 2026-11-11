import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from infercnvpy.io import genomic_position_from_gtf


@pytest.mark.parametrize("gtf_gene_id", ["gene_id", "gene_name"])
def test_gtf_duplicate_identifiers_are_dropped_by_the_matched_column(tmp_path, gtf_gene_id):
    """Duplicates have to be detected in the column used for matching (e.g. Ensembl ids when matching `gene_id`)."""
    rows = [
        # (chrom, start, end, gene_id, gene_name)
        ("1", 100, 200, "ENSG1", "A"),
        ("2", 300, 400, "ENSG2", "B"),
        ("3", 500, 600, "ENSG3", "C"),
        ("4", 700, 800, "ENSG4", "D"),
    ]
    # make one identifier of the *matched* column ambiguous; the other column stays unique
    if gtf_gene_id == "gene_id":
        rows.append(("X", 900, 1000, "ENSG4", "E"))
    else:
        rows.append(("X", 900, 1000, "ENSG5", "D"))
    gtf = tmp_path / "dup.gtf"
    gtf.write_text(
        "".join(
            f'{c}\ttest\tgene\t{s}\t{e}\t.\t+\t.\tgene_id "{gid}"; gene_name "{gn}";\n' for c, s, e, gid, gn in rows
        )
    )
    ids = {"gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4"], "gene_name": ["A", "B", "C", "D"]}[gtf_gene_id]
    adata = AnnData(X=np.ones((2, 4)), var=pd.DataFrame(index=ids))
    genomic_position_from_gtf(gtf, adata, gtf_gene_id=gtf_gene_id)
    assert adata.var["chromosome"].notna().tolist() == [True, True, True, False]  # the ambiguous one is skipped

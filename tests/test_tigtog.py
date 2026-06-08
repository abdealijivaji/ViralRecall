import random
from pathlib import Path

import pandas as pd

from viralrecall.tigtog import run_tigtog


def test_run_tigtog_with_test_data_directory():
    test_root = Path(__file__).resolve().parent / "test_data"
    data_dir = test_root
    table_file = test_root / "Chlamy_punui_contig_viralregions.annot.tsv"

    assert table_file.exists(), f"Missing fixture: {table_file}"

    table = pd.read_csv(table_file, sep="\t", header=0)
    vreg = "".join(random.Random(0).choices("ACGT", k=100))

    pred, prob = run_tigtog(vreg, table, data_dir)

    assert pred == "Imitervirales"
    assert prob == "98.71"

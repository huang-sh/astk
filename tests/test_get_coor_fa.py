
from pathlib import Path
import pandas as pd
from astk.utils import get_coor_fa


def test_get_coor_fa(shared_datadir):
    df = pd.read_csv(f"{shared_datadir}/coord.tsv", sep="\t", header=None)
    gfa = "/home/public/ref/genome/mm/release_M25/GRCm38.primary_assembly.genome.fa"
    if not Path(gfa).exists():
        return
    output = f"{shared_datadir}/coord.testout"
    get_coor_fa(df, gfa, output)
    with open(output, "r") as th, open(f"{shared_datadir}/coord.val") as vh:
        for l1, l2 in zip(th.readlines(), vh.readlines()):
            if l1.startswith(">"):
                continue
            assert l1 == l2

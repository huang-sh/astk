import sys
import subprocess
from pathlib import Path

import astk.utils.func  as ul
from astk.constant import BASE_DIR, RBP_sp_dic


def co_splice_net(output, fasta, psimeta, tqmeta, organism, database, gtf):
    sp = RBP_sp_dic.get(organism, None)
    if sp is None:
        print("--organism is not valid")
        sys.exit(1)

    rscript = BASE_DIR / "R" / "CoSpliceNet.R"
    csv_out = Path(output).with_suffix(".csv")
    html_out = Path(output).with_suffix(".html")
    param_dic = {
        "output": csv_out,
        "fasta": fasta, 
        "psiMeta": psimeta, 
        "tqMeta": tqmeta,
        "organism": sp,
        "database": database,
        "gtf": gtf
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])
    ul.coSpliceNet(csv_out, html_out)
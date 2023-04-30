import os
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Sequence

from astk.constant import *
from . import func as ulf
from astk.utils.meta_template import Template


def meta(
    output: str,
    ctrl_file: Sequence[str], 
    case_file: Sequence[str], 
    ctrl_rep: Sequence[int], 
    case_rep: Sequence[int], 
    groupname: Sequence[str], 
    **kwargs
):
    app = kwargs.pop("app")
    try:
        tp = Template(kwargs.pop("condition", None))
        tp.complete_df(groupname, ctrl_file, case_file, ctrl_rep, case_rep, **kwargs)
    except BaseException as e:
        print(e)
        sys.exit()

    if app == "SUPPA2":        
        tp.to_csv(output)
        tp.to_json(output)
    elif app == "rMATS":
        for gn, gdf in tp.df.groupby("group"):
            for cn, cdf in gdf.groupby("condition"):
                path_str = ",".join(cdf["path"].tolist())
                out = Path(output).with_suffix(f".{gn}.{cn}.txt")
                with open(out, "w") as f:
                    f.write(path_str)


def install(requirement, OrgDb, cran, bioconductor, java, mirror, conda):

    pdir = BASE_DIR
    rscript = pdir / "R" / "install.R"
    if conda:
        if cran:            
            cran_pkg = " ".join([f"r-{i.lower()}" for i in cran])
            os.system(f"conda install -y -c conda-forge -c r  {cran_pkg}")
        if bioconductor:
            bioc_pkg = " ".join([f"bioconductor-{i.lower()}" for i in bioconductor])
            os.system(f"conda install -y -c bioconda {bioc_pkg}")
    if java:
        print("install ChromHMM")
        ch.install(pdir)
        shutil.copyfile(pdir/"ChromHMM.jar", pdir/"ChromHMM/ChromHMM.jar")
    param_dic = {
        "requirement": requirement,
        "OrgDb": OrgDb or False,
        "CRAN": cran or False,
        "bioconductor": bioconductor or False,
        "mirror": mirror,
    }
    param_ls = ulf.parse_cmd_r(**param_dic)
    subprocess.run([ulf.Rscript_bin(), rscript, *param_ls])


def getmeme(motifid, meme, database, organism, output):

    sp = RBP_sp_dic.get(organism, None)
    if meme:
        meme_file = meme
    elif database and sp:
        meme_file = BASE_DIR  / f"data/motif/{database}/{sp}.meme"
    else:
        print("--meme or --db and -org muset be set")    
        sys.exit()

    try:
        ulf.get_meme(motifid, meme_file, output)
    except BaseException as e:
        print(e)


def list_(OrgDb, RBPSp):
    if OrgDb:
        for k, v in OrgDb_dic.items():
            print(f"{k}: {v}")      
    if RBPSp:
        for k, v in RBP_sp_dic.items():
            print(f"{k}: {v}")


def anchor(file, output, index, sideindex, offset5, offset3, strand_sp):
    try:
        ulf.gen_anchor_bed(file, output, index, sideindex, offset5, offset3, strand_sp)
    except BaseException as e:
        print(e)
 

def mkTxDb(gtf, organism, output):
    sp = RBP_sp_dic.get(organism, None)
    sp = sp.replace("_", " ")
    rscript = BASE_DIR / "R" / "makeTxDb.R"
    if not output:
        output = Path(gtf).with_suffix(".txdb")
    param_dic = {
        "organism": sp,
        "gtf": gtf,
        "output": output
    }
    param_ls = ulf.parse_cmd_r(**param_dic)
    subprocess.run([ulf.Rscript_bin(), rscript, *param_ls])


def getgene(file, output, unique):
    import pandas as pd
    from astk.event import SuppaEventID

    dpsi_df = pd.read_csv(file, sep="\t", index_col=0)
    gene_df = pd.DataFrame(
        {"geneID": [SuppaEventID(i).gene_id for i in dpsi_df.index]}
    )
    if unique:
        gene_df.drop_duplicates(inplace=True)
    if output is None:
        output = sys.stdout
    gene_df.to_csv(output, index=False, header=False, sep="\t")

    

from ast import Pass
import sys
import shutil
import subprocess
from pathlib import Path

import astk.ChromHMM as ch
from astk.constant import *
from . import func as ulf

from astk.utils.meta_template import Template

def meta(output, replicate, groupname, control, treatment, replicate1, replicate2, **kwargs):

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit() 
    if replicate:
        repN1 = [int(i) for i in replicate]
        repN2 = [int(i) for i in replicate]
    elif all([replicate1, replicate1]):
        repN1 = [int(i) for i in replicate1]
        repN2 = [int(i) for i in replicate2]
    else:
        print("repN1 and repN2 must be set!")    
        sys.exit()
    try:
        tp = Template(kwargs.pop("condition", None))
        tp.complete_df(groupname, control, treatment, repN1, repN2, **kwargs)
        tp.to_csv(output)
        tp.to_json(output)
    except BaseException as e:
        print(e)


def install(requirement, OrgDb, cran, bioconductor, java, mirror):

    pdir = BASE_DIR
    rscript = pdir / "R" / "install.R"
    
    if java:
        print("install ChromHMM")
        ch.install(pdir)
        shutil.copyfile(pdir/"ChromHMM.jar", pdir/"ChromHMM/ChromHMM.jar")
    param_dic = {
        "requirement": requirement,
        "OrgDb": OrgDb if OrgDb else False,
        "CRAN": cran if cran else False, 
        "bioconductor": bioconductor if bioconductor else False, 
        "mirror": mirror
    }
    param_ls = ulf.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


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
 

def getcoor(file, output, start, end, strand_sp, anchor, upstream_w, downstream_w, fasta):

    try:
        coor_df = ulf.get_coor_bed(file, start, end, strand_sp, anchor, upstream_w, downstream_w)
        coor_df.to_csv(output, index=False, header=False, sep="\t")
        if fasta:
            ulf.get_coor_fa(coor_df, fasta, Path(output).with_suffix(".fa"))
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
    subprocess.run(["Rscript", rscript, *param_ls])    

import sys
import subprocess
from pathlib import Path

from astk.constant import RBP_sp_dic, BASE_DIR
import astk.utils.func  as ul


def motif_enrich(tfasta, cfasta, outdir, database, organism):

    Path(outdir).mkdir(exist_ok=True)
    rscript = BASE_DIR / "R" / "motifEnrich.R"

    if cfasta and len(tfasta) != len(cfasta):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()
    elif cfasta is None:
        cfasta = ["0" for _ in tfasta]

    sp = RBP_sp_dic.get(organism, "0")

    param_dic = {
        "tfile": tfasta,
        "cfile": cfasta,
        "outdir": outdir,
        "database": database,
        "organism": sp
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


def motif_find(tfasta, cfasta, outdir, database, organism, pvalue, evalue, minw, maxw):

    Path(outdir).mkdir(exist_ok=True)
    rscript = BASE_DIR / "R" / "motifFind.R"

    sp = RBP_sp_dic.get(organism, "0")

    if cfasta and len(tfasta) != len(cfasta):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()
    elif cfasta is None:
        cfasta = ["0" for _ in tfasta]

    param_dic = {
        "outdir": outdir,
        "tfile": tfasta, 
        "cfile": cfasta, 
        "pvalue": pvalue,
        "minw": minw, 
        "maxw": maxw,
        "database": database,
        "organism": sp,
        "eval": evalue
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])

               
def motif_plot(motifid, database, organism, meme, output, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "motifPlot.R"
    sp = RBP_sp_dic.get(organism, None)
    if meme:
        meme_file = meme
    elif database and sp:
        meme_file = BASE_DIR  / f"data/motif/{database}/{sp}.meme"
    else:
        print("--meme or --db and -org muset be set")    
        sys.exit()
    param_dic = {
        "motifId": motifid,
        "meme": meme_file,
        "fmt": fmt,
        "width": width, 
        "height": height,
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


def mmap(fasta, name, center, meme, outdir, binsize, step, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "motifMap.R"
    Path(outdir).mkdir(exist_ok=True)

    param_dic = {
        "fasta": fasta,
        "outdir": outdir, 
        "meme": meme,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "bin": binsize,
        "step": step,
        "fmt": fmt,
        "center": center if len(center) ==len(fasta) else ["0"] * len(fasta),
        "seqid": name if len(name)==len(fasta) else [Path(i).stem for i in fasta]
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])
        
import sys
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile

import astk.utils as ul
from astk.constant import *


def motif_enrich(tevent, cevent, outdir, database, organism, fasta):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    case_fa_ls = []
    df_dic = ul.get_ss_bed(tevent, exon_width=150, intron_width=150)
    for ssi, df in df_dic.items():
        ssi_dir = outdir / ssi
        ssi_dir.mkdir(exist_ok=True)
        ssi_bed = ssi_dir / f"{ssi}.bed"
        ssi_fa = ssi_dir / f"{ssi}.fa"
        ## the complete event_id may cause some bug in motif analysis, eg fimo 
        df["name"] = [f"event{i}" for i in range(df.shape[0])]
        df.to_csv(ssi_bed, index=False, header=False, sep="\t")
        ul.get_coor_fa(df, fasta, ssi_fa)
        case_fa_ls.append(ssi_fa)
    
    if cevent:
        ctrl_fa_ls = []
        cdf_dic = ul.get_ss_bed(cevent, exon_width=150, intron_width=150)
        for cssi, cdf in cdf_dic.items():
            cssi_dir = outdir / cssi
            cssi_dir.mkdir(exist_ok=True)
            cssi_bed = cssi_dir / f"{cssi}_ctrl.bed"
            cssi_fa = cssi_dir / f"{cssi}_ctrl.fa"
            cdf["name"] = [f"event{i}" for i in range(cdf.shape[0])]
            cdf.to_csv(cssi_bed, index=False, header=False, sep="\t")
            ul.get_coor_fa(cdf, fasta, cssi_fa)
            ctrl_fa_ls.append(cssi_fa)
    else:
        ctrl_fa_ls = ["0" for _ in case_fa_ls]

    if case_fa_ls and len(case_fa_ls) != len(ctrl_fa_ls):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()

    sp = RBP_sp_dic.get(organism, "0")
    param_dic = {
        "tfile": case_fa_ls,
        "cfile": ctrl_fa_ls,
        "outdir": outdir,
        "database": database,
        "organism": sp,
        "meme_path": ul.get_meme_path()
    }
    rscript = BASE_DIR / "R" / "motifEnrich.R"
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def motif_find(tevent, cevent, outdir, database, organism, pvalue, evalue, minw, maxw, fasta):
    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)

    case_fa_ls = []
    df_dic = ul.get_ss_bed(tevent, exon_width=150, intron_width=150)
    for ssi, df in df_dic.items():
        ssi_dir = outdir / ssi
        ssi_dir.mkdir(exist_ok=True)
        ssi_bed = ssi_dir / f"{ssi}.bed"
        ssi_fa = ssi_dir / f"{ssi}.fa"
        ## ## the complete event_id may cause some bug in motif analysis, eg fimo 
        df["name"] = [f"event{i}" for i in range(df.shape[0])]
        df.to_csv(ssi_bed, index=False, header=False, sep="\t")
        ul.get_coor_fa(df, fasta, ssi_fa)
        case_fa_ls.append(ssi_fa)
    
    if cevent:
        ctrl_fa_ls = []
        cdf_dic = ul.get_ss_bed(cevent, exon_width=150, intron_width=150)
        for cssi, cdf in cdf_dic.items():
            cssi_dir = outdir / cssi
            cssi_dir.mkdir(exist_ok=True)
            cssi_bed = cssi_dir / f"{cssi}_ctrl.bed"
            cssi_fa = cssi_dir / f"{cssi}_ctrl.fa"
            cdf["name"] = [f"event{i}" for i in range(cdf.shape[0])]
            cdf.to_csv(cssi_bed, index=False, header=False, sep="\t")
            ul.get_coor_fa(cdf, fasta, cssi_fa)
            ctrl_fa_ls.append(cssi_fa)
    else:
        ctrl_fa_ls = ["0" for _ in case_fa_ls]

    rscript = BASE_DIR / "R" / "motifFind.R"

    sp = RBP_sp_dic.get(organism, "0")

    if case_fa_ls and len(case_fa_ls) != len(ctrl_fa_ls):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()

    param_dic = {
        "outdir": outdir,
        "tfile": case_fa_ls, 
        "cfile": ctrl_fa_ls, 
        "pvalue": pvalue,
        "minw": minw, 
        "maxw": maxw,
        "database": database,
        "organism": sp,
        "eval": evalue,
        "meme_path": ul.get_meme_path()
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])

               
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
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def mmap(event, fasta, name, center, meme, outdir, binsize, step, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "motifMap.R"
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    case_fa_ls = []
    df_dic = ul.get_ss_bed(event, exon_width=150, intron_width=150)
    for ssi, df in df_dic.items():
        ssi_dir = outdir / ssi
        ssi_dir.mkdir(exist_ok=True)
        ssi_bed = ssi_dir / f"{ssi}.bed"
        ssi_fa = ssi_dir / f"{ssi}.fa"
        ## ## the complete event_id may cause some bug in motif analysis, eg fimo 
        df["name"] = [f"event{i}" for i in range(df.shape[0])]
        df.to_csv(ssi_bed, index=False, header=False, sep="\t")
        ul.get_coor_fa(df, fasta, ssi_fa)
        case_fa_ls.append(ssi_fa)

    param_dic = {
        "fasta": case_fa_ls,
        "outdir": outdir, 
        "meme": meme,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "bin": binsize,
        "step": step,
        "fmt": fmt,
        "center": center if len(center) ==len(case_fa_ls) else ["0"] * len(case_fa_ls),
        "seqid": name if len(name)==len(case_fa_ls) else [Path(i).stem for i in case_fa_ls],
        "meme_path": ul.get_meme_path()
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])
        
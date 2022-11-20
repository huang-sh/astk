import sys
import subprocess
from pathlib import Path

import astk.utils as ul
from astk.constant import *


def _centrimo_enrich(outdir, tfa, cfa, motif, ethresh=1):
    param_dic = {
        "oc": outdir,
        "local": True,
        "norc": True,
        "ethresh": ethresh,
        "neg": cfa
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    exe_bin = "centrimo"
    subprocess.run([exe_bin, *param_ls, tfa, motif])


def _ame_enrich(outdir, tfa, cfa, motif):
    param_dic = {
        "oc": outdir,
        "control": cfa if cfa else "--shuffle--",
        "text": True
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    exe_bin = "ame"
    subprocess.run([exe_bin, *param_ls, tfa, motif])


def _streme_enrich(outdir, tfa, cfa, pvalue, evalue, maxw, minw):
    param_dic = {
        "oc": outdir,
        "rna": True,
        "thresh": pvalue,
        "evalue": evalue,
        "minw": minw,
        "maxw": maxw,
        "p": tfa
    }
    if cfa:
        param_dic["n"] = cfa
    param_ls = ul.parse_cmd_r(**param_dic)
    exe_bin = "streme"
    subprocess.run([exe_bin, *param_ls])


def _tomtom_enrich(outdir, motif, motif_db):
    param_dic = {
        "oc": outdir,
        "thresh": 5,
        "evalue": True,
        "norc": True,
        "min-overlap": 5,
        "no-ssc": True
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    exe_bin = "tomtom"
    subprocess.run([exe_bin, *param_ls, motif, motif_db])


def motif_enrich(tevent, cevent, outdir, database, organism, fasta, method):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    motif_path = ul.get_motif_path(database, organism)
    df_dic = ul.get_ss_bed(tevent, exon_width=150, intron_width=150)
    if cevent is not None:
        cdf_dic = ul.get_ss_bed(cevent, exon_width=150, intron_width=150)
    else:
        cdf_dic = {}
    for ssi, df in df_dic.items():
        ssi_dir = outdir / ssi
        ssi_dir.mkdir(exist_ok=True)
        ssi_bed = ssi_dir / f"{ssi}.bed"
        ssi_fa = ssi_dir / f"{ssi}.fa"
        ## the complete event_id may cause some bug in motif analysis, eg fimo 
        df.to_csv(ssi_bed, index=False, header=False, sep="\t")
        ul.get_coor_fa(df, fasta, ssi_fa, strandedness=True, rna=True)
        if (cdf := cdf_dic.get(ssi, None)) is None:
            cssi_fa = None
        else:
            cssi_bed = ssi_dir / f"{ssi}_ctrl.bed"
            cssi_fa = ssi_dir / f"{ssi}_ctrl.fa"
            cdf.to_csv(cssi_bed, index=False, header=False, sep="\t")
            ul.get_coor_fa(cdf, fasta, cssi_fa, strandedness=True, rna=True)
        if method == "centrimo":
            _centrimo_enrich(ssi_dir, ssi_fa, cssi_fa, motif_path)
        else:
            _ame_enrich(ssi_dir, ssi_fa, cssi_fa, motif_path)


def motif_find(tevent, cevent, outdir, **kwargs):
    fasta = kwargs.pop("fasta")
    motifcmp = kwargs.pop("motifcmp")
    motif_db = ul.get_motif_path(kwargs.pop("database"), kwargs.pop("organism"))
    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)
    df_dic = ul.get_ss_bed(tevent, exon_width=150, intron_width=150)
    if cevent is not None:
        cdf_dic = ul.get_ss_bed(cevent, exon_width=150, intron_width=150)
    else:
        cdf_dic = {}    
    for ssi, df in df_dic.items():
        ssi_dir = outdir / ssi
        ssi_dir.mkdir(exist_ok=True)
        ssi_bed = ssi_dir / f"{ssi}.bed"
        ssi_fa = ssi_dir / f"{ssi}.fa"
        df.to_csv(ssi_bed, index=False, header=False, sep="\t")
        ul.get_coor_fa(df, fasta, ssi_fa, strandedness=True, rna=True)
        if (cdf := cdf_dic.get(ssi, None)) is None:
            cssi_fa = None
        else:
            cssi_bed = ssi_dir / f"{ssi}_ctrl.bed"
            cssi_fa = ssi_dir / f"{ssi}_ctrl.fa"
            cdf.to_csv(cssi_bed, index=False, header=False, sep="\t")
            ul.get_coor_fa(cdf, fasta, cssi_fa, strandedness=True, rna=True)        
        _streme_enrich(ssi_dir, ssi_fa, cssi_fa, **kwargs)
        if motifcmp:
            _tomtom_enrich(ssi_dir, ssi_dir/"streme.txt", motif_db)


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
        ul.get_coor_fa(df, fasta, ssi_fa, strandedness=True, rna=True)
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
        
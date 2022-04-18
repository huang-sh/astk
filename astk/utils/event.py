from pathlib import Path

import astk.utils.func  as ul
import astk.utils.select as sl
from astk.constant import AS_type


def len_dist(infile, output, custom_len, cluster, width, len_weight, max_len):
    import pandas as pd

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   
    ioe_df = pd.read_csv(infile, sep="\t")
    if "event_id" not in ioe_df.columns:
        print("file does not support!")
        exit()

    info_df = ul.extract_info(ioe_df)
    info_df["event_id"] = ioe_df["event_id"]
    output = Path(output).with_suffix(".png")
    if custom_len:
        lens = [int(i) for i in custom_len]
        df = ul.custome_cluster_len(info_df, output, lens, width=width, max_len=max_len)
    else:
        df = ul.cluster_len(info_df, output, n_cls=cluster, max_len=max_len, len_weight=len_weight, width=width)
    df.to_csv(Path(output).parent / f"{Path(output).stem}_cls_info.csv", index=False)


def len_cluster(files, indir, outdir, lenrange):

    outdir = Path(outdir)
    outdir = Path(outdir).absolute()
    od_name = outdir.name

    lrs = list(map(int, lenrange))
    coor_ls = [(lrs[i], lrs[i+1]) for i in range(len(lrs)-1)]

    if indir:
        files = list(Path(indir).glob("*psi"))
    else:
        files = [Path(i) for i in files]
    
    if len(files) == 1:
        file = Path(files[0])
        outdir.mkdir(exist_ok=True)
        for s, e in coor_ls:
            outfile = outdir / f"{file.stem}_{s}-{e}{file.suffix}"
            ul.df_len_select(file, outfile, s, e)
    else:
        for s, e in coor_ls:
            s_outdir = outdir.with_name(f"{od_name}_{s}-{e}")
            s_outdir.mkdir(exist_ok=True)
            for file in files:
                outfile = s_outdir / Path(file).name
                ul.df_len_select(file, outfile, s, e+1)


def len_pick(infile, output, len_range):
    import astk.utils.event_id as ei
    import pandas as pd

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   
    AS_len = lambda x: ei.SuppaEventID(x).alter_element_len

    df = pd.read_csv(infile, sep="\t", index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    s, e = len_range
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(output, index=True, sep="\t", na_rep="nan", index_label=False)


def diff_splice(outdir, metadata, gtf, event_type, method, exon_len, poolgenes, tpm_col):

    from .func import DiffSplice
    if event_type == "all":
        event_types = AS_type
    else:
        event_types = [event_type]
    try:
        dsi = DiffSplice(outdir, metadata, gtf, event_types, exon_len, poolgenes)
        dsi.ds(method)
    except BaseException as e:
        print(e)


def sigfilter(files, outdir, dpsi, pval, abs_dpsi, psifile1, psifile2, fmt):

    for idx, dpsi_file in enumerate(files):
        if len(files) == len(psifile1) and len(files) == len(psifile2):
            psifiles = (psifile1[idx], psifile2[idx])
        else:
            psifiles = ()
            
        sf = sl.SigFilter(dpsi_file, outdir, dpsi, pval, abs_dpsi, psifiles, fmt)
        sf.run()


def psi_filter(file, output, psi, quantile):

    if file:
        pf = sl.PsiFilter(file, psi, quantile)
        pf.run(output)
  

def intersect(file_a, file_b, output):
    import pandas as pd

    fa = Path(file_a)
    fb = Path(file_b)
    dfa = pd.read_csv(fa, sep="\t", index_col=0)
    dfb = pd.read_csv(fb, sep="\t", index_col=0)

    share_id = list(set(dfa.index) & set(dfb.index))

    out_a = Path(f"{output}_a").with_suffix(fa.suffix)
    out_b = Path(f"{output}_b").with_suffix(fb.suffix)
    
    dfa.loc[share_id, :].to_csv(out_a, sep="\t")
    dfb.loc[share_id, :].to_csv(out_b, sep="\t")

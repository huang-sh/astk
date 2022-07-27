from pathlib import Path

import astk.utils.func  as ul
import astk.utils.select as sl
from astk.types import FilePath
from ._func import SuppaEventID


def len_dist(infile, output, custom_len, cluster, width, len_weight, max_len):
    import pandas as pd

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    edf = pd.read_csv(infile, sep="\t", index_col=0)
    ae_lens = [SuppaEventID(eid).alter_element_len for eid in edf.index]
    counts, bin_edges = ul.len_hist(ae_lens, width, max_len)
    cluster_ls = []
    clens = [0] + list(custom_len)
    print(clens)
    for i in range(1, len(clens)):
        print( clens[i-1], clens[i])
        subset = [l for l in  ae_lens if clens[i-1] < l <= clens[i]]
        cluster_ls.append(subset)
    else:
        subset = [l for l in  ae_lens if l > clens[-1]]
        cluster_ls.append(subset)
    output = Path(output).with_suffix(".png")
    ul.plot_hist_cluster(output, cluster_ls, bin_edges)


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
    import pandas as pd

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   
    AS_len = lambda x: SuppaEventID(x).alter_element_len

    df = pd.read_csv(infile, sep="\t", index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    s, e = len_range
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(output, index=True, sep="\t", na_rep="nan", index_label=False)


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
  

def intersect(
    file_a: FilePath,
    file_b: FilePath,
    output: FilePath, 
    ioeb: FilePath, 
    ignoreb: bool
    ) -> None:

    from pandas import read_csv

    outname = Path(output).stem
    outdir = Path(output).parent
    fa = Path(file_a)
    dfa = read_csv(fa, sep="\t", index_col=0)

    if file_b:
        fb = Path(file_b)
        dfb = read_csv(fb, sep="\t", index_col=0)
        fb_idx = dfb.index
        out_a = outdir / Path(f"{outname}_a{fa.suffix}")
    elif ioeb:
        dfb = read_csv(ioeb, sep="\t")
        fb_idx = dfb["event_id"]
        out_a = outdir / Path(f"{outname}{fa.suffix}")
    
    share_id = list(set(dfa.index) & set(fb_idx))   
    dfa.loc[share_id, :].to_csv(out_a, sep="\t")

    if file_b and (not ignoreb):
        out_b = outdir / Path(f"{outname}_b{fb.suffix}")
        dfb.loc[share_id, :].to_csv(out_b, sep="\t")

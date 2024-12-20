from pathlib import Path

from ..ctypes import FilePath
from .eid import SuppaEventID
from ..lazy_loader import LazyLoader
from .. import utils as ul
from ..constant import ALT_IDX


pd = LazyLoader("pd", globals(), "pandas")


def len_dist(infile, output, custom_len, cluster, width, len_weight, max_len, fmt):
    einfo = ul.detect_file_info(infile)
    app, etype = einfo["app"], einfo["etype"]
    coori = ul.get_event_coord(infile, app)
    df_ss = coori.df_ss
    a1, a2 = ALT_IDX[app][etype][:2]
    ae_lens = df_ss.iloc[:, a2] - df_ss.iloc[:, a1] + 1
    counts, bin_edges = ul.len_hist(ae_lens, width, max_len)
    cluster_ls = []
    clens = [0] + list(custom_len)
    for i in range(1, len(clens)):
        print( clens[i-1], clens[i])
        subset = [l for l in  ae_lens if clens[i-1] < l <= clens[i]]
        cluster_ls.append(subset)
    subset = [l for l in  ae_lens if l > clens[-1]]
    cluster_ls.append(subset)
    output = Path(output).with_suffix(f".{fmt}")
    ul.plot_hist_cluster(output, cluster_ls, bin_edges)


def len_cluster(files, indir, outdir, lenrange):

    outdir = Path(outdir)
    outdir = Path(outdir).absolute()
    outdir.mkdir(exist_ok=True)

    lrs = list(map(int, lenrange))
    coor_ls = [(lrs[i], lrs[i+1]) for i in range(len(lrs)-1)]

    files = list(Path(indir).glob("*psi")) if indir else [Path(i) for i in files]
    if len(files) == 1:
        file = Path(files[0])
        outdir.mkdir(exist_ok=True)
        for s, e in coor_ls:
            outfile = outdir / f"{file.stem}_{s}-{e}{file.suffix}"
            ul.df_len_select(file, outfile, s, e)
    else:
        for s, e in coor_ls:
            len_outdir = outdir / f"{s}-{e}"
            len_outdir.mkdir(exist_ok=True)
            for file in files:
                outfile = len_outdir / Path(file).name
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


# it is deprecated
def _sigfilter(files, outdir, dpsi, pval, abs_dpsi, psifile1, psifile2, fmt):

    for idx, dpsi_file in enumerate(files):
        if len(files) == len(psifile1) and len(files) == len(psifile2):
            psifiles = (psifile1[idx], psifile2[idx])
        else:
            psifiles = ()
            
        sf = ul.igFilter(dpsi_file, outdir, dpsi, pval, abs_dpsi, psifiles, fmt)
        sf.run()


def sigfilter(file, output, dpsi, pval, qval, abs_dpsi, sep, app):
    from pandas import read_csv

    if app == "auto":
        app = ul.detect_file_info(file)["app"]
    dpsi_df = read_csv(file, sep="\t", index_col=0).dropna()
    kwargs = {"dpsi":dpsi, "abs_dpsi": abs_dpsi, "pval": pval, "qval": qval,"app": app}
    if app == "SUPPA2":
        old_col = dpsi_df.columns
        dpsi_col = old_col[0]
        dpsi_df.columns = ["dpsi", "pval"]
        kwargs.pop("qval")
        df_fil = ul.sig_filter(dpsi_df, **kwargs)
        df_fil.columns = old_col
    elif app == "rMATS":
        dpsi_col = "IncLevelDifference"
        df_fil = ul.sig_filter(dpsi_df, **kwargs)
    if not output:
        output = Path(file).with_suffix(f".sig{Path(file).suffix}")
        if sep:
            pos_out = output = Path(file).with_suffix(f".sig+{Path(file).suffix}")
            neg_out = output = Path(file).with_suffix(f".sig-{Path(file).suffix}")
    elif sep:
        pos_out = Path(output).with_suffix(f".d+{Path(output).suffix}")
        neg_out = Path(output).with_suffix(f".d-{Path(output).suffix}")

    df_fil.to_csv(output, sep="\t")
    if sep and abs_dpsi:
        pos_df = df_fil.loc[df_fil[dpsi_col] > 0, ]
        neg_df = df_fil.loc[df_fil[dpsi_col] < 0, ]
        pos_df.to_csv(pos_out, sep="\t")
        neg_df.to_csv(neg_out, sep="\t")


def psi_filter(
    file: str, 
    output: str, 
    minv: float = 0, 
    maxv: float = 1, 
    minq: float = 0, 
    maxq: float = 1, 
    app: str = "auto"
):
    from pandas import read_csv

    def _rmats_mean_psi(vstr):
        strings = vstr.split(",")
        values = [float(i) for i in strings if i != "NA"]
        return sum(values) / len(values)
        
    dpsi_df = read_csv(file, sep="\t", index_col=0).dropna()
    if app == "auto":
        app = ul.detect_file_info(file)["app"]
    if app == "rMATS":        
        dpsi_df["PSI"] = dpsi_df["IncLevel1"].apply(_rmats_mean_psi)
    elif app == "SUPPA2":
        dpsi_df["PSI"] = dpsi_df.apply(lambda row: sum(row)/len(row), axis=1)

    min_psi = max([dpsi_df["PSI"].quantile(minq), minv])
    max_psi = min([dpsi_df["PSI"].quantile(maxq), maxv])
    keep = (dpsi_df["PSI"] >= min_psi) & (dpsi_df["PSI"] <= max_psi)
    df_fil = dpsi_df.loc[keep, ]
    del df_fil["PSI"]
    df_fil.to_csv(output, sep="\t")
    return df_fil


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

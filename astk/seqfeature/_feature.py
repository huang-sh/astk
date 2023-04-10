from pathlib import Path
import multiprocessing as mp

import astk.utils as ul
from astk.lazy_loader import LazyLoader

np = LazyLoader("np", globals(), "numpy")
pd = LazyLoader("pd", globals(), "pandas")
plt = LazyLoader("plt", globals(), "matplotlib.pyplot")


def _seq_gcc(seq):
    seq = seq.upper()
    GN = seq.count("G")
    GC = seq.count("C")
    return( GN + GC) / len(seq)


def _get_seq_gcc(seq, binsize):
    if len(seq) <= binsize or binsize < 10:
        gcc = [_seq_gcc(seq)]
    else:
        gcc = [_seq_gcc(seq[i:i+binsize]) for i in range(len(seq)-binsize-1)]
    return gcc


def get_seq_gcc(fasta, binsize, process=4):
    with open(fasta, "r") as rh:
        seqs = ul.read_fasta(rh)
        pool = mp.Pool(process)
        gccs = [pool.apply(_get_seq_gcc, args=(seq[1], binsize)) for seq in seqs]
        pool.close()
    
    return pd.DataFrame(gccs)


def get_gcc(coor_dic, outdir, gfasta, binsize):
    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)

    df_ls = []
    fig, axes = plt.subplots(1, len(coor_dic), figsize= (10, 5))
    for idx, (ssn, (df_up, df_dw)) in enumerate(coor_dic.items()):
        ss_dir = outdir / ssn
        ss_dir.mkdir(exist_ok=True)
        ul.get_coor_fa(df_up, gfasta, ss_dir / f"{ssn}_ups.fa", strandedness=True)
        ul.get_coor_fa(df_dw, gfasta, ss_dir / f"{ssn}_dws.fa", strandedness=True)

        df_up.to_csv(ss_dir / f"{ssn}_ups.bed", index=False, header=False, sep="\t")
        df_dw.to_csv(ss_dir / f"{ssn}_dws.bed", index=False, header=False, sep="\t")

        ups_gcc = get_seq_gcc(ss_dir / f"{ssn}_ups.fa", binsize)
        dws_gcc = get_seq_gcc(ss_dir / f"{ssn}_dws.fa", binsize)
        ups_gcc.to_csv(ss_dir / f"{ssn}_ups_gcc.csv")
        dws_gcc.to_csv(ss_dir / f"{ssn}_dws_gcc.csv")
        if "5SS" in ssn:
            if binsize == 0:
                ups_gcc.columns =  [f"{ssn}_exon"]
                dws_gcc.columns = [f"{ssn}_intron"]
            else:
                ups_gcc.columns = [f"{ssn}_exon_b{i}" for i in range(1-ups_gcc.shape[1], 1)]
                dws_gcc.columns = [f"{ssn}_intron_b{i}" for i in range(1, dws_gcc.shape[1]+ 1)]
        elif "3SS" in ssn:
            if binsize == 0:
                ups_gcc.columns = [f"{ssn}_intron"]
                dws_gcc.columns = [f"{ssn}_exon"]
            else:
                ups_gcc.columns = [f"{ssn}_intron_b{i}" for i in range(-ups_gcc.shape[1], 0)]
                dws_gcc.columns = [f"{ssn}_exon_b{i}" for i in range(dws_gcc.shape[1])]
        df = pd.concat([ups_gcc, dws_gcc], axis=1)
        # df.index = df_up.iloc[:, 3]
        df.index = df_up.index
        gcc_mean = df.mean()
        gcc_mean.index = range(-ups_gcc.shape[1], dws_gcc.shape[1])
        axes[idx].plot(gcc_mean)
        axes[idx].set_ylim([min(0.35, min(gcc_mean)), max(0.65, max(gcc_mean))])
        df_ls.append(df)
    dfs = pd.concat(df_ls, axis=1)
    dfs.to_csv(outdir / "gcc.csv")
    plt.tight_layout()
    plt.savefig(outdir / "gcc.png")


def get_elen(file, outdir, scale, mode, app):
    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)
    df_len = ul.get_all_ss_distance(file, mode, app)
    ylabel = "length"
    if scale == "log":
        ylabel = "log2(length)"
        df_len = np.log2(df_len)
    elif scale == "MinMaxScaler":
        from sklearn.preprocessing import MinMaxScaler
        df_len.iloc[:] = MinMaxScaler().fit_transform(np.log2(df_len))
        ylabel = "scale(length)"
    df_len.to_csv(outdir / "element_len.csv")
    fig, axes = plt.subplots(1, df_len.shape[1], sharey=True)
    for idx in range(df_len.shape[1]):
        data = df_len.iloc[:, idx]
        axes[idx].boxplot(data)
        axes[idx].set_ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outdir / "element_len.png")

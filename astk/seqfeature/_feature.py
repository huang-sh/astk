from pathlib import Path
import multiprocessing as mp

import astk.utils as ul


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
    from pandas import DataFrame

    with open(fasta, "r") as rh:
        seqs = ul.read_fasta(rh)
        pool = mp.Pool(process)
        gccs = [pool.apply(_get_seq_gcc, args=(seq[1], binsize)) for seq in seqs]
        pool.close()
    
    return DataFrame(gccs)


def get_gcc(file, outdir, gfasta, binsize, **kwargs):
    import matplotlib.pyplot as plt
    from pandas import concat
    
    if app := kwargs.get("app", "auto") == "auto":
        app = ul.detect_file_info(file)["app"]

    kwargs = {
             "split": True,
             "excludeSS": not kwargs["includess"],
             "exon_width": kwargs["exonflank"],
             "intron_width": kwargs["intronflank"]
             }
    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)
    coord_dic = ul.get_ss_bed(file, app=app, **kwargs)

    fig, axes = plt.subplots(1, len(coord_dic), figsize= (10, 5))
    for idx, (ssn, (df_up, df_dw)) in enumerate(coord_dic.items()):
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
        df = concat([ups_gcc, dws_gcc], axis=1)
        gcc_mean = df.mean()
        gcc_mean.index = range(-ups_gcc.shape[1], dws_gcc.shape[1])
        axes[idx].plot(gcc_mean)
        axes[idx].set_ylim([min(0.35, min(gcc_mean)), max(0.65, max(gcc_mean))])

    plt.tight_layout() 
    plt.savefig(outdir / "gcc.png")

def get_elen(file, outdir, app, log):
    import matplotlib.pyplot as plt
    import numpy as np

    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)

    df_len = ul.get_ss_range(file, app)
    df_len.to_csv(outdir / "element_len.csv")
    ylabel = "log2(length)" if log else "length"
    fig, axes = plt.subplots(1, df_len.shape[1], sharey=True)
    for idx in range(df_len.shape[1]):
        if log:
            data = np.log2(df_len.iloc[:, idx])
        else:
            data = df_len.iloc[:, idx]
        axes[idx].boxplot(data)
        axes[idx].set_ylabel(ylabel)
    plt.tight_layout() 
    plt.savefig(outdir / "element_len.png")

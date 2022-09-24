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

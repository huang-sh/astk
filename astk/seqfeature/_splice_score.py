"""
This is the Python implementation of MaxEntScan perl scripts
The perl scripts is from http://hollywood.mit.edu/burgelab/maxent/download/fordownload/
"""
import json
from math import log2
from pathlib import Path
import multiprocessing as mp

from astk.constant import BASE_DIR
from astk.utils._getfasta import get_coor_fa
from astk.utils import get_ss_bed, read_fasta, detect_file_info


def ss5_consensus_score(seq):
    bdg = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
    cons1 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
    cons2 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}
    addscore = cons1[seq[3]] * cons2[seq[4]] /( bdg[seq[3]] * bdg[seq[4]])
    return addscore


def ss5_score(seq):
    ss5_score_json = BASE_DIR / "data/maxent/ss5seq_score.json"
    with open(ss5_score_json, "r") as f:
        ss5_score_dic = json.load(f)    
    seq = seq.upper()[:9]
    score = log2(ss5_consensus_score(seq) * ss5_score_dic.get(seq[:3]+seq[5:], 1e-12))
    return score


def ss5_seqs_score(fasta, process=4):
    with open(fasta, "r") as rh:
        seqs = read_fasta(rh)
        pool = mp.Pool(process)
        scores = [pool.apply(ss5_score, args=(seq[1],)) for seq in seqs]
        pool.close()
    return scores


def hashseq(seq):
    seqn = {"A": 0, "C": 1, "G": 2, "T": 3}
    four = (1,4,16,64,256,1024,4096,16384)
    num = 0
    slen = len(seq)
    for i in range(0, slen):
        num += seqn[seq[i]] * four[slen-i-1]
    return num


def ss3_consensus_score(seq):
    bdg = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
    cons1 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
    cons2 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}    
    addscore = cons1[seq[18]] * cons2[seq[19]]/ (bdg[seq[18]] * bdg[seq[19]]) 
    return addscore


def ss3_maxentscore(seq):
    ss3_score_json = BASE_DIR / "data/maxent/ss3_score.json"
    with open(ss3_score_json, "r") as f:
        ss3_score_dic = json.load(f)
    sc = [0] * 9
    sc[0] = ss3_score_dic['0'][hashseq(seq[0:7])]
    sc[1] = ss3_score_dic['1'][hashseq(seq[7:7+7])]
    sc[2] = ss3_score_dic['2'][hashseq(seq[14:14+7])]
    sc[3] = ss3_score_dic['3'][hashseq(seq[4:4+7])]
    sc[4] = ss3_score_dic['4'][hashseq(seq[11:11+7])]
    sc[5] = ss3_score_dic['5'][hashseq(seq[4:4+3])]
    sc[6] = ss3_score_dic['6'][hashseq(seq[7:7+4])]
    sc[7] = ss3_score_dic['7'][hashseq(seq[11:11+3])]
    sc[8] = ss3_score_dic['8'][hashseq(seq[14:14+4])]
    finalscore = sc[0]*sc[1]*sc[2]*sc[3]*sc[4]/(sc[5]*sc[6]*sc[7]*sc[8])
    return finalscore


def ss3_score(seq):
    seq = seq.upper()[:23]
    score = log2(ss3_consensus_score(seq) * ss3_maxentscore(seq[:18]+seq[20:]))
    return score


def ss3_seqs_score(fasta, process=4):
    with open(fasta, "r") as rh:
        seqs = read_fasta(rh)
        # scores = [ss3_score(seq[1]) for seq in seqs]
        pool = mp.Pool(process)
        scores = [pool.apply(ss3_score, args=(seq[1],)) for seq in seqs]
        pool.close()
    return scores     


def splice_score(file, outdir, gfasta, app, process):
    from pandas import DataFrame

    if app == "auto":
        app = detect_file_info(file)["app"]

    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)
    coord_dic = get_ss_bed(file, sss=True, app=app)

    df_score = DataFrame()
    for idx, (ssn, df) in enumerate(coord_dic.items()):
        ss_dir = outdir / ssn
        ss_dir.mkdir(exist_ok=True)
        get_coor_fa(df, gfasta, ss_dir / f"{ssn}.fa", strandedness=True)
        df.to_csv(ss_dir / f"{ssn}.bed", index=False, header=False, sep="\t")
        if ssn.endswith("5SS"):
            df_score[ssn] = ss5_seqs_score(ss_dir / f"{ssn}.fa", process)
        elif ssn.endswith("3SS"):
            df_score[ssn] = ss3_seqs_score(ss_dir / f"{ssn}.fa", process)

    df_score.index = df.iloc[:, 3]
    df_score.to_csv(outdir / "splice_scores.csv")
    ax = df_score.plot.box()
    fig = ax.get_figure()
    fig.savefig(outdir / "splice_scores_box.png")

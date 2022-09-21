from pathlib import Path
from itertools import product
from functools import partial
from collections import Counter
from typing import Tuple, Union
import subprocess

import numpy as np

from astk.constant import BASE_DIR
import astk.utils.func  as ul


def seq_aac(seq, k=1, gap=0, lam=0, count=False):
    """ extract amino acid sequence composition feature
    :param seq: an amino acid seq 
    :param raa: representative aa, list
    :param k: 
    :param gap: 
    :param lam: 
    :param count:
    :return:
    """
    NT = ["A", "T", "G", "C"]
    def three_f(idx, k, gap, lam, seq_dic):
        a = (seq_dic.get(idx+idx*gap+(lam+1)*i, '') for i in range(k))
        return ''.join(a)
    seq_dic = {i: v for i,v in enumerate(seq)}
    nt = [''.join(n) for n in product(NT, repeat=k)]
    f3 = partial(three_f, gap=gap, k=k, lam=lam, seq_dic=seq_dic)
    nt_list = [f3(i) for i in range(len(seq))]
    nt_list = [i for i in nt_list if len(i) ==k]
    nt_dict = Counter(nt_list)
    all_count = len(nt_list)
    if count:
        all_count = 1 
    nt_fre = [nt_dict[i] / all_count for i in nt]
    return nt_fre


def read_fasta(seq):
    seq_ls = []
    for line in seq:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            if seq_ls:
                yield descr, ''.join(seq_ls)
                seq_ls.clear()
            descr = line
        else:
            seq_ls.append(line)
    else:
        yield descr, ''.join(seq_ls)


def extract_feature(file, k, gap, lam, count=False):
    with open(file, "r") as rh:
        seqs = read_fasta(rh)
        fea_func = partial(seq_aac, k=k, gap=gap, lam=lam, count=count)
        seq_vec = np.array([fea_func(sq[1]) for sq in seqs])
    return seq_vec


def write_array(file, *data, header=None):
    data = np.hstack(data)
    if header:
        np.savetxt(file, data, delimiter=",", fmt="%.6f", header=header, comments='')
    else:
        np.savetxt(file, data, delimiter=",", fmt="%.6f")
    
def seq_extract(fasta: Tuple[str],
                output: str,
                kmer: int,
                gap: int,
                lambda_: int,
                count: bool,
                featuremerge: bool
    ):
    iscount = count
    k, gap, lam = kmer, gap, lambda_

    xy_ls = []
    # aa_ls = [''.join(aa) for aa in product(raa, repeat=k)]
    for idx, file in enumerate(fasta):
        feature_file = Path(file)
        xy = extract_feature(feature_file, k, gap, lam, count=iscount)
        # new_aa_ls = aa_ls

        xy_ls.append(xy)
    # new_aa_ls.insert(0, 'label') if args.label_f else 0
    # header = ','.join(new_aa_ls)
    if featuremerge:
        write_array(output, *xy_ls)
    #     exit()
    # for idx, o in enumerate(args.output):
    #     write_array(Path(o), xy_ls[idx])


def seq_pcm(fasta: Union[str, Path], 
            output: Union[str, Path],
            fmt: str,
            width: int,
            height: int,
            resolution: int
    ):
    with open(fasta, "r") as rh:
        seqs = [seq for _, seq in read_fasta(rh)]
        seq_len = list(set(map(len, seqs)))
        if len(seq_len) > 1:
            print("fasta sequence lengths are not consistent")
            exit()
        nc_ls = []
        for pos in range(seq_len[0]):
            NC = Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
            nts = [seq[pos] for seq in seqs]
            NC.update(nts)
            nc_ls.append(NC)

    pwm = Path(output).with_suffix(".txt")
    with open(pwm, "w") as f:
        header = f">{Path(fasta).stem}\n"
        f.write(header)
        for pc in nc_ls:
           line = "\t".join([str((i+1)/(len(seqs)+1)) for i in pc.values()])
           f.write(line)
           f.write("\n")

    rscript = BASE_DIR / "R" / "seqLogo.R"

    param_dic = {
        "pwm":pwm,
        "fmt": fmt,
        "width": width, 
        "height": height,
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])

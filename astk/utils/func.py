import sys
import mmap
import json
from pathlib import Path
from functools import partial

from astk.types import FilePath
from astk.event import SuppaEventID
from astk.constant import OrgDb_dic, BASE_DIR, SSN, RBP_sp_dic


class RunConfigure:
    def __init__(self) -> None:
        self.path = Path.home() / ".astkrc"
        self.rc_dic = self.read()
    
    def read(self, path=None):
        if path is None:
            astkrc = self.path
        else:
            astkrc = path
            self.path = path
        if astkrc.exists():
            with open(astkrc, "r") as f:
                rc_dic = json.load(f)
        else:
            rc_dic = {}
        return rc_dic
    
    def update(self, **kwargs):
        self.rc_dic.update(**kwargs)

    def save(self):
        with open(self.path, "w") as f:
            json.dump(self.rc_dic, f, indent=4)
     

def Rscript_bin():
    rc = RunConfigure()
    path = rc.rc_dic.get("Rscript", "Rscript")
    return path


def sig_filter(df, dpsi=0, abs_dpsi=0, pval=0.05):
    dpsi_df = df
    pval_keep = dpsi_df.loc[:, "pval"] < pval
    if dpsi > 0:
        keep = (dpsi_df.loc[:, "dpsi"] > dpsi) & pval_keep
    elif dpsi < 0:
        keep = (dpsi_df.loc[:, "dpsi"] < dpsi) & pval_keep
    else:
        keep = pval_keep
    if abs_dpsi > 0:
        keep = (abs(dpsi_df.loc[:, "dpsi"]) > abs_dpsi) & keep
    fdf = dpsi_df.loc[keep, ]
    return fdf


def compute_feature_len(df):
    import pandas as pd

    AS_len = lambda x: SuppaEventID(x).alter_element_len
    lens = df["event_id"].apply(AS_len)
    len_df = pd.DataFrame(lens.tolist())
    return len_df


def extract_info(df):
    import pandas as pd

    AS_len = lambda x: SuppaEventID(x).alter_element_len
    chrs = df["event_id"].apply(lambda x: x.split(":")[1])
    genes = df["event_id"].apply(lambda x: x.split(";")[0].strip())
    lens = df["event_id"].apply(AS_len)
    len_df = pd.DataFrame(lens.tolist(), columns=["len"])
    data = {
        "geneId": genes,
        "chr": chrs,
    }
    new_df = pd.DataFrame(data)
    len_df.index = new_df.index
    new_df = pd.concat([new_df, len_df], axis=1)
    return new_df


def len_hist(len_counts, width, max_len):
    import numpy as np

    offset = (width - max_len % width) if max_len % width else 0
    bins =  (max_len + offset) // width
    counts, bin_edges = np.histogram(len_counts, bins=bins, range=(1, max_len+offset+1))
    return counts, bin_edges


def plot_hist_cluster(out, cluster_ls, bins):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for i, subset in enumerate(cluster_ls, 1):
        ax.hist(subset, bins=bins, alpha=0.5, label=f"Cluster {i}")
    ax.legend()
    plt.savefig(out)

                                                         
def custome_cluster_len(df, out, lens, width=10, max_len=500):
    import pandas as pd

    len_count = df["len"]
    counts, bin_edges = len_hist(len_count, width, max_len)
    cluster_ls = []
    start = 0
    for idx, item in enumerate(lens):
        tmp = len_count[start < len_count]
        subset = tmp[tmp <= item]
        cluster_ls.append(subset)

        df.loc[(start < len_count) & (len_count <= item), "cluster"] = idx+1
        start = item
    else:
        cluster_ls.append(len_count[len_count > item])
        df.loc[len_count >= item, "cluster"] = idx + 2

    ## cluster info saving
    # lps = [1, *lens, max(bin_edges)]
    # cn_ls = [f"cluster{idx+1}" for i in range(len(lens)+1)]
    # range_ls = [f"{lps[i]}-{lps[i+1]}" for i in range(len(lps)-1)]
    # count_ls = [len(i) for i in cluster_ls]
    # cluster_info = pd.DataFrame({"cluster": cn_ls, "range": range_ls, "counts": count_ls})
    # cluster_info.to_csv(Path(out).with_suffix(".cluster.csv"))

    plot_hist_cluster(out, cluster_ls, bin_edges)
    return df


def cluster_len(df, out, n_cls=5, width=10, max_len=500, len_weight=5):
    from sklearn.cluster import KMeans
    import pandas as pd
    import numpy as np

    len_count = df["len"]
    counts, bin_edges = len_hist(df["len"], width, max_len)
    lens = [(bin_edges[i] + bin_edges[i+1])/2-1  for i,v in enumerate(bin_edges[:-1])]
    lens = np.array(lens).reshape(-1,1) * len_weight
    counts = counts.reshape(-1,1)
    data = np.concatenate([counts, lens], axis = 1)

    cluster_id = KMeans(n_cls, random_state=2).fit_predict(data)
    new_cls_dic = {}

    range_ls = []
    for ci in np.unique(cluster_id):
        cls_idx = list(np.where(cluster_id == ci)[0])  
        start = bin_edges[cls_idx[0]]
        end = bin_edges[cls_idx[-1] + 1]
        tmp = len_count[start <= len_count]
        subset = tmp[tmp < end]
        new_cls_dic[start] = subset
        range_ls.append([start, end])

    cluster_ls = [new_cls_dic[i] for i in sorted(new_cls_dic.keys())]

    plot_hist_cluster(out, cluster_ls, bin_edges)

    range_ls = sorted(range_ls, key=lambda x: x[0])

    range_str_ls = [f"{i[0]}-{i[1]}" for i in range_ls]
    count_ls = [len(i) for i in cluster_ls]

    for idx, (s, e) in enumerate(range_ls, 1):
        print(s, e)
        df.loc[(s <= df["len"]) & ( df["len"] < e), "cluster"] = idx
        print(df.head())
    else:
         df.loc[df["len"] >= e, "cluster"] = idx

    cn_ls = [f"cluster{i+1}" for i in range(len(cluster_ls))]
    cluster_info = pd.DataFrame({"cluster": cn_ls, "range": range_str_ls, "counts": count_ls})
    cluster_info.to_csv(Path(out).with_suffix(".cls.csv") )
    print(df.head())
    return df


def select_OrgDb(org):
    return OrgDb_dic.get(org, None)


def check_kegg_RData(org):
    import subprocess
    rscript = BASE_DIR / "R" / "dl_keggdata.R"
    info = subprocess.Popen([Rscript_bin(), str(rscript), org])
    info.wait()
    # import datetime
    # date_str = datetime.datetime.now().strftime("%Y-%m-%d")
    # rpath = ".".join([org, date_str, "kegg.RData"])
    # p = Path(rpath)
    # if p.exists() and p.stat().st_size > 10000:
    #     return True
    # else:
    #     return False


def get_coor(event_id, start, end, strand_sp, anchor, upstream_w, downstream_w):
    eid = SuppaEventID(event_id)
    if eid.strand == "-" and strand_sp:
        coordinates = eid.coordinates[::-1]
    else:
        coordinates = eid.coordinates
    if not any([start, end]) and not anchor:
        s, e = eid.alter_element_coor
    elif all([start, end]):
        s, e = sorted([coordinates[start-1], coordinates[end-1]])
    elif anchor:
        anchor_coor = coordinates[anchor-1]
        if eid.strand == "-" and strand_sp:
            e = anchor_coor + upstream_w
            s = anchor_coor - downstream_w
        else:
            s = anchor_coor - upstream_w             
            e = anchor_coor + downstream_w
    else:
        print("param error")
        sys.exit()
    return eid.Chr, s, e, event_id, 0, eid.strand


#TODO re-code it
def get_coor_bed(dpsi_file, start, end, strand_sp, anchor, upstream_w, downstream_w):
    import pandas as pd

    wget_coor_coor = partial(get_coor, start=start, end=end, strand_sp=strand_sp,
                    anchor=anchor, upstream_w=upstream_w, downstream_w=downstream_w)
    dpsi_df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
    dpsi_df.drop_duplicates(inplace=True)
    dpsi_df["event_id"] = dpsi_df.index
    coors = dpsi_df["event_id"].apply(wget_coor_coor)
    coor_df = pd.DataFrame(coors.tolist())
    
    return coor_df


def get_anchor_coor(event_id, index, sideindex, offset5, offset3, strand_sp):
    eid = SuppaEventID(event_id)
    if eid.strand == "-" and strand_sp:
        coordinates = eid.coordinates[::-1]
    else:
        coordinates = eid.coordinates

    if index is None and sideindex is None:
        s, e = eid.alter_element_coor
        anchor =  s + round((e - s) / 2)
    elif index:
        anchor = coordinates[index-1]
    elif sideindex:
        s, e = sorted([coordinates[i-1] for i in sideindex])
        anchor =  s + round((e - s) / 2)

    if eid.strand == "-" and strand_sp:
        new_anchor = anchor + offset5 - offset3
    else:
        new_anchor = anchor - offset5 + offset3

    return eid.Chr, new_anchor, eid.strand


def gen_anchor_bed(dpsi_file, out, index, sideindex, offset5, offset3, strand_sp):
    import pandas as pd
    from functools import partial

    wget_anchor_coor = partial(get_anchor_coor, index=index, 
            sideindex=sideindex, offset5=offset5, offset3=offset3, 
            strand_sp=strand_sp)
    dpsi_df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
    dpsi_df.drop_duplicates(inplace=True)
    dpsi_df["event_id"] = dpsi_df.index
    coors = dpsi_df["event_id"].apply(wget_anchor_coor)
    coor_df = pd.DataFrame(coors.tolist())
    coor_df.to_csv(out, index=False, header=False, sep="\t")


def get_evnet_ss_bed(
    event_file: FilePath, 
    ups_width: int, 
    dws_width: int
    ):
    from pandas import concat
    dic = {}
    etype = sniff_AS_type(event_file)
    for i in range(1, SSN[etype]+1):
        ps_coor_df = get_coor_bed(event_file, None, None, False, 
                        i, ups_width, dws_width)
        ns_coor_df = get_coor_bed(event_file, None, None, False, 
                        SSN[etype]+1-i, dws_width, ups_width)
        ai_coor_df = concat([
                        ps_coor_df.loc[ps_coor_df[5] == "+", ],
                        ns_coor_df.loc[ns_coor_df[5] == "-", ]])
        dic[f"A{i}"] = ai_coor_df
    return dic


def parse_cmd_r(**param_dic):
    param_ls = []
    for k, v in param_dic.items():
        if isinstance(v, (tuple, list)) and bool(v):
            param_ls.append(f"--{k}")
            param_ls.extend([str(i) for i in v])
        elif isinstance(v, bool):
            if v: param_ls.append(f"--{k}")
        elif v is None:
            pass
        else:
            param_ls.append(f"--{k}")
            param_ls.append(str(v))
    return param_ls


def get_meme(motif_ids, meme, output):
    with open(meme, "r") as hi, open(output, "w") as ho:
        flag = True
        for line in hi:
            if line.startswith("MOTIF"):
                mid = line.split()[1]
                if mid in motif_ids:
                    ho.write(line)
                    flag = True 
                else:
                    flag = False 
            elif flag:
                ho.write(line)


def sep_name(name, sep, *idx):
    str_ls = name.split(sep)
    name_ls = [str_ls[int(i)-1] for i in idx]
    return ".".join(name_ls)


def df_len_select(infile, outfile, s, e):
    import pandas as pd

    AS_len = lambda x: SuppaEventID(x).alter_element_len
    sep = sniff_file_sep(infile)
    df = pd.read_csv(infile, sep=sep, index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(outfile, index=True, sep=sep, na_rep="nan", index_label=False)


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def sniff_AS_type(file):
    with open(file, "r") as f:
        next(f)
        lines = f.readlines()
        events = [SuppaEventID(line.split()[0]) for line in lines]
        if len(set([e.AS_type for e in events])) != 1:
            exit()
        
    return events[0].AS_type


def get_meme_path():
    import os

    meme_path = Path(os.popen("which meme").read().strip())
    return str(Path(meme_path).parent)


def sniff_fig_fmt(file, fmts=None):
    if fmts is None:
        fmts = ['png', 'pdf', 'pptx']
    fsuffix = Path(file).suffix[1:]
    if fsuffix in fmts:
        suffix = fsuffix
    else:
        suffix = "png"
    return suffix


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


def detect_file_info(file):
    info_dic = {}
    with open(file, "r") as f:
        line1_ls = f.readline().split()
        ## just simple judgment
        if line1_ls[:5] == ["ID","GeneID","geneSymbol","chr","strand"]:
            info_dic["app"] = "rMATS"

            if "riExonStart_0base" in line1_ls:
                info_dic["etype"] = "RI"
            elif "2ndExonStart_0base" in line1_ls:
                info_dic["etype"] = "MXE"
            elif "exonStart_0base" in line1_ls:
                info_dic["etype"] = "SE"
            else:
                line2_ls = f.readline().split()
                if line2_ls[4] == "+":
                    if line2_ls[5] == line2_ls[7]:
                        info_dic["etype"] = "A5SS"
                    else:
                        info_dic["etype"] = "A3SS"
                else:
                    if line2_ls[5] == line2_ls[7]:
                        info_dic["etype"] = "A3SS"
                    else:
                        info_dic["etype"] = "A5SS"
        else:
            info_dic["app"] = "SUPPA2"
            lines = f.readlines()
            events = [SuppaEventID(line.split()[0]) for line in lines]
            if len(set([e.AS_type for e in events])) != 1:
                raise ValueError("input SUPPA2 must contain one AS type!")
            info_dic["etype"] = events[0].AS_type

    return info_dic                        


def sniff_file_sep(file):
    """simply check the sep of file
    """
    sep_dic = {}
    with open(file, "r") as f:
        line = f.readline()        
        sep_dic[","] = line.split(",")
        sep_dic["\t"] = line.split("\t")
    sep = sorted(sep_dic, key=lambda x: len(sep_dic[x]), reverse=True)[0]
    return sep


def merge_files(files, output, axis, rmdup, rmna):
    from pandas import read_csv, concat

    sep = sniff_file_sep(files[0])
    df_ls = []
    for file in files:
        df = read_csv(file, sep=sep, index_col=0)
        if axis == 0:
            df.columns = [f"c{i}" for i in range(df.shape[1])]
        df_ls.append(df)
    df = concat(df_ls, axis=axis)
    if rmdup == "all":        
        df["event_id"] = df.index
        df = df.drop_duplicates()
        del df["event_id"]
    elif  rmdup == "index":
        df["event_id"] = df.index
        df = df.drop_duplicates(subset=["event_id"])
        del df["event_id"]
    elif  rmdup == "content":
        df = df.drop_duplicates()
    if rmna:
        df.dropna(inplace=True)
    df.to_csv(output, index=True, sep=sep, na_rep="nan")


def shift2nease(file):
    import pandas as pd

    get_geneid = lambda x: SuppaEventID(x).gene_id.split(".")[0]
    get_start = lambda x: SuppaEventID(x).alter_element_coor[0]
    get_end = lambda x: SuppaEventID(x).alter_element_coor[1]

    df = pd.read_csv(file, sep="\t", index_col=0)
    df["event_id"] = df.index
    
    nease_input = pd.DataFrame({
                    "Gene stable ID": df["event_id"].apply(get_geneid),
                    "new_start": df["event_id"].apply(get_start), 
                    "new_end": df["event_id"].apply(get_end),
                    "beta": df.iloc[:, 0].values})
    return nease_input


def get_motif_path(database, organism):
    org = RBP_sp_dic.get(organism)
    motif = BASE_DIR / "data/motif" / database / f"{org}.meme"
    return motif

import sys
import mmap
from pathlib import Path
from functools import partial

from astk.constant import OrgDb_dic, BASE_DIR
# from astk.suppa.lib.gtf_store import *
# from astk.suppa.lib.tools import *
from . import event_id as ei


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
    from .event_id import SuppaEventID

    AS_len = lambda x: SuppaEventID(x).alter_element_len
    lens = df["event_id"].apply(AS_len)
    len_df = pd.DataFrame(lens.tolist())
    return len_df

def extract_info(df):
    import pandas as pd
    from .event_id import SuppaEventID

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
    plt.savefig(Path(out).with_suffix(".png"))  

                                                         
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

    lps = [1, *lens, max(bin_edges)]
    cn_ls = [f"cluster{idx+1}" for i in range(len(lens)+1)]
    range_ls = [f"{lps[i]}-{lps[i+1]}" for i in range(len(lps)-1)]
    count_ls = [len(i) for i in cluster_ls]
    cluster_info = pd.DataFrame({"cluster": cn_ls, "range": range_ls, "counts": count_ls})
    cluster_info.to_csv(Path(out).with_suffix(".cluster.csv"))

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
    info = subprocess.Popen(["Rscript", str(rscript), org])
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
    eid = ei.SuppaEventID(event_id)
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


def get_coor_fa(df, fasta, out):
    import pybedtools

    lines = [f"{v[1]}\t{v[2]-1}\t{v[3]}\t{v[4]}" for v in df.itertuples()]
    bed = pybedtools.BedTool("\n".join(lines), from_string=True)
    bed.sequence(fi=fasta, fo=out, name=True)


def get_anchor_coor(event_id, index, sideindex, offset5, offset3, strand_sp):
    eid = ei.SuppaEventID(event_id)
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

    AS_len = lambda x: ei.SuppaEventID(x).alter_element_len
    df = pd.read_csv(infile, sep="\t", index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(outfile, index=True, sep="\t", na_rep="nan", index_label=False)


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

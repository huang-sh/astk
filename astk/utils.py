# -*- coding: utf-8 -*-
"""
astk.utils
~~~~~~~~~~~~~~~~~
This module provides some utility functions.
"""
from os import name
import sys
import json
import hashlib
import logging
from pathlib import Path
from functools import partial
from itertools import chain, repeat
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

from .suppa.lib.event import make_events
from .suppa.lib.gtf_store import *
from .suppa.lib.tools import *
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
    
    lens = df["event_id"].apply(fl.AS_len)
    len_df = pd.DataFrame(lens.tolist())
    return len_df

def extract_info(df):
    from . import feature_len as fl
    AS_len = lambda x: fl.EventID(x).alter_element_len
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

def check_gtf_used(gtf):
    gtf = Path(gtf)
    gtf_hash = hashlib.blake2b()
    gtf_size = gtf.stat().st_size
    gtf_hash.update(str(gtf_size).encode('utf-8'))
    f = gtf.open("rb")
    for _ in range(10):
        gtf_hash.update(f.readline().strip())
    hash_val = gtf_hash.hexdigest()
    return hash_val

def ioe_event(gtf, event_type, edge_exon_len, pool_genes=False):

    mode = "logging.INFO" 

    # Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    # Setting the level of the loggers in lib
    # setToolsLoggerLevel(mode)

    home = Path.home()
    laas_dir = home / f".astk"
    laas_dir.mkdir(exist_ok=True)

    ref_dir = laas_dir / "ref"
    ref_dir.mkdir(exist_ok=True)

    gtf_hash = check_gtf_used(gtf)
    output_dir = ref_dir / gtf_hash
    
    if not output_dir.exists():
        output_dir.mkdir(exist_ok=True)
        my_genome = Genome()
        logger.info("Reading input data.")
        try:
            fetched_exons = gtf_reader(gtf, logger)
        except UnicodeDecodeError as e:
            print(e)
            sys.exit(1)

        if len(fetched_exons) == 0:
            logger.info("No exons found. Check format and content of your GTF file.")
            sys.exit(1)

        for exon_meta in fetched_exons:
            my_genome.add_to_genes(exon_meta)
        
        my_genome.sort_transcripts()

        if pool_genes:
            my_genome.split_genes()
            logger.info("Pooling genes")
            my_genome.poll_genes()

        out_prefix = output_dir / "annotation"
        make_events(event_type, my_genome, gtf, out_prefix, edge_exon_len,
                logger, b_type="S", th=10)
    else:
         logger.info("Loading cache data.")
    return output_dir  

def ioe_psi(event_row, tpm):
    alternative_transcripts = event_row["alternative_transcripts"].split(",")
    total_transcripts = event_row["total_transcripts"].split(",")
    al_tpts_val = [i for i in alternative_transcripts if i in tpm.index]
    all_tpts_val = [i for i in total_transcripts if i in tpm.index]

    if len(all_tpts_val) == 0:
        psi_ls = ["nan" for _ in tpm.columns]
    elif len(al_tpts_val) == 0 and len(all_tpts_val) > 0:
        psi_ls = [0 for _ in tpm.columns]
    else:
        get_val = lambda tid, colid: tpm.loc[tid, colid]
        psi_ls = []
        for colid in tpm.columns:
            get_val = lambda tid: tpm.loc[tid, colid]
            al_tpts_abundance = sum(map(get_val, al_tpts_val))
            all_tpts_abundance = sum(map(get_val, all_tpts_val))
            if all_tpts_abundance < 0.001:
                psi_ls.append("nan")
                continue
            try:
                psi = al_tpts_abundance / all_tpts_abundance
            except ZeroDivisionError:
                psi = "nan"
            psi_ls.append(psi)
    return psi_ls


def len_hist(len_counts, width, max_len):
    offset = (width - max_len % width) if max_len % width else 0
    bins =  (max_len + offset) // width
    counts, bin_edges = np.histogram(len_counts, bins=bins, range=(1, max_len+offset+1))
    return counts, bin_edges


def plot_hist_cluster(out, cluster_ls, bins):
    fig, ax = plt.subplots()
    for i, subset in enumerate(cluster_ls, 1):
        ax.hist(subset, bins=bins, alpha=0.5, label=f"Cluster {i}")
    ax.legend()
    plt.savefig(Path(out).with_suffix(".png"))  

                                                         
def custome_cluster_len(df, out, lens, width=10, max_len=500):
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


class DiffSplice:
    def __init__(self, outdir, meta, gtf, as_types,
                 exon_len, pool_genes) -> None:
        self.psi = {}
        self.mk_outdir(outdir)
        self.ioe_dir = self.get_ioe_event(gtf, exon_len, pool_genes)
        self.parse_meta(meta)
        for at in as_types:
            self.get_psi(at)
        self.write_data()
   
    @staticmethod
    def get_ioe_event(gtf, exon_len, pool_genes):
        events = ['SE', "SS", "MX", "RI", 'FL']
        ioe_dir = ioe_event(gtf, events, exon_len, pool_genes=pool_genes)
        return ioe_dir
    
    def parse_meta(self, meta):
        self.tpm = {}
        tpm_col = 4
        self.meta_file = meta
        with open(meta, "r") as f:
            try:
                meta_dic = json.load(f)
            except json.decoder.JSONDecodeError:
                print("ERROR: JSON file is excepted!")
                sys.exit()
            self.meta = meta_dic
        read_tpm = partial(self.read_tpm, tpm_col=tpm_col)
        for gn, gdic in meta_dic.items():
            control_tpms = [read_tpm(sp["path"], sp["name"]) for sp in gdic["control"]["samples"]]
            treatment_tpms = [read_tpm(sp["path"], sp["name"]) for sp in gdic["treatment"]["samples"]]
            control_tpm = pd.concat(control_tpms, axis=1)
            treatment_tpm = pd.concat(treatment_tpms, axis=1)
            self.tpm[gn] = [control_tpm, treatment_tpm]

    def get_psi(self, as_type):
        self.psi.setdefault(as_type, {})

        ioe_file = Path(self.ioe_dir) / f"annotation_{as_type}_strict.ioe"
        ioe_df = pd.read_csv(ioe_file, sep="\t")        
        for gn, (tpm1, tpm2) in self.tpm.items():
            _get_psi1 = partial(self._get_psi, tpm=tpm1)
            _get_psi2 = partial(self._get_psi, tpm=tpm2)
            self._get_psi(ioe_df.iloc[1, :], tpm1)

            psi1 = ioe_df.apply(_get_psi1, axis=1)
            psi2 = ioe_df.apply(_get_psi2, axis=1)

            psi1_df = pd.DataFrame(psi1.tolist(), 
                        columns=tpm1.columns, index=ioe_df["event_id"])
            psi2_df = pd.DataFrame(psi2.tolist(), 
                        columns=tpm2.columns, index=ioe_df["event_id"])
            self.psi[as_type][gn] = [psi1_df, psi2_df]

    @staticmethod
    def _get_psi(event_row, tpm):
        alter_tpt = event_row["alternative_transcripts"].split(",")
        total_tpt = event_row["total_transcripts"].split(",")
        al_tpts_val = [i for i in alter_tpt if i in tpm.index]
        all_tpts_val = [i for i in total_tpt if i in tpm.index]

        if len(all_tpts_val) == 0:
            psi_ls = ["nan" for _ in tpm.columns]
        elif len(al_tpts_val) == 0 and len(all_tpts_val) > 0:
            psi_ls = [0 for _ in tpm.columns]
        else:
            get_val = lambda tid, colid: tpm.loc[tid, colid]
            psi_ls = []
            for colid in tpm.columns:
                get_val = lambda tid: tpm.loc[tid, colid]
                al_tpts_abundance = sum(map(get_val, al_tpts_val))
                all_tpts_abundance = sum(map(get_val, all_tpts_val))
                if all_tpts_abundance < 0.001:
                    psi_ls.append("nan")
                    continue
                try:
                    psi = al_tpts_abundance / all_tpts_abundance
                except ZeroDivisionError:
                    psi = "nan"
                psi_ls.append(psi)
        return psi_ls


    @staticmethod
    def read_tpm(file, colname, tpm_col):
        quant_df = pd.read_csv(file, sep = "\t")
        tpm_df = pd.DataFrame({colname: quant_df.iloc[:, tpm_col-1]})
        tpm_df.index = quant_df.iloc[:, 0]
        return tpm_df

    def mk_outdir(self, outdir):
        outdir = Path(outdir)
        outdir.mkdir(exist_ok=True)
        tpm_dir = outdir / "tpm"
        tpm_dir.mkdir(exist_ok=True)
        psi_dir = outdir / "psi"
        psi_dir.mkdir(exist_ok=True)
        dpsi_dir =outdir / "dpsi"
        dpsi_dir.mkdir(exist_ok=True)

        self.outdir = outdir
        self.tpm_dir = tpm_dir
        self.psi_dir = psi_dir
        self.dpsi_dir = dpsi_dir

    def update_meta(self):
        meta = self.meta
        for as_type, group_dic in self.psi_files.items():
           
            for gn, (psi1, psi2) in group_dic.items():
                meta[gn].setdefault("dpsi", {})
                meta[gn]["control"].setdefault("psi", {})
                meta[gn]["control"]["psi"][as_type] = str(psi1)
                meta[gn]["treatment"].setdefault("psi", {})
                meta[gn]["treatment"]["psi"][as_type] = str(psi2)
                meta[gn]["control"]["tpm"] = str(self.tpm_files[gn][0])
                meta[gn]["treatment"]["tpm"] = str(self.tpm_files[gn][1])
                meta[gn]["dpsi"][as_type] = str(self.dpsi_files[as_type][gn])

        with open(self.meta_file, "w") as f:
            json.dump(meta, f, indent=4)

    def write_data(self):
        tpm_dir = self.tpm_dir
        psi_dir = self.psi_dir
        self.tpm_files = {}
        self.psi_files = {}
        for gn, (tpm1, tpm2) in self.tpm.items():
            tpm1_file = tpm_dir / f"{gn}_c1.tpm"
            tpm2_file = tpm_dir / f"{gn}_c2.tpm"
            tpm1.to_csv(tpm1_file, sep="\t")
            tpm2.to_csv(tpm2_file, sep="\t")
            self.tpm_files[gn] = [tpm1_file, tpm2_file]

        for as_type, group_dic in self.psi.items():
            self.psi_files.setdefault(as_type, {})
            for gn, (psi1, psi2) in group_dic.items():
                psi1_file = psi_dir / f"{gn}_{as_type}_c1.psi"
                psi2_file = psi_dir / f"{gn}_{as_type}_c2.psi"
                psi1.to_csv(psi1_file, sep="\t")
                psi2.to_csv(psi2_file, sep="\t")
                self.psi_files[as_type][gn] = [psi1_file, psi2_file]

    #TODO 
    def ds(self, method):
        import shutil
        try:
            shutil.copytree(self.ioe_dir, self.outdir / "ref")
        except FileExistsError as e:
            print(e)
            # print(f"please rm {self.outdir}/ref directory")
        self.method = method
        self.dpsi_files = {}
        from .suppa.lib.diff_tools import multiple_conditions_analysis
        mca = multiple_conditions_analysis

        for as_type, group_dic in self.psi_files.items():
            self.dpsi_files.setdefault(as_type, {})
            for gn, psi_files in group_dic.items():
                expr_files = self.tpm_files[gn]
                dpis_file = self.dpsi_dir / f"{gn}_{as_type}"
                ioe_file = self.ioe_dir / f"annotation_{as_type}_strict.ioe"
                mca(method, psi_files, expr_files, ioe_file, 1000, 0, 
                    False, True, 0.05, True, False, False, 1, 0, str(dpis_file))

                self.dpsi_files[as_type][gn] = dpis_file.with_suffix(".dpsi")
        
        self.update_meta()


OrgDb_dic = {
     "hs": "org.Hs.eg.db", "mm": "org.Mm.eg.db", "rn": "org.Rn.eg.db",
     "dm": "org.Dm.eg.db", "at": "org.At.tair.db", "sc": "org.Sc.sgd.db",
     "dr": "org.Dr.eg.db", "ce": "org.Ce.eg.db", "bt": "org.Bt.eg.db",
     "ss": "org.Ss.eg.db", "mmu": "org.Mmu.eg.db", "gg": "org.Gg.eg.db", 
     "cf": "org.Cf.eg.db", "eck12": "org.EcK12.eg.db", "xl": "org.Xl.eg.db",
     "pt": "org.Pt.eg.db", "ag": "org.Ag.eg.db", "pf": "org.Pf.plasmo.db",
     "ecsakai": "org.EcSakai.eg.db", "mxanthus": "org.Mxanthus.db"
     }

def select_OrgDb(org):
    return OrgDb_dic.get(org, None)


def check_kegg_RData(org):
    import subprocess
    rscript = Path(__file__).parent / "R" / "dl_keggdata.R"
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
    return eid.Chr, s, e, event_id


def get_coor_bed(dpsi_file, start, end, strand_sp, anchor, upstream_w, downstream_w):
    wget_coor_coor = partial(get_coor, start=start, end=end, strand_sp=strand_sp,
                    anchor=anchor, upstream_w=upstream_w, downstream_w=downstream_w)
    dpsi_df = pd.read_csv(dpsi_file, sep="\t", index_col=0)
    dpsi_df["event_id"] = dpsi_df.index
    coors = dpsi_df["event_id"].apply(wget_coor_coor)
    coor_df = pd.DataFrame(coors.tolist())
    coor_df.drop_duplicates(inplace=True)
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
    dpsi_df["event_id"] = dpsi_df.index
    coors = dpsi_df["event_id"].apply(wget_anchor_coor)
    coor_df = pd.DataFrame(coors.tolist())
    coor_df.drop_duplicates(inplace=True)
    coor_df.to_csv(out, index=False, header=False, sep="\t")


def parse_cmd_r(**param_dic):
    param_ls = []
    for k, v in param_dic.items():
        if isinstance(v, (tuple, list)) and bool(v):
            param_ls.append(f"--{k}")
            param_ls.extend(v)
        elif isinstance(v, bool):
            if v: param_ls.append(f"--{k}")
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
    AS_len = lambda x: ei.SuppaEventID(x).alter_element_len
    df = pd.read_csv(infile, sep="\t", index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(outfile, index=True, sep="\t", na_rep="nan", index_label=False)


def coSpliceNet(file, output):
    from collections import Counter
    from pyecharts import options as opts
    from pyecharts.charts import Graph
    import pandas as pd

    link_df = pd.read_csv(file)
    link_sort = link_df.sort_values(by=["cor", "pval"], key=abs, ascending=False)
    link_fil = link_sort.drop_duplicates(subset=["rbp", "event"])

    rbp_symbol = []
    event_id = []
    links = []

    for row in link_fil.itertuples():
        rbp = row[-1].split(":")[0]
        event = row[-2]
        rbp_symbol.append(rbp)
        event_id.append(event)
        links.append({"source": rbp, "target": event})

    rbp_nodes = []
    for k, v in Counter(rbp_symbol).items():
        node = {"name": k, "symbolSize": v, "category": k,
                'label': {'normal': {'show': 'True'}}}
        rbp_nodes.append(node)

    category = [{"name": k} for k in set(rbp_symbol)]
    event_nodes = [{"name": k, "symbolSize": max([v,5])} for k, v in Counter(event_id).items()]
    nodes = rbp_nodes + event_nodes

    c = (
        Graph().add(
            "",
            nodes,
            links,
            category,
            repulsion=50,
            linestyle_opts=opts.LineStyleOpts(curve=0.2),
            label_opts=opts.LabelOpts(is_show=False),
            edge_symbol=[None, 'arrow'],
            edge_symbol_size=4
        )
        .set_global_opts(
            legend_opts=opts.LegendOpts(is_show=False),
            title_opts=opts.TitleOpts(title="Co-splicing network"),
        )
        .render(output)
    )
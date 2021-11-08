from functools import partial

import pandas as pd
import pysam


def site_flanking(chrN, site, sam, control_sam=None, window=150, bins=15):

    def signal_func(chrN, start, end, sam, control_sam):
        if control_sam:
            treat = sam.count(chrN, start, end)
            control = control_sam.count(chrN, start, end)
            signal = treat / control
        else:
            signal = 1e6*sam.count(chrN, start, end)/sam.mapped
        return signal
    
    psignal_func = partial(signal_func, sam=sam, control_sam=control_sam)

    up_bin_start_idx = sorted(list(range(site+1, site-window, -bins)))[:-1]
    down_bin_end_idx = list(range(site+1, site+window+2, bins))[:-1]

    up_signal = [psignal_func(chrN, i-1, i+bins-1) for i in up_bin_start_idx]
    down_signal = [psignal_func(chrN, i-1, i+bins-1) for i in down_bin_end_idx]
    df = pd.concat([
        pd.DataFrame({"signal": up_signal, "direction": "upstream"}),
        pd.DataFrame({"signal": down_signal, "direction": "downstream"})
    ])
    return df 


def epi_signal(out, achor_dic, bam_meta, width, binsize):
    df = pd.read_csv(bam_meta)
    cds = set(df.condition)
    df_ls = []

    if len(cds) == 1:
        tdf_iter = [sdf for _,sdf in df.iterrows()]
        cdf_iter = [None] * df.shape[0]
    elif cds == {'control', 'treatment'}:
        tdf = df.loc[df.condition == "treatment", ]
        cdf = df.loc[df.condition == "control", ]
        if tdf.shape != cdf.shape:
            raise ValueError("control input number dismatched treatment input")
        tdf_iter = [sdf for _,sdf in tdf.iterrows()]
        cdf_iter = [sdf for _,sdf in cdf.iterrows()]
    
    for c, t in zip(cdf_iter, tdf_iter):
        if c is None:
            csam = None
        else:
            csam = pysam.AlignmentFile(c["path"])

        tsam = pysam.AlignmentFile(t["path"])

        for an in achor_dic:
            anchor_df = pd.read_csv(achor_dic[an], sep="\t", header=None)
            for ai, site_row in anchor_df.iterrows():
                chrN, pos = site_row[0], int(site_row[1])
                sdf = site_flanking(chrN, pos, tsam, control_sam=csam, window=width, bins=binsize)
                sdf["event_idx"] = ai
                sdf["mark"] = t["name"]
                sdf["rep"] = t["replicate"]
                sdf["anchor"] = an
                df_ls.append(sdf)
    dfs = pd.concat(df_ls)
    dfs["bin"] = dfs.index
    dfs.to_csv(out, index=False)

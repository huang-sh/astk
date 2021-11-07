import pandas as pd
import pysam


def site_flanking(chr, site, sam, control_sam=None, window=150, bins=15):
    if control_sam:
        pass
    total_read =  sam.mapped

    up_bin_start_idx = sorted(list(range(site+1, site-window, -bins)))[:-1]
    down_bin_end_idx = list(range(site+1, site+window+2, bins))[:-1]
    up_signal = [1e6*sam.count(chr, i-1, i+bins-1)/total_read for i in up_bin_start_idx]
    down_signal = [1e6*sam.count(chr, i-1, i+bins-1)/total_read for i in down_bin_end_idx]
    df = pd.concat([
        pd.DataFrame({"signal": up_signal, "direction": "upstream"}),
        pd.DataFrame({"signal": down_signal, "direction": "downstream"})
    ])
    return df 


def epi_signal(out, achor_dic, bam_meta, width, binsize):
    bam_df = pd.read_csv(bam_meta)
    df_ls = []
    for bi, bam_row in bam_df.iterrows():
        sam = pysam.AlignmentFile(bam_row["path"])
        for an in achor_dic:
            anchor_df = pd.read_csv(achor_dic[an], sep="\t", header=None)
            for ai, site_row in anchor_df.iterrows():
                chrN, pos = site_row[0], int(site_row[1])
                sdf = site_flanking(chrN, pos, sam, window=width, bins=binsize)
                sdf["event_id"] = ai
                sdf["mark"] = bam_row["group"]
                sdf["anchor"] = an
                df_ls.append(sdf)
        sam.close()
    dfs = pd.concat(df_ls)
    dfs.to_csv(out)

from astk.utils import rMATSEventCoord 
from astk.constant import rMATS_POS_COLS
import pandas as pd

def test_sss_rMATSEventCoord():

    files = {
        "A3SS": "data/A3SS.MATS.JC.txt",
        "A5SS": "data/A5SS.MATS.JCEC.txt",
        "SE": "data/SE.MATS.JC.txt",
        "RI": "data/RI.MATS.JCEC.txt",
        "MXE":"data/MXE.MATS.JC.txt"
    }
    etype_ss_names = {
        "A3SS": ["A1_3SS", "A2_5SS", "A3_3SS", "A4_3SS", "A5_5SS"],
        "A5SS": ["A1_3SS", "A2_5SS", "A3_5SS", "A4_3SS", "A5_5SS"],
        "SE": ["A1_3SS", "A2_5SS", "A3_3SS", "A4_5SS", "A5_3SS", "A6_5SS"],
        "RI": ["A1_3SS", "A2_5SS", "A3_3SS", "A4_5SS"],
        "MXE": ["A1_3SS", "A2_5SS", "A3_3SS", "A4_5SS", "A5_3SS", "A6_5SS", "A7_3SS", "A8_5SS"]
    }
    for etype, file in files.items():
        dpsi_df = pd.read_csv(file, sep="\t", index_col=0)
        dpsi_df.drop_duplicates(inplace=True)
        dpsi_df.dropna(inplace=True)

        # tansform 1-based index for unified management
        if etype == "SE":
            end_col = ["upstreamEE", "exonEnd", "downstreamEE"]
        elif etype in ["A5SS", "A3SS"]:
            end_col = ["longExonEnd", "flankingEE", "shortEE"]
        elif etype == "MXE":
            end_col = ["upstreamEE", "1stExonEnd", "2ndExonEnd", "downstreamEE"]
        elif etype == "RI":
            end_col = ["upstreamEE", "downstreamEE"]
        dpsi_df.loc[:, end_col] -= 1

        ps_df = dpsi_df.loc[dpsi_df["strand"] == "+", ]
        ns_df = dpsi_df.loc[dpsi_df["strand"] == "-", ]

        eventcoor = rMATSEventCoord(file)
        df_dic = eventcoor.get_all_flank_bed(sss=True)
        assert list(df_dic.keys()) == etype_ss_names[etype]
        for idx,( key, df) in enumerate(df_dic.items()):
            if key.endswith("3SS"):
                assert set(df["end"] - df["start"]) == {23}
                assert all(ps_df.loc[:, rMATS_POS_COLS[etype]["+"][idx]] - 20 == df.loc[ps_df.index, "start"])
                assert all(ns_df.loc[:, rMATS_POS_COLS[etype]["-"][idx]] - 2 == df.loc[ns_df.index, "start"])                
                assert all(ps_df.loc[:, rMATS_POS_COLS[etype]["+"][idx]] + 2 == df.loc[ps_df.index, "end"]) -1
                assert all(ns_df.loc[:, rMATS_POS_COLS[etype]["-"][idx]] + 20 == df.loc[ns_df.index, "end"]) -1
            if key.endswith("5SS"):
                assert set(df["end"] - df["start"]) == {9}
                assert all(ps_df.loc[:, rMATS_POS_COLS[etype]["+"][idx]] - 2 == df.loc[ps_df.index, "start"])
                assert all(ns_df.loc[:, rMATS_POS_COLS[etype]["-"][idx]] - 6 == df.loc[ns_df.index, "start"]) 
                assert all(ps_df.loc[:, rMATS_POS_COLS[etype]["+"][idx]] + 6 == df.loc[ps_df.index, "end"]) -1
                assert all(ns_df.loc[:, rMATS_POS_COLS[etype]["-"][idx]] + 2 == df.loc[ns_df.index, "end"]) -1

from astk.utils import SuppaEventCoord 
from astk.event import SuppaEventID
from astk.constant import rMATS_POS_COLS

import pandas as pd

def test_sss_SuppaEventCoord():

    files = {
        "A3": "data/fb_e11_12_A3.dpsi",
        "A5": "data/fb_e11_12_A5.sig.dpsi",
        "SE": "data/fb_e11_12_SE.sig.dpsi",
        "RI": "data/fb_e11_12_RI.sig.dpsi",
        "MX":"data/fb_e11_12_MX.sig.dpsi",
        "AF":"data/fb_e11_12_AF.sig.dpsi",
        "AL":"data/fb_e11_12_AL.sig.dpsi",
    }
    etype_ss_names = {
        "A3": ["A0_5SS", "A1_3SS", "A2_3SS"],
        "A5": ["A0_5SS", "A1_5SS", "A2_3SS"],
        "SE": ["A0_5SS", "A1_3SS", "A2_5SS", "A3_3SS"],
        "RI": ["A0_3SS", "A1_5SS", "A2_3SS", "A3_5SS"],
        "MX": ["A0_5SS", "A1_3SS", "A2_5SS", "A3_3SS", "A4_5SS", "A5_3SS"],
        "AF": ["A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS"],
        "AL": ["A0_5SS", "A1_3SS", "A2_5SS", "A3_3SS"],
    }

    def _get_coor(event_id):
        ei = SuppaEventID(event_id)
        if ei.strand == "-":
            coords = reversed(ei.coordinates)
        else:
            coords = ei.coordinates
        row = ei.Chr, ei.gene_id, ei.strand, *coords
        return row


    for etype, file in files.items():
        dpsi_df = pd.read_csv(file, sep="\t", index_col=0)
        dpsi_df.drop_duplicates(inplace=True)
        dpsi_df.dropna(inplace=True)
        dpsi_df["event_id"] = dpsi_df.index
        event_df = pd.DataFrame(dpsi_df["event_id"].apply(_get_coor).tolist(),
                            index=dpsi_df.index)
        df_ss = event_df.iloc[:, 2:]

        ps_df = df_ss.loc[df_ss.iloc[:,0] == "+", ]
        ns_df = df_ss.loc[df_ss.iloc[:,0] == "-", ]
        del ps_df[2]
        del ns_df[2]

        eventcoor = SuppaEventCoord(file)
        df_dic = eventcoor.get_all_flank_bed(sss=True)
        assert list(df_dic.keys()) == etype_ss_names[etype]

        for idx,( key, df) in enumerate(df_dic.items()):
            if key.endswith("3SS"):
                assert set(df["end"] - df["start"]) == {23}
                if etype == "AF":                    
                    assert all(ps_df.iloc[:, idx+1] - 20 - 1 == df.loc[ps_df.index, "start"])
                    assert all(ns_df.iloc[:, idx+1] - 2 - 1 == df.loc[ns_df.index, "start"])
                    assert all(ps_df.iloc[:, idx+1] + 2 == df.loc[ps_df.index, "end"])
                    assert all(ns_df.iloc[:, idx+1] + 20 == df.loc[ns_df.index, "end"])
                else:
                    assert all(ps_df.iloc[:, idx] - 20 - 1 == df.loc[ps_df.index, "start"])
                    assert all(ns_df.iloc[:, idx] - 2 - 1 == df.loc[ns_df.index, "start"])
                    assert all(ps_df.iloc[:, idx] + 2 == df.loc[ps_df.index, "end"])
                    assert all(ns_df.iloc[:, idx] + 20 == df.loc[ns_df.index, "end"])                    
            if key.endswith("5SS"):
                assert set(df["end"] - df["start"]) == {9}
                if etype == "AF":
                    assert all(ps_df.iloc[:, idx+1] - 2 -1 == df.loc[ps_df.index, "start"])
                    assert all(ns_df.iloc[:, idx+1] - 6 -1 == df.loc[ns_df.index, "start"]) 
                    assert all(ps_df.iloc[:, idx+1] + 6 == df.loc[ps_df.index, "end"])
                    assert all(ns_df.iloc[:, idx+1] + 2 == df.loc[ns_df.index, "end"])
                else:
                    assert all(ps_df.iloc[:, idx] - 2 -1 == df.loc[ps_df.index, "start"])
                    assert all(ns_df.iloc[:, idx] - 6 -1 == df.loc[ns_df.index, "start"]) 
                    assert all(ps_df.iloc[:, idx] + 6 == df.loc[ps_df.index, "end"])
                    assert all(ns_df.iloc[:, idx] + 2 == df.loc[ns_df.index, "end"])

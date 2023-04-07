from astk.utils import SuppaEventCoord 
from astk.event import SuppaEventID
from astk.constant import rMATS_POS_COLS

import pandas as pd

def test_sss_SuppaEventCoord(shared_datadir):  # sourcery skip: assign-if-exp

    files = {
        "A3": f"{shared_datadir}/fb_e11_12_A3.dpsi",
        "A5": f"{shared_datadir}/fb_e11_12_A5.sig.dpsi",
        "SE": f"{shared_datadir}/fb_e11_12_SE.sig.dpsi",
        "RI": f"{shared_datadir}/fb_e11_12_RI.sig.dpsi",
        "MX":f"{shared_datadir}/fb_e11_12_MX.sig.dpsi",
        "AF":f"{shared_datadir}/fb_e11_12_AF.sig.dpsi",
        "AL":f"{shared_datadir}/fb_e11_12_AL.sig.dpsi",
    }
    etype_ss_names = {
        "A3": ["A1_5SS", "A2_3SS", "A3_3SS"],
        "A5": ["A1_5SS", "A2_5SS", "A3_3SS"],
        "SE": ["A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS"],
        "RI": ["A1_3SS", "A2_5SS", "A3_3SS", "A4_5SS"],
        "MX": ["A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS", "A5_5SS", "A6_3SS"],
        "AF": ['A1_pse_3SS', 'A2_5SS', 'A3_pse_3SS', 'A4_5SS', 'A5_3SS'],
        "AL": ['A1_5SS', 'A2_3SS', 'A3_pse_5SS', 'A4_3SS', 'A5_pse_5SS'],
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
        del df_ss[2]

        eventcoor = SuppaEventCoord(file)
        df_dic = eventcoor.get_all_flank_bed(sss=True)
        assert list(df_dic.keys()) == etype_ss_names[etype]

        for idx,( key, df) in enumerate(df_dic.items()):
            if key.endswith("3SS"):
                assert set(df["end"] - df["start"]) == {23}

                assert all(ps_df.iloc[:, idx] - 20 - 1 == df.loc[ps_df.index, "start"])
                assert all(ns_df.iloc[:, idx] - 2 - 1 == df.loc[ns_df.index, "start"])
                assert all(ps_df.iloc[:, idx] + 2 == df.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] + 20 == df.loc[ns_df.index, "end"])                    
            if key.endswith("5SS"):
                assert set(df["end"] - df["start"]) == {9}
                assert all(ps_df.iloc[:, idx] - 2 -1 == df.loc[ps_df.index, "start"])
                assert all(ns_df.iloc[:, idx] - 6 -1 == df.loc[ns_df.index, "start"]) 
                assert all(ps_df.iloc[:, idx] + 6 == df.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] + 2 == df.loc[ns_df.index, "end"])

        for idx, key in enumerate(etype_ss_names[etype]):
            df_up_bed, df_dw_bed = eventcoor.get_ss_flank_bed(idx, ups_w=100, dws_w=100, split=True, excludeSS=False)
            assert set(df_up_bed["end"] - df_up_bed["start"]) == {100} 
            assert {100} == set(df_dw_bed["end"] - df_dw_bed["start"])
            if key.endswith("3SS"):                    
                assert all(ps_df.iloc[:, idx] - 100 -1 == df_up_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] -1 == df_up_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx]  == df_up_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] + 100 == df_up_bed.loc[ns_df.index, "end"])

                assert all(ps_df.iloc[:, idx] - 1 == df_dw_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] + 100 -1 == df_dw_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 100  == df_dw_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] == df_dw_bed.loc[ns_df.index, "end"])
            if key.endswith("5SS"):
                assert all(ps_df.iloc[:, idx] - 100 == df_up_bed.loc[ps_df.index, "start"])                
                assert all(ps_df.iloc[:, idx] == df_up_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 1  == df_up_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] + 100 - 1  == df_up_bed.loc[ns_df.index, "end"])

                assert all(ps_df.iloc[:, idx]  == df_dw_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] +100 == df_dw_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 100 - 1  == df_dw_bed.loc[ns_df.index, "start"])        
                assert all(ns_df.iloc[:, idx] - 1  == df_dw_bed.loc[ns_df.index, "end"])

            df_up_bed, df_dw_bed = eventcoor.get_ss_flank_bed(idx, exon_width=100, intron_width=200, split=True, excludeSS=False)
            if key.endswith("3SS"):
                assert all(ps_df.iloc[:, idx] - 200 -1 == df_up_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] -1 == df_up_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx]  == df_up_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] + 200 == df_up_bed.loc[ns_df.index, "end"])

                assert all(ps_df.iloc[:, idx] - 1 == df_dw_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] + 100 -1 == df_dw_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 100  == df_dw_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] == df_dw_bed.loc[ns_df.index, "end"])
            if key.endswith("5SS"):
                assert all(ps_df.iloc[:, idx] - 100 == df_up_bed.loc[ps_df.index, "start"])                
                assert all(ps_df.iloc[:, idx] == df_up_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 1  == df_up_bed.loc[ns_df.index, "start"])
                assert all(ns_df.iloc[:, idx] + 100 - 1  == df_up_bed.loc[ns_df.index, "end"])

                assert all(ps_df.iloc[:, idx]  == df_dw_bed.loc[ps_df.index, "start"])
                assert all(ps_df.iloc[:, idx] + 200 == df_dw_bed.loc[ps_df.index, "end"])
                assert all(ns_df.iloc[:, idx] - 200 - 1  == df_dw_bed.loc[ns_df.index, "start"])        
                assert all(ns_df.iloc[:, idx] - 1  == df_dw_bed.loc[ns_df.index, "end"])
        
        df_dic = eventcoor.get_all_flank_bed(sss=True, ss_idx=[0, 1])
        assert len(df_dic.keys()) == 2

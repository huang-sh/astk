from astk.utils import rMATSEventCoord 


def test_sss_rMATSEventCoord():
    A3SS = rMATSEventCoord("data/A3SS.MATS.JC.txt")
    A3SS_ss_dic = A3SS.get_all_flank_bed(sss=True)
    assert list(A3SS_ss_dic.keys()) == ["A0_3SS", "A1_5SS", "A2_3SS", "A3_3SS", "A4_5SS"]
    for key, df in A3SS_ss_dic.items():
        if key.endswith("3SS"):
            assert set(df["end"] - df["start"]) == {23}
        if key.endswith("5SS"):
            assert set(df["end"] - df["start"]) == {9}

    A5SS = rMATSEventCoord("data/A5SS.MATS.JCEC.txt")
    A5SS_ss_dic = A5SS.get_all_flank_bed(sss=True)
    assert list(A5SS_ss_dic.keys()) == ["A0_3SS", "A1_5SS", "A2_5SS", "A3_3SS", "A4_5SS"]
    for key, df in A5SS_ss_dic.items():
        if key.endswith("3SS"):
            assert set(df["end"] - df["start"]) == {23}
        if key.endswith("5SS"):
            assert set(df["end"] - df["start"]) == {9}

    SE = rMATSEventCoord("data/SE.MATS.JC.txt")
    SE_ss_dic = SE.get_all_flank_bed(sss=True)
    assert list(SE_ss_dic.keys()) == ["A0_3SS", "A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS", "A5_5SS"]
    for key, df in SE_ss_dic.items():
        if key.endswith("3SS"):
            assert set(df["end"] - df["start"]) == {23}
        if key.endswith("5SS"):
            assert set(df["end"] - df["start"]) == {9}

    RI = rMATSEventCoord("data/RI.MATS.JCEC.txt")
    RI_ss_dic = RI.get_all_flank_bed(sss=True)
    assert list(RI_ss_dic.keys()) == ["A0_3SS", "A1_5SS", "A2_3SS", "A3_5SS"]
    for key, df in RI_ss_dic.items():
        if key.endswith("3SS"):
            assert set(df["end"] - df["start"]) == {23}
        if key.endswith("5SS"):
            assert set(df["end"] - df["start"]) == {9}

    MXE = rMATSEventCoord("data/MXE.MATS.JC.txt")
    MXE_ss_dic = MXE.get_all_flank_bed(sss=True)
    assert list(MXE_ss_dic.keys()) == ["A0_3SS", "A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS", "A5_5SS", "A6_3SS", "A7_5SS"]
    for key, df in MXE_ss_dic.items():
        if key.endswith("3SS"):
            assert set(df["end"] - df["start"]) == {23}
        if key.endswith("5SS"):
            assert set(df["end"] - df["start"]) == {9}

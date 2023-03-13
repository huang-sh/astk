
from pathlib import Path
import pandas as pd
from astk.utils import get_ss_bed


def test_SUPPA2_get_ss_bed():
    file = "data/fb_e11_12_SE.sig.dpsi"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2], app="SUPPA2")
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_3SS", "A3_5SS"]

    coor_dic1 = get_ss_bed(file, ss_idx=[], app="SUPPA2")
    assert len(coor_dic1) == 4
    assert list(coor_dic1.keys()) == ["A1_5SS", "A2_3SS", "A3_5SS", "A4_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[1, 2, 3])
    assert len(coor_dic) == 3
    assert list(coor_dic.keys()) == ["A2_3SS", "A3_5SS", "A4_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[0, 3])
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A1_5SS", "A4_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[1, 2, 3], altidx=True)
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A2_3SS", "A3_5SS"]

    file = "data/fb_e11_12_A3.dpsi"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2])
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_3SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[0, 1, 2])
    assert len(coor_dic) == 3
    assert list(coor_dic.keys()) == ["A1_5SS", "A2_3SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, altidx=True)
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A2_3SS", "A3_3SS"]

    file = "data/fb_e11_12_A5.sig.dpsi"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2])
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_5SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[0, 2])
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A1_5SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, altidx=True)
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A1_5SS", "A2_5SS"]

    file = "data/fb_e11_12_AF.sig.dpsi"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2])
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_5SS", "A3_pse_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[0, 2, 4])
    assert len(coor_dic) == 3
    assert list(coor_dic.keys()) == ["A1_pse_3SS", "A3_pse_3SS", "A5_3SS"]

    coor_dic = get_ss_bed(file, altidx=True, app="SUPPA2")
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A2_5SS", "A4_5SS"]


def test_rMATS_get_ss_bed():
    file = "data/A3SS.MATS.JC.txt"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2], app="rMATS")
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_5SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[1, 2, 3], altidx=True)
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A3_3SS", "A4_3SS"]
 
    file = "data/MXE.MATS.JC.txt"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 6])
    assert len(coor_dic1) == 2
    assert list(coor_dic1.keys()) == ["A2_5SS", "A7_3SS"]

    coor_dic = get_ss_bed(file, ss_idx=[0, 1, 2])
    assert len(coor_dic) == 3
    assert list(coor_dic.keys()) == ["A1_3SS", "A2_5SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, altidx=True)
    assert len(coor_dic) == 4
    assert list(coor_dic.keys()) == ["A3_3SS", "A4_5SS", "A5_3SS", "A6_5SS"]

    file = "data/RI.MATS.JCEC.txt"
    coor_dic1 = get_ss_bed(file, ss_idx=[0, 1, 2])
    assert len(coor_dic1) == 3
    assert list(coor_dic1.keys()) == ["A1_3SS", "A2_5SS", "A3_3SS"]

    coor_dic = get_ss_bed(file, altidx=True)
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A2_5SS", "A3_3SS"]

    file = "data/SE.MATS.JC.txt"
    coor_dic1 = get_ss_bed(file, ss_idx=[1, 2, 3])
    assert len(coor_dic1) == 3
    assert list(coor_dic1.keys()) == ["A2_5SS", "A3_3SS", "A4_5SS"]

    coor_dic = get_ss_bed(file, altidx=True, app="rMATS")
    assert len(coor_dic) == 2
    assert list(coor_dic.keys()) == ["A3_3SS", "A4_5SS"]

    coor_dic1 = get_ss_bed(file, ss_idx=())
    assert len(coor_dic1) == 6


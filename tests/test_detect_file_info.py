from astk.utils import detect_file_info




def test_detect_file_info(shared_datadir):
    A3SS_test = detect_file_info(f"{shared_datadir}/A3SS.MATS.JC.txt")
    A5SS_test = detect_file_info(f"{shared_datadir}/A5SS.MATS.JCEC.txt")
    MXE_test = detect_file_info(f"{shared_datadir}/MXE.MATS.JC.txt")
    RI_test = detect_file_info(f"{shared_datadir}/RI.MATS.JCEC.txt")
    SE_test = detect_file_info(f"{shared_datadir}/SE.MATS.JC.txt")
    assert A3SS_test == {"app": "rMATS", "etype": "A3SS"}
    assert A5SS_test == {"app": "rMATS", "etype": "A5SS"}
    assert MXE_test == {"app": "rMATS", "etype": "MXE"}
    assert RI_test == {"app": "rMATS", "etype": "RI"}
    assert SE_test == {"app": "rMATS", "etype": "SE"}

    A3_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_A3.dpsi")
    A5_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_A5.sig.dpsi")
    MX_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_MX.sig.dpsi")
    RI_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_RI.sig.dpsi")
    SE_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_SE.sig.dpsi")
    AF_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_AF.sig.dpsi")
    AL_suapp2 = detect_file_info(f"{shared_datadir}/fb_e11_12_AL.sig.dpsi")
    assert A3_suapp2 == {"app": "SUPPA2", "etype": "A3"}
    assert A5_suapp2 == {"app": "SUPPA2", "etype": "A5"}
    assert MX_suapp2 == {"app": "SUPPA2", "etype": "MX"}
    assert RI_suapp2 == {"app": "SUPPA2", "etype": "RI"}
    assert SE_suapp2 == {"app": "SUPPA2", "etype": "SE"}
    assert AF_suapp2 == {"app": "SUPPA2", "etype": "AF"}
    assert AL_suapp2 == {"app": "SUPPA2", "etype": "AL"}


    neg_ep_csv = detect_file_info(f"{shared_datadir}/se.event.sig.neg.csv")
    neg_ep_tsv = detect_file_info(f"{shared_datadir}/se.event.sig.neg.tsv")

    assert neg_ep_csv == {"app": "EventPointer", "etype": "Cassette Exon"}
    assert neg_ep_tsv == {"app": "EventPointer", "etype": "Cassette Exon"}

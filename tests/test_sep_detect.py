from astk.utils import sniff_file_sep


def test_sniff_file_sep(shared_datadir):
    A3SS_test = sniff_file_sep(f"{shared_datadir}/A3SS.MATS.JC.txt")
    A5SS_suppa2 = sniff_file_sep(f"{shared_datadir}/fb_e11_12_A5.sig.dpsi")
    SE_ep = sniff_file_sep(f"{shared_datadir}/se.event.sig.neg.csv")
    SE_ep_tsv = sniff_file_sep(f"{shared_datadir}/se.event.sig.neg.tsv")
    SE_ep_1tsv = sniff_file_sep(f"{shared_datadir}/se.event.sig.pos1.tsv")
    assert A3SS_test == "\t"
    assert A5SS_suppa2 == "\t"
    assert SE_ep == ","
    assert SE_ep_tsv == "\t"
    assert SE_ep_1tsv == " "

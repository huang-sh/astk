from astk.seqfeature._splice_score import ss5_score, ss3_score


def test_ss5_score():
    assert round(ss5_score("cagGTAAGT"), 2) == 10.86
    assert round(ss5_score("gagGTAAGT"), 2) == 11.08
    assert round(ss5_score("taaATAAGT"), 2) == -0.12


def test_ss3_score():
    assert round(ss3_score("ttccaaacgaacttttgtAGgga"), 2) == 2.89
    assert round(ss3_score("tgtctttttctgtgtggcAGtgg"), 2) == 8.19
    assert round(ss3_score("ttctctcttcagacttatAGcaa"), 2) == -0.08

from astk.suppa import read_tpm




def test_read_tpm(shared_datadir):
    file = shared_datadir / "quant.sf"
    salmon_df = read_tpm(file, "salmon", 4)
    assert salmon_df.columns == ["salmon"]


from pandas import read_csv
import astk.event as et


def test_rmats_psi_filter():
    rmats_file = "data/test_psi_filter_rmats.txt"

    dpsi_df = read_csv(rmats_file, sep="\t", index_col=0).dropna()

    def _rmats_mean_psi(vstr):
        strings = vstr.split(",")
        values = [float(i) for i in strings if i != "NA"]
        return sum(values) / len(values)
    
    dpsi_df["PSI"] = dpsi_df["IncLevel1"].apply(_rmats_mean_psi)

    rmats_df_min = et.psi_filter(rmats_file, "/tmp/tpm11111", minv=0.2)
    df_fil_min = dpsi_df.loc[dpsi_df["PSI"]>=0.2, dpsi_df.columns[:-1]]
    assert rmats_df_min.equals(df_fil_min)

    rmats_df_max = et.psi_filter(rmats_file, "/tmp/tpm11111", maxv=0.8)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]<=0.8, dpsi_df.columns[:-1]]
    assert rmats_df_max.equals(df_fil_max)

    minqv = dpsi_df["PSI"].quantile(0.2)
    rmats_df_max = et.psi_filter(rmats_file, "/tmp/tpm11111", minq = 0.2)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]>=minqv, dpsi_df.columns[:-1]]
    assert rmats_df_max.equals(df_fil_max)

    maxqv = dpsi_df["PSI"].quantile(0.8)
    rmats_df_max = et.psi_filter(rmats_file, "/tmp/tpm11111", maxq=0.8)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]<=maxqv, dpsi_df.columns[:-1]]
    assert rmats_df_max.equals(df_fil_max)


def test_suppa2_psi_filter():
    suppa_file = "data/test_psi_filter_suppa2.psi"
        
    dpsi_df = read_csv(suppa_file, sep="\t", index_col=0).dropna()

    dpsi_df["PSI"] = dpsi_df.apply(lambda row: sum(row)/len(row), axis=1)

    suppa_df_min = et.psi_filter(suppa_file, "/tmp/tpm11111", minv=0.2)
    df_fil_min = dpsi_df.loc[dpsi_df["PSI"]>=0.2, dpsi_df.columns[:-1]]
    assert suppa_df_min.equals(df_fil_min)

    suppa_df_max = et.psi_filter(suppa_file, "/tmp/tpm11111", maxv=0.8)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]<=0.8, dpsi_df.columns[:-1]]
    assert suppa_df_max.equals(df_fil_max)

    minqv = dpsi_df["PSI"].quantile(0.2)
    suppa_df_max = et.psi_filter(suppa_file, "/tmp/tpm11111", minq = 0.2)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]>=minqv, dpsi_df.columns[:-1]]
    assert suppa_df_max.equals(df_fil_max)

    maxqv = dpsi_df["PSI"].quantile(0.8)
    suppa_df_max = et.psi_filter(suppa_file, "/tmp/tpm11111", maxq=0.8)
    df_fil_max = dpsi_df.loc[dpsi_df["PSI"]<=maxqv, dpsi_df.columns[:-1]]
    assert suppa_df_max.equals(df_fil_max)

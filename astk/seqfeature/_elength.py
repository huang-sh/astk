from pathlib import Path

from astk.utils._coord import  get_ss_range


def get_elen(file, outdir, app):
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    Path(outdir).mkdir(exist_ok=True)

    df_len = get_ss_range(file, app)    
    df_len.to_csv(outdir / "element_len.csv")

    fig, axes = plt.subplots(1, df_len.shape[1]) 
    for idx in range(df_len.shape[1]):
       axes[idx].boxplot(df_len.iloc[:, idx], showfliers=False, sym="")
    plt.tight_layout() 
    plt.savefig(outdir / "element_len.png")

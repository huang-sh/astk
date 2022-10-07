from pathlib import Path

from astk.utils import sniff_fig_fmt


def cmp_sss(files, output, test, **kwargs):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statannotations.Annotator import Annotator

    df_ls = [pd.read_csv(file, index_col=0) for file in files]
    gns = kwargs.get("groupnames")
    for i, df in enumerate(df_ls):
        df["condition"] = gns[i]

    df = pd.concat(df_ls)    
    dft = df.melt(id_vars=["condition"], var_name="splice_site", value_name="score")
    sns.set_theme()
    pairs = [[(i, gn) for gn in gns] for i in df.columns[:-1]]
    fig_type = kwargs.get("figtype", "point")
    if fig_type == "strip":
        kwargs = {"split": True}
    else:
        kwargs = {}
    g = sns.catplot(
        data=dft, x="splice_site", y="score", 
        hue="condition", kind=fig_type, legend=False, **kwargs
    )
    annotator = Annotator(g.ax, pairs, data=dft, x="splice_site", y="score", hue="condition")
    annotator.configure(test=test, text_format="star", show_test_name=False)
    annotator.apply_test()
    annotator.annotate()
    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    fig_fmt = kwargs.get("figformat", "auto")
    if fig_fmt == "auto":
        fig_fmt = sniff_fig_fmt(output, fmts=['png', 'pdf', 'tiff', 'jpeg'])

    output = Path(output).with_suffix(f".{fig_fmt}")    
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor="w")  

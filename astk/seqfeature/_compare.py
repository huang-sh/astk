from pathlib import Path
from itertools import combinations

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
        fig_kwargs = {"split": True}
    else:
        fig_kwargs = {}
    if all([kwargs["width"], kwargs["height"]]):
        fig_kwargs["height"] = kwargs["height"]
        fig_kwargs["aspect"] = kwargs["width"] / kwargs["height"]
    g = sns.catplot(
        data=dft, x="splice_site", y="score", 
        hue="condition", kind=fig_type, legend=False, **fig_kwargs
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


def cmp_value(files, output, test, **kwargs):
    import numpy as np
    import pandas as pd    
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statannotations.Annotator import Annotator

    df_ls = [pd.read_csv(file, index_col=0) for file in files]
    gns = kwargs.get("groupnames")
    origin_col = df_ls[0].columns
    for i, df in enumerate(df_ls):
        if kwargs.get("log", False):
            df.iloc[:, :df.shape[1]] = np.log2(df)
        xlabel = kwargs["xlabel"]
        if (xlabel is not None) and (len(xlabel) == df.shape[1]):
            if len(xlabel) == len(set(xlabel)):
                columns = [(i, j) for i,j in zip(xlabel, xlabel)]
            else:
                columns = [(i, j) for i,j in zip(origin_col, xlabel)]
        else:            
            columns = [(col, "_".join(col.split("_")[:2])) for col in origin_col]
        columns = pd.MultiIndex.from_tuples(columns)
        df.columns = columns
        df["condition"] = gns[i]
    df = pd.concat(df_ls)
    ss_cols = df.columns[:-1]
    xn, yn = kwargs["xtitle"], kwargs["ytitle"]    
    dft = df.melt(id_vars=["condition"]) # var_name=xn, value_name=yn
    if kwargs.pop("merge_ss"):
        ss_df = dft["variable_0"].str.split("_", n=1, expand=True)
        dft["variable_0"] = ss_df[1]
        dft["variable_1"] = ss_df[1]
        ss_cols = [["5SS", "5SS"], ["3SS","3SS"]]
    dft.rename(columns={"variable_0": xn, "variable_1": "item", "value": yn}, inplace=True)
    sns.set_theme()
    pairs = []
    for col in ss_cols:
        for gn_c in combinations(gns, 2):
            pairs.append(((col[0], gn_c[0]), (col[0], gn_c[1])))
    fig_kwargs = {}
    if (fig_type := kwargs["figtype"]) == "strip":
        fig_kwargs["split"] = True
    if kwargs.get("facet", False):
        fig_kwargs["sharex"] = False
        fig_kwargs["sharey"] = False
        fig_kwargs["col"] = "item"
    if all([kwargs["width"], kwargs["height"]]):
        fig_kwargs["height"] = kwargs["height"]
        fig_kwargs["aspect"] = kwargs["width"] / kwargs["height"]
    g = sns.catplot(
        data=dft, kind=fig_type, x=xn, y=yn, 
        hue="condition", legend=False, **fig_kwargs)
    test_kwargs = {
        "comparisons_correction": kwargs["multicorrect"], 
        "show_test_name": False,
        "text_format": kwargs["pvaltext"]
    }        
    if kwargs.get("facet", False):
        l_col_dic = {}
        for s_col, l_col in ss_cols:
            l_col_dic.setdefault(l_col, [])
            l_col_dic[l_col].append(s_col)
        s_col_dic = dict(ss_cols)
        l_col_ls = [lcol for lcol in l_col_dic.keys()]
        for i, pair in enumerate(pairs):            
            col_g1 = pair[0][0]
            ax = g.axes[0][l_col_ls.index(s_col_dic[col_g1])]
            annotator = Annotator(ax, [pair],
                    data=dft.loc[dft["item"]==s_col_dic[col_g1],],  x=xn, y=yn, 
                    hue="condition",order=l_col_dic[s_col_dic[col_g1]])
            annotator.configure(test=test, **test_kwargs)
            annotator.apply_test()
            annotator.annotate()
    else:
        annotator = Annotator(
                g.ax, pairs,  data=dft, 
                x=xn, y=yn, hue="condition")
        annotator.configure(test=test, **test_kwargs)
        annotator.apply_test()
        annotator.annotate() 
    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    fig_fmt = kwargs.get("figformat", "auto")
    if fig_fmt == "auto":
        fig_fmt = sniff_fig_fmt(output, fmts=['png', 'pdf', 'tiff', 'jpeg'])
    output = Path(output).with_suffix(f".{fig_fmt}")
    plt.xticks(rotation=kwargs["xrotation"])
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor="w")

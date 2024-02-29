import subprocess
from pathlib import Path

from ..constant import *
from ..utils import func as ul
from ..ctypes import *
from ..lazy_loader import LazyLoader

np = LazyLoader("np", globals(), "numpy")
pd = LazyLoader("pd", globals(), "pandas")
sns = LazyLoader("sns", globals(), "seaborn")
plt = LazyLoader("plt", globals(), "matplotlib.pyplot")


def gseplot(termid, output, rdata, fmt, width, height, resolution):
    
    rscript = BASE_DIR / "R" / "gseaPlot.R"
    param_dic = {
        "termid": termid,
        "fmt": fmt, 
        "RData": rdata,
        "output": output,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def upset(files, output, xlabel, dg, fmt, width, height, resolution):

    rscript = BASE_DIR / "R" / "upset.R"

    param_dic = {
        "file": files,
        "dg": dg,
        "fmt": fmt,
        "width": width,
        "height": height,
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}"),
        "name": xlabel or [str(i) for i in range(len(files))],
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def volcano(file, output, adpsi, pvalue, figfmt, width, height, *args, **kwargs):
    df = pd.read_csv(file, sep=kwargs["sep"], index_col=0).dropna()
    dpsi_col, pval_col = kwargs["dpsi_col"], kwargs["pval_col"]
    if {dpsi_col, pval_col}  == {"1", "2"}:
        df.columns = ["FC", "pval"]
    else:
        df_cols = list(df.columns)        
        df_cols[df_cols.index(dpsi_col)] = "FC"
        df_cols[df_cols.index(pval_col)] = "pval"
        df.columns = df_cols

    non_sig = df['pval'] > pvalue
    pos_sig = (df['FC'] > 0) & (df['pval'] < pvalue)
    neg_sig = (df['FC'] < 0) & (df['pval'] < pvalue)

    fig, ax = plt.subplots(figsize=(width, height))
    plt.scatter(df.loc[non_sig, 'FC'], -np.log10(df.loc[non_sig, 'pval']), s=5, color='#bebebe', label="Stable")
    plt.scatter(df.loc[pos_sig, 'FC'], -np.log10(df.loc[pos_sig, 'pval']), s=5, color='#f8766d', label="Up")
    plt.scatter(df.loc[neg_sig, 'FC'], -np.log10(df.loc[neg_sig, 'pval']), s=5, color='#619cff', label="Down") 
    plt.xlim(-1, 1)

    line_kwargs = {"color": 'black', "linestyle": '--', "linewidth": 1, "alpha": 0.5}
    plt.axvline(x=adpsi, **line_kwargs)
    plt.axvline(x=-adpsi, **line_kwargs)
    plt.axhline(y=-np.log10(pvalue), **line_kwargs)
    plt.legend(loc="lower right")

    ## plot scatter dot label
    pos_sig_df = df.loc[pos_sig, ].sort_values(by=["FC", "pval"], ascending=[False, True])
    neg_sig_df = df.loc[neg_sig, ].sort_values(by=["FC", "pval"], ascending=[True, True])
    for i in range(kwargs["top_label"]):
        pid = pos_sig_df.iloc[i, ]
        nid = neg_sig_df.iloc[i, ]
        plabel, nlabel = pid.name.split(";")[0], nid.name.split(";")[0]
        plt.annotate(plabel, (pid["FC"], -np.log10(pid["pval"])), xytext=(-60, 3), textcoords='offset points')
        plt.annotate(nlabel, (nid["FC"], -np.log10(nid["pval"])), xytext=(-60, 3), textcoords='offset points')

    # Add axis labels and title
    plt.xlabel('dPSI')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.savefig(output)


def heatmap(files, output, figfmt, colormap, width, height):
    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    df_ls = [pd.read_csv(file, sep="\t", index_col=0) for file in files]
    dfm = pd.concat(df_ls, axis=1, join="inner")
    sns.set(rc={'figure.facecolor':'white'})
    figsize = (width, height)
    cbar_pos = (1, .35, 0.05, 0.3)
    fig = sns.clustermap(dfm, cmap=colormap, cbar_pos=cbar_pos,
                         yticklabels=False, col_cluster=False, figsize=figsize)
    out = Path(output).with_suffix(f".{figfmt}")
    fig.savefig(out, bbox_inches="tight")
    

def barplot(output, files, xlabel, dg, figfmt, width, height, *args, **kwargs):
    df_ls = []
    for idx, file in enumerate(files):
        file = Path(file)
        df = pd.read_csv(file, sep=kwargs["sep"], index_col=0).dropna()
        dpsi_col, pval_col = kwargs["dpsi_col"], kwargs["pval_col"]
        if {dpsi_col, pval_col}  == {"1", "2"}:
            df.columns = ["dPSI", "pval"]
        else:
            df_cols = list(df.columns)        
            df_cols[df_cols.index(dpsi_col)] = "dPSI"
            df_cols[df_cols.index(pval_col)] = "pval"
            df.columns = df_cols

        if xlabel is None:
            if kwargs["fn_split"] is not None:
                split_param = kwargs["fn_split"]
                ssym, sidx = split_param[0], split_param[1:]
                fn = file.stem
                df["name"] = "_".join([fn.split(ssym)[int(i)] for i in sidx])
            else:
                df["name"] = file.stem
        else:
            df["name"] = xlabel[idx]
        df["type"] = "u"
        df.loc[df["dPSI"] < 0, "type"] = "-"
        df.loc[df["dPSI"] > 0, "type"] = "+"
        if df.shape[0] > kwargs["count_cutoff"]:
            df_ls.append(df)
    dfm = pd.concat(df_ls, axis=0)    
    fig, ax = plt.subplots(figsize=(width, height))
    palette = sns.color_palette("Set2")
    if dg:
        sns.countplot(data=dfm, x="name", hue="type", hue_order=["+", "-"], palette=palette)
    else:
        sns.countplot(x=dfm["name"], palette=palette)
    plt.xticks(rotation=kwargs["x_rotation"], fontsize=kwargs["x_fontsize"])
    plt.tight_layout()
    plt.savefig(output)


def plot_signal_heatmap(
    files: Sequence[FilePath],
    output: FilePath, 
    stype: str,   
    label: Sequence[str],
    width: int,
    height: int,
    colormap: str,
    figfmt: str,
    **kwargs
):
    df_ls = [pd.read_csv(file, index_col=0) for file in files]
    if len({df.shape[1] for df in df_ls}) > 1:
        print("Feature files have different columns!")
        exit()
    if label is None:
        label = [Path(file).stem for file in files]
    cols = {col[:2] for col in df_ls[0].columns}
    fig, axs = plt.subplots(
            len(df_ls)+1, len(cols), 
            figsize=(width, height),
            constrained_layout=True,
            height_ratios = [1, *[3 for _ in files]])
    nbins = df_ls[0].shape[1] // len(cols)
    zmaxs, zmins = [], []
    for df in df_ls:
        zmins.append(np.percentile(df, 1.0))
        zmaxs.append(np.percentile(df, 98.0))

    for di, df in enumerate(df_ls, 1):
        interpolation = 'bilinear' if df.shape[0] >= 1000 else 'nearest'
        for si in range(len(cols)):
            sdf = df.iloc[:, si*nbins:(si+1)*nbins].copy()
            # sort rows by row mean value
            sdf["mean"] = sdf.mean(axis=1)
            sdf = sdf.sort_values('mean',ascending=False)
            del sdf["mean"]
            # plot heatmap
            ax = axs[di, si]
            im = ax.imshow(
                    np.clip(sdf, min(zmins), max(zmaxs)), 
                    aspect='auto', 
                    cmap=colormap, 
                    vmax=max(zmaxs), 
                    vmin=min(zmins),
                    interpolation=interpolation)
            ax.set_yticks([])
            ax.set_xticks([sdf.shape[1]//2])
            ss_label = "_".join(sdf.columns[0].split("_")[1:-1])
            ax.set_xticklabels([ss_label], fontdict={"fontsize": 8})
            ## plot summary profile line
            mvalue = sdf.agg(stype)
            axs[0, si].plot(range(len(mvalue)), mvalue, alpha=0.9, label=label[di-1])
            axs[0, si].set_xticks([])
            if si != 0:
                axs[0, si].set_yticks([])
    axs[0, len(cols)-1].legend(fontsize=6, frameon=False)
    fig.colorbar(im, ax=axs[1:, :])
    plt.savefig(output, bbox_inches='tight', pad_inches=0.1, format=figfmt)
    plt.close()


def plot_cor_heatmap(**kwargs):
    if (glabels := kwargs["grouplabel"]) is not None:
        df_ls = []
        for idx, file in enumerate(kwargs["files"]):
            df = pd.read_csv(file, index_col=0)
            df.columns = [f"{glabels[idx]}_{col}" for col in df.columns]
            df_ls.append(df)        
    else:
        df_ls = [pd.read_csv(file, index_col=0) for file in kwargs["files"]]
    dfm = pd.concat(df_ls, axis=1)
    cor = dfm.corr(method=kwargs["method"])
    fig, ax = plt.subplots(figsize=(kwargs["width"], kwargs["height"]))
    plt.title(kwargs["title"])
    sns.heatmap(cor, mask=np.zeros_like(cor, dtype=np.bool), cmap=kwargs["colormap"],
                square=True, ax=ax)
    plt.savefig(kwargs["output"], bbox_inches="tight")

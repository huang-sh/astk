import subprocess
from pathlib import Path

from astk.constant import *
import astk.utils.func as ul
from astk.ctypes import *


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


def volcano(file, output, fmt, width, height, resolution):
    rscript = BASE_DIR / "R" / "volcano.R"
    param_dic = {
        "file": file,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def pca(files, output, fmt, width, height, resolution, groupname):

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = BASE_DIR / "R" / "pca.R"
    param_dic = {
        "file": files,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}"),
        "groupname": groupname
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def heatmap(files, output, fmt, colormap, width, height):
    import seaborn as sns
    from pandas import read_csv, concat

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    df_ls = [read_csv(file, sep="\t", index_col=0) for file in files]
    dfm = concat(df_ls, axis=1, join="inner")
    sns.set(rc={'figure.facecolor':'white'})
    figsize = (width, height)
    cbar_pos = (1, .35, 0.05, 0.3)
    fig = sns.clustermap(dfm, cmap=colormap, cbar_pos=cbar_pos,
                         yticklabels=False, col_cluster=False, figsize=figsize)
    out = Path(output).with_suffix(f".{fmt}")
    fig.savefig(out, bbox_inches="tight")
    

def barplot(output, files, xlabel, dg, fmt, width, height, resolution):
    rscript = BASE_DIR / "R" / "barplot.R"
    param_dic = {
        "file": files,
        "name": xlabel,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output,
        "dg": dg
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def plot_signal_heatmap(
    files: Sequence[FilePath],
    output: FilePath, 
    stype: str,   
    label: Sequence[str],
    width: int,
    height: int,
    colormap: str,
    fmt: str
):
    import matplotlib.pyplot as plt
    from numpy import percentile, clip
    from pandas import read_csv
    

    df_ls = [read_csv(file, index_col=0) for file in files]

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
        zmins.append(percentile(df, 1.0))
        zmaxs.append(percentile(df, 98.0))

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
                    clip(sdf, min(zmins), max(zmaxs)), 
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
    axs[0, 0].get_shared_y_axes().join(*axs[0, :])
    fig.colorbar(im, ax=axs[1:, :])
    plt.savefig(output, bbox_inches='tight', pad_inches=0.1, format=fmt)
    plt.close()


def plot_cor_heatmap(**kwargs):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

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

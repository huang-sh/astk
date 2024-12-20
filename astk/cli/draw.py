from .config import *
from astk import draw
from astk.utils import sniff_fig_fmt


@cli_fun.command(help="Gene Set Enrichment Analysis ploting")
@click.option('-id', '--id', "termid", cls=MultiOption, type=str, 
                required=True, help="term id")
@click.option('-o', '--output', help="output figure path")
@click.option('-rd', '--RData', help="output figure path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-fw', '--width', default=6, help="fig width, default=6 inches")
@click.option('-fh', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def gseplot(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])
    draw.gseplot(*args, **kwargs)


@cli_fun.command(help="draw UpSet plots for AS events")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")           
@click.option('-o', '--output', required=True, help="output figure path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str,  help="x xlabel")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("AS events can be divided into two groups based on dPSI values \
                   (group +: dPSI > 0, group -: dPSI < 0)"))                 
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-fw', '--width', default=6, help="fig width, default=6 inches")
@click.option('-fh', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def upset(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])    
    draw.upset(*args, **kwargs)


@cli_fun.command(name="volcano", help="Volcano plot analysis for dPSI; short alias: vol")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                required=True, help="input dpsi file")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-adpsi', '--adpsi', default=0.1, show_default=True, help="absolute dpsi cut-off value")
@click.option('-pval', '--pvalue', default=0.05, show_default=True, help="p-value cut-off value")
# @click.option('--dpsi-col', default=1, type=int, show_default=True, help="dpsi column index")
@click.option('--dpsi-col', default="1", type=click.Choice(["1", "IncLevelDifference", "Delta PSI"]), 
                show_default=True, help="dpsi column index, 1 for SUPPA2")
@click.option('--pval-col', default="2", type=click.Choice(["2", "PValue", "FDR", "Pvalue"]), 
                show_default=True, help="p-value column, 2 for SUPPA2")
@click.option('--top-label', default=0, type=int, show_default=True, help="show the top number of AS events' labels")
@click.option('-sep', '--sep', type=click.Choice([",", r"\t"]), default=r"\t", show_default=True,
                help="separator of file content")
@fig_common_options()
def sc_volcano(*args, **kwargs):
    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = sniff_fig_fmt(kwargs["output"])
    if kwargs["dpsi_col"] == "1":
        if kwargs["pval_col"] != "2":
            raise UsageError("You need to set both --dpsi-col 1 and --pval-col 2 at the same time")
    if kwargs["dpsi_col"] == "IncLevelDifference":
        if  kwargs["pval_col"] not in ["PValue", "FDR"]:
            msg = "You need to set both --dpsi-col IncLevelDifference and --pval-col PValue or FDR at the same time"
            raise UsageError(msg)
    if kwargs["dpsi_col"] == "Delta PSI":
        if kwargs["pval_col"] != "Pvalue":
            raise UsageError("You need to set both --dpsi-col 'Delta PSI' and --pval-col Pvalue at the same time")
    kwargs["sep"] = "\t" if kwargs["sep"] in ("\\t", "t") else kwargs["sep"]
    draw.volcano(*args, **kwargs)


@cli_fun.command(name="pca", help="PCA analysis for PSI")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="input psi files")
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-gl', '--groupLabel', cls=MultiOption, type=str, help="file label")
@click.option('-sep', '--sep', type=click.Choice([",", r"\t"]), 
                help="separator of file content")
@click.option('-gb', '--groupBy', type=click.Choice(["row", "col"]), required=True, 
                help=("this option is used to choose whether the sample information is stored in rows or columns;"
                       "AS event PSI use col, feature values use row"))
@fig_common_options()
def sc_pca(*args, **kwargs):
    import pandas as pd
    from .ml import plot_pca

    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = sniff_fig_fmt(kwargs["output"])
    if gn := kwargs.get("grouplabel", []):
        if len(gn) != len(kwargs["files"]):
            raise UsageError("-gl/--groupLabel number must be same as -i/--input")
    else:
        gn = [idx for idx, _ in enumerate(kwargs["files"])]
    sep = "\t" if kwargs["sep"] in ("\\t", "t") else kwargs["sep"]
    axis = 1 if kwargs["groupby"] == "col" else 0
    dfs, labels = [], []
    for idx, file in enumerate(kwargs["files"]):
        df = pd.read_csv(file, sep=sep, index_col=0)
        labels.extend([gn[idx]] * df.shape[axis])
        dfs.append(df)
    dfm = pd.concat(dfs, axis=axis, join='inner').dropna()
    dfm = dfm.T if axis == 1 else dfm
    plot_pca(kwargs["output"], dfm, labels)


@cli_fun.command(help="Heatmap plot for PSI; short alias: hm")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="input psi files")
@click.option('-o', '--output', type=click.Path(), required=True, help="figure output path")
@click.option('-cmap', '--colormap', default="crest", help="matplotlib colormap name")
@fig_common_options()
def heatmap(*args, **kwargs):
    from matplotlib.pyplot import colormaps

    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = sniff_fig_fmt(kwargs["output"])
    if kwargs["colormap"] not in colormaps()+["crest"]:
        msg  = f"'{kwargs['colormap']}' is not a valid value for colormap name"
        msg += f"; supported values are {', '.join(map(repr, colormaps()))}"        
        raise BadParameter(msg)
    draw.heatmap(*args, **kwargs)


@cli_fun.command(help="barplot")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="input psi/dpsi files")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str,
             help="x tick labels for files")
@click.option('-dg', '--dg', is_flag=True, default=False,
              help=("AS events can be divided into two groups based on dPSI values \
                   (group +: dPSI > 0, group -: dPSI < 0)"))
@click.option('-app','--app', required=True, type=click.Choice(["SUPPA2"]), 
                default="SUPPA2", help="the program that generates event file")
@click.option('--dpsi-col', default="1", type=click.Choice(["1", "IncLevelDifference", "Delta PSI"]), 
                show_default=True, help="dpsi column index, 1 for SUPPA2")
@click.option('--pval-col', default="2", type=click.Choice(["2", "PValue", "FDR", "Pvalue"]), 
                show_default=True, help="p-value column, 2 for SUPPA2")
@click.option('-sep', '--sep', type=click.Choice([",", r"\t"]), default=r"\t", show_default=True,
                help="separator of file content")
@click.option('--fn-split', cls=MultiOption, type=str, help="file name split symbol and index")
@click.option('--count-cutoff', default=0, type=int, show_default=True,
                help="Choose values that exceed the cutoff to display")
@fig_common_options()
def barplot(*args, **kwargs):
    nlabel = len(kwargs["files"]) if kwargs["xlabel"] is None else len(kwargs["xlabel"])
    if nlabel != len(kwargs["files"]):
        raise UsageError("-i/--input values number must be same as -xl/--xlabel")
    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = sniff_fig_fmt(kwargs["output"])        
    kwargs["sep"] = "\t" if kwargs["sep"] in ("\\t", "t") else kwargs["sep"]
    draw.barplot(*args, **kwargs)

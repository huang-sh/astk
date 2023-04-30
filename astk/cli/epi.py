from pathlib import Path

import astk.utils as ul
from .config import *
from astk import epi


@cli_fun.command(name="signalExtract", help="Extract bigwig signal of splice sites; short alias: se")
@click.option('-e', '--event', type=click.Path(exists=True), required=True, 
                help="AS event file")
@click.option('-bw', '--bigwig', type=click.Path(exists=True), help="bigwig file")
@click.option('-st', '--signalType', "stype", default="max", show_default=True,
                type=click.Choice(["mean", "min", "max", "coverage", "std"]), 
                help="bigwig bin signal extraction type")
@click.option('-si', '--siteIndex', "sites", cls=MultiOption, type=int, 
                default=[], help="splice site index, 1-index")
@click.option('-alt', '--altIdx', is_flag=True, default=False, show_default=True,
                help="get alternative splicing sites index")
@click.option('-ew', '--exonWidth', "exon_width", type=int, default=150, 
                show_default=True, help="exon flank window width")
@click.option('-iw', '--intronWidth', "intron_width", type=int, default=150, 
                show_default=True,help="intron flank window width")
@click.option('-bs', '--binSize', type=int, default=5, show_default=True,
                help="bin size")
@click.option('-o', '--output', required=True, help="output name")                                 
def sc_extract_signal(*args, **kwargs):
    coor_dic = ul.get_ss_bed(
        kwargs["event"], 
        exon_width=kwargs["exon_width"], 
        intron_width=kwargs["intron_width"],
        ss_idx=[i+1 for i in kwargs["sites"]],
        altidx=kwargs["altidx"]
    )
    nbins = (kwargs["exon_width"] + kwargs["intron_width"]) // kwargs["binsize"]
    epi.extract_signal(
        kwargs["output"], 
        coor_dic, 
        kwargs["bigwig"], 
        kwargs["stype"],
        nbins
    )


@cli_fun.command(name="signalHeatmap", help="Plot signal heatmap; short alias: shm")
@click.option('-s', '--score', "files", cls=MultiOption, type=click.Path(exists=True),  
                required=True, help="AS event file")
@click.option('-o', '--output', required=True, help="output name")
@click.option('-st', '--summaryType', "stype", default="mean", 
                type=click.Choice(["mean","median","min","max","std","sum"]), 
                help="""Define the type of statistic that should be plotted in the summary\n\b
                        image above the heatmap, default=mean""")
@click.option('--label', cls=MultiOption, help="AS event file labels")
@click.option('-fw', '--width', default=6, help="fig width, default=6 inches")
@click.option('-fh', '--height', default=6, help="fig height, default=6 inches")
@click.option('-cmap', '--colormap', default="RdYlBu", help="heatmap color")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf']),
                default="auto", help="output figure format")                                 
def sc_plot_heatmap(*args, **kwargs):
    from matplotlib.pyplot import colormaps
    from astk import draw

    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = ul.sniff_fig_fmt(kwargs["output"])
    if kwargs["colormap"] not in colormaps():
        msg  = f"'{kwargs['colormap']}' is not a valid value for colormap name"
        msg += f"; supported values are {', '.join(map(repr, colormaps()))}"        
        BadParameter(msg)
    draw.plot_signal_heatmap(*args, **kwargs)


@click.command(help = "generate ChromHMM mark file")
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-ct', '--cellType', cls=MultiOption, type=str, required=True, help="cell types")
@click.option('-bed', '--bed', cls=MultiOption, type=click.Path(exists=True), help="bed files")
@click.option('-mn', '--markNum', cls=MultiOption, type=int, required=True,
                help="mark count of every cell types") 
@click.option('-sep', help="split symbol for splitting bed file names") 
@click.option('-mi', "--markIndex", type=int, 
            help="the mark index when bed files are splited by -sep symbol")
@click.option('--markName', cls=MultiOption, type=str, help="mark names")
@click.option('--stacked', is_flag=True, 
            help="This flag replaces the mark entry with an entry of the form cell_mark.")
def mark(*args, **kwargs):
    pass


CHROMSIZES_dir = Path(__file__).parent / "ChromHMM/CHROMSIZES"
genomes = [i.stem for i in CHROMSIZES_dir.glob("*.txt")]

@click.command(["learnState", "lcs"], 
    help = "learns a chromatin state model, Wrapper of ChromHMM BinarizeBed and LearnModel")
@click.option('-n', "--numstates", default=2, type=int, help="states number")
@click.option('-mark', '--markfile', type=click.Path(exists=True), required=True, help="cell mark file")
@click.option('-dir', '--directory', type=click.Path(exists=True), 
                default=".", help="work directory")
@click.option('-bd', '--binaryDir', required=True, help="binary output directory")
@click.option('-od', '--outdir', required=True,  help="result output directory")
@click.option('-b', '--binsize', type=int, default=150, help="binsize")
@click.option('-g', '--genome', required=True, type=click.Choice(genomes), help="genome assembly")
@click.option('-mx', default="2400M", help="memory for java running")
@click.option('-p', "--processor", default=2, type=int, help="maxprocessors")
@click.option('-anchor', "--anchorDir", type=click.Path(exists=True), help="anchor files directory")
@click.option('-coor', "--coorDir", type=click.Path(exists=True), help="coordinate files directory")
@click.option('-nb', "--no_binary", is_flag=True, default=False, help="not run BinarizeBed")
@click.option('--name', default="", help="file name")
@click.option('--stacked', is_flag=True, 
            help="This flag replaces the mark entry with an entry of the form cell_mark.")
@click.option('--nostrand', is_flag=True, 
            help="This flag is present then strand information is not used in the enrichment calculations.")
@click.option('-DC', '--defaultCoor', is_flag=True, 
            help="This flag is present then default is used in the enrichment calculations.")            
def LearnState(*args, **kwargs):
    pass
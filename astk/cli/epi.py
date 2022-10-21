from pathlib import Path

from astk.utils import sniff_fig_fmt
from .config import *
from astk import epi

plotType_help_text = """
lines" will plot the profile line based on the average type selected. "fill" fills the region between zero and the profile curve. The fill in color is semi transparent to
distinguish different profiles. "se" and "std" color the region between the profile and the standard error or standard deviation of the data. (default: lines)")
"""

@click.command(help="signal profile withing splicing sites")
@click.option('-o', '--output', required=True, type=click.Path(),
                help="output path")
@click.option('-e', '--eventFile', "event_file", type=click.Path(exists=True),
                help="AS event file that including AS event id")
@click.option('-bw', '--bwFile', "bw_files", cls=MultiOption, type=click.Path(exists=True),
                help="AS event files that including AS event id")
@click.option('-bs', '--binSize', "bin_size", default=5, type=int, 
                help="bin size, default=5")
@click.option('-uw', '--upStreamWidth', "ups_width", default=150, type=int, 
                help="upstream width, default=150")
@click.option('-dw', '--downStreamWidth', "dws_width", default=150, type=int, 
                help="downstream width, default=150")
@click.option('-p', '--process', "threads", default=4, type=int, 
                help="running processors, default=4")
@click.option('--plotType', 'plot_type', type=click.Choice(['lines',"fill","se",'std']), 
                default="lines", help=plotType_help_text)
@click.option('--colorMap', 'color_map', default="RdYlBu", help="Color map to use for the heatmap.")
@click.option('-fmt', '--figureFormat', "plot_format", type=click.Choice(['auto', 'png',"pdf","svg",'plotly']), 
                default="auto", help="Image format type")
@click.option('--height', 'fig_height', default=28, type=int, 
                help="plot height in cm, default=28")
@click.option('--width', 'fig_width', default=4, type=int, 
                help="plot width in cm, default=4")
@click.option('-ssl', '--splicingSiteLabel', "regionsLabel", cls=MultiOption, 
                help="splicing site labels")
@click.option('--samplesLabel', 'samplesLabel', cls=MultiOption, 
                help="bw samples labels")
def signal_profile(*args, **kwargs):
    if kwargs["plot_format"] == "auto":
        kwargs["plot_format"] = sniff_fig_fmt(kwargs["output"], fmts=['png',"pdf","svg",'plotly'])       
    epi.signal_heatmap(*args, **kwargs)


@click.command(help="signal profile comparision")
@click.option('-o', '--output', required=True, type=click.Path(),
                help="output path")
@click.option('-mat', 'mat_file', type=click.Path(exists=True), cls=MultiOption,
                help="matrix file from signalProfile")
@click.option('-n', '--name', default="signal", help="name")
@click.option('-gn', '--groupName', cls=MultiOption, type=str, help="event group name")
@click.option('--height', 'fig_height', default=8, type=int, 
                help="plot height in cm, default=28")
@click.option('--width', 'fig_width', default=8, type=int, 
                help="plot width in cm, default=4")
@click.option('-fmt', '--figureFormat', "fig_format", type=click.Choice(['auto', 'png',"pdf"]), 
                default="auto", help="Image format type")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")                
def signal_profile2(*args, **kwargs):
    if kwargs["fig_format"] == "auto":
        kwargs["fig_format"] = sniff_fig_fmt(kwargs["output"], fmts=['png',"pdf"])       
    epi.signal_metaplot(*args, **kwargs)


@click.command(help="epi feature extract")
@click.option('-e', '--eventFile', "event_file", cls=MultiOption, type=click.Path(exists=True),
                help="AS event files that including AS event id")
@click.option('-el', '--eventLabel', "event_label", cls=MultiOption, type=str, 
                help="AS event files label")
@click.option('-bam', '--bamFile', "bam_file", cls=MultiOption, type=click.Path(exists=True),
                help="AS event files that including AS event id")
@click.option('-bl', '--bamLabel', "bam_label", cls=MultiOption, type=str,
                help="bam files label")
@click.option('-w', '--width',cls=MultiOption, type=str,
                help="window width, default=300")
@click.option('-bs', '--binSize', "bin_size", default=5, type=int, help="bin size, default=5")
@click.option('-o', '--output', required=True, help="output name")
@click.option('-nm', '--normalMethod', default="count", type=click.Choice(['count', 'CPM']))
@click.option('-pe', '--pairedEnd', "paired_end", is_flag=True, default=False)
@click.option('--title', default="AS sites signal Profile", help="plot title")
@click.option('--bamMerge', "bam_merge", is_flag=True, default=False, help="bam feature merge")
@click.option('-et', "--eventType", "as_type", help="bam feature merge")
def epi_sc(*args, **kwargs):
    epi.epi_sc(*args, **kwargs)

@click.command(help="epi signal profile")
@click.argument('file', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('--title', default="AS sites signal Profile", help="plot title")
@click.option('--ylim', type=(float, float), help="y limit")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def epi_profile_sc(*args, **kwargs):
    epi.epi_profile(*args, **kwargs)

@click.command(help="epi signal heatmap")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def epihm(*args, **kwargs):
    epi.epihm(*args, **kwargs)


@click.command(help="epi signal compare")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def sigcmp(*args, **kwargs):
    epi.sigcmp(*args, **kwargs)


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
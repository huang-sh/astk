from email.policy import default
from pathlib import Path

from .config import *
from astk import epi


@click.command(help="epi feature extract")
@click.option('-e', '--eventFile', "event_file", cls=MultiOption, type=tuple, 
                default=(), help="AS event files that including AS event id")
@click.option('-el', '--eventLabel', "event_label", cls=MultiOption, type=tuple, 
                default=(), help="AS event files label")
@click.option('-bam', '--bamFile', "bam_file", cls=MultiOption, type=tuple, default=(),
                help="AS event files that including AS event id")
@click.option('-bl', '--bamLabel', "bam_label", cls=MultiOption, type=tuple, default=(),
                help="bam files label")
@click.option('-w', '--width',cls=MultiOption, type=tuple, default=(),
                help="window width, default=300")
@click.option('-bs', '--binSize', "bin_size", default=5, type=int, help="bin size, default=5")
@click.option('-o', '--output', required=True, help="output name")
@click.option('-nm', '--normalMethod', default="count", type=click.Choice(['count', 'CPM']))
@click.option('-pe', '--pairedEnd', "paired_end", is_flag=True, default=False)
@click.option('--title', default="AS sites signal Profile", help="plot title")
@click.option('--bamMerge', "bam_merge", is_flag=True, default=False, help="bam feature merge")
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
@click.option('-ct', '--cellType', cls=MultiOption, type=tuple, required=True, help="cell types")
@click.option('-bed', '--bed', cls=MultiOption, type=tuple, help="bed files")
@click.option('-mn', '--markNum', cls=MultiOption, type=tuple, required=True,
                 help="mark count of every cell types") 
@click.option('-sep', help="split symbol for splitting bed file names") 
@click.option('-mi', "--markIndex", type=int, 
            help="the mark index when bed files are splited by -sep symbol")
@click.option('--markName', cls=MultiOption, type=tuple, help="mark names")
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
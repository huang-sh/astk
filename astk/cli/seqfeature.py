"""
astk.cli.seqfeature
~~~~~~~~~~~~~~~~~
This module provide sequence feature extraction cli api
"""

from .config import *
import astk.seqfeature.feature as sf
from astk.seqfeature import splice_score, get_elen, get_gcc
from astk.seqfeature import cmp_sss, cmp_value


@click.command(help = "extract DNA sequence feature ")
@click.option('-fa', '--fasta', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="DNA sequence fasta file")
@click.option('-o', '--output', type=click.Path(), 
                help="output path")                
@click.option('-k', '--kmer', type=int, default=1, help="K-tuple or k-mer value, default=1")
@click.option('-g', '--gap', type=int, default=0, help="gap value, default=0")
@click.option('-l', '--lambda', "lambda_", type=int, default=0, help="lambda-correlation value, default=0")
@click.option('--count', is_flag=True, default=False, help="feature count")
@click.option('-FM', "--featureMerge", is_flag=True, default=False, help="feature count")
def sc_extract(*args, **kwargs):
    sf.seq_extract(*args, **kwargs)


@click.command(help = "draw seqLogo")
@click.option('-fa', '--fasta', type=click.Path(exists=True), 
                required=True, help="DNA sequence fasta file")
@click.option('-o', '--output', type=click.Path(), 
                help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")                
def sc_seqlogo(*args, **kwargs):
    sf.seq_pcm(*args, **kwargs)


@click.command(["spliceScore", "ss"], help="Compute 5/3 Splice site strength")
@click.option('-e', "--event", 'file', type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-fi', 'gfasta', type=click.Path(exists=True), help="genome fasta")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
@click.option('-p', '--process', default=4, help="process number, default=4")
def sc_splice_score(*args, **kwargs):
    splice_score(*args, **kwargs)


@click.command(["getlen"], help="Compute element length")
@click.option('-e', "--event", 'file', type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-log', "--log", is_flag=True, default=False, help="log2 transformation")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
def sc_get_elen(*args, **kwargs):
    get_elen(*args, **kwargs)


@click.command(["gcc"], help="Compute GC content")
@click.option('-e', "--event", 'file', type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-fi', 'gfasta', type=click.Path(exists=True), help="genome fasta")
@click.option('-bs', '--binsize', default=0, help="use bin size or slide window to compute exon/intron GC content")
@click.option('-ef', '--exonFlank', default=150, help="the exon flank width")
@click.option('-if', '--intronFlank', default=150, help="the intron flank width")
@click.option('--includeSS', is_flag=True, default=False, help="include splice site flank region")
# @click.option('-ele', '--element', is_flag=True, default=False, help="compute the whole exon/intron GC content")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
@click.option('-p', '--process', default=4, help="process number, default=4")
def sc_get_gcc(*args, **kwargs):
    get_gcc(*args, **kwargs)


@click.command(["ssscmp"], help="Compare Splice site strength between two conditions")
@click.option('-e', "--events", 'files', type=click.Path(exists=True), cls=MultiOption, 
                required=True, help="event files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-test', '--test', default="Mann-Whitney", 
                type=click.Choice(['Mann-Whitney', 't-test_ind', 't-test_welch', "Wilcoxon"]), 
                help=" statistical test method, default='Mann-Whitney'")
@click.option('-gn', '--groupNames', cls=MultiOption, type=str, 
                help="group names, default= g1 g2 ")
@click.option('-ft', '--figType',  default="box", help="figure display type",
                type=click.Choice(["point", 'strip', 'box', 'boxen', 'violin', 'bar']))
@click.option('-ff', '--figFormat', type=click.Choice(['auto', 'png', 'pdf', 'tiff', 'jpeg']), 
                default="auto", help="output figure format")       
def sc_ssscmp(*args, **kwargs):
    fn = len(kwargs["files"])
    if fn != 2:
        raise UsageError("only support two files for -e/--events")
    if gn := kwargs.get("groupnames"):
        if len(gn) != fn:
            raise UsageError("-gn/--groupNames parameter number must be same as -e/--events")
    else:
        kwargs["groupnames"] = [f"g{i}" for i in range(1, fn+1)]
    cmp_sss(*args, **kwargs)


@click.command(["vcmp"], help="Compare sequence feature value between two conditions")
@click.option('-e', "--events", 'files', type=click.Path(exists=True), cls=MultiOption, 
                required=True, help="event files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-test', '--test', default="Mann-Whitney", 
                type=click.Choice(['Mann-Whitney', 't-test_ind', 't-test_welch', "Wilcoxon"]), 
                help=" statistical test method, default='Mann-Whitney'") 
@click.option('-facet', "--facet", is_flag=True, default=False, 
                help="If true, the facets will not x axes across rows.")
@click.option('-log', "--log", is_flag=True, default=False, help="log2 transformation")
@click.option('-gn', '--groupNames', cls=MultiOption, type=str, 
                help="group names, default= g1 g2 ")
@click.option("--merge-ss", is_flag=True, default=False, help="merge 5'/3' splice sites")      
@click.option('-mc', '--multiCorrect', type=click.Choice(['bonf', 'HB', 'holm', 'BH' ,'BY']), 
                help="multiple test correction method")
@click.option('--pvalText', type=click.Choice(['star', 'simple']), default="star",
                help="p-value display format, default='star'")
@click.option('--xtitle', default="item", help="x title, default='item'")
@click.option('--ytitle', default="value", help="y title, default='value'")
@click.option('--xlabel', cls=MultiOption, type= str, help="x labels, default is input file colnames")
@click.option('--xrotation', type=int, default=0, help="x tick labels rotation")
@click.option('-fs', '--figSize', type=(float, float), help="figure size")
@click.option('-ft', '--figType',  default="box", help="figure display type",
                type=click.Choice(["point", 'strip', 'box', 'boxen', 'violin', 'bar']))
@click.option('-ff', '--figFormat', type=click.Choice(['auto', 'png', 'pdf', 'tiff', 'jpeg']),
                default="auto", help="output figure format")
@click.option('-fw', '--width', type=float, help="figure width")
@click.option('-fh', '--height', type=float, help="figure height")  
def sc_cmp_value(*args, **kwargs):
    fn = len(kwargs["files"])
    if fn < 2:
        raise UsageError("-e/--events input must be greater than one file")
    if gn := kwargs.get("groupnames"):
        if len(gn) != fn:
            raise UsageError("-gn/--groupNames parameter number must be same as -e/--events")
    else:
        kwargs["groupnames"] = [f"g{i}" for i in range(1, fn+1)]
    cmp_value(*args, **kwargs)

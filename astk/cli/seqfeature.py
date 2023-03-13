"""
astk.cli.seqfeature
~~~~~~~~~~~~~~~~~
This module provide sequence feature extraction cli api
"""

from .config import *
import astk.seqfeature.feature as sf
import astk.utils as ul
from astk.seqfeature import splice_score, get_elen, get_gcc
from astk.seqfeature import cmp_value


@cli_fun.command(name="seqfeature", help = "extract DNA sequence feature; short alias: seqf")
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


@cli_fun.command(name="seqlogo", help = "draw seqLogo")
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


@cli_fun.command(name="spliceScore", help="Compute 5/3 Splice site strength; short alias: sss")
@click.option('-e', "--event", type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-fi', "--fasta",'gfasta', type=click.Path(exists=True), help="genome fasta")
@click.option('-si', '--siteIndex', "sites", cls=MultiOption, type=int,
                default=[], help="splice site index, 0-index")
@click.option('-alt', '--altIdx', is_flag=True, default=False, show_default=True,
                help="get alternative splicing sites index")
@click.option('-app','--app', type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", show_default=True, help="the program that generates event file")
@click.option('-p', '--process', type=int, default=4, show_default=True, help="process number")
def sc_splice_score(*args, **kwargs):
    coor_dic = ul.get_ss_bed(
        kwargs["event"], 
        sss=True,
        app=kwargs["app"],
        ss_idx=[i+1 for i in kwargs["sites"]],
        altidx=kwargs["altidx"]        
    )
    splice_score(
        coor_dic,
        kwargs["outdir"],
        kwargs["gfasta"],
        kwargs["process"]
    )


@cli_fun.command(name="getlen", help="Compute element length")
@click.option('-e', "--event", 'file', type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-log', "--log", is_flag=True, default=False, help="log2 transformation")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file")
def sc_get_elen(*args, **kwargs):
    get_elen(*args, **kwargs)


@cli_fun.command(name="gcc", help="Compute GC content")
@click.option('-e', "--event", type=click.Path(exists=True), required=True,
                help="event file")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-fi', '--fasta', 'gfasta', type=click.Path(exists=True), help="genome fasta")
@click.option('-si', '--siteIndex', "sites", cls=MultiOption, type=int,
                default=[], help="splice site index, 1-index")
@click.option('-alt', '--altIdx', is_flag=True, default=False, show_default=True,
                help="get alternative splicing sites index")
@click.option('-bs', '--binsize', default=15, 
                help="use bin size or slide window to compute exon/intron GC content, default=15")
@click.option('-ef', '--exonFlank', default=150, help="the exon flank width, default=150")
@click.option('-if', '--intronFlank', default=150, help="the intron flank width, default=150")
@click.option('--includeSS', is_flag=True, default=False, help="include splice site flank region")
# @click.option('-ele', '--element', is_flag=True, default=False, help="compute the whole exon/intron GC content")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]), 
                default="auto", help="the program that generates event file, default='auto'")
@click.option('-p', '--process', default=4, help="process number, default=4")
def sc_get_gcc(*args, **kwargs):  
    coor_dic = ul.get_ss_bed(
        kwargs["event"], 
        split=True,
        app=kwargs["app"],
        excludeSS=(not kwargs["includess"]),
        exon_width=kwargs["exonflank"],
        intron_width=kwargs["intronflank"],
        ss_idx=[i+1 for i in kwargs["sites"]],
        altidx=kwargs["altidx"]          
    )    
    get_gcc(
        coor_dic,
        kwargs["outdir"],
        kwargs["gfasta"],
        kwargs["binsize"]
    )


@cli_fun.command(name="vcmp", help="Compare sequence feature value among multiple conditions")
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
                help="group names, default= g1 g2... ")
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

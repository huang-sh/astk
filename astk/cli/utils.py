"""
astk.cli.utils
~~~~~~~~~~~~~~~~~
This module provide some utility function
"""

from email.policy import default
from .config import *
import astk.utils._cli_func as ul


@click.command(help="generate metadata template file")
@click.option('-c1', '--ctrl', cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 1")
@click.option('-c2', '--case', cls=MultiOption, type=click.Path(exists=True),
                required=True, help="file path for condtion 2")
@click.option('-gn', '--groupName', cls=MultiOption, type=str,
                default=(), help="group name")                
@click.option('-repN', '--replicate', cls=MultiOption, type=int,
                help="replicate, number")
@click.option('-o', '--output', required=True, type=click.Path(),
                help='metadata output path')                
@click.option('-repN1', '--replicate1', cls=MultiOption, type=int,
                help="replicate1, number")
@click.option('-repN2', '--replicate2', cls=MultiOption, type=int, 
                help="replicate2, number")
@click.option('--condition', cls=MultiOption, type=str, help="condition name")
@click.option('-fn', '--filename', 'is_fn', is_flag=True, help="file name")
@click.option('--split', cls=MultiOption, type=str, help="name split symbol and index")    
def meta(*args, **kwargs):
    ul.meta(*args, **kwargs)


@click.command(help = "install R packages")
@click.option('-r', '--requirement', is_flag=True, default=False,
                help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", cls=MultiOption, type=str, 
                default=(), help="install Genome wide annotation package")
@click.option('-cran', '--cran', cls=MultiOption, type=str, 
                default=(), help="install CRAN package")
@click.option('-bioc', '--bioconductor', cls=MultiOption, type=str, 
                default=(),help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
@click.option('-m', '--mirror',  is_flag=True, default=False,
                help="use tsinghua mirrors.")
@click.option('--conda',  is_flag=True, default=False,
                help="install packages via conda")
def install(*args, **kwargs):
    ul.install(*args, **kwargs)


@click.command(help = "get motif from meme file")
@click.argument('motifId', nargs=-1, required=True)   
@click.option('-mm', '--meme', type=click.Path(exists=True), help="meme motif file")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                help="RBP motif database")
@click.option('-org', '--organism', help="RBP organism")
@click.option('-o', '--output', required=True, help="output path")
def getmeme(*args, **kwargs):

    ul.getmeme(*args, **kwargs)

@click.command(help = "list OrgDb")
@click.option('-orgdb', '--OrgDb', "OrgDb", is_flag=True, help="list OrgDb")
@click.option('-rbpsp', '--RBPSp', "RBPSp", is_flag=True, help="RNA binding protein  ")
def list_(*args, **kwargs):
    ul.list_(*args, **kwargs)


@click.command(help = "generate ChromHMM anchor file")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-idx', '--index', type=int, help="element index")
@click.option('-si', '--sideIndex', type=(int, int), help="the center of two side index")
@click.option('-u', '--upstreamOffset', "offset5", type=int, default=0, help="upstream offset")
@click.option('-d', '--downstreamOffset', "offset3", type=int, default=0, help="downstream offset")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
def anchor(*args, **kwargs):
    ul.anchor(*args, **kwargs)
 

@click.command(help = "generate bed file according to selected coordinates")
@click.option('-i', '--input', 'file', type=click.Path(exists=True),
                required=True,  help='AS event file')
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-a', '--anchor', cls=MultiOption, type=int, 
                help="splice site index. if not set, it will use all splice sites. 1-based")
@click.option('-uw', '--upstreamWidth', "upstream_w", type=int, default=50,
                help="flank width of splice site upstream")
@click.option('-dw', '--downstreamWidth', "downstream_w", type=int, default=50,
                help="flank width of splice site downstream")
@click.option('--interval', type=(int, int), default=(None, None),
                help="splice site start and end index, 1-based")
@click.option('--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
@click.option('-fi', 'fasta', type=click.Path(exists=True), 
                help="Input FASTA file. if set, the fasta sequence will be extracted")
def getcoor(*args, **kwargs):
    ul.getcoor(*args, **kwargs)


@click.command(help = "Make the TxDb object")
@click.argument('gtf', type=click.Path(exists=True), required=True)
@click.option('-org', '--organism', required=True, help="organism")
@click.option('-o', '--output', help="file output path")
def mktxdb(*args, **kwargs):
    ul.mkTxDb(*args, **kwargs)


@click.command(help = "get gene ID from AS event file")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), help="file output path")
@click.option('-u', '--unique', is_flag=True, help="only save unique gene ID")
def getgene(*args, **kwargs):
    ul.getgene(*args, **kwargs)


@click.command(["merge"], help="merge file")
@click.option('-i','--input', 'files' ,cls=MultiOption, type=click.Path(exists=True),
                required=True , help="input files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-axis', '--axis', type=click.IntRange(min=0, max=1), default=0,
                help="merge direction, 0 for row merge and 1 for column merge")
@click.option('-rmdup', '--rmdup', type=click.Choice(["index", 'all', 'content']), 
                help="remove duplicate rows")
@click.option('-rmna', '--rmna', is_flag=True, help="remove NA data")                
def sub_merge_files(*args, **kwargs):
    from astk.utils.func import merge_files
    merge_files(*args, **kwargs)

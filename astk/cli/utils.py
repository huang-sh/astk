from .config import *
import astk.utils as ul


@click.command(help="generate metadata template file")
@click.option('-p1', '--control', cls=MultiOption,type=tuple, default=(),
                required=True, help="file path for condtion 1")
@click.option('-p2', '--treatment', cls=MultiOption,type=tuple, default=(),
                required=True, help="file path for condtion 2")
@click.option('-gn', '--groupName', cls=MultiOption, type=tuple,
                default=(), help="group name")                
@click.option('-repN', '--replicate', cls=MultiOption,type=tuple, 
                help="replicate, number")
@click.option('-o', '--output', required=True, type=click.Path(),
                help='metadata output path')                
@click.option('-repN1', '--replicate1', cls=MultiOption,type=tuple, 
                help="replicate1, number")
@click.option('-repN2', '--replicate2', cls=MultiOption,type=tuple, 
                help="replicate2, number")
@click.option('-fn', '--filename', 'is_fn', is_flag=True, help="file name")
@click.option('--split', cls=MultiOption, type=tuple, help="name split symbol and index")    
def meta(*args, **kwargs):
    ul.meta(*args, **kwargs)



@click.command(help = "install R packages")
@click.option('-r', '--requirement', is_flag=True, default=False,
                help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", cls=MultiOption, type=tuple, 
                default=(), help="install Genome wide annotation package")
@click.option('-cran', '--cran', cls=MultiOption, type=tuple, 
                default=(), help="install CRAN package")
@click.option('-bioc', '--bioconductor', cls=MultiOption, type=tuple, 
                default=(),help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
@click.option('-m', '--mirror',  is_flag=True, default=False,
                help="use tsinghua mirrors.")
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
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-s', '--start', type=int, help="start index")
@click.option('-e', '--end', type=int, help="end index")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
@click.option('-a', '--anchor', type=int, help="element index")
@click.option('-u', '--upstreamWidth', "upstream_w", type=int, default=150, help="width of right flank")
@click.option('-d', '--downstreamWidth', "downstream_w", type=int, default=150, help="width of left flank")
@click.option('-fi', 'fasta', type=click.Path(exists=True), 
            help="Input FASTA file. if set, the fasta sequence will be extracted")
def getcoor(*args, **kwargs):
    ul.getcoor(*args, **kwargs)
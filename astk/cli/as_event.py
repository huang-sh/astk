from .config import *
from astk.event import _cli_func as et


@click.command(help="length distribution")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True),
                required=True,  help='AS ioe file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-cl', '--custom_len', 'custom_len', cls=MultiOption, type=int, help="custom length")
@click.option('-nc', '--cluster', type=int, default=4, help="number of cluster")
@click.option('-bw', '--width', type=int, default=3, help="bin width")
@click.option('-lw', '--len_weight', type=float, default=2, help="length weight")
@click.option('--max_len', type=int, default=500, help="the max length of exon in clustering")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'jpeg', 'pdf', 'tiff']), 
                default="png", help="output figure format")
def len_dist(*args, **kwargs):
    et.len_dist(*args, **kwargs)


@click.command(help="length cluster")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-id', '--inDir', type=click.Path(exists=True), help='input direcory')
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-lr', '--lenRange', cls=MultiOption, type=int,
                default=(1, ), help="custom length")
def len_cluster(*args, **kwargs):
    et.len_cluster(*args, **kwargs)


@click.command(help="pick AS event with specific length")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True), help='AS ioe file')
@click.option('-o', '--output', help="output path")
@click.option('-rg', '--range', "len_range", type=(int, int), required=True, help="length range")
def len_pick(*args, **kwargs):
    et.len_pick(*args, **kwargs)


# it is deprecated
@click.command(help="filter significant result")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outDir', required=True, help="output directory")
@click.option('-dpsi', '--dpsi', type=float, default=0, help="dpsi threshold value")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value")
@click.option('-adpsi', '--abs_dpsi', type=float, default=0, help="absulte dpsi threshold value")
@click.option('-pf1', '--psiFile1', cls=MultiOption, type=click.Path(exists=True),
                default=(), help="psi files of condtion 1")
@click.option('-pf2', '--psiFile2', cls=MultiOption, type=click.Path(exists=True),
                default=(), help="psi files of condtion 2")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
                default="tsv", help="out  file format ")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]),
                default="auto", help="the program that generates event file")
def _sigfilter(*args, **kwargs):
    et.sigfilter(*args, **kwargs)


@click.command(help="filter psi result")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input psi file")
#@click.option('-md', '--metadata', type=click.Path(exists=True), help="metadata file")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-psi', '--psi', type=float, default=0, help="psi threshold value, defualt=0")
@click.option('-qt', '--quantile', type=float, default=0, help="quantile threshold value, defualt=0")
# @click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
#                 default="tsv", help="out  file format ")
def psi_filter(*args, **kwargs):
    et.psi_filter(*args, **kwargs)

 
@click.command(help="intersect AS event")
@click.option('-a', "file_a", type=click.Path(exists=True), required=True, help="dpsi or psi file A")
@click.option('-b', "file_b", type=click.Path(exists=True), help="dpsi or psi file B")
@click.option('-o', "--output", required=True, help="output suffix")
@click.option('-ioeb', "--ioeB", type=click.Path(exists=True), help="ioe file file B")
@click.option('-ib', "--ignoreB", type=bool, default=False, help="ignore file B intersection output")
def intersect(*args, **kwargs):
    if any([kwargs.get("file_b"), kwargs.get("ioeb")]):
        et.intersect(*args, **kwargs)


@click.command(help="filter significant result")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input dpsi file")
@click.option('-o', '--output', help="output directory")
@click.option('-dpsi', '--dpsi', type=float, default=0, help="dpsi threshold value, defualt=0")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value, defualt=0.05")
@click.option('-q', '--qval', type=float, default=1, help="qval threshold value, defualt=1")
@click.option('-adpsi', '--abs_dpsi', type=float, default=0, 
                help="absulte dpsi threshold value, defualt=0")
@click.option('-sep', '--sep', "sep", is_flag=True, default=False, 
                help="split file into two files according to dpsi > 0 and dpsi < 0")
@click.option('-app','--app', required=True, type=click.Choice(["auto", "SUPPA2", "rMATS"]),
                default="auto", help="the program that generates event file, defualt='auto'")
def sigfilter(*args, **kwargs):
    et.sigfilter(*args, **kwargs)
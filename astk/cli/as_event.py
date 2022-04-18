from .config import *
from astk.constant import AS_type
from astk.utils import event as et


@click.command(help="length distribution")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True),
                required=True,  help='AS ioe file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-cl', '--custom_len', 'custom_len', cls=MultiOption, type=tuple, help="custom length")
@click.option('-nc', '--cluster', type=int, default=4, help="number of cluster")
@click.option('-bw', '--width', type=int, default=3, help="bin width")
@click.option('-lw', '--len_weight', type=float, default=2, help="length weight")
@click.option('--max_len', type=int, default=500, help="the max length of exon in clustering")
def len_dist(*args, **kwargs):
    et.len_dist(*args, **kwargs)


@click.command(help="length cluster")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-id', '--inDir', type=click.Path(exists=True), help='input direcory')
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-lr', '--lenRange', cls=MultiOption, type=tuple,
              default=(1, ), help="custom length")
def len_cluster(*args, **kwargs):
    et.len_cluster(*args, **kwargs)


@click.command()
@click.option('-i', '--input', 'infile', type=click.Path(exists=True), help='AS ioe file')
@click.option('-o', '--output', help="output path")
@click.option('-rg', '--range', "len_range", type=(int, int), required=True, help="length range")
def len_pick(*args, **kwargs):
    et.len_pick(*args, **kwargs)


@click.command(help="differential splicing analysis")
@click.option('-od', '--outdir', help="output directory")
@click.option('-md', '--metadata', type=click.Path(exists=True),
             help="contrast group metadata, generated by meta")
@click.option('-gtf', '--gtf', type=click.Path(exists=True), help="gene annotation gtf file")
@click.option('-t', '--event_type', type=click.Choice(['all']+AS_type),
             default="all", help="gene annotation gtf file")
@click.option('-m', '--method', type=click.Choice(['empirical', 'classical']),
             default="empirical", help="gene annotation gtf file")
@click.option('--exon_len', type=int, default=100,
             help="Defines the number of nucleotides to display in the output GTF. (Default: 100 nt)")             
@click.option('-pg', '--poolGenes', default=False, is_flag=True, 
             help="pool together overlapping genes")
@click.option('--tpm_col', type=int, default=4, help="TPM columns index")
def diff_splice(*args, **kwargs):
    et.diff_splice(*args, **kwargs)


@click.command(help="filter significant result")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)  
@click.option('-od', '--outDir', required=True, help="output directory")
@click.option('-dpsi', '--dpsi', type=float, default=0, help="dpsi threshold value")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value")
@click.option('-adpsi', '--abs_dpsi', type=float, default=0, help="absulte dpsi threshold value")
@click.option('-pf1', '--psiFile1', cls=MultiOption,type=tuple, default=(),
             help="psi files of condtion 1")
@click.option('-pf2', '--psiFile2', cls=MultiOption,type=tuple, default=(),
             help="psi files of condtion 2")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
                default="tsv", help="out  file format ")
def sigfilter(*args, **kwargs):
    et.sigfilter(*args, **kwargs)


@click.command(help="filter psi result")
@click.argument('file', type=click.Path(exists=True), required=True)
#@click.option('-md', '--metadata', type=click.Path(exists=True), help="metadata file")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-psi', '--psi', type=float, default=0, help="psi threshold value")
@click.option('-qt', '--quantile', type=float, default=0, help="quantile threshold value")
# @click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
#                 default="tsv", help="out  file format ")
def psi_filter(*args, **kwargs):
    et.psi_filter(*args, **kwargs)


@click.command(help="intersect AS event")
@click.option('-a', "file_a", type=click.Path(exists=True), required=True, help="file a")
@click.option('-b', "file_b", type=click.Path(exists=True), required=True, help="file b")
@click.option('-o', "--output", required=True, help="output suffix")
def intersect(*args, **kwargs):
    et.intersect(*args, **kwargs)

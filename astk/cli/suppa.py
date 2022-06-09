"""
astk.cli.suppa
~~~~~~~~~~~~~~~~~
This module provide cli api for suppa2 function
"""

from .config import *
from astk.suppa import *
from astk.constant import AS_type


@click.command(help="generate the events from the GTF file")
@click.option('-gtf', '--gtf', type=click.Path(exists=True), help="a GTF format file")
@click.option('-et', '--eventType', "event_types", type=click.Choice(['ALL', 'SE', "RI"]), 
                default="ALL", help="AS event Type")
@click.option('-o', '--output', required=True, 
                help="name of the output file without any extension")         
@click.option('-ps', '--promoterSplit', "promoterSplit", is_flag=True, 
                help="save AS event that within promoter to a new ioe file")
def generateEvents(*args, **kwargs):
    generate_events(*args, **kwargs)


@click.command(help="calculates AS events PSI values")
@click.option('-o', '--output', type=click.Path(), help="output file path")
@click.option('-ioe', '--ioe', type=click.Path(exists=True), help="ioe file path")
@click.option('-qf', '--quantifyFile', "tpm_files", cls=MultiOption, type=click.Path(exists=True),
                help="transcript quantification files from salmon")
@click.option('--tpmThreshold', "tpm_th", type=float, default=0, 
                help="minimum transcript TPM value that using for calculates PSI")                
@click.option('--tpmCol', "tpm_col", type=int, default=2, help="TPM columns index,0-based")
@click.option('--txCol', "tx_col", type=int, default=0, 
                help="transcript ID columns index, 0-based")
def psiPerEvent(*args, **kwargs):
    output = kwargs.pop("output")
    psi_df = calculate_psi(*args, **kwargs)
    psi_df.to_csv(output, sep="\t")


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
def diffSplice(*args, **kwargs):
    et.diff_splice(*args, **kwargs)

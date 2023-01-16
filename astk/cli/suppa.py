"""
astk.cli.suppa
~~~~~~~~~~~~~~~~~
This module provide cli api for suppa2 function
"""

from .config import *
from astk.suppa import *
from astk.constant import AS_TYPE


@cli_fun.command(name="generateEvents", help="generate the events from the GTF file; short alias:ge")
@click.option('-gtf', '--gtf', type=click.Path(exists=True), help="a GTF format file")
@click.option('-et', '--eventType', "event_types", type=click.Choice(['ALL'] + AS_TYPE), 
                default="ALL", help="AS event Type")
@click.option('-o', '--output', required=True, 
                help="name of the output file without any extension")
@click.option('--idType', default="SUPPA2", type=click.Choice(["SUPPA2"]), 
                help="output event ID type, default='SUPPA2'") 
@click.option('-ep','--event-pos', type=click.Choice(['body', 'FT', 'LT']), 
                help="AS event exon that overlapping the transcript first or last \
                        terminal exon startCodon and stopCodon will save separately")
def generateEvents(*args, **kwargs):
    generate_events(*args, **kwargs)


@cli_fun.command(name="generatePsi", help="calculates AS events PSI values; alias: ge, psiPerEvent")
@click.option('-o', '--output', type=click.Path(), required=True, help="output file path")
@click.option('-ioe', '--ioe', type=click.Path(exists=True), help="ioe file path")
@click.option('-qf', '--quantifyFile', "tpm_files", cls=MultiOption, type=click.Path(exists=True),
                help="transcript quantification files from salmon")
@click.option('--tpmThreshold', "tpm_th", type=float, default=0, 
                help="minimum transcript TPM value that using for calculates PSI, default=0")                
@click.option('--tpmCol', "tpm_col", type=int, default=4, help="TPM columns index,0-based, default=4")
@click.option('--txCol', "tx_col", type=int, default=0, 
                help="transcript ID columns index, 0-based")
def generatePsi(*args, **kwargs):
    from pathlib import Path
    
    output = kwargs.pop("output")
    psi_df, tpm_df = calculate_psi(*args, **kwargs)
    psi_df.to_csv(output, sep="\t")
    tpm_df.to_csv(Path(output).with_suffix(".tpm"), sep="\t", index_label=False)


@cli_fun.command(name="diffSplice", help="differential splicing analysis; short alias: ds")
@click.option('-psi', '--psiFile', "psi_files", cls=MultiOption, required=True, 
                type=click.Path(exists=True),  help="AS event PSI files")
@click.option('-exp', '--expressFile', "exp_files", cls=MultiOption, required=True, 
                type=click.Path(exists=True), help="transcript quantification files from salmon")
@click.option('-ref', '--reference', type=click.Path(exists=True),
                help="ioe reference file")
@click.option('-o', '--output', help="output directory")
@click.option('-m', '--method', type=click.Choice(['empirical', 'classical']), default="empirical",
                help="The method to calculate the significance, default=empirical")
def diffSplice(*args, **kwargs):
    diff_splice(*args, **kwargs)


@cli_fun.command(help="differential splicing analysis workflow")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-md', '--metadata', type=click.Path(exists=True),
             help="contrast group metadata, generated by meta")
@click.option('-gtf', '--gtf', type=click.Path(exists=True), help="gene annotation gtf file")
@click.option('-et', '--event-type', type=click.Choice(['ALL']+AS_TYPE), cls=MultiOption,
                default="ALL", help="gene annotation gtf file")
@click.option('-m', '--method', type=click.Choice(['empirical', 'classical']), default="empirical",
                help="The method to calculate the significance, default=empirical")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value, default=0.05")
@click.option('-adpsi', '--abs-dpsi', type=float, default=0.1, help="absulte dpsi threshold value, default=0.1")
@click.option('--id-type', default="SUPPA2", type=click.Choice(["SUPPA2"]), 
                help="output event ID type, default='SUPPA2'")          
@click.option('--exon-len', type=int, default=100,
             help="Defines the number of nucleotides to display in the output GTF. (Default: 100 nt)")             
@click.option('-pg', '--pool-genes', default=False, is_flag=True, 
             help="pool together overlapping genes")
@click.option('--tpm-col', type=int, help="TPM columns index, default=4")
@click.option('--tpm-threshold', type=float, default=0, 
                help="Minimum TPM value to be included in the analysis. default=0")
@click.option('-ep','--event-pos', type=click.Choice(['body', 'FT', 'LT']), 
                help="AS event exon that overlapping the transcript first or last \
                        terminal exon startCodon and stopCodon will save separately")                
def dsflow(*args, **kwargs):
    if "ALL" in kwargs.get("event_type"):
        etypes = AS_TYPE
    else:
        etypes = kwargs.get("event_type")
    ds_flow(kwargs.get("metadata"),
        kwargs.get("gtf"),
        etypes,
        kwargs.get("outdir"),
        kwargs.get("method"),
        kwargs.get("id_type"),
        kwargs.get("pval"),
        kwargs.get("abs_dpsi"),
        kwargs.get("tpm_threshold"),
        kwargs.get("event_pos")
    )

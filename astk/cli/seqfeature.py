"""
astk.cli.seqfeature
~~~~~~~~~~~~~~~~~
This module provide sequence feature extraction cli api
"""

from .config import *
import astk.seqfeature.feature as sf



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
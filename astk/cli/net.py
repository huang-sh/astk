import click
from astk.network import co_splice_net

@click.command(help="Co-Splice Network")
@click.option('-o', '--output', required=True,
                help="AS sites flank fasta")
@click.option('-fa', '--fasta', required=True, type=click.Path(exists=True),
                help="AS sites flank fasta")
@click.option('-psi', '--psiMeta', required=True, type=click.Path(exists=True),
                help="psi meta file")
@click.option('-tq', '--tqMeta', required=True, type=click.Path(exists=True),
                help="transcript quantification meta file")
@click.option('-org', '--organism', required=True, help="organism")
@click.option('-db', '--database', type=click.Choice(['ATtRACT']),default='ATtRACT',
                help="RBP motif database")
@click.option('-txdb', '--txdb', type=click.Path(exists=True),
                help="TxDb file")
def sc_coSpliceNet(*args, **kwargs):
    co_splice_net(*args, **kwargs)
from .config import *
import astk.motif as mo


@click.command(["motifEnrich", "me"], help = "Motif Enrichment")
@click.option('-tf', "--tfasta", cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="input fasta files")
@click.option('-cf', "--cfasta", cls=MultiOption, type=click.Path(exists=True), 
                help="control fasta files")                
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                 default="CISBP-RNA", help="RBP motif database default=CISBP-RNA")
@click.option('-org', '--organism', help="RBP organism")
# @click.option('-mm', '--meme', type=click.Path(exists=True), 
#                 required=True, help="path to .meme format file")
def motif_enrich(*args, **kwargs):
    mo.motif_enrich(*args, **kwargs)

@click.command(help = "Motif Discovery and similarity comparision")
@click.option('-tf', "--tfasta", cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="fasta file")
@click.option('-cf', "--cfasta", cls=MultiOption, type=click.Path(exists=True), 
                help="control fasta files")                     
@click.option('-od', '--outdir', type=click.Path(), default=".",
                 help="output directory")
# @click.option('-mcmp', '--motifcmp', is_flag=True, 
#                 default = False, help="motif similarity comparision")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                 default="CISBP-RNA", help="RBP motif database default=CISBP-RNA")
@click.option('-org', '--organism', help="RBP organism")                              
@click.option('-pval', '--pvalue', type=float, default=0.05,
                help="Discovery pvalue cutoff, default=0.05")
@click.option('-eval', '--evalue', type=float, default=0.5,
                help="motif comparison pvalue cutoff, default=0.5")
@click.option('-minw', '--minw', type=int, default=5,
                 help="minimal motifs width,default=5")                    
@click.option('-maxw', '--maxw', type=int, default=15,
                 help="maximal motifs width, default=15")               
def motif_find(*args, **kwargs):
    mo.motif_find(*args, **kwargs)


@click.command(help = "Motif plot")
@click.option('-mi', "--motifId", "motifid", cls=MultiOption, type=str, 
                required=True, help="motif id")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                help="RBP motif database")
@click.option('-org', '--organism', help="RBP organism")                
@click.option('-mm', "--meme", type=click.Path(exists=True), 
                help="meme motif file")
@click.option('-o', '--output', type=click.Path(), required=True,
                 help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")                 
def motif_plot(*args, **kwargs):
    mo.motif_plot(*args, **kwargs)


@click.command(help = "generate motif map")
@click.option('-fa', '--fasta', cls=MultiOption, type=click.Path(exists=True), 
                required=True, help="fasta files")
@click.option('-n', '--name', cls=MultiOption, type=str,
                help="fasta file names")
@click.option('-c', '--center', cls=MultiOption, type=str,
                help="fasta files names")                
@click.option('-mm', '--meme', required=True, type=click.Path(exists=True), 
                help="meme motif file")
@click.option('-od', '--outdir', default=".", type=click.Path(), 
                help="meme motif file")
@click.option('-b', '--binsize', default=20, 
                help="the window width for scanning motif, default=20")
@click.option('-s', '--step', default=10, 
                help="the slide window size, default=10")                
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def mmap(*args, **kwargs):
    mo.mmap(*args, **kwargs)


@click.command(help="Eukaryotic Linear Motif Searching")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-g', '--genome', type=click.Choice(['mm10', 'hg38']),
                  required=True, help="genome assembly")
def elms(*args, **kwargs):
    mo.search_elm(*args, **kwargs)

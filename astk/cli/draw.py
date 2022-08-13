from .config import *
from astk import draw


@click.command(help="Gene Set Enrichment Analysis ploting")
@click.option('-id', '--id', "termid", cls=MultiOption, type=str, 
                help="term id")
@click.option('-o', '--output', help="output figure path")
@click.option('-rd', '--RData', help="output figure path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def gseplot(*args, **kwargs):
    draw.gseplot(*args, **kwargs)

@click.command(help="draw UpSet plots for AS events")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")           
@click.option('-o', '--output', required=True, help="output figure path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str, 
                help="x xlabel")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("This flag is present then a dpsi file will divide "
                  "two part according to |dpsi| > 0 and |dpsi| < 0"))                     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def upset(*args, **kwargs):
    draw.upset(*args, **kwargs)


@click.command(help="Volcano plot analysis for dPSI")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input psi files")    
@click.option('-o', '--output', help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def volcano(*args, **kwargs):
    draw.volcano(*args, **kwargs)


@click.command(help="PCA analysis for PSI")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def pca(*args, **kwargs):
    draw.pca(*args, **kwargs)


@click.command(help="Heatmap plot for PSI")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def heatmap(*args, **kwargs):
    draw.heatmap(*args, **kwargs)


@click.command(help="barplot ")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str,
             help="input dpsi names")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("This flag is present then a dpsi file will divide "
                  "two part according to |dpsi| > 0 and |dpsi| < 0"))             
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def barplot(*args, **kwargs):
    draw.barplot(*args, **kwargs)

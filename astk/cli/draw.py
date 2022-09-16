from .config import *
from astk import draw
from astk.utils import sniff_fig_fmt


@click.command(help="Gene Set Enrichment Analysis ploting")
@click.option('-id', '--id', "termid", cls=MultiOption, type=str, 
                help="term id")
@click.option('-o', '--output', help="output figure path")
@click.option('-rd', '--RData', help="output figure path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def gseplot(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])
    draw.gseplot(*args, **kwargs)

@click.command(help="draw UpSet plots for AS events")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")           
@click.option('-o', '--output', required=True, help="output figure path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str,  help="x xlabel")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("AS events can be divided into two groups based on dPSI values \
                   (group +: dPSI > 0, group -: dPSI < 0)"))                 
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def upset(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])    
    draw.upset(*args, **kwargs)


@click.command(help="Volcano plot analysis for dPSI")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input psi files")    
@click.option('-o', '--output', help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def volcano(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])       
    draw.volcano(*args, **kwargs)


@click.command(help="PCA analysis for PSI")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
@click.option('-gn', '--groupName', cls=MultiOption, type=str, help="group names")
def pca(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])
    if gn := kwargs.get("groupname", False):
        if len(gn) != len(kwargs["files"]):
            print("-gn/--groupName values number must be same as -i/--input")
            exit()
    draw.pca(*args, **kwargs)


@click.command(help="Heatmap plot for PSI")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                 default="auto", help="output figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def heatmap(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])       
    draw.heatmap(*args, **kwargs)


@click.command(help="barplot ")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input psi files")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str,
             help="input dpsi names")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("AS events can be divided into two groups based on dPSI values \
                   (group +: dPSI > 0, group -: dPSI < 0)"))   
@click.option('-fmt', '--format', "fmt", type=click.Choice(['auto', 'png', 'pdf', 'pptx']),
                default="auto", help="output figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def barplot(*args, **kwargs):
    if kwargs["fmt"] == "auto":
        kwargs["fmt"] = sniff_fig_fmt(kwargs["output"])
    draw.barplot(*args, **kwargs)

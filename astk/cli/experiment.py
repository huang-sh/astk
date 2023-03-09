"""
astk.cli.experiment
~~~~~~~~~~~~~~~~~
This module provide some experimental function
"""

from .config import *
import astk.utils._cli_func as ul
from astk import draw


@cli_fun.command(name="vcor", help="correlation heatmap")
@click.option('-i','--input', 'files', cls=MultiOption, type=click.Path(exists=True),
                required=True , help="input files")
@click.option('-o', '--output', type=click.Path(), help="output path")
@click.option('-m', '--method', type=click.Choice(['pearson', 'kendall', 'spearman']), 
                default="spearman", show_default=True, help="method of correlation")
@click.option('-fw', '--width', default=6, show_default=True, help="figure width")
@click.option('-fh', '--height', default=6, show_default=True, help="figure height")
@click.option('-cmap', '--colormap', default="RdYlBu", show_default=True, help="matplotlib colormap name")
@click.option('-title', '--title', default="Correlation Plot", help="figure title")
@click.option('-gl', '--groupLabel', cls=MultiOption, help="label prefix for multiple files")
def sub_vcor(*args, **kwargs):
    from matplotlib.pyplot import colormaps

    if kwargs["colormap"] not in colormaps():
        msg  = f"'{kwargs['colormap']}' is not a valid value for colormap name"
        msg += f"; supported values are {', '.join(map(repr, colormaps()))}"        
        BadParameter(msg)
    
    draw.plot_cor_heatmap(**kwargs)

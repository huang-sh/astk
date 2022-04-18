
from .utils import *
from .net import *
from .motif import *
from .gsea import *
from .as_event import *
from .epi import *
from .draw import *
from .config import CustomMultiCommand


@click.group(cls=CustomMultiCommand, 
        context_settings=dict(help_option_names=['-h', '--help']))
def cli_fun():
     pass


## AS event processing module
cli_fun.add_command(len_dist, name=["lenDist", "ld"])
cli_fun.add_command(len_cluster, name=["lenCluster", "lc"])
cli_fun.add_command(len_pick, name=["lenPick", "lp"])
cli_fun.add_command(diff_splice, name=["diffSplice", "ds"])
cli_fun.add_command(sigfilter, name=["sigfilter", "sf"])
cli_fun.add_command(psi_filter, name=["psiFilter", "pf"])
cli_fun.add_command(intersect, name=["intersect"])

# co-splicing network module
cli_fun.add_command(coSpliceNet, name=["coSpliceNet", "csnet"])

# motif analysis module
cli_fun.add_command(motif_enrich, name=["motifEnrich", "me"])
cli_fun.add_command(motif_find, name=["motifFind", "mf"])
cli_fun.add_command(motif_plot, name=["motifPlot", "mp"])
cli_fun.add_command(mmap)
cli_fun.add_command(elms)


## gene set enrichment analysis module
cli_fun.add_command(gsea_fun, name=["gsea"])
cli_fun.add_command(enrich)
cli_fun.add_command(enrich_cmp, name=["enrichCompare", "ecmp"])
cli_fun.add_command(nease_sc, name=["nease"])
# cli_fun.add_command(enrich_lc, help=["enrichLenCluster", "elc"])

## epigenetic analysis module
cli_fun.add_command(epi_sc, name=["epi"])
cli_fun.add_command(epihm)
cli_fun.add_command(epiline)

## drawing module
cli_fun.add_command(gseplot)
cli_fun.add_command(upset)
cli_fun.add_command(volcano, name=["volcano", "vol"])
cli_fun.add_command(pca)
cli_fun.add_command(heatmap, name=["heatmap", "hm"])
cli_fun.add_command(barplot)

## utils module
cli_fun.add_command(meta)
cli_fun.add_command(anchor)
cli_fun.add_command(install)
cli_fun.add_command(getmeme)
cli_fun.add_command(getcoor)
cli_fun.add_command(list_, name=["list", "ls"])
cli_fun.add_command(mktxdb)

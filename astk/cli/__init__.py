
from .utils import *
from .net import *
from .motif import *
from .gsea import *
from .as_event import *
from .epi import *
from .draw import *
from .suppa import *
from .seqfeature import *
from .config import CustomMultiCommand, sc_setting


@click.group(cls=CustomMultiCommand, 
        context_settings=dict(help_option_names=['-h', '--help']))
def cli_fun():
     pass

## astk configure setting
cli_fun.add_command(sc_setting, name=["config"])

## suppa2 module
cli_fun.add_command(generateEvents, name=["generateEvent", "ge"])
cli_fun.add_command(compute_psi, name=["generatePsi", "gp", "psiPerEvent"])
cli_fun.add_command(diffSplice, name=["diffSplice", "ds"])
cli_fun.add_command(dsflow, name=["dsflow"])

## AS event processing module
cli_fun.add_command(len_dist, name=["lenDist", "ld"])
cli_fun.add_command(len_cluster, name=["lenCluster", "lc"])
cli_fun.add_command(len_pick, name=["lenPick", "lp"])
cli_fun.add_command(sigfilter, name=["sigFilter", "sf"])
cli_fun.add_command(psi_filter, name=["psiFilter", "pf"])
cli_fun.add_command(intersect, name=["intersect"])

# co-splicing network module
cli_fun.add_command(sc_coSpliceNet, name=["coSpliceNet", "csnet"])

# motif analysis module
cli_fun.add_command(motif_enrich, name=["motifEnrich", "me"])
cli_fun.add_command(motif_find, name=["motifFind", "mf"])
cli_fun.add_command(motif_plot, name=["motifPlot", "mp"])
# cli_fun.add_command(mmap)
cli_fun.add_command(elms)


## gene set enrichment analysis module
cli_fun.add_command(gsea_fun, name=["gsea"])
cli_fun.add_command(enrich)
cli_fun.add_command(enrich_cmp, name=["enrichCompare", "ecmp"])
cli_fun.add_command(nease_sc, name=["nease"])
cli_fun.add_command(sc_neasecmp, name=["neaseCompare", "necmp"])

## epigenetic analysis module
#cli_fun.add_command(epi_sc, name=["epi"])
# cli_fun.add_command(epihm)
# cli_fun.add_command(epi_profile_sc, name=["epiProfile", "ep"])
cli_fun.add_command(signal_profile, name=["signalProfile", "sp"])
cli_fun.add_command(signal_profile2, name=["sp2"])

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
cli_fun.add_command(getgene)
cli_fun.add_command(sub_merge_files, name=["merge"])


## seq feature
cli_fun.add_command(sc_extract, name=["seqfeature", "seqf"])
cli_fun.add_command(sc_seqlogo, name=["seqlogo"])
cli_fun.add_command(sc_splice_score, name=["sss", "spliceScore"])
cli_fun.add_command(sc_get_elen, name=["elen"])
cli_fun.add_command(sc_get_gcc, name=["gcc"])
cli_fun.add_command(sc_ssscmp, name=["ssscmp"])
cli_fun.add_command(sc_cmp_value, name=["vcmp"])


from click_option_group import optgroup

from .config import *
from astk.constant import NEASE_DATABASE
from astk import gsea
from astk.utils import detect_file_info


@cli_fun.command(name="gsea", help="Gene Set Enrichment Analysis")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                required=True, help="input dpsi files")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-n', '--name', default="GSEA", help="output name prefix")
@click.option('-pval', '--pvalue', type=float, default=0.2, show_default=True, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO']), show_default=True, 
                default="GO", help="enrich database")
@click.option('-gt', '--gene-id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", show_default=True, help="gene ID type, defualt='ENSEMBL'") 
@click.option('-ont', '--ontology', type=click.Choice(['BP', 'MF', 'CC']), show_default=True, 
                default="BP", help="one of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three")  
@click.option('-org', '--organism', type=click.Choice(['hs', 'mm']), required=True,  help="organism")
@click.option('-app','--app', required=True, type=click.Choice(["SUPPA2", "rMATS", "EventPointer"]), 
                help="the software that generates event file")                 
def gsea_fun(*args, **kwargs):
    if kwargs["app"] is None: # it will not run
        kwargs["app"] = detect_file_info(kwargs["file"])["app"] 
    gsea.gsea_fun(*args, **kwargs)


@cli_fun.command(help="Over representation enrichment analysis")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                required=True, help="input dpsi files")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-pval', '--pvalue', type=float, default=0.1, show_default=True, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, show_default=True, help="qvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG', 'Reactome']), 
                default="GO", show_default=True, help="enrich database")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP", show_default=True,
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three")                
@click.option('-gi', '--gene-id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", show_default=True, help="gene ID type")
@click.option('-org', '--organism', type=click.Choice(['hs', 'mm']), required=True, 
                show_default=True, help="organism")
@click.option('--simple', is_flag=True, default=False, show_default=True, help="simplify GO enrichment")
@click.option('-app','--app', type=click.Choice(["auto", "SUPPA2", "rMATS", "EventPointer"]), 
                default="auto", show_default=True, help="the software that generates event file")
@fig_common_options()          
def enrich(*args, **kwargs):
    if kwargs["app"] == "auto":
        kwargs["app"] = detect_file_info(kwargs["file"])["app"]
    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = "png"
    gsea.enrich(*args, **kwargs)


@cli_fun.command(name="enrichCompare", help="function enrichment comparison; short alias: ecmp")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                required=True, help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")          
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG', 'Reactome']), 
                default="GO", show_default=True, help="enrich database")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP", show_default=True,
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three.")                  
@click.option('-pval', '--pvalue', type=float, default=0.1, show_default=True, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, show_default=True, help="qvalue cutoff")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str, help="xlabel")
@click.option('-gi', '--gene-id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", show_default=True, help="gene ID type")
@click.option('-org', '--organism', type=click.Choice(['hs', 'mm']), required=True, 
                show_default=True, help="organism")
@click.option('-app','--app', type=click.Choice(["auto", "SUPPA2", "rMATS", "EventPointer"]), 
                default="auto", show_default=True, help="the software that generates event file")
@fig_common_options()
def enrich_cmp(*args, **kwargs):
    if kwargs["app"] == "auto": 
        # Assuming that the input files are all in the same format
        kwargs["app"] = detect_file_info(kwargs["files"][0])["app"]
    if kwargs["figfmt"] == "auto":
        kwargs["figfmt"] = "png"
    gsea.enrich_cmp(*args, **kwargs)


@cli_fun.command(name="nease", help="Functional enrichment with NEASE")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                required=True, help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-db', '--database', cls=MultiOption, type=click.Choice(NEASE_DATABASE), 
                default=["Reactome"], help="nease enrich database, \
                [PharmGKB|HumanCyc|Wikipathways|Reactome|KEGG|SMPDB|Signalink|NetPath|EHMN|INOH|BioCarta|PID]")
@click.option('-org', '--organism', default='Human', type=click.Choice(['Human']),
                help="organism")
def nease_sc(*args, **kwargs):
    gsea.nease_sc(*args, **kwargs)


@cli_fun.command(name="neaseCompare", 
                help="Functional enrichment comparision with NEASE; short alias: ncmp")
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-n', '--showNumber', 'num', type=int, default=15, help="qvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="qvalue cutoff")
@click.option('-db', '--database', cls=MultiOption, type=click.Choice(NEASE_DATABASE), 
                default=["Reactome"], help="nease enrich database, \
                [PharmGKB|HumanCyc|Wikipathways|Reactome|KEGG|SMPDB|Signalink|NetPath|EHMN|INOH|BioCarta|PID]")
@click.option('-org', '--organism', default='Human', type=click.Choice(['Human']),
                help="organism")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str, help="xlabel")
@click.option('-ff', '--figFormat', type=click.Choice(['png', 'pdf', 'tiff', 'jpeg']), 
                default="png", help="output figure format")
@click.option('-fw', '--width', type=float, default=6, help="figure width, default=6 inches")
@click.option('-fh', '--height', type=float, default=6, help="figure height, default=6 inches")        
def sc_neasecmp(*args, **kwargs):
    fn = len(kwargs["files"])
    if kwargs["xlabel"] is None:
        kwargs["xlabel"] = [f"c{i+1}" for i in range(fn)]
    if fn < 2:
        raise UsageError("for -i/--input should be more than two files")
    gsea.neasecmp(*args, **kwargs)

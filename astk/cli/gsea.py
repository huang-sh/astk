from .config import *
from astk.constant import NEASE_DATABASE
from astk import gsea


@click.command(help="Gene Set Enrichment Analysis")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-n', '--name', default="GSEA", help="output name prefix")
@click.option('-pval', '--pvalue', type=float, default=0.2, help="pvalue cutoff, defualt=0.2")
@click.option('-db', '--database', type=click.Choice(['GO']), 
                default="GO", help="enrich database, defualt='GO'")
@click.option('-gt', '--geneId', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type, defualt='ENSEMBL'")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -org' to view more ")
@click.option('-ont', type=click.Choice(['BP', 'MF', 'CC']), 
                default="BP", help="one of 'BP', 'MF', and 'CC' subontologies.")  
@click.option('-org', '--keggOrganism', "organism", default = "",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")                            
def gsea_fun(*args, **kwargs):
    gsea.gsea_fun(*args, **kwargs)


@click.command(help="Over representation enrichment analysis")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff, default=0.1")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff, default=0.1")
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG', 'Reactome']), 
                default="GO", help="enrich database GO|KEGG|Reactome, default=GO")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP",
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. default=BP")                
@click.option('-gene_id', '-gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")
@click.option('-org', '--organism', type=click.Choice(['hs', 'mm']), 
                required=True, help="organism")
@click.option('--simple', is_flag=True, help="simplify GO enrichment")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="pdf", help="output figure format") 
@click.option('-w', '--width', default=10, help="fig width, default=10 inches")
@click.option('--height', default=12, help="fig height, default=12 inches")    
def enrich(*args, **kwargs):
    gsea.enrich(*args, **kwargs)


@click.command()
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")          
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG', 'Reactome']), 
                default="GO", help="enrich database, default='GO'")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP",
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. default=BP")                  
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-xl', '--xlabel', cls=MultiOption, type=str, help="xlabel")
@click.option('-gene_id', '--gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-org', '--organism', type=click.Choice(['hs', 'mm']), required=True,
                help="organism: hs|mm")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                default="pdf", help="output figure format")
@click.option('-fw', '--width', type=float, help="figure width")
@click.option('-fh', '--height', type=float, help="figure height")
def enrich_cmp(*args, **kwargs):
    gsea.enrich_cmp(*args, **kwargs)


@click.command()
@click.option('-i', '--input', "files", cls=MultiOption, type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                required=True, help="cluster information file")
@click.option('-m', '--merge', is_flag=True, help="enrich with multiple merged files")              
@click.option('-db', '--database', type=click.Choice(['GO']), 
                default="GO", help="enrich database")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-gene_id', '--geneId', "gene_id", type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -orgdb' to view more ")
@click.option('-org', '--keggOrganism', "kegg_organism",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html") 
def enrich_lc(*args, **kwargs):
    pass


@click.command(help="Functional enrichment with NEASE")
@click.option('-i', '--input', "file", type=click.Path(exists=True),
                help="input dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-db', '--database', cls=MultiOption, type=click.Choice(NEASE_DATABASE), 
                default=["Reactome"], help="nease enrich database, \
                [PharmGKB|HumanCyc|Wikipathways|Reactome|KEGG|SMPDB|Signalink|NetPath|EHMN|INOH|BioCarta|PID]")
@click.option('-org', '--organism', default='Human', type=click.Choice(['Human']),
                help="organism")
def nease_sc(*args, **kwargs):
    gsea.nease_sc(*args, **kwargs)


@click.command(help="Functional enrichment comparision with NEASE")
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
@click.option('-fw', '--width', type=float, help="figure width")
@click.option('-fh', '--height', type=float, help="figure height")           
def sc_neasecmp(*args, **kwargs):
    fn = len(kwargs["files"])
    if kwargs["xlabel"] is None:
        kwargs["xlabel"] = [f"c{i+1}" for i in range(fn)]
    if fn < 2:
        raise UsageError("for -i/--input should be more than two files")
    gsea.neasecmp(*args, **kwargs)

# -*- coding: utf-8 -*-
"""
astk.astk
~~~~~~~~~~~~~~~~~
This module provides program command line api.
"""

import sys
import shutil
import subprocess
from pathlib import Path

from .cli_config import *
from .constant import *

@click.group(cls=CustomMultiCommand, 
        context_settings=dict(help_option_names=['-h', '--help']))
def cli():
     pass


@cli.command(help="generate metadata template file")
@click.option('-p1', '--control', cls=MultiOption,type=tuple, default=(),
                required=True, help="file path for condtion 1")
@click.option('-p2', '--treatment', cls=MultiOption,type=tuple, default=(),
                required=True, help="file path for condtion 2")
@click.option('-gn', '--groupName', cls=MultiOption, type=tuple,
                default=(), help="group name")                
@click.option('-repN', '--replicate', cls=MultiOption,type=tuple, 
                help="replicate, number")
@click.option('-o', '--output', required=True, type=click.Path(),
                help='metadata output path')                
@click.option('-repN1', '--replicate1', cls=MultiOption,type=tuple, 
                help="replicate1, number")
@click.option('-repN2', '--replicate2', cls=MultiOption,type=tuple, 
                help="replicate2, number")
@click.option('-fn', '--filename', 'is_fn', is_flag=True, help="file name")
@click.option('--split', cls=MultiOption, type=tuple, help="name split symbol and index")    
def meta(output, replicate, groupname, control, treatment, replicate1, replicate2, **kwargs):
    from .meta_template import Template

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit() 
    if replicate:
        repN1 = [int(i) for i in replicate]
        repN2 = [int(i) for i in replicate]
    elif all([replicate1, replicate1]):
        repN1 = [int(i) for i in replicate1]
        repN2 = [int(i) for i in replicate2]
    else:
        print("repN1 and repN2 must be set!")    
        sys.exit()
    try:
        tp = Template()
        tp.complete_df(groupname, control, treatment, repN1, repN2, **kwargs)
        tp.to_csv(output)
        tp.to_json(output)
    except BaseException as e:
        print(e)

@cli.command(["lenDist", "ld"], help="length distribution")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True),
                required=True,  help='AS ioe file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-cl', '--custom_len', 'custom_len', cls=MultiOption, type=tuple, help="custom length")
@click.option('-nc', '--cluster', type=int, default=4, help="number of cluster")
@click.option('-bw', '--width', type=int, default=3, help="bin width")
@click.option('-lw', '--len_weight', type=float, default=2, help="length weight")
@click.option('--max_len', type=int, default=500, help="the max length of exon in clustering")
def len_dist(infile, output, custom_len, cluster, width, len_weight, max_len):
    import pandas as pd
    from . import utils  as ul

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   
    ioe_df = pd.read_csv(infile, sep="\t")
    if "event_id" not in ioe_df.columns:
        print("file does not support!")
        exit()

    info_df = ul.extract_info(ioe_df)
    info_df["event_id"] = ioe_df["event_id"]
    output = Path(output).with_suffix(".png")
    if custom_len:
        lens = [int(i) for i in custom_len]
        df = ul.custome_cluster_len(info_df, output, lens, width=width, max_len=max_len)
    else:
        df = ul.cluster_len(info_df, output, n_cls=cluster, max_len=max_len, len_weight=len_weight, width=width)
    df.to_csv(Path(output).parent / f"{Path(output).stem}_cls_info.csv", index=False)

@cli.command(["lenCluster", "lc"], help="length cluster")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-id', '--inDir', type=click.Path(exists=True), help='input direcory')
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-lr', '--lenRange', cls=MultiOption, type=tuple,
              default=(1, ), help="custom length")
def len_cluster(files, indir, outdir, lenrange):
    from . import utils  as ul

    outdir = Path(outdir)
    
    outdir = Path(outdir).absolute()
    od_name = outdir.name

    lrs = list(map(int, lenrange))
    coor_ls = [(lrs[i], lrs[i+1]) for i in range(len(lrs)-1)]

    if indir:
        files = list(Path(indir).glob("*psi"))
    else:
        files = [Path(i) for i in files]
    
    if len(files) == 1:
        file = Path(files[0])
        outdir.mkdir(exist_ok=True)
        for s, e in coor_ls:
            outfile = outdir / f"{file.stem}_{s}-{e}{file.suffix}"
            ul.df_len_select(file, outfile, s, e)
    else:
        for s, e in coor_ls:
            s_outdir = outdir.with_name(f"{od_name}_{s}-{e}")
            s_outdir.mkdir(exist_ok=True)
            for file in files:
                outfile = s_outdir / Path(file).name
                ul.df_len_select(file, outfile, s, e+1)

@cli.command(["lenPick", "lp"])
@click.option('-i', '--input', 'infile', type=click.Path(exists=True), help='AS ioe file')
@click.option('-o', '--output', help="output path")
@click.option('-rg', '--range', "len_range", type=(int, int), required=True, help="length range")
def len_pick(infile, output, len_range):
    from . import event_id as ei
    import pandas as pd

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   
    AS_len = lambda x: ei.SuppaEventID(x).alter_element_len

    df = pd.read_csv(infile, sep="\t", index_col=0)
    cols = df.columns
    df["event_id"] = df.index
    df["len"] = df["event_id"].apply(AS_len)
    s, e = len_range
    pdf = df.loc[(s <= df["len"]) & ( df["len"] < e), cols]
    pdf.to_csv(output, index=True, sep="\t", na_rep="nan", index_label=False)


AS_type = ['SE', "A5", "A3", "MX", "RI", 'AF', 'AL']

@cli.command(["diffSplice", "ds"], help="differential splicing analysis")
@click.option('-od', '--outdir', help="output directory")
@click.option('-md', '--metadata', type=click.Path(exists=True),
             help="contrast group metadata, generated by meta")
@click.option('-gtf', '--gtf', type=click.Path(exists=True), help="gene annotation gtf file")
@click.option('-t', '--event_type', type=click.Choice(['all']+AS_type),
             default="all", help="gene annotation gtf file")
@click.option('-m', '--method', type=click.Choice(['empirical', 'classical']),
             default="empirical", help="gene annotation gtf file")
@click.option('--exon_len', type=int, default=100,
             help="Defines the number of nucleotides to display in the output GTF. (Default: 100 nt)")             
@click.option('-pg', '--poolGenes', default=False, is_flag=True, 
             help="pool together overlapping genes")
@click.option('--tpm_col', type=int, default=4, help="TPM columns index")
def diff_splice(outdir, metadata, gtf, event_type, method, exon_len, poolgenes, tpm_col):
    from . import utils  as ul

    if event_type == "all":
        event_types = AS_type
    else:
        event_types = [event_type]
    try:
        dsi = ul.DiffSplice(outdir, metadata, gtf, event_types, exon_len, poolgenes)
        dsi.ds(method)
    except BaseException as e:
        print(e)


@cli.command(["sigfilter", "sf"], help="filter significant result")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)  
@click.option('-od', '--outDir', required=True, help="output directory")
@click.option('-dpsi', '--dpsi', type=float, default=0, help="dpsi threshold value")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value")
@click.option('-adpsi', '--abs_dpsi', type=float, default=0, help="absulte dpsi threshold value")
@click.option('-pf1', '--psiFile1', cls=MultiOption,type=tuple, default=(),
             help="psi files of condtion 1")
@click.option('-pf2', '--psiFile2', cls=MultiOption,type=tuple, default=(),
             help="psi files of condtion 2")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
                default="tsv", help="out  file format ")
def sigfilter(files, outdir, dpsi, pval, abs_dpsi, psifile1, psifile2, fmt):
    from . import select as sl

    for idx, dpsi_file in enumerate(files):
        if len(files) == len(psifile1) and len(files) == len(psifile2):
            psifiles = (psifile1[idx], psifile2[idx])
        else:
            psifiles = ()
            
        sf = sl.SigFilter(dpsi_file, outdir, dpsi, pval, abs_dpsi, psifiles, fmt)
        sf.run()


@cli.command(["psiFilter", "pf"], help="filter psi result")
@click.argument('file', type=click.Path(exists=True), required=True)
#@click.option('-md', '--metadata', type=click.Path(exists=True), help="metadata file")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-psi', '--psi', type=float, default=0, help="psi threshold value")
@click.option('-qt', '--quantile', type=float, default=0, help="quantile threshold value")
# @click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
#                 default="tsv", help="out  file format ")
def psi_filter(file, output, psi, quantile):
    from . import select as sl

    if file:
        pf = sl.PsiFilter(file, psi, quantile)
        pf.run(output)
  

@cli.command(help="Over representation enrichment analysis")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO']), 
                default="GO", help="enrich database")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP",
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. default=BP")                
@click.option('-gene_id', '-gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -org' to view more ")
@click.option('-ko', '--keggOrganism', "kegg_organism",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="pdf", help="out figure format") 
@click.option('-w', '--width', default=10, help="fig width, default=10 inches")
@click.option('-h', '--height', default=12, help="fig height, default=12 inches")    
def enrich(file, outdir, pvalue, qvalue, database, ontology,  
            gene_id, orgdb, kegg_organism, fmt, width, height):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "enrich.R"
    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -org' to view more")
    if database == "KEGG":
        if not kegg_organism:
            print("Error: --kegg_organism is required!")
            exit()
    else:
        kegg_organism = "0"

    Path(outdir).mkdir(exist_ok=True)
    if database == "KEGG":
        ul.check_kegg_RData(kegg_organism)

    param_dic = {
        "file": file,
        "outdir": outdir, 
        "orgdb": org_db,
        "database": database,
        "pval": pvalue,
        "qval": qvalue,
        "keggorganism": kegg_organism,
        "genetype": gene_id,
        "ontology": ontology,
        "width": width,
        "height": height,
        "fmt": fmt
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(["enrichCompare", "ecmp"])
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)  
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")            
@click.option('-db', '--database', type=click.Choice(['GO']), 
                default="GO", help="enrich database")
@click.option('-ont', '--ontology', type=click.Choice(['ALL', 'BP', 'CC','MF']), default="BP",
                help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. default=BP")                  
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-xl', '--xlabel', cls=MultiOption, type=tuple, 
                help="xlabel")
@click.option('-gene_id', '--gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -orgdb' to view more ")
@click.option('-ko', '--keggOrganism', "kegg_organism",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="pdf", help="out figure format") 
@click.option('-w', '--width', default=10, help="fig width, default=6 inches")
@click.option('-h', '--height', default=12, help="fig height, default=6 inches")                
def enrich_cmp(files, outdir, cluster, database, ontology, pvalue, qvalue, 
                xlabel, gene_id, orgdb, kegg_organism, fmt, width, height):
    from . import utils  as ul

    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -orgdb' to view more")
    if database == "KEGG":
        ul.check_kegg_RData(kegg_organism)
        if not kegg_organism:
            print("Error: --kegg_organism is required!")
            exit()
    else:
        kegg_organism = "0"
    rscript = Path(__file__).parent / "R" / "enrichCompare.R"
    Path(outdir).mkdir(exist_ok=True)
    cluster = cluster if cluster else "0"
    param_dic = {
        "files": files,
        "outdir": outdir, 
        "clusterfile": cluster, 
        "orgdb": org_db,
        "database": database,
        "pval": pvalue,
        "qval": qvalue,
        "xlabel": xlabel, 
        "keggorganism": kegg_organism,
        "genetype": gene_id,
        "ontology": ontology,
        "width": width,
        "height": height,
        "fmt": fmt
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(["enrichLenCluster", "elc"])
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
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
def enrich_lc(files, outdir, cluster, merge, database, pvalue, qvalue,
              gene_id, orgdb, kegg_organism):
    from . import utils  as ul

    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -orgdb' to view more")                
    if database == "KEGG":
        ul.check_kegg_RData(kegg_organism)
        if not kegg_organism:
            print("Error: --kegg_organism is required!")
            exit()
    else:
        kegg_organism = "0"
    rscript = Path(__file__).parent / "R" / "enrichLenCluster.R"
    Path(outdir).mkdir(exist_ok=True)
    merge = "1" if merge else "0"
    for file in files:
        params = [outdir, str(pvalue), str(qvalue), database, cluster,
                 org_db, gene_id, kegg_organism, file]
        info = subprocess.Popen(["Rscript", str(rscript), *params])
        if database == "KEGG" and not ul.check_kegg_RData(kegg_organism):
             info.wait() 
    else:
        if merge:
            params = [outdir, str(pvalue), str(qvalue), database, cluster,
                 org_db, gene_id, kegg_organism, *files]
            merge_info = subprocess.Popen(["Rscript", str(rscript), *params])
            merge_info.wait()
        else:
            info.wait()

@cli.command(["nease"], help="Functional enrichment with NEASE")
@click.argument('file', type=click.Path(exists=True))
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(NEASE_DATABASE), 
                default="Reactome", help="nease enrich database")
@click.option('-org', '--organism', default='Human', type=click.Choice(['Human']),
                help="organism") 
def nease_sc(file, outdir, pvalue, database, organism):
    import pandas as pd
    from .enrich import nease_enrich
    from . import event_id as ei

    get_geneid = lambda x: ei.SuppaEventID(x).gene_id.split(".")[0]
    get_start = lambda x: ei.SuppaEventID(x).alter_element_coor[0]
    get_end = lambda x: ei.SuppaEventID(x).alter_element_coor[1]

    df = pd.read_csv(file, sep="\t", index_col=0)
    df["event_id"] = df.index
    
    nease_input = pd.DataFrame({
                    "Gene stable ID": df["event_id"].apply(get_geneid),
                    "new_start": df["event_id"].apply(get_start), 
                    "new_end": df["event_id"].apply(get_end),
                    "beta": df[df.columns[0]].values})

    nease_enrich(nease_input, outdir, database=[database], organism=organism, cutoff=pvalue)


@cli.command(["necmp"], help="Functional enrichment comparision with NEASE")
def neasecmp_sc():
    pass

@cli.command(["volcano", "vol"], help="Volcano plot analysis for dPSI")
@click.argument('file', type=click.Path(exists=True), required=True)     
@click.option('-o', '--output', help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def volcano(file, output, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "volcano.R"
    param_dic = {
        "file": file,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help="PCA analysis for PSI")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)     
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def pca(files, output, fmt, width, height, resolution):
    from . import utils  as ul

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = Path(__file__).parent / "R" / "pca.R"
    param_dic = {
        "file": files,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(["heatmap", "hm"], help="Heatmap plot for PSI")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="figure output path")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def heatmap(files, output, cluster, fmt, width, height, resolution):
    from . import utils  as ul

    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = Path(__file__).parent / "R" / "heatmap.R"

    param_dic = {
        "file": files,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "clusterinfo": cluster if cluster else "0",
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help="barplot ")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=tuple, default=(),
             help="input dpsi names")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("This flag is present then a dpsi file will divide "
                  "two part according to |dpsi| > 0 and |dpsi| < 0"))             
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def barplot(output, files, xlabel, dg, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "barplot.R"
    param_dic = {
        "file": files,
        "name": xlabel,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output,
        "dg": dg
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help = "install R packages")
@click.option('-r', '--requirement', is_flag=True, default=False,
                help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", cls=MultiOption, type=tuple, 
                default=(), help="install Genome wide annotation package")
@click.option('-cran', '--cran', cls=MultiOption, type=tuple, 
                default=(), help="install CRAN package")
@click.option('-bioc', '--bioconductor', cls=MultiOption, type=tuple, 
                default=(),help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
@click.option('-m', '--mirror',  is_flag=True, default=False,
                help="use tsinghua mirrors.")
def install(requirement, OrgDb, cran, bioconductor, java, mirror):
    from . import ChromHMM as ch
    from . import utils  as ul

    pdir = Path(__file__).parent
    rscript = pdir / "R" / "install.R"
    
    if java:
        print("install ChromHMM")
        ch.install(pdir)
        shutil.copyfile(pdir/"ChromHMM.jar", pdir/"ChromHMM/ChromHMM.jar")
    param_dic = {
        "requirement": requirement,
        "OrgDb": OrgDb if OrgDb else False,
        "CRAN": cran if cran else False, 
        "bioconductor": bioconductor if bioconductor else False, 
        "mirror": mirror
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])




@cli.command(["list", "ls"], help = "list OrgDb")
@click.option('-orgdb', '--OrgDb', "OrgDb", is_flag=True, help="list OrgDb")
@click.option('-rbpsp', '--RBPSp', "RBPSp", is_flag=True, help="RNA binding protein  ")
def list_(OrgDb, RBPSp):
    from .constant import OrgDb_dic, RBP_sp_dic
    if OrgDb:
        for k, v in OrgDb_dic.items():
            print(f"{k}: {v}")      
    if RBPSp:
        for k, v in RBP_sp_dic.items():
            print(f"{k}: {v}")
   

@cli.command(help = "generate ChromHMM mark file")
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-ct', '--cellType', cls=MultiOption, type=tuple, required=True, help="cell types")
@click.option('-bed', '--bed', cls=MultiOption, type=tuple, help="bed files")
@click.option('-mn', '--markNum', cls=MultiOption, type=tuple, required=True,
                 help="mark count of every cell types") 
@click.option('-sep', help="split symbol for splitting bed file names") 
@click.option('-mi', "--markIndex", type=int, 
            help="the mark index when bed files are splited by -sep symbol")
@click.option('--markName', cls=MultiOption, type=tuple, help="mark names")
@click.option('--stacked', is_flag=True, 
            help="This flag replaces the mark entry with an entry of the form cell_mark.")
def mark(output, celltype, bed, marknum, sep, markindex, markname, stacked):
    from itertools import chain
    import pandas as pd

    if all([sep, markindex]):
        marks = [Path(i).stem.split(sep)[markindex-1] for i in bed]
    else:
        marks = markname

    if len(marks) != len(bed):
        print("mark number dismatchs with bed files ")
        sys.exit()
    if len(celltype) == len(marknum):
        cells = [i for i in chain(*[[celltype[i]]*int(marknum[i]) for i in range(len(celltype))])]
    elif len(celltype) == 1:
        cells = [celltype[0] for _ in bed]
    else:
        print("ERROR: -ct/--cellType and -mn/--markNum is not dismatched")
        sys.exit()
     
    if len(cells) != len(marks):
        print("-mn/--markNum is wrong")
        sys.exit()

    if stacked:
        marks = [f"{cells[idx]}_{marks[idx]}" for idx in range(len(marks))]
        cells = [f"genome" for _ in range(len(marks))]
    df = pd.DataFrame({
        "cell": cells,
        "mark": marks,
        "bed": bed
    })

    df.to_csv(output, sep="\t", index=False, header=False)


@cli.command(help = "generate ChromHMM anchor file")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-idx', '--index', type=int, help="element index")
@click.option('-si', '--sideIndex', type=(int, int), help="the center of two side index")
@click.option('-u', '--upstreamOffset', "offset5", type=int, default=0, help="upstream offset")
@click.option('-d', '--downstreamOffset', "offset3", type=int, default=0, help="downstream offset")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
def anchor(file, output, index, sideindex, offset5, offset3, strand_sp):
    from . import utils  as ul

    try:
        ul.gen_anchor_bed(file, output, index, sideindex, offset5, offset3, strand_sp)
    except BaseException as e:
        print(e)
 
@cli.command(help = "generate bed file according to selected coordinates")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-s', '--start', type=int, help="start index")
@click.option('-e', '--end', type=int, help="end index")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
@click.option('-a', '--anchor', type=int, help="element index")
@click.option('-u', '--upstreamWidth', "upstream_w", type=int, default=150, help="width of right flank")
@click.option('-d', '--downstreamWidth', "downstream_w", type=int, default=150, help="width of left flank")
@click.option('-fi', 'fasta', type=click.Path(exists=True), 
            help="Input FASTA file. if set, the fasta sequence will be extracted")
def getcoor(file, output, start, end, strand_sp, anchor, upstream_w, downstream_w, fasta):
    from . import utils  as ul

    try:
        coor_df = ul.get_coor_bed(file, start, end, strand_sp, anchor, upstream_w, downstream_w)
        coor_df.to_csv(output, index=False, header=False, sep="\t")
        if fasta:
            ul.get_coor_fa(coor_df, fasta, Path(output).with_suffix(".fa"))
    except BaseException as e:
        print(e)


CHROMSIZES_dir = Path(__file__).parent / "ChromHMM/CHROMSIZES"
genomes = [i.stem for i in CHROMSIZES_dir.glob("*.txt")]

@cli.command(["learnState", "lcs"], 
    help = "learns a chromatin state model, Wrapper of ChromHMM BinarizeBed and LearnModel")
@click.option('-n', "--numstates", default=2, type=int, help="states number")
@click.option('-mark', '--markfile', type=click.Path(exists=True), required=True, help="cell mark file")
@click.option('-dir', '--directory', type=click.Path(exists=True), 
                default=".", help="work directory")
@click.option('-bd', '--binaryDir', required=True, help="binary output directory")
@click.option('-od', '--outdir', required=True,  help="result output directory")
@click.option('-b', '--binsize', type=int, default=150, help="binsize")
@click.option('-g', '--genome', required=True, type=click.Choice(genomes), help="genome assembly")
@click.option('-mx', default="2400M", help="memory for java running")
@click.option('-p', "--processor", default=2, type=int, help="maxprocessors")
@click.option('-anchor', "--anchorDir", type=click.Path(exists=True), help="anchor files directory")
@click.option('-coor', "--coorDir", type=click.Path(exists=True), help="coordinate files directory")
@click.option('-nb', "--no_binary", is_flag=True, default=False, help="not run BinarizeBed")
@click.option('--name', default="", help="file name")
@click.option('--stacked', is_flag=True, 
            help="This flag replaces the mark entry with an entry of the form cell_mark.")
@click.option('--nostrand', is_flag=True, 
            help="This flag is present then strand information is not used in the enrichment calculations.")
@click.option('-DC', '--defaultCoor', is_flag=True, 
            help="This flag is present then default is used in the enrichment calculations.")            
def LearnState(numstates, markfile, directory, binarydir , outdir, binsize, genome, 
                mx, processor, anchordir, coordir, no_binary, name, stacked, nostrand,defaultcoor):
    ChromHMM_dir = Path(__file__).parent / "ChromHMM"
    ChromHMM_jar = f"java -mx{mx} -jar {ChromHMM_dir/'ChromHMM.jar'}"
    
    uv_params = ""
    if coordir:
        p_coordir = Path(coordir)
        uv_params += f"-u {p_coordir.absolute()} "
        (p_coordir / f"{genome}").mkdir(exist_ok=True)
        if defaultcoor:
            shutil.copytree(ChromHMM_dir/f"COORDS/{genome}", p_coordir/f"{genome}", dirs_exist_ok=True)
        for file in p_coordir.glob("*"):
            if file.is_file(): shutil.copy(file, p_coordir/f"{genome}")
    if anchordir:
        p_anchordir = Path(anchordir)
        uv_params += f"-v {p_anchordir.absolute()}"
        (p_anchordir / f"{genome}").mkdir(exist_ok=True)
        if defaultcoor:
            shutil.copytree(ChromHMM_dir/f"ANCHORFILES/{genome}", p_anchordir/f"{genome}", dirs_exist_ok=True)
        for file in p_anchordir.glob("*"):
            if file.is_file(): shutil.copy(file, p_anchordir/f"{genome}")
    
    stacked_param = "-stacked" if stacked else ""
    chrom_len = ChromHMM_dir / f"CHROMSIZES/{genome}.txt"
    strand_param = "-nostrand" if nostrand else ""
    name_param = f"-i {name}"

    ChromHMM_bin = f"BinarizeBed {stacked_param} -peaks -b {binsize} {chrom_len} {directory} {markfile} {binarydir}"
    ChromHMM_lm = f"LearnModel {uv_params} {strand_param} {name_param} -p {processor} \
                    -b {binsize} {binarydir} {outdir} {numstates} {genome}"
    if not no_binary:
        print(ChromHMM_bin)
        info = subprocess.Popen([*ChromHMM_jar.split(), *ChromHMM_bin.split()])
        info.wait()
    print(ChromHMM_lm)
    info = subprocess.Popen([*ChromHMM_jar.split(), *ChromHMM_lm.split()])
    info.wait()

    rscript = Path(__file__).parent / "R" / "ChromHMM_hm.R"
    for of in Path(outdir).glob("*_overlap.txt"):
        info = subprocess.Popen(["Rscript", str(rscript), of])
    info.wait()


@cli.command(["motifEnrich", "me"], help = "Motif Enrichment")
@click.option('-tf', "--tfasta", cls=MultiOption, type=tuple, 
                required=True, help="input fasta files")
@click.option('-cf', "--cfasta", cls=MultiOption, type=tuple, 
                help="control fasta files")                
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                 default="CISBP-RNA", help="RBP motif database default=CISBP-RNA")
@click.option('-org', '--organism', help="RBP organism")
# @click.option('-mm', '--meme', type=click.Path(exists=True), 
#                 required=True, help="path to .meme format file")
def motif_enrich(tfasta, cfasta, outdir, database, organism):
    from .constant import RBP_sp_dic
    from . import utils  as ul

    Path(outdir).mkdir(exist_ok=True)
    rscript = Path(__file__).parent / "R" / "motifEnrich.R"

    if cfasta and len(tfasta) != len(cfasta):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()
    elif cfasta is None:
        cfasta = ["0" for _ in tfasta]

    sp = RBP_sp_dic.get(organism, "0")

    param_dic = {
        "tfile": tfasta,
        "cfile": cfasta,
        "outdir": outdir,
        "database": database,
        "organism": sp
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])

@cli.command(["motifFind", "mf"], help = "Motif Discovery and similarity comparision")
@click.option('-tf', "--tfasta", cls=MultiOption, type=tuple, 
                required=True, help="fasta file")
@click.option('-cf', "--cfasta", cls=MultiOption, type=tuple, 
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
def motif_find(tfasta, cfasta, outdir, database, organism, pvalue, evalue, minw, maxw):
    from .constant import RBP_sp_dic
    from . import utils  as ul

    Path(outdir).mkdir(exist_ok=True)
    rscript = Path(__file__).parent / "R" / "motifFind.R"

    sp = RBP_sp_dic.get(organism, "0")

    if cfasta and len(tfasta) != len(cfasta):
        print("-tf/tfasta number must be same as -cf/cfasta")
        sys.exit()
    elif cfasta is None:
        cfasta = ["0" for _ in tfasta]

    param_dic = {
        "outdir": outdir,
        "tfile": tfasta, 
        "cfile": cfasta, 
        "pvalue": pvalue,
        "minw": minw, 
        "maxw": maxw,
        "database": database,
        "organism": sp,
        "eval": evalue
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(["motifPlot", "mp"], help = "Motif plot")
@click.option('-mi', "--motifId", "motifid", cls=MultiOption, type=tuple, 
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
def motif_plot(motifid, database, organism, meme, output, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "motifPlot.R"
    sp = RBP_sp_dic.get(organism, None)
    if meme:
        meme_file = meme
    elif database and sp:
        meme_file = Path(__file__).parent  / f"data/motif/{database}/{sp}.meme"
    else:
        print("--meme or --db and -org muset be set")    
        sys.exit()
    param_dic = {
        "motifId": motifid,
        "meme": meme_file,
        "fmt": fmt,
        "width": width, 
        "height": height,
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help="Gene Set Enrichment Analysis")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-n', '--name', default="GSEA", help="output name prefix")
@click.option('-pval', '--pvalue', type=float, default=0.2, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO']), 
                default="GO", help="enrich database")
@click.option('-gt', '--geneId', type=click.Choice(['ENSEMBL', 'ENTREZID', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -org' to view more ")
@click.option('-ont', type=click.Choice(['BP', 'MF', 'CC']), 
                default="BP", help="one of 'BP', 'MF', and 'CC' subontologies.")  
@click.option('-org', '--keggOrganism', "organism", default = "",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")                            
def gsea(file, outdir, name, pvalue, database, geneid, orgdb, ont, organism):
    from . import utils  as ul

    Path(outdir).mkdir(exist_ok=True)
    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -org' to view more")

    rscript = Path(__file__).parent / "R" / "gsea.R"
    params = [outdir, str(pvalue), org_db, ont, geneid, name, database, organism, file]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


@cli.command(help="Gene Set Enrichment Analysis ploting")
@click.option('-id', '--id', "termid", cls=MultiOption, type=tuple, 
                required=True, help="term id")
@click.option('-o', '--output', help="output figure path")
@click.option('-rd', '--RData', help="output figure path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def gseplot(termid, output, rdata, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "gseaPlot.R"
    param_dic = {
        "termid": termid,
        "fmt": fmt, 
        "RData": rdata,
        "output": output,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}")
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])

@cli.command(help="draw UpSet plots for AS events")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)           
@click.option('-o', '--output', required=True, help="output figure path")
@click.option('-xl', '--xlabel', cls=MultiOption, type=tuple, 
                help="x xlabel")
@click.option('-dg', '--dg', is_flag=True, default = False,
              help=("This flag is present then a dpsi file will divide "
                  "two part according to |dpsi| > 0 and |dpsi| < 0"))                     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format") 
@click.option('-w', '--width', default=6, help="fig width, default=6 inches")
@click.option('--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def upset(files, output, xlabel, dg, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "upset.R"

    param_dic = {
        "file": files,
        "dg": dg,
        "fmt": fmt, 
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "output": Path(output).with_suffix(f".{fmt}"),
        "name": xlabel if xlabel else [str(i) for i  in range(len(files))]
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help = "get motif from meme file")
@click.argument('motifId', nargs=-1, required=True)   
@click.option('-mm', '--meme', type=click.Path(exists=True), help="meme motif file")
@click.option('-db', '--database', type=click.Choice(['ATtRACT', 'CISBP-RNA']),
                help="RBP motif database")
@click.option('-org', '--organism', help="RBP organism")
@click.option('-o', '--output', required=True, help="output path")
def getmeme(motifid, meme, database, organism, output):
    from . import utils  as ul

    sp = RBP_sp_dic.get(organism, None)
    if meme:
        meme_file = meme
    elif database and sp:
        meme_file = Path(__file__).parent  / f"data/motif/{database}/{sp}.meme"
    else:
        print("--meme or --db and -org muset be set")    
        sys.exit()

    try:
        ul.get_meme(motifid, meme_file, output)
    except BaseException as e:
        print(e)


@cli.command(help = "generate motif map")
@click.option('-fa', '--fasta', required=True, cls=MultiOption, type=tuple, 
                help="fasta files")
@click.option('-n', '--name', cls=MultiOption, type=tuple, default=(),
                help="fasta file names")
@click.option('-c', '--center', cls=MultiOption, type=tuple, default=(),
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
def mmap(fasta, name, center, meme, outdir, binsize, step, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "motifMap.R"
    Path(outdir).mkdir(exist_ok=True)

    param_dic = {
        "fasta": fasta,
        "outdir": outdir, 
        "meme": meme,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "bin": binsize,
        "step": step,
        "fmt": fmt,
        "center": center if len(center) ==len(fasta) else ["0"] * len(fasta),
        "seqid": name if len(name)==len(fasta) else [Path(i).stem for i in fasta]
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])
    

@cli.command(help="epi signal")
@click.option('-o', '--output', required=True, help="output path")
@click.option('-md', '--metadata', type=click.Path(exists=True),
             help="contrast group metadata, generated by meta")
@click.option('-anchor', '--anchor', cls=MultiOption, type=tuple, default=(),
                help="site anchor files")
@click.option('-n', '--name', cls=MultiOption, type=tuple, default=(),
                help="site anchor names")
@click.option('-w', '--width', default=150, type=int, help="flank window width, default=150")
@click.option('-bs', '--binSize', default=15, type=int,help="bin size, default=15")
def epi(output, metadata, anchor, name, width, binsize):
    from . import epi

    names = name if len(name)==len(anchor) else list(range(1, len(anchor)+1))

    anchor_dic = dict(zip(names, anchor))
    epi.epi_signal(output, anchor_dic, metadata, width, binsize)
    

@cli.command(help="epi signal heatmap")
@click.argument('files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def epihm(output, files, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "signalHeatmap.R"
    param_dic = {
        "file": files,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help="epi signal line")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=4, help="fig height, default=4 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def epiline(output, file, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "signalProfile.R"
    param_dic = {
        "file": file,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


@cli.command(help="epi signal compare")
@click.argument('file', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', required=True, help="output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf', 'pptx']),
                 default="png", help="out figure format")
@click.option('-w', '--width', default=8, help="fig width, default=8 inches")
@click.option('-h', '--height', default=6, help="fig height, default=6 inches")
@click.option('-res', '--resolution', default=72, help="resolution, default=72 ppi")
def sigcmp(output, file, fmt, width, height, resolution):
    from . import utils  as ul

    rscript = Path(__file__).parent / "R" / "signalCompare.R"
    param_dic = {
        "file": file,
        "width": width, 
        "height": height, 
        "resolution": resolution,
        "fmt": fmt,
        "output": output
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])


if __name__ == '__main__':
    cli()

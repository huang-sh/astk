from collections import defaultdict
import os
import sys
import json
import shutil
import subprocess
from pathlib import Path

import click
import pandas as pd 

from .cli_config import *
from . import utils  as ul
from . import ChromHMM as ch



@click.group(cls=CustomMultiCommand, 
        context_settings=dict(help_option_names=['-h', '--help']))
def cli():
     pass


@cli.command(help="generate metadata template file")
@click.option('-o', '--output', required=True, help='metadata output path')
@click.option('-g', '--group', type=int, required=True, help="contrast groups, number")
@click.option('-repN', '--replicate', cls=OptionEatAll,type=tuple, help="replicate, number")
@click.option('-gn', '--group_name', cls=OptionEatAll, type=tuple, help="group name")
@click.option('-p1', cls=OptionEatAll,type=tuple, help="file path for condtion 1")
@click.option('-p2', cls=OptionEatAll,type=tuple, help="file path for condtion 2")
@click.option('-repN1', '--replicate1', cls=OptionEatAll,type=tuple, help="replicate1, number")
@click.option('-repN2', '--replicate2', cls=OptionEatAll,type=tuple, help="replicate2, number")
def meta(output, group, replicate, group_name, p1, p2, replicate1, replicate2):
    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    output = Path(output).with_suffix(".csv")
    ul.meta_template(output, group, replicate, group_name,
                    p1, p2, repN1=replicate1, repN2=replicate2)


@cli.command(["lenCluster", "lc"], help="length cluster")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True),
                required=True,  help='AS ioe file')
@click.option('-o', '--output', required=True, help="output path")
@click.option('-cl', '--custom_len', 'custom_len', cls=OptionEatAll, type=tuple, help="custom length")
@click.option('-nc', '--cluster', type=int, default=4, help="number of cluster")
@click.option('-bw', '--width', type=int, default=3, help="bin width")
@click.option('-lw', '--len_weight', type=float, default=2, help="length weight")
@click.option('--max_len', type=int, default=500, help="the max length of exon in clustering")
def len_cluster(infile, output, custom_len, cluster, width, len_weight, max_len):
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


@cli.command(["lenPick", "lp"])
@click.option('-i', '--input', 'infile', type=click.Path(exists=True), help='AS ioe file')
@click.option('-o', '--output', help="output path")
@click.option('-rg', '--range', "len_range", type=(int, int), required=True, help="length range")
def len_pick(infile, output, len_range):
    if not (pdir:= Path(output).parent).exists():
        print(f"{pdir} doest not exist")
        exit()                   

    ioe_df = pd.read_csv(infile, sep="\t")
    info_df = ul.extract_info(ioe_df)
    info_df["event_id"] = ioe_df["event_id"]
    s, e = len_range
    pdf = info_df.loc[(s <= info_df["len"]) & ( info_df["len"] < e), :]
    pdf.to_csv(Path(output).parent / f"{Path(output).stem}_{s}-{e}", index=False)


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
@click.option('--tpm_col', type=int, default=4, help="TPM columns index")
def diff_splice(outdir, metadata, gtf, event_type, exon_len, tpm_col, method):
    if event_type == "all":
        event_types = AS_type
    else:
        event_types = [event_type]
    dsi = ul.DiffSplice(outdir, metadata, gtf, event_types, exon_len)
    dsi.ds(method)


@cli.command(["sigfilter", "sf"], help="filter significant result")
@click.option('-i', '--input', 'infile', type=click.Path(exists=True), help="dpsi file")
@click.option('-md', '--metadata', type=click.Path(exists=True), help="metadata file")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-dpsi', '--dpsi', type=float, default=0, help="dpsi threshold value")
@click.option('-p', '--pval', type=float, default=0.05, help="pval threshold value")
@click.option('-adpsi', '--abs_dpsi', type=float, default=0, help="absulte dpsi threshold value")
@click.option('-pf', '--psiFile', cls=OptionEatAll,type=tuple, help="psi files")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['csv', 'tsv']), 
                default="tsv", help="out  file format ")
def sigfilter(infile, metadata, outdir, dpsi, pval, abs_dpsi, psifile, fmt):
    if infile:
        sf = ul.SigFilter(infile, outdir, dpsi, pval, abs_dpsi, psifile, fmt)
        sf.run()
    if metadata:
        with open(metadata, "r") as f:
            meta_dic = json.load(f)
        
        for i in meta_dic:
            for at in meta_dic[i]["dpsi"].keys():
                infile = meta_dic[i]["dpsi"][at]
                psi_file1 = meta_dic[i]["control"]["psi"][at]
                psi_file2 = meta_dic[i]["treatment"]["psi"][at]
                psi_file = [psi_file1, psi_file2]
                sf = ul.SigFilter(infile, outdir, dpsi, pval, abs_dpsi, psi_file, fmt)
                sf.run()

@cli.command()
@click.option('-i', '--input', 'infiles',  cls=OptionEatAll, type=tuple, required=True, help="dpsi files")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG']), 
                default="GO", help="enrich database")
@click.option('-gene_id', '-gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'REFSEQ', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -org' to view more ")
@click.option('-org', '--keggOrganism', "kegg_organism",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")          
def enrich(infiles, outdir, pvalue, qvalue, database, gene_id, orgdb, kegg_organism):
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

    for file in infiles:
        out = Path(outdir) / Path(file).stem
        out.mkdir(exist_ok=True)
        params = [str(out), str(pvalue), str(qvalue), database,
                 gene_id, org_db , kegg_organism, file]
        info = subprocess.Popen(["Rscript", str(rscript), *params])

    else:
        info.wait()


@cli.command(["enrichCompare", "ecmp"])
@click.option('-i', '--input', 'infiles',  cls=OptionEatAll, 
                required=True, type=tuple, help="dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")            
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG']), 
                default="GO", help="enrich database")
@click.option('-pval', '--pvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-qval', '--qvalue', type=float, default=0.1, help="pvalue cutoff")
@click.option('-gene_id', '-gene_id', type=click.Choice(['ENSEMBL', 'ENTREZID', 'REFSEQ', 'SYMBOL']), 
                default="ENSEMBL", help="gene ID type")                      
@click.option('-orgdb', '--orgdb', required=True,
                help="OrgDb for GO annotation, such as: hs for Human, mm for Mouse. \
                    run 'astk ls -orgdb' to view more ")
@click.option('-org', '--keggOrganism', "kegg_organism",
                help="KEGG organism short alias.This is required if -db is KEGG.\
                    Organism list in http://www.genome.jp/kegg/catalog/org_list.html")   
def enrich_cmp(infiles, outdir, cluster, database,
                pvalue, qvalue, gene_id, orgdb, kegg_organism):
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
    params = [outdir, str(pvalue), str(qvalue), database,
              cluster, org_db, gene_id, kegg_organism, *infiles]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()
    print(info)


@cli.command(["enrichLenCluster", "elc"])
@click.option('-i', '--input', 'infiles',  cls=OptionEatAll, 
                required=True, type=tuple, help="dpsi files")
@click.option('-od', '--outdir', required=True, help="output directory")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                required=True, help="cluster information file")
@click.option('-m', '--merge', is_flag=True, help="enrich with multiple merged files")              
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG']), 
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
def enrich_lc(infiles, outdir, cluster, merge, database, pvalue, qvalue,
              gene_id, orgdb, kegg_organism):
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
    for file in infiles:
        params = [outdir, str(pvalue), str(qvalue), database, cluster,
                 org_db, gene_id, kegg_organism, file]
        info = subprocess.Popen(["Rscript", str(rscript), *params])
        if database == "KEGG" and not ul.check_kegg_RData(kegg_organism):
             info.wait() 
    else:
        if merge:
            params = [outdir, str(pvalue), str(qvalue), database, cluster,
                 org_db, gene_id, kegg_organism, *infiles]
            merge_info = subprocess.Popen(["Rscript", str(rscript), *params])
            merge_info.wait()
        else:
            info.wait()


@cli.command(["volcano", "vol"], help="Volcano plot analysis for dPSI")
@click.option('-i', '--input', "dpsi", cls=OptionEatAll, type=tuple, help="dpsi files")
@click.option('-od', '--outdir', type=click.Path(exists=True), help="output directory")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf']),
                 default="png", help="out figure format ")
def volcano(dpsi, outdir, fmt):
    rscript = Path(__file__).parent / "R" / "volcano.R"
    for f in dpsi:
        out = Path(outdir) / Path(f).with_suffix(f".volcano.{fmt}").name
        info = subprocess.Popen(["Rscript", str(rscript), f, out])
        info.wait()


@cli.command(help="PCA analysis for PSI")
@click.option('-i', '--input', 'infiles',  cls=OptionEatAll, 
                required=True, type=tuple, help="psi files")
@click.option('-o', '--output', 'outpath', required=True, help="figure output path")
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf']),
                 default="png", help="out figure format ")
def pca(infiles, outpath, fmt):
    if not (pdir:= Path(outpath).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = Path(__file__).parent / "R" / "pca.R"
    outpath = Path(outpath).with_suffix(f".{fmt}")
    params = [str(outpath), *infiles]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


@cli.command(["heatmap", "hm"], help="Heatmap plot for PSI")
@click.option('-i', '--input', 'infiles',  cls=OptionEatAll, 
                required=True, type=tuple, help="psi files")
@click.option('-o', '--output', 'outpath', required=True, help="figure output path")
@click.option('-cls', '--cluster', type=click.Path(exists=True),
                help="cluster information file")     
@click.option('-fmt', '--format', "fmt", type=click.Choice(['png', 'pdf']),
                 default="png", help="out figure format ")
def heatmap(infiles, outpath, cluster, fmt):
    if not (pdir:= Path(outpath).parent).exists():
        print(f"{pdir} doest not exist")
        exit()
    rscript = Path(__file__).parent / "R" / "heatmap.R"
    cluster = cluster if cluster else "0"
    outpath = Path(outpath).with_suffix(f".{fmt}")
    params = [str(outpath), cluster, *infiles]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


@cli.command(help = "install R packages")
@click.option('-r', '--requirement', is_flag=True, help="install astk requirement R packages")
@click.option('-OrgDb', '--OrgDb', "OrgDb", help="install Genome wide annotation package")
@click.option('-cran', '--cran', help="install CRAN package")
@click.option('-bio', '--bioconductor', help="install Bioconductor package")
@click.option('-j', '--java',  is_flag=True, help="install java software")
def install(requirement, OrgDb, cran, bioconductor, java):
    pdir = Path(__file__).parent
    rscript = pdir / "R" / "install.R"
    rq = "1" if requirement else "0"
    OrgDb = OrgDb if OrgDb else "0"
    cran = cran if cran else "0"
    bioconductor = bioconductor if bioconductor else "0"
    params = [rq, OrgDb, cran, bioconductor]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()
    if java:
        print("install ChromHMM")
        ch.install(pdir)
        shutil.copyfile(pdir/"ChromHMM.jar", pdir/"ChromHMM/ChromHMM.jar")


@cli.command(["list", "ls"], help = "list OrgDb")
@click.option('-orgdb', '--OrgDb', "OrgDb",
                is_flag=True, help="list OrgDb")
def list(OrgDb):
    if OrgDb:
        for k, v in ul.OrgDb_dic.items():
            print(f"{k}: {v}")

@cli.command(help = "generate ChromHMM mark file")
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-ct', '--cellType', cls=OptionEatAll, type=tuple, required=True, help="cell types")
@click.option('-bed', '--bed', cls=OptionEatAll, type=tuple, help="bed files")
@click.option('-mn', '--markNum', cls=OptionEatAll, type=tuple, required=True,
                 help="mark count of every cell types") 
@click.option('-sep', help="split symbol for splitting bed file names") 
@click.option('-mi', "--markIndex", type=int, 
            help="the mark index when bed files are splited by -sep symbol")
@click.option('--markName', cls=OptionEatAll, type=tuple, help="mark names")
@click.option('--stacked', is_flag=True, 
            help="This flag replaces the mark entry with an entry of the form cell_mark.")
def mark(output, celltype, bed, marknum, sep, markindex, markname, stacked):
    from itertools import chain

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
@click.option('-i', '--input', 'infile', required=True, type=click.Path(exists=True), help="dpsi file")
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-idx', '--index', type=int, help="element index")
@click.option('-si', '--sideIndex', type=(int, int), help="the center of two side index")
@click.option('-u', '--upstreamOffset', "offset5", type=int, default=0, help="upstream offset")
@click.option('-d', '--downstreamOffset', "offset3", type=int, default=0, help="downstream offset")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
def anchor(infile, output, index, sideindex, offset5, offset3, strand_sp):
    try:
        ul.gen_anchor_bed(infile, output, index, sideindex, offset5, offset3, strand_sp)
    except BaseException as e:
        print(e)
 
@cli.command(help = "generate bed file according to selected coordinates")
@click.option('-i', '--input', 'infile', required=True, type=click.Path(exists=True), help="dpsi file")
@click.option('-o', '--output', required=True, help="file output path")
@click.option('-s', '--start', type=int, help="start index")
@click.option('-e', '--end', type=int, help="end index")
@click.option('-ss', '--strandSpecifc', "strand_sp", is_flag=True, help="strand specifc")
@click.option('-anchor', '--anchor', type=int, help="element index")
@click.option('-u', '--upstreamWidth', "upstream_w", type=int, default=150, help="width of right flank")
@click.option('-d', '--downstreamWidth', "downstream_w", type=int, default=150, help="width of left flank")
def getcoor(infile, output, start, end, strand_sp, anchor, upstream_w, downstream_w):
    try:
        ul.get_coor_bed(infile, output, start, end, strand_sp, anchor, upstream_w, downstream_w)
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


@cli.command(["motifEnrich", "me"], help = "motif enrichment")
@click.option('-fa', "--fasta", cls=OptionEatAll, type=tuple, required=True, help="fasta files")
@click.option('-od', '--outdir', type=click.Path(), default=".", help="output directory")
@click.option('-mm', '--meme', type=click.Path(exists=True), required=True, help="path to .meme format file")
def motif_enrich(fasta, outdir, meme):
    Path(outdir).mkdir(exist_ok=True)

    rscript = Path(__file__).parent / "R" / "motifEnrich.R"
    params = [outdir, meme, *fasta]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


@cli.command(help="Gene Set Enrichment Analysis")
@click.option('-i', '--input', 'infile',  type=click.Path(exists=True), 
                required=True, help="input dpsi file")
@click.option('-od', '--outdir', default=".", help="outdir")
@click.option('-n', '--name', default="GSEA", help="output name prefix")
@click.option('-pval', '--pvalue', type=float, default=0.2, help="pvalue cutoff")
@click.option('-db', '--database', type=click.Choice(['GO', 'KEGG']), 
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
def gsea(infile, outdir, name, pvalue, database, geneid, orgdb, ont, organism):
    Path(outdir).mkdir(exist_ok=True)
    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -org' to view more")

    rscript = Path(__file__).parent / "R" / "gsea.R"
    params = [outdir, str(pvalue), org_db, ont, geneid, name, database, organism, infile]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


@cli.command(help="Gene Set Enrichment Analysis ploting")
@click.option('-id', '--id', "termid", cls=OptionEatAll, type=tuple, 
                required=True, help="term id")
@click.option('-o', '--output', help="output figure path")
@click.option('-rd', '--RData', help="output figure path")     
def gseplot(termid, output, rdata):
    rscript = Path(__file__).parent / "R" / "gsea_plot.R"
    params = [output, rdata, *termid]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()


if __name__ == '__main__':
    cli()
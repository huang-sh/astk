import subprocess
from pathlib import Path

from astk.constant import BASE_DIR, PATHWAY_DB_ORG
import astk.utils.func  as ul
from astk.types import FilePath
from astk.event import SuppaEventID
                
def gsea_fun(file, outdir, name, pvalue, database, geneid, orgdb, ont, organism):
    Path(outdir).mkdir(exist_ok=True)
    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -org' to view more")

    rscript = BASE_DIR / "R" / "gsea.R"
    params = [outdir, str(pvalue), org_db, ont, geneid, name, database, organism, file]
    info = subprocess.Popen(["Rscript", str(rscript), *params])
    info.wait()

   
def enrich(
    file: FilePath,
    outdir: FilePath, 
    pvalue: float, 
    qvalue: float, 
    database: str, 
    ontology: str,  
    gene_id: str, 
    organism: str, 
    fmt: str, 
    width: float, 
    height: float
) -> None:
    rscript = BASE_DIR / "R" / "enrich.R"
    if not (org_db := ul.select_OrgDb(organism)):
        print(f"{organism} is wrong! Please run 'astk ls -org' to view more")
    if database in ["KEGG", "Reactome"]:
        organism = PATHWAY_DB_ORG[database].get(organism, None)
    Path(outdir).mkdir(exist_ok=True)
    param_dic = {
        "file": file,
        "outdir": outdir, 
        "orgdb": org_db,
        "database": database,
        "pval": pvalue,
        "qval": qvalue,
        "organism": organism,
        "genetype": gene_id,
        "ontology": ontology,
        "width": width,
        "height": height,
        "fmt": fmt
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run(["Rscript", rscript, *param_ls])

             
def enrich_cmp(files, outdir, cluster, database, ontology, pvalue, qvalue, 
                xlabel, gene_id, orgdb, kegg_organism, fmt, width, height):

    if not (org_db := ul.select_OrgDb(orgdb)):
        print(f"{orgdb} is wrong! Please run 'astk ls -orgdb' to view more")
    if database == "KEGG":
        ul.check_kegg_RData(kegg_organism)
        if not kegg_organism:
            print("Error: --kegg_organism is required!")
            exit()
    else:
        kegg_organism = "0"
    rscript = BASE_DIR / "R" / "enrichCompare.R"
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


def enrich_lc(files, outdir, cluster, merge, database, pvalue, qvalue,
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
    rscript = BASE_DIR / "R" / "enrichLenCluster.R"
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


def nease_enrich(nease_input, outdir, n=15, database=['Reactome'], organism='Human', cutoff=0.05):
    import nease
    import numpy as np
    import matplotlib.pyplot as plt

    events = nease.run(nease_input, organism=organism, p_value_cutoff=cutoff)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    stats_file = outdir / "affected_stats_fig.pdf"
    enrich_csv = outdir / "nease_enrichment.csv"
    affected_domain_csv = outdir / "affected_domain.csv"
    affected_elm_csv = outdir / "affected_elm.csv"
    affected_pdb_csv = outdir / "affected_pdb.csv"
    enrich_bar_pdf = outdir / "nease_enrichment_bar.pdf"
    
    events.get_stats(file_path=stats_file)
    events.get_domains().to_csv(affected_domain_csv)
    events.get_elm().to_csv(affected_elm_csv)
    events.get_pdb().to_csv(affected_pdb_csv)
    
    nease_enr = events.enrich(database=database)
    if nease_enr is not None:
        nease_enr.to_csv(enrich_csv)
        top_terms = min([nease_enr.shape[0], n])
        nease_enr = nease_enr.sort_values(by='adj p_value')
        Term = nease_enr['Pathway name'][:top_terms]
        Pvalues = nease_enr['adj p_value'][:top_terms]
        Pvalues = [ -np.log10(x) for x in Pvalues]
        Term = [x.split('Homo')[0] for x in Term]

        plt.figure()
        plt.barh(Term[::-1],Pvalues[::-1], color="#83cdfd")
        plt.title('NEASE enrichment')
        plt.ylabel('Terms')
        plt.xlabel('-log10(Adjusted P-value)')
        plt.savefig(enrich_bar_pdf, format='pdf',bbox_inches='tight')


def nease_sc(file, outdir, pvalue, database, organism):
    import pandas as pd

    get_geneid = lambda x: SuppaEventID(x).gene_id.split(".")[0]
    get_start = lambda x: SuppaEventID(x).alter_element_coor[0]
    get_end = lambda x: SuppaEventID(x).alter_element_coor[1]

    df = pd.read_csv(file, sep="\t", index_col=0)
    df["event_id"] = df.index
    
    nease_input = pd.DataFrame({
                    "Gene stable ID": df["event_id"].apply(get_geneid),
                    "new_start": df["event_id"].apply(get_start), 
                    "new_end": df["event_id"].apply(get_end),
                    "beta": df[df.columns[0]].values})
    nease_enrich(nease_input, outdir, database=database, organism=organism, cutoff=pvalue)


def neasecmp_sc():
    pass

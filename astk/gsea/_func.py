import subprocess
from pathlib import Path
from typing import Sequence

from astk.constant import BASE_DIR, PATHWAY_DB_ORG
import astk.utils.func  as ul
from astk.ctypes import FilePath
from astk.event import SuppaEventID


def gsea_fun(
    file: FilePath, 
    outdir: FilePath, 
    name: str, 
    pvalue: float, 
    database: str, 
    gene_id: str, 
    ontology: str, 
    organism: str,
    app: str
):
    Path(outdir).mkdir(exist_ok=True)
    if not (org_db := ul.select_OrgDb(organism)):
        print(f"{organism} is wrong! Please run 'astk ls -org' to view more")
    if database in {"KEGG", "Reactome"}:
        organism = PATHWAY_DB_ORG[database].get(organism, None)
    rscript = BASE_DIR / "R" / "gsea.R"
    params = [outdir, str(pvalue), org_db, ontology, 
                gene_id, name, database, organism, file, app]
    info = subprocess.Popen([ul.Rscript_bin(), str(rscript), *params])
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
    simple: bool,
    figfmt: str, 
    width: float, 
    height: float,
    app: str
) -> None:
    rscript = BASE_DIR / "R" / "enrich.R"
    if not (org_db := ul.select_OrgDb(organism)):
        print(f"{organism} is wrong! Please run 'astk ls -org' to view more")
    if database in {"KEGG", "Reactome"}:
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
        "simple": simple,
        "width": width,
        "height": height,
        "fmt": figfmt,
        "app": app
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def enrich_cmp(
    files: Sequence[FilePath],
    outdir: FilePath, 
    database: str, 
    ontology: str,
    pvalue: float,
    qvalue: float, 
    xlabel: Sequence[str],
    gene_id: str,
    organism: str,
    figfmt: str,
    width: float,
    height: float,
    app: str
) -> None:

    if not (org_db := ul.select_OrgDb(organism)):
        print(f"{organism} is wrong! Please run 'astk ls -orgdb' to view more")
    if database in {"KEGG", "Reactome"}:
        organism = PATHWAY_DB_ORG[database].get(organism, None)
    rscript = BASE_DIR / "R" / "enrichCompare.R"
    Path(outdir).mkdir(exist_ok=True)
    param_dic = {
        "files": files,
        "outdir": outdir, 
        "orgdb": org_db,
        "database": database,
        "pval": pvalue,
        "qval": qvalue,
        "xlabel": xlabel, 
        "organism": organism,
        "genetype": gene_id,
        "ontology": ontology,
        "width": width,
        "height": height,
        "fmt": figfmt,
        "app": app
    }
    param_ls = ul.parse_cmd_r(**param_dic)
    subprocess.run([ul.Rscript_bin(), rscript, *param_ls])


def nease_enrich(nease_input, outdir, n=15, database=None, organism='Human', cutoff=0.05):
    if database is None:
        database = ['Reactome']
    import nease
    import numpy as np
    import matplotlib.pyplot as plt

    events = nease.run(nease_input, organism=organism, p_value_cutoff=cutoff)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    enrich_csv = outdir / "nease_enrichment.csv"
    enrich_bar_pdf = outdir / "nease_enrichment_bar.pdf"
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
    nease_input = ul.shift2nease(file)
    nease_enrich(nease_input, outdir, database=database, organism=organism, cutoff=pvalue)


def neasecmp(files, outdir, num, qvalue, database, organism, xlabel, figformat, width, height):
    import nease
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    def _nease_func(file):
        df = ul.shift2nease(file)
        nr = nease.run(df, organism=organism, p_value_cutoff=0.2)
        res = nr.enrich(database=database)
        return res
    res_ls = list(map(_nease_func, files))
    pathways = []
    top_res_ls = []
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    for label, res in zip(xlabel, res_ls):
        if res is None:
            continue
        output = outdir / f"{label}_enrichment.csv"
        res.to_csv(output)
        res["cluster"] = label
        res.index = res["Pathway ID"]
        res.drop(index=res[res["adj p_value"]>qvalue].index, inplace=True)
        pathways.append(set(res["Pathway ID"]))
        top_res_ls.append(res.head(num))
    # share_pathways = set.intersection(*pathways)
    # df_share = pd.concat([df.loc[share_pathways, :] for df in res_ls if df is not None])
    df_top = pd.concat(top_res_ls)
    # dfl = pd.concat([df_top, df_share])
    dfl = df_top
    dfl["Pathway name"] = dfl["Pathway name"].apply(lambda x: x.split("- Homo")[0])
    if all([width, height]):
        fig, ax = plt.subplots(figsize=(width, height))
    else:
        ax = None
    plt.grid(linestyle='dotted', axis='y', color='#808080', dashes=(1, 6), zorder=1)
    palette = sns.color_palette("crest", as_cmap=True)
    ax = sns.scatterplot(data=dfl, x="cluster", y="Pathway name", palette=palette, zorder=2,
                         ax=ax, size="Nease score", hue="adj p_value", sizes=(50, 250))
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax.set_xlim(-0.5, len(top_res_ls)-0.5)    
    plt.savefig(outdir / f"cmp_enrichment.{figformat}", bbox_inches='tight')

# -*- coding: utf-8 -*-

"""
astk.enrich
~~~~~~~~~~~~~~~~~
This module provides a gene sets function enrichment enrich analysis.
"""

from pathlib import Path

import nease
import numpy as np
import matplotlib.pyplot as plt


def nease_enrich(nease_input, outdir, n=15, database=['Reactome'], organism='Human', cutoff=0.05):
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

    nease_enr = events.enrich(database=database)
    nease_enr.to_csv(enrich_csv)

    events.get_domains().to_csv(affected_domain_csv)
    events.get_elm().to_csv(affected_elm_csv)
    events.get_pdb().to_csv(affected_pdb_csv)


    top_terms = min([nease_enr.shape[0], n])
    nease_enr = nease_enr.sort_values(by='adj p_value')
    Term = nease_enr['Pathway name'][:top_terms]
    Pvalues = nease_enr['adj p_value'][:top_terms]
    Pvalues = [ -np.log10(x) for x in Pvalues]
    Term = [x.split('Homo')[0] for x in Term]

    plt.figure()
    plt.barh(Term[::-1],Pvalues[::-1] )
    plt.title('NEASE enrichment')
    plt.ylabel('Terms')
    plt.xlabel('-log10(Adjusted P-value)')
    plt.savefig(enrich_bar_pdf, format='pdf',bbox_inches='tight')

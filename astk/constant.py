# -*- coding: utf-8 -*-

"""
astk.constant
~~~~~~~~~~~~~~~~~
This module stores constants.
"""
from pathlib import Path


OrgDb_dic = {
     "hs": "org.Hs.eg.db", "mm": "org.Mm.eg.db", "rn": "org.Rn.eg.db",
     "dm": "org.Dm.eg.db", "at": "org.At.tair.db", "sc": "org.Sc.sgd.db",
     "dr": "org.Dr.eg.db", "ce": "org.Ce.eg.db", "bt": "org.Bt.eg.db",
     "ss": "org.Ss.eg.db", "mmu": "org.Mmu.eg.db", "gg": "org.Gg.eg.db", 
     "cf": "org.Cf.eg.db", "eck12": "org.EcK12.eg.db", "xl": "org.Xl.eg.db",
     "pt": "org.Pt.eg.db", "ag": "org.Ag.eg.db", "pf": "org.Pf.plasmo.db",
     "ecsakai": "org.EcSakai.eg.db", "mxanthus": "org.Mxanthus.db"
     }


RBP_sp_dic = {
    "cg": 'Cricetulus_griseus', "ng": 'Naegleria_gruberi', "tp": 'Thalassiosira_pseudonana',
    "mm": 'Mus_musculus', "oc": 'Oryctolagus_cuniculus', "bd": 'Brachypodium_distachyon', 
    "zm": 'Zea_mays', "sm": 'Schistosoma_mansoni', "tn": 'Tetraodon_nigroviridis', 
    "pp": 'Physcomitrella_patens', "hs" :'Homo_sapiens', "ma": 'Mesocricetus_auratus', 
    "nc": 'Neurospora_crassa', "lm": 'Leishmania_major', "ce": 'Caenorhabditis_elegans', 
    "xt": 'Xenopus_tropicalis', "dm": 'Drosophila_melanogaster', "pr": 'Phytophthora_ramorum', 
    "tb": 'Trypanosoma_brucei', "an": 'Aspergillus_nidulans',  "dr": 'Danio_rerio', 
    "tv": 'Trichomonas_vaginalis', "gg": 'Gallus_gallus', "ro": 'Rhizopus_oryzae', 
    "vp": 'Vanderwaltozyma_polyspora', "at": 'Arabidopsis_thaliana', "bm": 'Bombyx_mori', 
    "sf" :'Spodoptera_frugiperda', "nv": 'Nematostella_vectensis', "ot": 'Ostreococcus_tauri', 
    "bt" :'Bos_taurus', "sc": 'Saccharomyces_cerevisiae', "ct": 'Chaetomium_thermophilum', 
    "pf": 'Plasmodium_falciparum', "scs": 'Saccharomyces_cerevisiae_s288c', "rn": 'Rattus_norvegicus', 
    "xl": 'Xenopus_laevis', "ol": 'Oryzias_latipes'
    }

NEASE_DATABASE = [
     'PharmGKB', 'HumanCyc', 'Wikipathways', 'Reactome', 'KEGG','SMPDB',
     'Signalink', 'NetPath', 'EHMN', 'INOH', 'BioCarta', 'PID'
]


AS_TYPE = ['SE', "A5", "A3", "MX", "RI", 'AF', 'AL']

BASE_DIR = Path(__file__).parent



GTF_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

SSN = {
    "SE": 4,
    "MX": 6,
    "A5": 3,
    "A3": 3,
    "RI": 4,   # In fact, it only has two
    "AF": 5,   # In fact, it only has four
    "AL": 5    # In fact, it only has four
}

PATHWAY_DB_ORG = {
    "KEGG": {"hs": "hsa", "mm": "mmu"},
    "Reactome": {"hs": "human", "mm": "mouse"}
}

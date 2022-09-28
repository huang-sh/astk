suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--pval", type="double", help="pval")
parser$add_argument("--qval", type="double", help="qval")
parser$add_argument("--database", help="database")
parser$add_argument("--ontology", help="ontology")
parser$add_argument("--orgdb", help="orgdb")
parser$add_argument("--genetype", help="gene type")
parser$add_argument("--organism", help="kegg organism")
parser$add_argument("--file", help="dpsi files")
parser$add_argument("--simple", action='store_true', help="simple")

args <- parser$parse_args()

out.dir <- args$outdir
pval  <- args$pval
qval  <- args$qval
db <- args$database
ontology <- args$ontology
gene_type <- args$genetype
org_db <- args$orgdb
organism <- args$organism
dpsi_file <- args$file

suppressMessages(library(org_db, character.only = T))


dpsi <- read_tsv(dpsi_file, 
                skip      = 1, 
                col_names = F,
                col_types = cols("c", "d", "d"))

colnames(dpsi) = c("event_id", "dpsi", "pval")

genes <- gsub("\\..*", "",  dpsi$event_id) 

if (db == "GO"){
    p <- enrichGOSep(genes, 
                out.dir,
                org_db, 
                name    = "GO",
                ont     = ontology,
                keyType = gene_type,
                pval    = pval, 
                qval    = qval,
                simple  = args$simple,
                width   = args$width, 
                height  = args$height,
                format  = args$fmt)
} else if (db == "KEGG") {
    enrichKEGGSep(genes, 
                  out.dir, 
                  org_db, 
                  name    = "KEGG",
                  org     = organism, 
                  keytype = gene_type,
                  pval    = pval, 
                  qval    = qval,
                  width   = args$width, 
                  height  = args$height,
                  format  = args$fmt)
} else if (db == "Reactome") {
    enrichReactomeSep(genes, 
                  out.dir, 
                  org_db, 
                  name    = "Reactome",
                  org     = organism, 
                  keytype = gene_type,
                  pval    = pval, 
                  qval    = qval,
                  width   = args$width, 
                  height  = args$height,
                  format  = args$fmt)
}


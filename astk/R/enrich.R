suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))


args <- commandArgs()

out.dir <- args[6]
pval  <- as.numeric(args[7])
qval  <- as.numeric(args[8])
db <- args[9]
gene_type <- args[10]
org_db <- args[11]
kegg_org <- args[12]
dpsi_file <- args[13]

suppressMessages(library(org_db, character.only = T))

script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


dpsi <- read_tsv(dpsi_file, 
                skip      = 1, 
                col_names = F,
                col_types = cols("c", "d", "d"))

colnames(dpsi) = c("event_id", "dpsi", "pval")

genes <- gsub("\\..*", "",  dpsi$event_id) 

if (db == "GO"){
    enrichGOSep(genes, 
                out.dir,
                org_db, 
                name    = "GO",
                keyType = gene_type,
                pval    = pval, 
                qval    = qval)
} else if (db == "KEGG") {
    enrichKEGGSep(genes, 
                  out.dir, 
                  org_db, 
                  name    = "KEGG",
                  org     = kegg_org, 
                  keytype = gene_type,
                  pval    = pval, 
                  qval    = qval)
}

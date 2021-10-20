suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))


args <- commandArgs()

out.dir <- args[6]
pval  <- as.numeric(args[7])
qval  <- as.numeric(args[8])
db <- args[9]
cluster <- args[10] 
org_db <- args[11]
gene_id <- args[12]
kegg_organism <- args[13]
dpsi_files <- args[14:length(args)]

suppressMessages(library(org_db, character.only = T))

script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


cls_ext <- tools::file_ext(cluster)

if (cls_ext == "tsv"){
    exon_cluster <- read_tsv(cluster)
}else if (cls_ext == "csv") {
    exon_cluster <- read_csv(cluster)
}

exon_cluster$geneId <- gsub("\\..*", "",  exon_cluster$geneId) 

unique_exon_cls <- exon_cluster %>% 
    distinct(cluster, geneId) 

if (length(dpsi_files) == 1){
    file  <-  dpsi_files[1]

    dpsi <- read_tsv(file, skip =1, col_names = F,
                           col_types = cols("c", "d", "d"))
    fn  <-  tools::file_path_sans_ext(basename(file))

    genes <- gsub("\\..*", "",  dpsi$X1)                             
    
    sig_genes <- unique_exon_cls %>% 
        filter(geneId %in% genes)

    clusters <- sort(unique(sig_genes$cluster))
    gene_ls <- lapply(clusters, function(x){
        cg <- sig_genes %>% 
            filter(cluster == x)
        cg$geneId
    })
    names(gene_ls) <- paste("cluster", clusters, sep = "")
    if (db == "GO"){
        compareClusterSep(gene_ls, out.dir, org_db, name=fn,
                            keyType = gene_id,
                            pval = pval, qval = qval)
    }else if (db == "KEGG") {
        compareKEGGCluster(gene_ls, out.dir, org_db, 
                        kegg_organism, name=paste("KEGG", fn, sep="_"), 
                        keyType = gene_id, pval = pval, qval = qval)
    }


} else if (length(dpsi_files) > 1){
    dpsi_ls <- lapply(dpsi_files, function(file){
        dpsi <- read_tsv(file, skip =1, col_names = F,
                               col_types = cols("c", "d", "d"))
        colnames(dpsi) = c("event_id", "dpsi", "pval")
        return(dpsi)
    })
    dpsi_merge <- Reduce(rbind, dpsi_ls)
    genes <-gsub("\\..*", "",  dpsi_merge$event_id) 

    sig_genes <- unique_exon_cls %>% 
        filter(geneId %in% genes)

    clusters <- sort(unique(sig_genes$cluster))
    gene_ls <- lapply(clusters, function(x){
        cg <- sig_genes %>% 
            filter(cluster == x)
        cg$geneId
    })

    names(gene_ls) <- paste("cluster", clusters, sep = "")

    if (db == "GO"){
        compareClusterSep(gene_ls, out.dir, org_db, name="merge",
                        keyType = gene_id,
                        pval = pval, qval = qval)
    }else if (db == "KEGG") {
        compareKEGGCluster(gene_ls, out.dir, org_db, 
                        kegg_organism, name="KEGG_merge", 
                        keyType = gene_id, pval = pval, qval = qval)
    }
}

print("Done")
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))


args <- commandArgs()

out.dir <- args[6]
pval  <- as.numeric(args[7])
qval  <- as.numeric(args[8])
db <- args[9]
cluster_file <- args[10]
orgdb <- args[11]
gene_type <- args[12]
kegg_organism <- args[13]

dpsi_files <- args[14:length(args)]

suppressMessages(library(orgdb, character.only = T))

script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))



filenames <- unlist(lapply(dpsi_files, function(file){
                            tools::file_path_sans_ext(basename(file))
                        }))   

if (file.exists(cluster_file)){

    cls_ext <- tools::file_ext(cluster_file)

    if (cls_ext == "tsv"){
        exon_cluster <- read_tsv(cluster_file)
    }else if (cls_ext == "csv") {
        exon_cluster <- read_csv(cluster_file)
    }

    exon_cluster$geneId <- gsub("\\..*", "",  exon_cluster$geneId) 

    unique_exon_cls <- exon_cluster %>% 
        distinct(cluster, geneId) 


    clusters <- sort(unique(unique_exon_cls$cluster))
    lapply(clusters, function(c){
        cluster_gene_ls <- lapply(dpsi_files, function(file){
            dpsi <- read_tsv(file, skip =1, col_names = F,
                                col_types = cols("c", "d", "d"))
            colnames(dpsi) = c("event_id", "dpsi", "pval")

            genes <- unique(gsub("\\..*", "",  dpsi$event_id)) 
            cluster_sig_genes <- unique_exon_cls %>% 
                filter(geneId %in% genes) %>% 
                filter(cluster == c)
            cluster_sig_genes$geneId
        })      

        names(cluster_gene_ls) <- filenames  

        tryCatch({
            if (db == "GO"){
                print(gene_id)
                compareClusterSep(cluster_gene_ls, out.dir, orgdb,
                                    name = paste("cluster", c, sep=""),
                                    keyType = gene_type,
                                    pval = pval, qval = qval)                
            }else if (db == "KEGG") {
               compareKEGGCluster(cluster_gene_ls, out.dir, orgdb, 
                               kegg_organism, name=paste("KEGG_cluster", c, sep=""), 
                               keyType = gene_type, pval = pval, qval = qval)
            }

            }, error = function(e) {
                print(e)
                next()
            }
          )
    })

} else {
    gene_ls <- lapply(dpsi_files, function(file){
        dpsi <- read_tsv(file, skip =1, col_names = F,
                            col_types = cols("c", "d", "d"))
        genes <- unique(gsub("\\..*", "",  dpsi$X1))   
        return(genes)
   })
    names(gene_ls) <- filenames
    if (db == "GO"){
        compareClusterSep(gene_ls, out.dir, orgdb,
                            name="GO", keyType = gene_type,
                            pval = pval, qval = qval)
    }else if (db == "KEGG") {
        compareKEGGCluster(gene_ls, out.dir, orgdb, kegg_organism,
                        name="KEGG", keyType = gene_type,
                        pval = pval, qval = qval)
    }
    

}






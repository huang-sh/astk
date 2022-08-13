suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))


args <- commandArgs()

out.dir <- args[6]
pval  <- as.numeric(args[7])
org_db <- args[8]
ont <- args[9]
gene_type <- args[10]
name <- args[11]
database <- args[12]
organism <- args[13]
dpsi_file <- args[14]

script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

OrgDb <- load_OrgDb(org_db)

dpsi_df <- read_tsv(dpsi_file, 
                    skip      = 1, 
                    col_names = F,
                    col_types = cols("c", "d", "d")) %>% 
                    drop_na()

colnames(dpsi_df) = c("event_id", "dpsi", "pval")

if (gene_type == "ENSEMBL"){
    dpsi_df$gene <- gsub("\\..*", "",  dpsi_df$event_id) 
} else {
    dpsi_df$gene <- gsub(";.*", "",  dpsi_df$event_id)
}


dpis_gene <- sapply(unique(dpsi_df$gene), function(x){
    scores <- dpsi_df[dpsi_df$gene == x, ]$dpsi
    pos_score <- sum(scores[scores > 0])
    neg_score <- sum(scores[scores < 0])
    if (pos_score > abs(neg_score)){
        return(pos_score)
    } else{
        return(neg_score)
    }
    
})


dpis_gene.sort <- sort(dpis_gene, decreasing=TRUE)

if (database == "GO"){
    gse <- gseGO(geneList    = dpis_gene.sort,
                OrgDb        = OrgDb,
                ont          = ont,
                pvalueCutoff = pval,
                verbose      = FALSE,
                keyType      = gene_type)

} else if (database == "KEGG") {
    if (gene_type == 'ENTREZID'){
        dpis_gene.sort <- dpis_gene.sort
    } else {
        genes <- names(dpis_gene.sort)
        gene_df <- bitr(genes, fromType = gene_type, toType = "ENTREZID", OrgDb = OrgDb)
        new.gene.sort <- dpis_gene.sort[gene_df[, gene_type]]
        names(new.gene.sort) <- gene_df$ENTREZID
    }    
    gse <- gseKEGG(geneList    = new.gene.sort,
                  organism     = organism,
                  pvalueCutoff = pval,
                  verbose      = FALSE)
}


RData.out <- file.path(out.dir, sprintf("%s.%s", name, "RData"))
csv.out <- file.path(out.dir, sprintf("%s.%s", name, "csv"))

if (dim(gse)[1] == 0){
    write_csv(as.data.frame(gse), file = csv.out)
    save(gse, file = RData.out) 

    # # return()
    # stop("No enrichment found in input...")
} else {
    keyType <- ifelse(database == "KEGG", "ENTREZID", gene_type)
    ego <- setReadable(gse, OrgDb = OrgDb, keyType = keyType)
    save(gse, file = RData.out)
    write_csv(as.data.frame(gse), file = csv.out)        
}

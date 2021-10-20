suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))


args <- commandArgs()

out.dir <- args[6]
pval  <- as.numeric(args[7])
org_db <- args[8]
ont <- args[9]
gene_type <- args[10]
name <- args[11]
dpsi_file <- args[12]

script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

OrgDb <- load_OrgDb(org_db)


dpsi_df <- read_tsv(dpsi_file, skip =1, col_names = F,
                col_types = cols("c", "d", "d")) %>% drop_na()

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
    }else{
        return(neg_score)
    }
    
})


dpis_gene.sort <- sort(dpis_gene, decreasing=TRUE)


ego <- gseGO(geneList     = dpis_gene.sort,
              OrgDb        = OrgDb,
              ont          = ont,
              pvalueCutoff = pval,
              verbose      = FALSE,
              keyType = gene_type)


RData.out <- file.path(out.dir, sprintf("%s.%s", name, "RData"))
csv.out <- file.path(out.dir, sprintf("%s.%s", name, "csv"))

save(ego, file = RData.out)

ego <- setReadable(ego, OrgDb = OrgDb, keyType = gene_type)

write_csv(as.data.frame(ego), file = csv.out)

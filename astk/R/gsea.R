suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))



dpsi_file <- "/home/huangshenghui/project/asla/output/fb_e10/sig_015/6_SE.sig.dpsi"

dpsi_df <- read_tsv(dpsi_file, skip =1, col_names = F,
                col_types = cols("c", "d", "d")) %>% drop_na()

colnames(dpsi_df) = c("event_id", "dpsi", "pval")


dpsi_df$gene <- gsub("\\..*", "",  dpsi_df$event_id) 


dpis_gene <- sapply(unique(dpsi_df$gene), function(x){
    scores <- dpsi_df[dpsi_df$gene == x, ]$dpsi
    pos_score <- sum(scores[scores > 0])
    neg_score <- sum(scores[scores < 0])
    if (pos_score > abs(neg_score)){
        return(-pos_score)
    }else{
        return(-neg_score)
    }
    
})


dpis_gene.sort <- sort(dpis_gene, decreasing=TRUE)

gmt <- read.gmt("/home/huangshenghui/project/asla/geneset.gmt")


gene_list <- merge.dpsi_df  %>% 
    select(gene, dpsi) %>% 
    arrange(desc(dpsi))


geneList <- merge.dpsi_df$sum_dpsi
names(geneList) <- merge.dpsi_df$gene




geneList[c("ENSMUSG00000026527", "ENSMUSG00000062372")]

library(org.Mm.eg.db)

eg = bitr(names(geneList), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")


geneList[]

names(geneList) 
new_eg <- eg[eg$ENSEMBL %in% names(geneList), ]

new_gene_list <- geneList[new_eg$ENSEMBL]

names(new_gene_list) <- new_eg$SYMBOL




egmt2 <- GSEA(new_gene_list, TERM2GENE=gmt, verbose=FALSE, pvalueCutoff =0.5)

eg

length(geneList)
dim(dpsi_df)
length(geneList[eg$ENSEMBL])

length(unique(names(geneList)))
geneList["ENSMUSG00000021820"]

table(names(geneList))


gmt$gene <- stringr::str_to_title(gmt$gene)


intersect(names(new_gene_list), gmt$gene)

dpis_gene.sort

ego4 <- gseGO(geneList     = dpis_gene.sort,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              keyType = "ENSEMBL")


p1 <- gseaplot(ego4, geneSetID = "GO:0007399", by = "all", title = "GSEA")

p1


View(as.data.frame(ego3))
library(httpgd)
hgd()



R.version

ego3


gseGO
clusterProfiler:::GSEA_internal

enrichplot::gseaplot2(ego3, geneSetID = "GO:0007399", title = ego3$Description[1])



table(dpsi_df$gene)



dpis_gene["ENSMUSG00000005656"]

sort(dpis_gene, decreasing=TRUE)
dpsi_df[dpsi_df$gene == "ENSMUSG00000025871", ]

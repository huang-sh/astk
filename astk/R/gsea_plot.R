args <- commandArgs()

out.file <- args[6]
RData.file <- args[7]
term_ID <- args[8:length(args)]

load(RData.file)

if (length(term_ID) == 1){
    p <- enrichplot::gseaplot(gse, 
                  geneSetID = term_ID,
                  by        = "all", 
                  title     = gse[term_ID, "Description"])

}else {
    p <- enrichplot::gseaplot2(
                gse, 
                geneSetID = term_ID,
                title     = "GSEA")
}

ggplot2::ggsave(out.file, plot = p)

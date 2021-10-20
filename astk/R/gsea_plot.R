suppressMessages(library(enrichplot))



args <- commandArgs()

out.file <- args[6]
RData.file <- args[7]
term_ID <- args[8:length(args)]

print(term_ID)

load(RData.file)

if (length(term_ID) == 1){
    p <- gseaplot(ego, geneSetID = term_ID,
        by = "all", title = ego[term_ID, "Description"])

}else {
    p <- gseaplot2(ego, geneSetID = term_ID,
        title = "GSEA")
}


ggplot2::ggsave(out.file, plot = p)



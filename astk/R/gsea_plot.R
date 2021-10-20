suppressMessages(library(enrichplot))



args <- commandArgs()

term_ID <- args[6]
out.file <- args[7]
RData.file <- args[8]


load(RData.file)


p <- gseaplot(ego, geneSetID = term_ID,
     by = "all", title = ego[term_ID, "Description"])


ggplot2::ggsave(out.file, plot = p)

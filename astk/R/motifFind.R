suppressMessages(library(memes))
suppressMessages(library(tidyverse))
suppressMessages(library(universalmotif))


args <- commandArgs()

outdir <- args[6]
fa_file <- args[7]


streme_out <- runStreme(fa_file, "shuffle", outdir = outdir, minw = 5)


p <- streme_out %>% 
  to_list() %>% 
  head() %>%
  view_motifs()

ggsave(file.path(outdir, "motif_n10.pdf"), plot = p)

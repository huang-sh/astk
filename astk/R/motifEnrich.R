library(memes)
library(tidyverse)



args <- commandArgs()

outdir <- args[6]
meme_db_path <- args[7]

files <- args[8:length(args)]



res_ls <- lapply(files, function(file){
    out.dir <- file.path(outdir, tools::file_path_sans_ext(basename(file)))
    runAme(file , database=meme_db_path, outdir=out.dir)    
})

file_names <- sapply(files, function(file)tools::file_path_sans_ext(basename(file)))
names(res_ls) <- file_names


data <- res_ls %>% 
    dplyr::bind_rows(.id = "condition") %>% 
    dplyr::group_by(condition) %>% 
    dplyr:: slice_head(n = 25) %>%   ## 选取了前25个motif
    ungroup()

plot <- plot_ame_heatmap(data, group = condition)

ggsave(file.path(outdir, "heatmap.pdf"), plot = plot)

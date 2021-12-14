library(memes)
library(magrittr)


script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- argparse::ArgumentParser()
parser$add_argument("--tfile",  nargs='+', help="treatment file path")
parser$add_argument("--cfile",  nargs='+', help="control file path")
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--database", help="database")
parser$add_argument("--organism", help="organism")

args <- parser$parse_args()

if (args$database == "CISBP"){
   motif_path <- file.path(dirname(dirname(script_path)), "data", "CIS_RBP_motif.txt")
    
} else {
   motif_path <- file.path(dirname(dirname(script_path)), "data", "ATtRACT_RBP_motif.txt")
}

motif_ls <- universalmotif::read_cisbp(motif_path)

if (args$organism != "0"){
    print(args$organism)
    mf <- universalmotif::filter_motifs(motif_ls, organism = args$organism)
} else {
    mf <- motif_ls
}

readFasta <- function(file){
    dna <- Biostrings::readDNAStringSet(file)
    rna <- Biostrings::BStringSet(Biostrings::RNAStringSet(dna))
    return(rna)
}

file_ls  <-  list()

for (i in seq(1, length(args$tfile))){
    file_ls[[i]] <- c(args$tfile[i], args$cfile[i])
}

res_ls <- lapply(file_ls, function(files){
    tfa <- readFasta(files[1])
    
    if (files[2] == "0"){
        cfa <- "shuffle"
    } else {
       cfa <- readFasta(files[2])
    }
    
    out.dir <- file.path(args$outdir, tools::file_path_sans_ext(basename(files[1])))
    runAme(tfa , database = mf, outdir=out.dir, control = cfa)    
})

# file_names <- sapply(files, function(file)tools::file_path_sans_ext(basename(file)))
# names(res_ls) <- file_names


data <- res_ls %>% 
    dplyr::bind_rows(.id = "condition") %>% 
    dplyr::group_by(condition) %>% 
    dplyr:: slice_head(n = 20) %>%   ## 选取了前20个motif
    dplyr::ungroup()

plot <- plot_ame_heatmap(data, group = condition)

ggplot2::ggsave(file.path(args$outdir, "heatmap.pdf"), plot = plot, width = 10)

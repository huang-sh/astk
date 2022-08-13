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
parser$add_argument("--meme_path ", help="meme_path ")

args <- parser$parse_args()

motif_dir <- file.path(dirname(dirname(script_path)), "data/motif")
if (args$database == "CISBP-RNA"){
   motif_path <- file.path(motif_dir, "CISBP-RNA", sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
    
} else {
   motif_path <- file.path(motif_dir, "ATtRACT",  sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
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
    runAme(tfa , database = mf, outdir=out.dir, control = cfa, meme_path =args$meme_path )    
})

# file_names <- sapply(files, function(file)tools::file_path_sans_ext(basename(file)))
# names(res_ls) <- file_names

df = tryCatch({
        data <- res_ls %>% 
            dplyr::bind_rows(.id = "condition") %>% 
            dplyr::group_by(condition) %>% 
            dplyr:: slice_head(n = 20) %>%   ## 选取了前20个motif
            dplyr::ungroup()

        plot <- plot_ame_heatmap(data, group = condition)

        ggplot2::ggsave(file.path(args$outdir, "heatmap.pdf"), plot = plot, width = 10)

    }, warning = function(w) {
        print(w)
    }, error = function(e) {    
    
    }
    )
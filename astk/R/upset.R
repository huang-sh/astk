
parser <- argparse::ArgumentParser()

parser$add_argument("--file", nargs='+', help="file names")
parser$add_argument("--name", nargs='+', help="file path")
parser$add_argument("--output", help="file path")
parser$add_argument("--fmt", help="figure format")
parser$add_argument("--resolution", help="resolution")
parser$add_argument("--width", type="double", help="figure width")
parser$add_argument("--height", type="double", help="figure height")
parser$add_argument("--dg", action="store_true", default=F, help="dpsi group")


script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

args <- parser$parse_args()

files <- args$file
names(files) <- args$name

if (args$dg){
    rep_num <- 2
    dpsi_ls_n <- lapply(files, function(x){
        dpsi <- readr::read_tsv(x, skip = 1, 
                    col_names = c("event_id", "dPSI", "pval"),
                    col_types = readr::cols("c", "d", "d"))
        pos <- dpsi[dpsi$dPSI>0, ]$event_id
        neg <- dpsi[dpsi$dPSI<0, ]$event_id
        list(pos=pos, neg=neg)
    })

    dpsi_ls <- unlist(dpsi_ls_n, recursive = F)

} else {
    rep_num <- 1
    dpsi_ls <- lapply(files, function(x){
        dpsi <- readr::read_tsv(x, skip = 1, 
                    col_names = c("event_id", "dPSI", "pval"),
                    col_types = readr::cols("c", "d", "d"))
        dpsi$event_id
    })
}


pre_colors <- c("#a1d9f9", "#d69bc5", "#ff9999", "#ffc799", "#99d1b8")

colors <- pre_colors[1:length(files)]
colors[is.na(colors)] <- "black"


p <- UpSetR::upset(UpSetR::fromList(dpsi_ls), 
                   nsets          = length(files), 
                   keep.order     = TRUE, 
                   sets           = names(dpsi_ls), 
                   matrix.color   = "red", 
                   main.bar.color = "black",
                   sets.bar.color = rep(colors, each=rep_num))


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+', help="file names")
parser$add_argument("--name", nargs='+', help="file path")


args <- parser$parse_args()

files <- args$file
names(files) <- args$name

evnet_id.ls <- lapply(files, function(x){
    df <- readr::read_tsv(x, skip = 1)
    colnames(df)[1] <- "event_id"
    df$event_id
})

pre_colors <- c("#a1d9f9", "#d69bc5", "#ff9999", "#ffc799", "#99d1b8")

colors <- pre_colors[1:length(files)]
colors[is.na(colors)] <- "black"


p <- UpSetR::upset(UpSetR::fromList(evnet_id.ls), 
                   nsets          = length(files), 
                   keep.order     = TRUE, 
                   sets           = names(evnet_id.ls), 
                   matrix.color   = "red", 
                   main.bar.color = "black",
                   sets.bar.color = colors)


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

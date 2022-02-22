
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyverse))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+', help="file names")


args <- parser$parse_args()

plot_ls <- lapply(args$file, function(file){
    signal_df <- read_csv(file)

    filter.signal <- signal_df %>% 
        group_by(event_idx, mark, direction, anchor) %>% 
        summarise(value = mean(signal)) %>% 
        group_by(mark, direction, anchor) %>% 
        summarise(value = mean(value)) %>% 
        arrange(mark)
            
    anchors <- unique(filter.signal$anchor)    

    signal_mtx <- sapply(anchors, function(x){
        up_sig <- filter.signal %>% filter(anchor == x & direction == "upstream")
        down_sig <- filter.signal %>% filter(anchor == x & direction == "downstream")
        signals <- c(up_sig$value, down_sig$value)
        return(signals)
    })

    signal_mtx <- matrix(signal_mtx, nrow = length(unique(filter.signal$mark)))

    rownames(signal_mtx) <- unique(filter.signal$mark)

    ud <- rep(c("up", "down"), times = length(anchors))
    anchors <- rep(anchors, each = 2)
    colnames(signal_mtx) <- paste(anchors, ud, sep="_")
    filename  <-  tools::file_path_sans_ext(basename(file))
    p <- Heatmap(signal_mtx, name = filename, cluster_columns = F, cluster_rows = F)
    return(p)
})


ps <- Reduce("+", plot_ls) 

save_fig(ps, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

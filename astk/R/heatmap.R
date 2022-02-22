suppressMessages(library(ComplexHeatmap))


script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs="+", help="psi file path")
parser$add_argument("--clusterinfo", help="cluster info file path")


args <- parser$parse_args()


psi_df_ls <- lapply(args$file, function(file){
    psi_df <- read.delim(file, row.names = 1) # col_types = cols("c", "d", "d")
    return(psi_df)
})


my_inner_join <- function(x, y){
    rows <- intersect(rownames(x), rownames(y))
    cbind(x[rows, ], y[rows, ])
}

dat <- Reduce(my_inner_join, psi_df_ls)
dat <- na.omit(dat)

if (file.exists(args$clusterinfo)){

    cls_ext <- tools::file_ext(args$clusterinfo)

    if (cls_ext == "tsv"){
        exon_cluster <- read.delim(args$clusterinfo, row.names = 1)
    } else if (cls_ext == "csv") {
        exon_cluster <- read.delim(args$clusterinfo, row.names = 1)
    }

    exon_cluster <- as.data.frame(exon_cluster)
    rownames(exon_cluster) <- exon_cluster$event_id

    exon_cluster[rownames(dat), ]$cluster
    split <- paste("Cluster", exon_cluster[rownames(dat), ]$cluster, sep = "")
} else {
   split <- NULL
}


col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 9), 
                c('#0077B6', '#00B4D8', '#90E0EF',
                  '#CAF0F8', '#FAE0E4', '#F9BEC7',
                  '#FF99AC', '#FF7096', '#FF477E'))

p <- Heatmap(as.matrix(dat), 
             name            = "Heatmap", 
             show_row_dend   = T, 
             border          = TRUE, 
             col             = col_fun,
             cluster_columns = F, 
             show_row_names  = F,
             row_split       = split, 
             row_gap         = unit(c(3), "mm"))

save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

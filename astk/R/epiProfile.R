suppressMessages(library(metagene2))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+')
parser$add_argument("--title")
# parser$add_argument("--ylim", nargs=2, type="double")

args <- parser$parse_args()

files <- args$file
title <- args$title


df_ls <- lapply(files, function(file) {
    read.csv(file)
})

cols <- as.vector(sapply(df_ls, colnames))

marks <- unique(sapply(str_split(cols, "_"), function(x) {x[1]}))
regions <- unique(sapply(str_split(cols, "_"), function(x) {x[2]}))
anchors <- unique(sapply(str_split(cols, "_"), function(x) {x[3]}))

names(marks) <- marks
print(anchors)
print(regions)
print(marks)

split_coverage <- lapply(marks, function(m) {
  
    ras <- paste(m, rep(regions, each = length(anchors)), anchors, sep = "_")
    names(ras) <- paste(rep(regions, each = length(anchors)), anchors, sep = "_")
    lapply(ras, function(ra){
        for (df in df_ls){
            df_cols <- colnames(df)
            ncol <- df_cols[str_detect(df_cols, ra)]
            if (length(ncol)){
                if (length(ncol) == 1){
                    sdf <- data.frame("nol" = df[, ncol])
                    colnames(sdf) <- ncol
                }else{
                    sdf <- df[, ncol]
                }
                
                break()
            }
        }
        sdf
    })
})
print(length(split_coverage[[1]]))

ci_df <- metagene2:::calculate_matrices_ci(
            split_coverage,
            1000, 
            0.05,
            "bin",
            metagene2:::Parallel_Job$new(4))   # metagene2:::Parallel_Job$new(1)


ci_df <- as.data.frame(ci_df)
ci_df$group <- factor(sapply(str_split(ci_df$region, "_"), function(x) {x[2]}))
p <- plot_metagene(ci_df, facet_by=~group, group_by="region") + ggplot2::ggtitle(title)


# if (! is.null(args$ylim)){
#     p <- p + ggplot2:: coord_cartesian(ylim = args$ylim)
# }

save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

suppressMessages(library(tidyverse))

parser <- fig_cmd_parser()


parser$add_argument("--file", nargs='+', help="psi file path")
parser$add_argument("--groupname", nargs='+', help="group names")



args <- parser$parse_args()

psi_df_ls <- lapply(args$file, function(file){
    psi_df <- read.delim(file, row.names = 1)
    return(psi_df)
})

my_inner_join <- function(x, y){
    rows <- intersect(rownames(x), rownames(y))
    cbind(x[rows, ], y[rows, ])
}

dat <- Reduce(my_inner_join, psi_df_ls)
dat <- na.omit(dat)


dat <- na.omit(dat)
print(dim(dat))
df_pca <- prcomp(t(dat))
df <- as.data.frame(df_pca$x)

if (is.null(args$groupname)){
    names <- lapply(args$file, function(file){
        name <- tools::file_path_sans_ext(basename(file))
        psi_df <- read.delim(file, row.names = 1) # col_types = cols("c", "d", "d")
        rep(name, dim(psi_df)[2])
    })    
} else{
    names <- lapply(1:length(args$file), function(idx){
        psi_df <- read.delim(args$file[idx], row.names = 1) # col_types = cols("c", "d", "d")
        rep(args$groupname[idx], dim(psi_df)[2])
    })
}
df$group <- unlist(names)

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)

percentage <- paste(colnames(df), "(", paste( as.character(percentage), "%", ")", sep="") )

p <- ggplot(df, aes(x=PC1,y=PC2, color=group)) +
    geom_point() +
    # stat_ellipse(level = 0.95, show.legend = F) +
    xlab(percentage[1]) + 
    ylab(percentage[2])


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

library(stringr)
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(readr))


args <- commandArgs()

out <- args[6]
psi_files <- args[7:length(args)]


psi_df_ls <- lapply(psi_files, function(file){
    psi_df <- read_tsv(file, col_types = cols("c", "d", "d"))
    return(psi_df)
})

names <- lapply(psi_files, function(file){
    name <- tools::file_path_sans_ext(basename(file))
    psi_df <- read_tsv(file, col_types = cols("c", "d", "d"))
    rep(name, dim(psi_df)[2]-1)
})

my_inner_join <- function(x, y){
    inner_join(x, y, by = "event_id")
}

dat <- Reduce(my_inner_join, psi_df_ls)


dat <- na.omit(dat)
df_pca <- prcomp(t(dat[, -1]))

df <- as.data.frame(df_pca$x)

df$group <- unlist(names)



percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)

percentage <- paste(colnames(df), "(", paste( as.character(percentage), "%", ")", sep="") )

p <- ggplot(df, aes(x=PC1,y=PC2, color=group)) +
    geom_point()+xlab(percentage[1]) + 
    ylab(percentage[2])


ggsave(out, plot = p)

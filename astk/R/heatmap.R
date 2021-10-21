library(stringr)
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(ComplexHeatmap))

args <- commandArgs()

out <- args[6]
cluster_file <- args[7]
psi_files <- args[8:length(args)]


psi_df_ls <- lapply(psi_files, function(file){
    psi_df <- read_tsv(file, col_types = cols("c", "d", "d"))
    return(psi_df)
})


my_inner_join <- function(x, y){
    inner_join(x, y, by = "event_id")
}

dat <- Reduce(my_inner_join, psi_df_ls)
dat <- na.omit(dat)


if (file.exists(cluster_file)){

    cls_ext <- tools::file_ext(cluster_file)

    if (cls_ext == "tsv"){
        exon_cluster <- read_tsv(cluster_file, show_col_types = FALSE)
    } else if (cls_ext == "csv") {
        exon_cluster <- read_csv(cluster_file, show_col_types = FALSE)
    }

    exon_cluster <- as.data.frame(exon_cluster)
    rownames(exon_cluster) <- exon_cluster$event_id

    exon_cluster[dat$event_id, ]$cluster
    split <- paste("Cluster", exon_cluster[dat$event_id, ]$cluster, sep = "")
} else {
   split <- NULL
}



col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 9), 
            c('#0077B6', '#00B4D8', '#90E0EF', '#CAF0F8',
               '#FAE0E4', '#F9BEC7', '#FF99AC', '#FF7096', '#FF477E'))


p <- Heatmap(dat[, -1], 
             name = "Heatmap", 
             show_row_dend  = T, 
             border = TRUE, 
             col = col_fun,
             cluster_columns = F, 
             row_split  = split, 
             row_gap = unit(c(3), "mm"))


if (str_detect(out, "pdf")){
    pdf(out)
    print(p)
    dev.off()    
}else if (str_detect(out, "png")) {
    png(out)
    print(p)
    dev.off()   
}




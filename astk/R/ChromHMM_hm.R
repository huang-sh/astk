suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(ComplexHeatmap))

args <- commandArgs()

file <- args[6]


if (str_detect(file, "_overlap")){

    data <- read_tsv(file)
    mat <- data[-7, -1]
    colnames(mat) <- sapply(str_split(colnames(data)[-c(1)], "\\."), function(x)x[1])

    pmat <- scale(mat)
    name <- str_split(basename(file), "_overlap")[[1]][1]

    p <- Heatmap(pmat, name=name, cluster_columns = F, cluster_rows = F,
               column_title = "Category", column_title_side = "bottom",
               row_title = "State (Emission order)"
                )
    png(str_replace(file, ".txt", ".png"),
        width = dim(mat)[2] * 30, height = dim(mat)[1] * 50)
    print(p)
    dev.off()
}

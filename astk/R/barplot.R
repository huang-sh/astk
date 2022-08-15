suppressMessages(library(ggplot2))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+', help="file path")
parser$add_argument("--name", nargs='+', help="file names")
parser$add_argument("--dg", action='store_true', help="dg")


args <- parser$parse_args()

files <- args$file


filenames <- sapply(stringr::str_split(basename(files), "\\."), function(x) x[1])

if (length(filenames) == length(args$name)){
    filenames <- args$name
} 

names(files) <- filenames

if (args$dg){
    count.ls <- lapply(files, function(x){
        df <- readr::read_tsv(x, skip = 1, show_col_types = FALSE)
        colnames(df)[1:3] <- c("event_id", "dPSI", "pval")
        c(dim(df[df$dPSI>0, ])[1], dim(df[df$dPSI<0, ])[1])
    })

    filenames <- rep(filenames, each=2)
    types <- rep(c("+", "-"), time=length(files)) 

} else {
    count.ls <- lapply(files, function(x){
        df <- readr::read_tsv(x, skip = 1, show_col_types = FALSE)
        colnames(df)[1] <- "event_id"
        dim(df)[1]
    })   
   filenames <- filenames
   types <- "u"    
}


data <- data.frame(
    count = unlist(count.ls),
    name = filenames,
    type = types
)

data$name <- factor(data$name)

p <- ggplot(data) +
    theme_bw() +
    theme(panel.grid=element_blank()) + 
    geom_col(mapping = aes(x = name, y = count, fill=type), 
             position = "dodge") +
    scale_fill_manual(values = c("#83cdfd", "#ffee6b")) + 
    scale_y_continuous(expand = expansion(mult = c(0, .1)))


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        res    = args$resolution,
        units  = "in")

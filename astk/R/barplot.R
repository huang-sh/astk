suppressMessages(library(ggplot2))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+', help="file names")
parser$add_argument("--name", nargs='+', help="file path")
parser$add_argument("--dg", action='store_true', help="dg")


args <- parser$parse_args()

files <- args$file
names(files) <- args$name

# files <- list.files("~/project/astk/paper/fb_e11_ct/sig01",
#             pattern = "*_SE.sig.dpsi", full.names = T)



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

p <- ggplot(data) +
    geom_col(mapping = aes(x = name, y = count, fill=type), 
             position = "dodge")

save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        res    = args$resolution,
        units  = "in")

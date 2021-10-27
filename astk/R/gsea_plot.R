
script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--RData", help="RData file path")
parser$add_argument("--termid", nargs='+', help="term id")


args <- parser$parse_args()


load(args$RData)

term_ID <- args$termid

if (length(term_ID) == 1){
    p <- enrichplot::gseaplot(gse, 
                  geneSetID = term_ID,
                  by        = "all", 
                  title     = gse[term_ID, "Description"])

}else {
    p <- enrichplot::gseaplot2(
                gse, 
                geneSetID = term_ID,
                title     = "GSEA")
}


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)
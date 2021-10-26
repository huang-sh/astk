
parser <- argparse::ArgumentParser()

parser$add_argument("--output", help="output file path")
parser$add_argument("--RData", help="RData file path")
parser$add_argument("--termid", nargs='+', help="term id")
parser$add_argument("--fmt", help="figure format")
parser$add_argument("--resolution", type="integer", help="resolution")
parser$add_argument("--width", type="double", help="figure width")
parser$add_argument("--height", type="double", help="figure height")


args <- parser$parse_args()

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

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
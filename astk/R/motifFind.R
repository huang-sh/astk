
script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- argparse::ArgumentParser()

parser$add_argument("--outdir", help="ouput dirctory")
parser$add_argument("--fasta",  help="fasta file path")
parser$add_argument("--pvalue",  help="pvalue")
parser$add_argument("--minw",  help=" minimal motifs width")
parser$add_argument("--maxw",  help=" maximal motifs width")


args <- parser$parse_args()


streme_out <- memes::runStreme(args$fasta,
                     "shuffle",
                      outdir = args$outdir,
                      minw   = args$minw,
                      maxw   = args$maxw,
                      pvt = args$thresh)


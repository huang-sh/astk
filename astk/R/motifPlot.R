script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()

parser$add_argument("--meme", help="meme file path")
parser$add_argument("--motifId", nargs="+", help="meme file path")

args <- parser$parse_args()

meme_ls <- universalmotif::read_meme(args$meme)

mf <- universalmotif::filter_motifs(meme_ls, name = args$motifId)

p <- universalmotif::view_motifs(mf)


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

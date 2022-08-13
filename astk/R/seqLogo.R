script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()

parser$add_argument("--pwm", help="pwm file path")


args <- parser$parse_args()

pwm <- universalmotif::read_matrix(args$pwm, headers = ">", positions = "rows")


p <- universalmotif::view_motifs(pwm, use.type = "PPM")
# motif <- universalmotif::convert_type(pwm, "PPM")
# ## Only need the matrix itself
# motif <- motif["motif"]
# p <- seqLogo::seqLogo(motif)

save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

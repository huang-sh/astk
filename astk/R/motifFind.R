
script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- argparse::ArgumentParser()

parser$add_argument("--tfile", nargs='+', help="treatment file path")
parser$add_argument("--cfile", nargs='+', help="control file path")
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--pvalue",  help="pvalue")
parser$add_argument("--minw",  help=" minimal motifs width")
parser$add_argument("--maxw",  help=" maximal motifs width")
parser$add_argument("--database", help="database")
parser$add_argument("--organism", help="organism")
parser$add_argument("--eval", help="eval")
parser$add_argument("--meme_path ", help="meme_path ")

args <- parser$parse_args()


file_ls  <-  list()

for (i in seq(1, length(args$tfile))){
    file_ls[[i]] <- c(args$tfile[i], args$cfile[i])
}

motif_dir <- file.path(dirname(dirname(script_path)), "data/motif")

if (args$database == "CISBP-RNA"){
   motif_path <- file.path(motif_dir, "CISBP-RNA", sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
    
} else {
   motif_path <- file.path(motif_dir, "ATtRACT",  sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
}


res_ls <- lapply(file_ls, function(files){
    tfa <- readFasta(files[1])

    if (files[2] == "0"){
        cfa <- "shuffle"
    } else {
       cfa <- readFasta(files[2])
    }
    
    out.dir <- file.path(args$outdir, tools::file_path_sans_ext(basename(files[1])))
    streme_out <- memes::runStreme(tfa,
                    cfa,
                    outdir = out.dir,
                    minw   = args$minw,
                    maxw   = args$maxw,
                    pvt    = args$thresh,
                    alph   = "rna",
                    meme_path = args$meme_path
                    )
    memes::runTomTom(
        streme_out,
        database = mf,
        outdir = out.dir,
        thresh = args$eval,
        min_overlap = 5,
        evalue = TRUE,
        silent = TRUE,
         meme_path = args$meme_path
        )

})

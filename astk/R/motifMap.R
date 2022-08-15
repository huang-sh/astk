suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--outdir", help="output directory")
parser$add_argument("--meme", help="meme motif")
parser$add_argument("--fasta", nargs="+", help="fasta files")
parser$add_argument("--step", type="integer", help="slide window step size")
parser$add_argument("--bin", type="integer", help="bin size")
parser$add_argument("--seqid", nargs="+", help="fasta files id")
parser$add_argument("--center", nargs="+", type="integer", help="fasta files")
parser$add_argument("--meme_path ", help="meme_path ")

args <- parser$parse_args()

step <- args$step
win <- args$bin
meme <- args$meme

seq_files <- args$fasta
seq_names <- args$seqid
centers <- args$cente
meme_ls <- universalmotif::read_meme(meme)

if (class(meme_ls) != "list"){
    meme_ls  <-  list(meme_ls)

}

motif_names <- sapply(meme_ls, function(x) x@name)
motif_lens  <-  sapply(meme_ls, function(x) dim(x@motif)[2])
max_motif_len <- max(motif_lens)

names(seq_files)  <- seq_names
names(centers) <- seq_names

df_ls <- lapply(seq_names, function(n){
    faset <- readBStringSet(seq_files[n])
    seq_len <- width(faset)[1]
    starts <- seq(1, seq_len, by = step) 
    win_starts <- starts[(seq_len-starts+1) >= max_motif_len] 
    res_ls <- lapply(win_starts, function(s){
        if (s == win_starts[length(win_starts)]){
            subfa <- subseq(faset, start=s, end=min(width(faset))) 
        }
        else {
           subfa <- subseq(faset, start=s, end=NA, width=win) 
        }
        res <- memes::runFimo(subfa, meme_ls, meme_path = args$meme_path)
        return(res)
    })

    motif_count <- lapply(res_ls, function(x) {
        counts <- sapply(motif_names, function(n){
            xdf <- as.data.frame(x)
            dim(xdf[xdf$motif_id == n,])[1]
        })
        names(counts) <- motif_names
        counts
    })

    count_df <-  data.frame(t(sapply(motif_count,c)))
    len_factor <- motif_lens / max_motif_len
    count_df <- as.data.frame(mapply('/', count_df, len_factor))
    colnames(count_df) <- motif_names
    new_motif_names <- colnames(count_df)

    count_df$pos <- win_starts - centers[n]
    count_df$seqId <- n

    long_data <- count_df %>% 
        gather(all_of(new_motif_names), key="motif_id", value ="count") %>% 
        mutate(count = count / length(faset))
    long_data
})

all_df <- Reduce(rbind, df_ls)

r <- lapply(unique(all_df$motif_id), function(x){
    p <- ggplot(data = all_df[all_df$motif_id ==x, ]) +
        theme_bw() +
        geom_line(mapping = aes(x = pos, y=count, group=motif_id), color="#83cdfd") + 
        scale_fill_manual(values = c('pink', "#83cdfd")) + 
        facet_wrap(~seqId, ncol = length(seq_names), scales = "free_x") +
        labs(x = "Relative Position", y="Motif Coverage") + 
        ggtitle("RNA map")
    output <- file.path(args$outdir, sprintf("%s.%s", x, args$fmt))

    save_fig(p, 
            output, 
            format = args$fmt,
            width  = args$width, 
            height = args$height, 
            units  = "in",
            res    = args$resolution)
})

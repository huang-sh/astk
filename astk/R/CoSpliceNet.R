suppressMessages(library(tidyverse))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicFeatures))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- argparse::ArgumentParser()

parser$add_argument("--psiMeta",  help="psi meta file path")
parser$add_argument("--tqMeta", help="quant.sh meta file path")
parser$add_argument("--organism", help="organism")
parser$add_argument("--database", help="database")
parser$add_argument("--txdb", help="txdb")
parser$add_argument("--fasta", help="fasta")
parser$add_argument("--output", help="output")


args <- parser$parse_args()

motif_dir <- file.path(dirname(dirname(script_path)), "data/motif")
if (args$database == "CISBP-RNA"){
   motif_path <- file.path(motif_dir, "CISBP-RNA", sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
    
} else {
   motif_path <- file.path(motif_dir, "ATtRACT",  sprintf("%s.meme", args$organism))
   mf <- universalmotif::read_meme(motif_path)
}


psi_meta_file <- args$psiMeta
sf_meta_file <- args$tqMeta

psi_meta <- read.csv(psi_meta_file)
sf_meta <- read.csv(sf_meta_file)

mean_psi_ls <- lapply(psi_meta$path, function(file){
    sub_psi_df  <- read.delim(file, sep = "\t")
    mean_psi <- rowMeans(sub_psi_df[, -1])
    names(mean_psi) <- sub_psi_df[, 1]
    mean_psi
})

psi_df <- Reduce(cbind, mean_psi_ls)

colnames(psi_df) <- psi_meta$group


txdb <- AnnotationDbi::loadDb(args$txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

txi <- tximport(sf_meta$path, type = "salmon", tx2gene = tx2gene)

sampleTable <- data.frame(condition = sf_meta$group)
rownames(sampleTable) <- sf_meta$name

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

keep <- rowSums(counts(dds)) >= dim(sf_meta)[1] * 10
dds <- dds[keep,]

count.sf <-  estimateSizeFactors(dds)
normal.mtx <- counts(count.sf, normalized=T)


mean_exp_ls <- lapply(unique(sf_meta$group), function(g){
    gnames  <- sf_meta[sf_meta$group == g, "name"]
    rowMeans(normal.mtx[, gnames])
})

exp_df <- Reduce(cbind, mean_exp_ls)

colnames(exp_df) <- unique(sf_meta$group)
rownames(exp_df) <- sapply(strsplit(rownames(normal.mtx), "\\."), function(x)x[1])


fa <- unique(readFasta(args$fasta))

as_event_id <- sapply(strsplit(names(fa), "::"), function(x)x[1])

fa_pos <- sapply(strsplit(names(fa), "::"), function(x)x[2])


names(fa) <- paste(gsub(":", "_", as_event_id), fa_pos, sep="__")
fimo.out <-  memes::runFimo(fa, mf)  


fimo_df <- unique(as.data.frame(fimo.out))


motif_as_event_id <- sapply(strsplit(as.character(fimo_df$seqnames), "__"), function(x)x[1])
motif_as_event_id <- gsub("_", ":", motif_as_event_id)


rbp_id <- sapply(strsplit(fimo_df$motif_alt_id, ":"), function(x)x[2])

motif_link_df <- data.frame(
    "rbp" = fimo_df$motif_alt_id,
    "rbp_id" = rbp_id,
    "AS_event" = motif_as_event_id
)

cor_ls <- apply(motif_link_df, 1, function(x){
    rbp_exp <- as.numeric(exp_df[x["rbp_id"], ])
    as_event_psi <- as.numeric(psi_df[x["AS_event"], ])
    coor <- cor.test(rbp_exp, as_event_psi, method = "spearman")
    return(c(cor=coor$estimate, pval=coor$p.value, rbp=x["rbp_id"], AS_event=x["AS_event"]))

})

gene_cor <- as.data.frame(t(as.data.frame(cor_ls)))
colnames(gene_cor) <- c("cor", "pval", "rbp", "event")

gene_cor$motif_alt_id <- motif_link_df$rbp

high_link <- gene_cor[abs(as.numeric(gene_cor$cor)) > 0.8 & gene_cor$pval < 0.05, ]

write.csv(high_link, file = args$output, row.names=F)

suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--pval", type="double", help="pval")
parser$add_argument("--qval", type="double", help="qval")
parser$add_argument("--database", help="database")
parser$add_argument("--ontology", help="ontology")
parser$add_argument("--orgdb", help="orgdb")
parser$add_argument("--genetype", help="gene type")
parser$add_argument("--organism", help="organism")
parser$add_argument("--files",  nargs='+', help="dpsi files")
parser$add_argument("--xlabel",  nargs='+', help="dpsi files")

args <- parser$parse_args()

out.dir <- args$outdir
pval  <- args$pval
qval  <- args$qval
database <- args$database
orgdb <- args$orgdb
gene_type <- args$genetype
organism <- args$organism
dpsi_files <- args$files
xlabel <- args$xlabel
ontology <- args$ontology

suppressMessages(library(orgdb, character.only = T))


filenames <- unlist(lapply(dpsi_files, function(file){
                            tools::file_path_sans_ext(basename(file))
                        }))   

if (length(unique(xlabel)) == length(dpsi_files)){
    filenames <- xlabel
} else if (length(filenames) != length(unique(filenames))) {
   filenames <- paste(filenames, seq(1, length(filenames)), sep=".") 
} else {
   filenames <- filenames
}


gene_ls <- lapply(dpsi_files, function(file){
    dpsi <- read_tsv(file, 
                    skip      = 1, 
                    col_names = F,
                    col_types = cols("c", "d", "d"))
    genes <- unique(gsub("\\..*", "",  dpsi$X1))   
    return(genes)
})

names(gene_ls) <- filenames

if (database == "GO"){
    compareClusterSep(gene_ls, 
                        out.dir, 
                        orgdb,
                        ont     = ontology,
                        name    = "GO", 
                        keyType = gene_type,
                        pval    = pval, 
                        qval    = qval,
                        width   = args$width, 
                        height  = args$height,
                        format  = args$fmt)
}else if (database %in% c("KEGG", "Reactome")) {
    comparePathwayCluster(gene_ls, 
                        out.dir, 
                        orgdb, 
                        organism,
                        database,
                        name    = database, 
                        keyType = gene_type,
                        pval    = pval, 
                        qval    = qval,
                        width   = args$width, 
                        height  = args$height,
                        format  = args$fmt)
}

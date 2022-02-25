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
parser$add_argument("--clusterfile", help="cluster file")
parser$add_argument("--orgdb", help="orgdb")
parser$add_argument("--genetype", help="gene type")
parser$add_argument("--keggorganism", help="kegg organism")
parser$add_argument("--files",  nargs='+', help="dpsi files")
parser$add_argument("--xlabel",  nargs='+', help="dpsi files")

args <- parser$parse_args()

out.dir <- args$outdir
pval  <- args$pval
qval  <- args$qval
db <- args$database
cluster_file <- args$clusterfile
orgdb <- args$orgdb
gene_type <- args$genetype
kegg_organism <- args$keggorganism
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


if (file.exists(cluster_file)){

    cls_ext <- tools::file_ext(cluster_file)

    if (cls_ext == "tsv"){
        exon_cluster <- read_tsv(cluster_file)
    } else if (cls_ext == "csv") {
        exon_cluster <- read_csv(cluster_file)
    }

    exon_cluster$geneId <- gsub("\\..*", "",  exon_cluster$geneId) 

    unique_exon_cls <- exon_cluster %>% 
        distinct(cluster, geneId) 


    clusters <- sort(unique(unique_exon_cls$cluster))
    lapply(clusters, function(c){
        cluster_gene_ls <- lapply(dpsi_files, function(file){
            dpsi <- read_tsv(file, 
                             skip      = 1, 
                             col_names = F,
                             col_types = cols("c", "d", "d"))
            colnames(dpsi) = c("event_id", "dpsi", "pval")

            genes <- unique(gsub("\\..*", "",  dpsi$event_id)) 
            cluster_sig_genes <- unique_exon_cls %>% 
                filter(geneId %in% genes) %>% 
                filter(cluster == c)
            cluster_sig_genes$geneId
        })      

        names(cluster_gene_ls) <- filenames  

        tryCatch({
            if (db == "GO"){
                compareClusterSep(cluster_gene_ls, 
                                  out.dir, 
                                  orgdb,
                                  name    = paste("cluster", c, sep=""),
                                  keyType = gene_type,
                                  pval    = pval, 
                                  qval    = qval)                
            }else if (db == "KEGG") {
               compareKEGGCluster(cluster_gene_ls, 
                                  out.dir, 
                                  orgdb, 
                                  kegg_organism, 
                                  name=paste("KEGG_cluster", c, sep=""), 
                                  keyType = gene_type, 
                                  pval = pval, 
                                  qval = qval)
            }

            }, error = function(e) {
                print(e)
                next()
            }
          )
    })

} else {
    gene_ls <- lapply(dpsi_files, function(file){
        dpsi <- read_tsv(file, 
                        skip      = 1, 
                        col_names = F,
                        col_types = cols("c", "d", "d"))
        genes <- unique(gsub("\\..*", "",  dpsi$X1))   
        return(genes)
   })
    names(gene_ls) <- filenames
    if (db == "GO"){
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
    }else if (db == "KEGG") {
        compareKEGGCluster(gene_ls, 
                           out.dir, 
                           orgdb, 
                           kegg_organism,
                           name    = "KEGG", 
                           keyType = gene_type,
                           pval    = pval, 
                           qval    = qval,
                           width   = args$width, 
                           height  = args$height,
                           format  = args$fmt)
    }
}






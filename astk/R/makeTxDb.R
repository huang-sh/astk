
parser <- argparse::ArgumentParser()

parser$add_argument("--gtf", help="gtf")
parser$add_argument("--organism", help="organism")
parser$add_argument("--output", help="output")

args <- parser$parse_args()


txdb <- GenomicFeatures::makeTxDbFromGFF(
                        args$gtf,
                        format="gtf",
                        organism=args$organism,
                        dbxrefTag = "gene_name")


AnnotationDbi::saveDb(txdb, args$output)
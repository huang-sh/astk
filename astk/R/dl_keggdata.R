args <- commandArgs()

organism <- args[6]


script_path  <- str_split(args[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

prepare_KEGG(organism)

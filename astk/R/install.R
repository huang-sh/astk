
CRAN_url <- c("https://cloud.r-project.org/", 
        "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if (!("argparse" %in% rownames(installed.packages()))){
    options("repos" = CRAN_url)
    install.packages("argparse")
} 

parser <- argparse::ArgumentParser()
parser$add_argument("--requirement", action='store_true')
parser$add_argument("--OrgDb", nargs='+')
parser$add_argument("--CRAN", nargs='+')
parser$add_argument("--bioconductor", nargs='+')
parser$add_argument("--mirror", action='store_true')

args <- parser$parse_args()

if (args$mirror){
    options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
    options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
} else {
   options("repos" = CRAN_url)
}

if (args$requirement){
    ## Install packages from Cran

    cran.packages <- c("tidyverse", "ggplot2", "argparse", 
                       "UpSetR", "eoffice", "ggnewscale")

    for(i in cran.packages){
    if(!(i %in% rownames(installed.packages()))){
        message('Installing package: ', i)
        install.packages(i, dependencies = T)
    } else 
        next
    }
    ## Install packages from Bioconductor

    bioconductor.packages <- c('ComplexHeatmap', 'clusterProfiler', 'org.Mm.eg.db',
             'org.Hs.eg.db', 'simplifyEnrichment', 'universalmotif', 'memes',"tximport",
             "DESeq2", "metagene2")

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    for(i in bioconductor.packages){
        if(!(i %in% rownames(installed.packages()))){
            message('Installing package: ', i)
            BiocManager::install(i, dependencies = T)
    } else 
        next
    }    
}

if (! is.null(args$OrgDb)){

    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }else {
        for (i in args$OrgDb){
            if(!(i %in% rownames(installed.packages()))){
                message('Installing package: ', i)
                BiocManager::install(i)
            } else 
                next
        }
    }
}

if (! is.null(args$CRAN)){
    for (i in args$CRAN){
        if(!(i %in% rownames(installed.packages()))){
            message('Installing package: ', i)
            install.packages(i)
        } else 
            next
    }
    
}

if (! is.null(args$bioconductor)){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }
    for (i in args$bioconductor){
        if(!(i %in% rownames(installed.packages()))){
            message('Installing package: ', i)
            BiocManager::install(i)
        } else 
            next
    }    
}

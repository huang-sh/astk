args <- commandArgs()

requirement <- args[6]
OrgDb <- args[7]
cran <- args[8]
bio <- args[8]



if (requirement == "1"){
    ## Install packages of dependency
    ##----->> Install packages from Cran
    options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 

    cran.package.list <- c("tidyverse", "ggplot2")

    for(i in cran.package.list){
    if(!(i %in% rownames(installed.packages()))){
        message('Installing package: ', i)
        install.packages(i, dependencies = T)
    } else 
        next
    }
    ##----->> Install packages from Bioconductor
    options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

    bioconductor.package.list <- c('ComplexHeatmap', 'clusterProfiler',
            'org.Mm.eg.db', 'org.Hs.eg.db', 'simplifyEnrichment', 'memes')

    for(i in bioconductor.package.list){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        if(!(i %in% rownames(installed.packages()))){
            message('Installing package: ', i)
            BiocManager::install(i, dependencies = T)
    } else 
        next
    }    
}

if (OrgDb != "0"){
    options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }else {
       BiocManager::install(OrgDb)
    }
}

if (cran != "0"){
    install.packages(cran)
}

if (bio != "0"){
    options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }else {
       BiocManager::install(bio)
    }
    
}

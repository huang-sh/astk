suppressMessages(library(metagene2))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

parser <- argparse::ArgumentParser()
# parser$add_argument("--anchor", type="integer", nargs='+')
parser$add_argument("--binsize", type="integer")
parser$add_argument("--width", type="integer", nargs='+')
parser$add_argument("--bam", nargs='+')
parser$add_argument("--bamlabel", nargs='+')
parser$add_argument("--region", nargs='+')
parser$add_argument("--regionlabel", nargs='+')
parser$add_argument("--title")
parser$add_argument("--out")
parser$add_argument("--markmerge", action='store_true')
parser$add_argument("--normalmethod")
parser$add_argument("--pairedend", action='store_true')
parser$add_argument("--ASType")



args <- parser$parse_args()

if (args$ASType == "SE"){
    anchors <- seq(1, 4)
    names(anchors) <- paste("a", anchors, sep="")    
} else if (args$ASType == "RI") {
    anchors <- seq(2, 3)
    names(anchors) <- paste("a", anchors, sep="")    
}


binsize <- args$binsize
width <- args$width
if (length(width) == 1){
    width <- rep(width, length(anchors))
}
flank_width <- as.integer(width / 2)
bin_num <- seq(0, width[1], by=binsize)

bam_files  <-  args$bam
mark <- args$bamlabel
names(bam_files) <- mark
title <- args$title
outname <- args$out
markMerge <- args$markmerge
region_files <- args$region
region_name <- args$regionlabel
names(region_files) <- region_name

normalmethod <- args$normalmethod
peParam <-  args$pairedend
# anchors <- seq(1, 4)
# names(anchors) <- paste("a", anchors, sep="")
# binsize <- 5
# width <- c(300, 300, 300, 300)
# flank_width <- as.integer(width / 2)
# bin_num <- seq(0, width[1], by=binsize)
# bam_files  <-  c("/disk/public/hsh/fb_chip/H3K9ac.e16.5.fb_m.bam")
# mark <- "H3K9ac"
# names(bam_files) <- mark
# title <- "AS sites profile"
# outname <- "test"
# markMerge <- T


# region_files <- c("/home/huangshenghui/project/astk/demo3/diff/fb_e11_based/psi/fb_e11_16_SE_c2.0.8.psi",
#            "/home/huangshenghui/project/astk/demo3/diff/fb_e11_based/psi/fb_e11_16_SE_c2.0.2.psi")



# region_name <- c("high", "low")

names(region_files) <- region_name
coor_ls <- lapply(region_files, function(file) {
    df <- read.table(file, skip = 1)
    as_type  <- unique(sapply(str_split( df$V1, ";"), function(subs) { substr(subs[2], 1, 2) }))
    strand <- sapply(str_split( df$V1, ":"), function(subs) { subs[length(subs)] })
    chr <- sapply(str_split( df$V1, ":"), function(subs) { subs[2] })    
    if (as_type == "SE"){
        a12s <- lapply(str_split( df$V1, ":"), function(subs) { str_split(subs[3], "-")[[1]] })
        a1 <- sapply(a12s, function(x) {x[1]})
        a2 <- sapply(a12s, function(x) {x[2]})

        a34s <- lapply(str_split( df$V1, ":"), function(subs) { str_split(subs[4], "-")[[1]] })
        a3 <- sapply(a34s, function(x) {x[1]})
        a4 <- sapply(a34s, function(x) {x[2]})
        cdf <- data.frame(
            event_id = df$V1,
            chr = chr,
            a1 = a1,
            a2 = a2,
            a3 = a3,
            a4 = a4,
            strand = strand
        )
    }else if (as_type == "RI") {
        a2 <- sapply(str_split( df$V1, ":"), function(subs) { str_split(subs[4], "-")[[1]][1] })
        a3 <- sapply(str_split( df$V1, ":"), function(subs) { str_split(subs[4], "-")[[1]][2] })
        cdf <- data.frame(
            event_id = df$V1,
            chr = chr,
            a2 = a2,
            a3 = a3,
            strand = strand
        )
    }
    cdf
})

# data_ls <- lapply(anchors, function(i) {
data_ls <- lapply(seq(1, length(anchors)), function(i) {

    lapply(seq(1, length(coor_ls)) , function(dfi) {
        df <- coor_ls[[dfi]]
        anchor_pos <- as.numeric(df[, 2 + i])
        start <- anchor_pos - flank_width[i]
        end <- anchor_pos + flank_width[i]

        saf <- tibble(
            GeneID = df$event_id,
            Chr = df$chr,
            Start = start,
            End = end,
            Strand = df$strand
        )
        ps_saf <- saf[saf$Strand == "+", ]
        ns_saf <- saf[saf$Strand == "-", ]
        bin_num <- seq(0, width[i], by=binsize)
        ps_saf_ls <- lapply(bin_num[1:length(bin_num)-1], function(x) {          
                    saf_ <- ps_saf
                    bin_id <- paste("bin", x, sep = "")
                    saf_$GeneID <- paste(region_name[dfi], bin_id, "+", saf_$GeneID, sep="_")
                    saf_$Start <- saf_$Start + x
                    saf_$End <- saf_$Start + binsize -1
                    saf_
                })
        ns_saf_ls <- lapply(bin_num[1:length(bin_num)-1], function(x) {
                    saf_ <- ns_saf  
                    bin_id <- paste("bin", x, sep = "")
                    saf_$GeneID <- paste(region_name[dfi], bin_id, "-", saf_$GeneID, sep="_")
                    saf_$Start <- saf_$Start + x
                    saf_$End <- saf_$Start + binsize -1
                    saf_
                })

        ps_saf_df <- dplyr::bind_rows(ps_saf_ls)
        ns_saf_df <- dplyr::bind_rows(ns_saf_ls)        
        list(
            "+" = ps_saf_df,
            "-" = ns_saf_df
        )
    })
})

data_ls <- lapply(coor_ls, function(df) {
    # anchors
    lapply(seq(1, length(anchors)) , function(i) {
        anchor_pos <- as.numeric(df[, 2 + i])
        start <- anchor_pos - flank_width[i]
        end <- anchor_pos + flank_width[i]

        saf <- tibble(
            GeneID = df$event_id,
            Chr = df$chr,
            Start = start,
            End = end,
            Strand = df$strand
        )
        ps_saf <- saf[saf$Strand == "+", ]
        ns_saf <- saf[saf$Strand == "-", ]
        bin_num <- seq(0, width[i], by=binsize)
        ps_saf_ls <- lapply(bin_num[1:length(bin_num)-1], function(x) {          
                    saf_ <- ps_saf
                    bin_id <- paste("bin", x, sep = "")
                    saf_$GeneID <- paste(bin_id, saf_$GeneID, sep="_")
                    saf_$Start <- saf_$Start + x
                    saf_$End <- saf_$Start + binsize -1
                    saf_
                })

        ns_saf_ls <- lapply(bin_num[1:length(bin_num)-1], function(x) {
                    saf_ <- ns_saf  
                    bin_id <- paste("bin", x, sep = "")
                    saf_$GeneID <- paste(bin_id, saf_$GeneID, sep="_")
                    saf_$Start <- saf_$Start + x
                    saf_$End <- saf_$Start + binsize -1
                    saf_
                })

        ps_saf_df <- dplyr::bind_rows(ps_saf_ls)
        ns_saf_df <- dplyr::bind_rows(ns_saf_ls)        
        list(
            "+" = ps_saf_df,
            "-" = ns_saf_df
        )
    })
})


ps_count_ls <- lapply(data_ls, function(cdi) {
    
    lapply(cdi, function(df) {
        res <- Rsubread::featureCounts(
            bam_files,
            annot.ext = df[["+"]],
            isPairedEnd = peParam,
            nthreads = 4,
            countMultiMappingReads=T,
            allowMultiOverlap=T
        ) 
        lapply(seq(1, dim(res$counts)[2]), function(m) {
            mtx <- matrix(res$counts[, m], ncol=length(bin_num)-1)
            if (normalmethod == "CPM"){
                reads_sum <- sum(res$stat[, m+1])
                mtx <- mtx * 1e6 / reads_sum                
            }
            mtx

        })        
    })
})

ns_count_ls <- lapply(data_ls, function(cdi) {
    lapply(cdi, function(df) {
        res <- Rsubread::featureCounts(
            bam_files,
            annot.ext = df[["-"]],
            isPairedEnd = peParam,
            nthreads = 4,
            countMultiMappingReads=T,
            allowMultiOverlap=T
        )
        lapply(seq(1, dim(res$counts)[2]), function(m) {
            mtx <- matrix(res$counts[, m], ncol=length(bin_num)-1)
            if (normalmethod == "CPM"){
                reads_sum <- sum(res$stat[, m+1])
                mtx <- mtx * 1e6 / reads_sum                
            }
            mtx
        })
    })
})

mark_idx <- seq(1, length(bam_files))
names(mark_idx) <- mark

save.image("/home/huangshenghui/project/astk/dev/tmp/epi.RData")

count_ls <- lapply(mark_idx, function(m) {
    region_count_ls <- lapply(seq(1, length(ps_count_ls)), function(ri) {
        # anchors
        lapply(seq(1, length(anchors)), function(ai) {
            ps_mtx <- ps_count_ls[[ri]][[ai]][[m]]
            ns_mtx <- ns_count_ls[[ri]][[length(anchors)-ai+1]][[m]][, seq(length(bin_num)-1, 1)]
            print(dim(ps_mtx))
            print(dim(ns_mtx))
            mtx <- rbind(ps_mtx, ns_mtx)
            # if (is.null(ns_mtx)){
            #     mtx <- rbind(ps_mtx, matrix(ns_mtx))
            # } else{
            #     mtx <- rbind(ps_mtx, ns_mtx)
            # }
            
            bins <- paste("bin", seq(1, dim(mtx)[2]), sep = "")
            colnames(mtx) <- paste(mark[m], region_name[ri], names(anchors)[ai], bins, sep = "_")
            mtx
        })
    })
    names(region_count_ls) <- names(ps_count_ls)
    region_count_ls
    # lapply(rapply(count_ls, enquote, how="unlist"), eval)
})


if (markMerge){

    count_ls1 <- unlist(unlist(count_ls, recursive = F), recursive = F)
    over <- lapply(seq(1, length(count_ls[[1]])), function(x) {
            region_count_ls <- lapply(count_ls, function(mcl) {
                mcl[[x]]
            })
            mark_count_ls <- unlist(region_count_ls, recursive = F)
            data <- Reduce(cbind, mark_count_ls)
            outname <- sprintf("%s_%s.csv", outname, names(count_ls[[1]])[x])
            write.csv(data, file=outname, row.names=F)
    })     
}else {
    over <- lapply(seq(1, length(count_ls)), function(mi) {
       lapply(seq(1, length(count_ls[[mi]])), function(ri){
           df <- Reduce(cbind, count_ls[[mi]][[ri]])
           outname <- sprintf("%s_%s_%s.csv", outname, names(count_ls)[mi], names(count_ls[[mi]])[ri])
           write.csv(df, file=outname, row.names=F)

       })
   })
}

# save.image("/home/huangshenghui/project/astk/dev/diff/test.RData")

# load("/home/huangshenghui/project/astk/dev/diff/test.RData")

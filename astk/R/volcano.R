suppressMessages(library(stringr))
suppressMessages(library(ggplot2))


script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))

parser <- fig_cmd_parser()

parser$add_argument("--file", help="meme file path")

args <- parser$parse_args()


dpsi <- read.table(args$file,  header=T, row.names=1, sep="\t")
dat <- na.omit(dpsi)

dat$type=str_split(rownames(dat),'[;:]',simplify = T)[,2]
event = dat$type[1]
dat$threshold = factor(ifelse(dat[,2] < 0.05 & abs(dat[,1]) >= 0,
                              ifelse(dat[,1] >= 0  ,'Up','Down'),'stable'),
                             levels=c('Up','Down','stable'))

colnames(dat)[1:2] <- c('dPSI', 'p.val')
dat[,'p.val'] = -log10(dat[, 'p.val'])

P_dpsi <- ggplot(dat,aes(dPSI, p.val, color = threshold))+
        geom_point() +
        theme_bw() +
        ylab(paste0(event,' pval')) + #修改y轴名称
        xlab(paste0(event,' dPSI')) + #修改x轴名称
        geom_vline(xintercept = c(-0.1,0.1), lty=3, col="black", lwd=0.5) + #添加横线|FoldChange|>2
        geom_hline(yintercept = 0.05, lty=5, col="black", lwd=0.5) 


save_fig(P_dpsi, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)
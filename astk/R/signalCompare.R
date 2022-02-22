
suppressMessages(library(tidyverse))
suppressMessages(library(ggsignif))
suppressMessages(library(ggplot2))


script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", help="file path")

args <- parser$parse_args()

signal_df <- read_csv(args$file)


filter.signal <- signal_df %>% 
    group_by(event_idx, mark, direction, anchor) %>% 
    summarise(value = mean(signal)+0.1) %>% 
    mutate(anchor_d = paste(anchor, direction, sep = "_")) %>% 
    mutate(anchor_s = paste("a", anchor, sep = ""))
    
print(unique(filter.signal$anchor_s))
p <- ggplot(filter.signal, aes(x=anchor_s, y=value, group=anchor_s)) + 
  # scale_y_continuous(trans='log10') +
  geom_boxplot(aes(fill=anchor_s)) +
  geom_signif(test="wilcox.test", step_increase = 0.2,
              comparisons = combn(unique(filter.signal$anchor_s),2, simplify = F)) +
  facet_grid(~mark+direction) 


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)

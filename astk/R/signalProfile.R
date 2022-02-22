suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyverse))

script_path  <- stringr::str_split(commandArgs()[4], "=")[[1]][2]
source(file.path(dirname(script_path), "utils.R"))


parser <- fig_cmd_parser()

parser$add_argument("--file", nargs='+', help="file names")


args <- parser$parse_args()

signal_df <- read_csv(args$file, show_col_types = FALSE)
print(head(signal_df))


max_up_bin <- max(signal_df[signal_df$direction == "upstream", "bin"])

data_new <- signal_df %>%                               # Replacing values
  mutate(bin = replace(bin, direction == "upstream", bin - max_up_bin-1)) 
  

f.data <- data_new %>% 
    group_by(event_idx, mark, direction, anchor, bin) %>% 
    summarise(value = mean(signal)) %>% 
    group_by(mark, direction, anchor, bin) %>% 
    summarise(value = mean(value)) %>% 
    arrange(mark)

p <- ggplot(data = f.data) +
        geom_line(mapping = aes(x = bin, y = value, color = mark))+
        facet_wrap(~ anchor, nrow = 1)


save_fig(p, 
        args$output, 
        format = args$fmt,
        width  = args$width, 
        height = args$height, 
        units  = "in",
        res    = args$resolution)


    

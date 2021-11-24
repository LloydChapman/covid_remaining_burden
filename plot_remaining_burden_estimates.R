library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)

source("./R/plot_remaining_burden.R")

# Date of fitting - change according to run date
date_fitting = as.Date("2021-11-03")

# Set output directory
dir_out = paste0("./output/",date_fitting,"/")

# Set figure directory
dir_fig = "./figs/"

# Read in maximum remaining burden estimates
rem_burden_dt = fread(paste0(dir_out,"rem_burden_output.csv"))
ovrl_rem_burden_dt = fread(paste0(dir_out,"ovrl_rem_burden_output.csv"))

# Set plot theme
theme_set(theme_cowplot(font_size = 12) + theme(
    strip.background = element_blank(),
    plot.background = element_rect(fill="white"),
    legend.background = element_rect(fill="white"),
    panel.background = element_rect(fill="white")))

# Plot maximum remaining burden estimates
plot_remaining_burden(rem_burden_dt,ovrl_rem_burden_dt,dir_fig)

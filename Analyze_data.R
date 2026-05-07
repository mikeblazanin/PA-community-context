## Setup ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gcplyr)
library(GGally)
library(lme4)
library(ggrepel)

dir.create("plots", showWarnings = FALSE)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

dup_rows <- function(mydf, mut_col, from_vals, to_vals) {
  #this function takes in a data.frame
  #finds all rows where df[, mut_col] == from_vals
  #deletes those rows
  #and adds rows that replicate the original but with to_vals in the
  # mut_col instead
  
  mydf <- as.data.frame(mydf)
  
  #Generate vector of rows to duplicate
  rowsrep <- ifelse(mydf[, mut_col] %in% from_vals,
                    length(to_vals), 1)
  idx <- rep(1:nrow(mydf), rowsrep)
  
  #Duplicate rows
  newdf <- mydf[idx, ]
  
  #replace from_vals with to_vals in duplicated rows
  newdf[which(idx %in% which(rowsrep > 1)), mut_col] <- to_vals
  
  return(newdf)
}

mynormalize <- function(x, a_input, b_input, a_output, b_output) {
  #normalize a vector of data, x, by scaling values such that 
  #a_input becomes a_output and b_input becomes b_output
  stopifnot(b_input > a_input, b_output > a_output)
  
  return((b_output-a_output)*(x-a_input)/(b_input - a_input) + a_output)
}

## Prelim experiment: population data from run without competitors ----
prelimdens_NC_sum <- read.csv("./clean_data/prelimdens_NC_sum.csv",
                               stringsAsFactors = FALSE)

#Summarize pops for treat-level means
prelimdens_NC_sum_sum <- 
  summarize(group_by(prelimdens_NC_sum, Comm_treat, Phage_treat, 
                     Timepoint, Strain),
            dens_ml = 10**mean(log10(dens_ml+1)), freq_L = mean(freq_L))

#Plots
png("plots/prelimdens_NC.png", width = 5, height = 3.5, units = "in", res = 150)
ggplot(data = prelimdens_NC_sum,
       aes(x = Timepoint, y = dens_ml+1, color = Strain)) +
  #geom_point() + 
  geom_line(aes(group = paste(Population, Strain)), alpha = 0.5) +
  #geom_line(data = prelimdens_NC_sum_sum, lwd = 2, alpha = 0.8) +
  geom_hline(data = data.frame(Comm_treat = c("-", "-"),
                               Phage_treat = c("+", "-", "+", "-"),
                               limit = c(1000/5, 1000/5)),
             aes(yintercept = limit), lty = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("-" = "Only Pseudomonas"))) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlim(0, NA) +
  labs(y = "Density (cfu/mL or pfu/mL)") +
  scale_color_manual(limits = c("PA", "H6"),
                     labels = c("P. aeruginosa", "Phage H6"),
                     values = my_cols[c(3, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black"),
        legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(lwd = 2, alpha = 1))) +
  NULL
dev.off()

png("plots/prelimdens_NC_cfumorphfreq.png", width = 5, height = 2.5, units = "in", res = 150)
ggplot(data = prelimdens_NC_sum[prelimdens_NC_sum$Strain == "PA", ],
       aes(x = Timepoint, y = freq_L)) +
  geom_line(aes(group = paste(Population, Strain)), alpha = 0.5, lwd = 0.8) +
  #geom_line(data = prelimdens_NC_sum_sum[prelimdens_NC_sum_sum$Strain == "PA", ], lwd = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("-" = "Only Pseudomonas"))) +
  xlim(0, NA) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Frequency of Large\nColony Morphology") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text =  element_text(size = 12.5)) +
  NULL
dev.off()

## Prelim experiment: population data from run with competitors ----
prelimdens_C <- read.csv("./clean_data/prelimdens_C.csv", stringsAsFactors = FALSE)

#Summarize pops for treat-level means
prelimdens_C_sum <- 
  summarize(group_by(prelimdens_C, Comm_treat, Phage_treat, 
                     Timepoint, Strain),
            dens_ml = 10**mean(log10(dens_ml+1)), freq_L = mean(freq_L))

#Plots
png("plots/prelimdens_C.png", width = 5, height = 3.5, units = "in", res = 150)
ggplot(data = prelimdens_C,
       aes(x = Timepoint, y = dens_ml+1, color = Strain)) +
  #geom_point() + 
  geom_line(aes(group = paste(Population, Strain)), alpha = 0.5) +
  #geom_line(data = prelimdens_C_sum, lwd = 2, alpha = 0.8) +
  geom_hline(data = data.frame(Comm_treat = c("+", "+"),
                               Phage_treat = c("+", "-", "+", "-"),
                               limit = c(1000/2.5, 1000/2.5)),
             aes(yintercept = limit), lty = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "All Bacteria"))) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlim(0, NA) +
  labs(y = "Density (cfu/mL or pfu/mL)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black"),
        legend.position = "top") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(lwd = 2, alpha = 1))) +
  NULL
dev.off()

png("plots/prelimdens_C_cfumorphfreq.png", width = 5, height = 2.5, units = "in", res = 150)
ggplot(data = prelimdens_C[prelimdens_C$Strain == "PA", ],
       aes(x = Timepoint, y = freq_L)) +
  geom_line(aes(group = paste(Population, Strain)), alpha = 0.5, lwd = 0.8) +
  #geom_line(data = prelimdens_C_sum[prelimdens_C_sum$Strain == "PA", ], lwd = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "All Bacteria"))) +
  xlim(0, NA) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Frequency of Large\nColony Morphology") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text =  element_text(size = 12.5)) +
  NULL
dev.off()

## Population data ----
density <- read.csv("clean_data/density.csv", stringsAsFactors = FALSE)

#Change factor orders
density <- mutate(density,
                  treat_C = factor(treat_C, levels = c("NC", "C")),
                  treat_P = factor(treat_P, levels = c("NP", "P")))

#Write to file for other scripts to use
saveRDS(select(density,
               Timepoint, Pop, Strain, treat_P, treat_C,
               cfu_ml, bd, log_cfuml),
        file = "density.rds")

density_sum <- summarize(group_by(density,
                              Timepoint, treat_P, treat_C, Strain),
                     cfu_ml = 10**mean(log10(cfu_ml), na.rm = TRUE),
                     L_freq = mean(L_freq, na.rm = TRUE))

p1 <- ggplot(
  data = filter(density, Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
       aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                             "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                             "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black"))

png("plots/fig1_density_legendside.png", width = 7, height = 5, units = "in",
    res = 150)
p1  +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1))) +
  # geom_text(data = filter(density, Timepoint == 9, Strain == "H6"),
  #           aes(y = cfu_ml, label = Pop)) +
  NULL
dev.off()

png("plots/density.png", width = 5, height = 5, units = "in",
    res = 150)
p1 +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1),
                              nrow = 2, byrow = TRUE))
dev.off()

png("plots/density_bact.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain != "H6", Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                                "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/density_PA.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain == "PA", Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                             "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                             "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/density_otherbact.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain %in% c("AX", "EF", "SA"), 
                Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                             "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                             "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(4, 6, 8),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/density_phage.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain == "H6", Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = as.factor(Pop))) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                                "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  # scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
  #                    labels = c("P. aeruginosa", "S. aureus",
  #                               "E. faecalis", "A. xylosoxidans",
  #                               "Phage H6"),
  #                    values = my_cols[c(3, 1, 2, 7, 8)]) +
  scale_color_manual(values = my_cols) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/density_phage_Ponly.png", width = 3.5, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain == "H6", Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP"), treat_P == "P"),
  aes(x = Timepoint, y = cfu_ml)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                             "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                             "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Phage Density (pfu/mL)", x = "Timepoint (days)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/density_phageemphasized.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Pop %in% 3:10, !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain), lwd = Strain), alpha = 0.6) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria Added",
                                                "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                     labels = c("P. aeruginosa", "S. aureus",
                                "E. faecalis", "A. xylosoxidans",
                                "Phage H6"),
                     values = my_cols[c(3, 1, 2, 7, 8)]) +
  scale_linewidth_manual(breaks = c("AX", "EF", "PA", "SA", "H6"),
                         values = c(.5, .5, .5, .5, 1.5)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)),
         linewidth = "none")
dev.off()

png("plots/density_NC_treats.png", width = 7, height = 3, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain %in% c("PA", "H6"), Pop %in% 3:10, 
                treat_C == "NC", !is.na(cfu_ml),
                !(Pop == 9 & treat_P == "NP")),
  aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  #geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(treat_P ~ .,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                                "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = "Density (cfu/mL or pfu/mL)", x = "Timepoint (days)") +
  scale_color_manual(limits = c("PA", "H6"),
                     labels = c("P. aeruginosa", "Phage H6"),
                     values = my_cols[c(3, 8)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 12.5),
        #legend.position = c(0.3, 0.6),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(lwd = 2, alpha = 1)))
dev.off()

png("plots/cfumorphfreq.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = filter(density, Pop %in% 3:10, Strain == "PA",
                     !(Pop == 9 & treat_P == "NP")),
       aes(x = Timepoint, y = L_freq)) +
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5, lwd = 0.8) +
  #geom_line(data = density_sum[density_sum$Strain == "PA", ], lwd = 2) +
  facet_grid(treat_C ~ treat_P,
             labeller = labeller(treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                 treat_C = c("C" = "All Bacteria",
                                                "NC" = "Only Pseudomonas"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                     limits = c(0, NA)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Frequency of Large Colony Morphology") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text =  element_text(size = 12.5)) +
  NULL
dev.off()

# Stats
ggplot(data = filter(density, Timepoint > 1, Strain == "H6",
                     Pop %in% 3:10, !is.na(cfu_ml)),
       aes(x = Timepoint, y = cfu_ml, color = treat_C, 
           group = paste(treat_P, treat_C, Pop, Strain))) +
  geom_line() +
  scale_y_log10()

phagedens_lm <- lm(log_cfuml ~ as.factor(Timepoint) + as.factor(treat_C),
                   data = filter(density, Timepoint > 1, Strain == "H6",
                                 Pop %in% 3:10, !is.na(cfu_ml)))
summary(phagedens_lm)


# Resistance data ----
resis_isols <- read.csv("./clean_data/resis_isols.csv", stringsAsFactors = FALSE)

resis_isols <- mutate(as_tibble(resis_isols),
                      Bact_treat_C = as.factor(Bact_treat_C),
                      Bact_treat_P = as.factor(Bact_treat_P),
                      Bact_treat = as.factor(Bact_treat))

#At this point, all the NA EOP2_plot values are for one of:
# Phages that don't exist bc extinction
# Bact Pop 3 T0 #3 which we couldn't recover from frozen
# NP treatments

#Calculate average for ea bact isol across phages
resis_isols_bactsum <- summarize(
  group_by(resis_isols,
           Bact_treat, Bact_treat_C, Bact_treat_P, Bact_timepoint,
           Bact_pop, Bact_isol),
  clearing_avg = mean(clearing2, na.rm = TRUE),
  EOP_avg = mean(EOP2_plot, na.rm = TRUE), #(the only NA's are when phage are extinct)
  Plaque_morphology = names(sort(table(Plaque_morphology), decreasing = TRUE)[1]))
resis_isols_bactsum <- mutate(
  group_by(resis_isols_bactsum,
           Bact_treat, Bact_treat_C, Bact_treat_P, Bact_timepoint, Bact_pop),
  Bact_isol_forplotting = rank(EOP_avg, ties.method = "first", na.last = "keep"),
  Bact_isol_forplotting = 
    if(any(!is.na(Bact_isol_forplotting))) {
      1 + max(Bact_isol_forplotting, na.rm = TRUE) - Bact_isol_forplotting
      } else {NA})
  
resis_isols <-
  left_join(resis_isols,
            select(resis_isols_bactsum,
                   Bact_treat, Bact_treat_C, Bact_treat_P, Bact_timepoint,
                   Bact_pop, Bact_isol,
                   Bact_isol_forplotting))

#The columns that contain useful data in resis_isols/resis_isols_bactsum:
# EOP2_plot/EOP_avg 
#                     EOP data with bd data (NA's in EOP2) filled in as the lowest
#                     EOP in the dataset
# clearing2/clearing_avg
#                     binary of whether EOP is above detection or not
# Plaque_morphology
#                     categorical that works well for resistance categories
#                     'Clear' are susceptible, 'Clear?' are NP treats we never
#                     went back to check, 'Cloudy' or 'Tiny' are part-resis
#                     'None' are resis
#                     Note: this works even better than EOP, because EOP for
#                     clear vs cloudy strongly overlap, and even average EOP
#                     for ea bact isol overlaps slightly between clear and
#                     cloudy/tiny

resis_isols_bactsum_undup <- summarize(
  group_by(mutate(resis_isols_bactsum, 
         Bact_treat = ifelse(Bact_timepoint == 0, 
                             "Anc", 
                             as.character(Bact_treat)),
         Bact_treat_C = ifelse(Bact_timepoint == 0, 
                               "Anc", 
                               as.character(Bact_treat_C)),
         Bact_treat_P = ifelse(Bact_timepoint == 0, 
                               "Anc", 
                               as.character(Bact_treat_P))),
         Bact_treat, Bact_treat_C, Bact_treat_P, Bact_timepoint,
         Bact_pop, Bact_isol),
  clearing_avg = mean(clearing_avg), 
  EOP_avg = mean(EOP_avg), 
  Plaque_morphology = names(sort(table(Plaque_morphology), decreasing = TRUE)[1]))


#Calculate average across pops
resis_sum_treats <- summarize(group_by(resis_isols,
                                Bact_treat, Bact_treat_C,
                                Bact_treat_P, Bact_timepoint,
                                Phage_treat, Phage_treat_C,
                                Phage_treat_P, Phage_timepoint),
                       clearing_avg = mean(clearing, na.rm = TRUE))

#Calculate average for each pop
resis_sum_pops <- summarize(group_by(resis_isols,
                                     Bact_pop, Phage_pop,
                                     Bact_treat, Bact_treat_C,
                                     Bact_treat_P, Bact_timepoint,
                                     Phage_treat, Phage_treat_C,
                                     Phage_treat_P),
                            clearing_avg = mean(clearing2, na.rm = TRUE))
resis_sum_pops_sum <- summarize(group_by(resis_sum_pops,
                                         Bact_treat, Bact_treat_C,
                                         Bact_treat_P, Bact_timepoint),
                                clearing_avg = mean(clearing_avg))


#Run stats
resfreqlm <- lm(data = filter(resis_sum_pops, Bact_timepoint != 0, Bact_treat_P == "P"),
   clearing_avg ~ Bact_treat_C + Bact_timepoint:Bact_treat_C)
summary(resfreqlm)
anova(resfreqlm)

ggplot(resis_sum_pops,
       aes(x = Bact_timepoint, y = clearing_avg, color = Bact_treat)) +
  geom_point()

#Plots of EOP by plaque morphology
ggplot(data = mutate(resis_isols,
                     Plaq_morph = ifelse(Bact_timepoint == 0, "Anc", Plaque_morphology)), 
       aes(x = Plaq_morph, y = EOP2_plot)) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) + 
  scale_y_log10()

ggplot(data = mutate(resis_isols_bactsum_undup,
                     Plaq_morph = ifelse(Bact_timepoint == 0, "Anc", Plaque_morphology)), 
       aes(x = Plaq_morph, y = EOP_avg)) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) + 
  scale_y_log10() +
  geom_hline(yintercept = 0.5, lty = 2)

summarize(
  group_by(
    mutate(resis_isols_bactsum_undup, 
           Plaq_morph = ifelse(Bact_timepoint == 0, "Anc", Plaque_morphology)), 
    Plaq_morph), 
  max_eop = max(EOP_avg, na.rm = TRUE), 
  min_eop = min(EOP_avg, na.rm = TRUE), 
  avg_eop = mean(EOP_avg, na.rm = TRUE))

#Check which combos were done in which runs (spots, cross-streaks, EOP, EOP2)
agplots <- list()
for(runcol in grep("run", colnames(resis_isols), value = TRUE)) {
  agplots[[length(agplots)+1]] <-
    ggplot(data = resis_isols,
           aes(x = Bact_isol, 
               y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
    geom_tile(aes(fill = as.factor(.data[[runcol]]))) +
    scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
    # scale_x_discrete(labels = rep(c(0, "", 1, "", 2, "", 4, "", 8, ""),
    #                               times = c(1, 5, 1, 5, 1, 5, 1, 5, 1, 5))) +
    theme_bw() +
    theme(strip.text = element_text(size = 8),
          axis.text.y = element_text(vjust = -0.1),
          #axis.text.x = element_text(hjust = 0.2),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          panel.spacing.x = unit(.05, "lines")) +
    scale_fill_discrete(breaks = c(TRUE, NA)) +
    ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                        labeller = labeller(
                          Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                          Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"))) +
    labs(x = "Bacterial Isolate", y = "Phage timepoint",
         subtitle = "Population\nBacterial timepoint") +
    NULL
}

png("plots/resis_assay_grid.png", width = 8, height = 20, units = "in", res = 150)
cowplot::plot_grid(
  plotlist = agplots,
  nrow = 4)
dev.off()


png("plots/resis_grid.png", width = 8, height = 5, units = "in", res = 150)
ggplot(data = resis_isols,
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing2))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  # scale_x_discrete(labels = rep(c(0, "", 1, "", 2, "", 4, "", 8, ""),
  #                               times = c(1, 5, 1, 5, 1, 5, 1, 5, 1, 5))) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_manual(breaks = c(1, 0, NA),
                    labels = c("Susceptible", "Resistant", "Phage Extinct"),
                    values = c("#132B43", "#56B1F7", "gray"),
                    name = "") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Population\nBacterial timepoint") +
  NULL
dev.off()

png("plots/resis_grid_phagetreats.png", width = 10, height = 3.5, 
    units = "in", res = 150)
ggplot(data = filter(resis_isols, Bact_treat_P == "P"),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing2))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  scale_x_continuous(breaks = 1:6) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_manual(breaks = c(1, 0, NA),
                    labels = c("Susceptible", "Resistant", "Phage Extinct"),
                    values = c("#132B43", "#56B1F7", "gray"),
                    name = "") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Replicate population") +
  NULL
dev.off()

png("plots/resis_freq_time.png", width = 8, height = 3.5, units = "in", res = 150)
ggplot(data = resis_sum_pops,
       aes(x = as.factor(Bact_timepoint), y = 1-clearing_avg,
           color = Bact_treat_P)) +
  geom_point(position = position_jitter(width = 0.1, height = 0),
             size = 3, alpha = 0.5) +
  geom_line(data = resis_sum_pops_sum,
            aes(group = Bact_treat), linewidth = 2) +
  scale_color_manual(name = "Phage Treatment",
                     breaks = c("NP", "P"),
                     labels = c("No Phage", "Phage Added"),
                     values = my_cols[c(2, 6)]) +
  facet_grid(~ Bact_treat_C,
             labeller = labeller(
               Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
               Bact_treat_C = c(NC = "Only Pseudomonas", C = "All Bacteria"))) +
  labs(x = "Bacterial timepoint", y = "Frequency of\nPhage Resistance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 11.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

png("plots/resis_freq_time_NConly.png", width = 7, height = 4, units = "in", res = 150)
ggplot(data = filter(resis_sum_pops, Bact_treat_C == "NC"),
       aes(x = as.factor(Bact_timepoint), y = 1-clearing_avg,
           color = Bact_treat_P)) +
  geom_point(position = position_jitter(width = 0.1),
             size = 3, alpha = 0.5) +
  geom_line(data = filter(resis_sum_pops_sum, Bact_treat_C == "NC"),
            aes(group = Bact_treat), linewidth = 2) +
  scale_color_manual(name = "Phage Treatment",
                     breaks = c("NP", "P"),
                     labels = c("No Phage", "Phage Added"),
                     values = my_cols[c(2, 6)]) +
  labs(x = "Bacterial timepoint", y = "Bacterial Resistance\nto Phage") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 11.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

png("plots/EOP_grid_phagetreats.png", width = 10, height = 3.5, 
    units = "in", res = 300)
ggplot(data = filter(resis_isols, Bact_treat_P == "P"),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = log10(EOP2_plot))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  scale_x_continuous(breaks = 1:6) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_continuous(name = "log10(Susceptibility)",
                        type = "viridis") +
                        #high = "#132B43", low = "#56B1F7") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Replicate population") +
  NULL
dev.off()

png("plots/fig2_EOP_grid_phagetreats_reordered.png", width = 10, height = 3.5, 
    units = "in", res = 300)
ggplot(data = mutate(filter(resis_isols, Bact_treat_P == "P"),
                     Phage_pop = LETTERS[as.numeric(Phage_pop)-2]),
       aes(x = Bact_isol_forplotting, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = log10(EOP2_plot))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  scale_x_continuous(breaks = 1:6) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_continuous(name = "log10(Susceptibility)",
                        type = "viridis") +
  #high = "#132B43", low = "#56B1F7") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint") +
  NULL
dev.off()

mytimepoints <- list(c(0, 1), c(2, 4), c(8))
for(i in 1:length(mytimepoints)) {
  png(paste0("plots/EOP_grid_phagetreats2_reordered_split", i, ".png"), 
      width = ifelse(i == 3, 3, 5.3), height = 3, units = "in", res = 300)
  print(
    ggplot(data = 
             mutate(filter(resis_isols, Bact_treat_P == "P", 
                           Bact_timepoint %in% mytimepoints[[i]]),
                    Phage_pop = LETTERS[as.numeric(Phage_pop)-2]),
           aes(x = Bact_isol_forplotting, 
               y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
      geom_tile(aes(fill = log10(EOP2_plot))) +
      scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
      scale_x_continuous(breaks = 1:6) +
      theme_bw() +
      theme(strip.text = element_text(size = 8),
            axis.text.y = element_text(vjust = -0.1),
            #axis.text.x = element_text(hjust = 0.2),
            axis.text.x = element_blank(),
            legend.position = "bottom",
            plot.subtitle = element_text(size = 9),
            panel.spacing.x = unit(.05, "lines")) +
      scale_fill_continuous(name = "log10(Susceptibility)",
                            type = "viridis") +
      #high = "#132B43", low = "#56B1F7") +
      ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                          labeller = labeller(
                            Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                            Bact_treat_C = c(NC = "Only\nPseud.", C = "All Bacteria"),
                            Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
      labs(x = "Bacterial Isolate", y = "Phage timepoint") +
      NULL)
  dev.off()
}

png("plots/Plaquemorph_grid_phagetreats.png", width = 10, height = 3.5, 
    units = "in", res = 150)
ggplot(data = mutate(
  filter(resis_isols, Bact_treat_P == "P"),
  Plaque_morphology = ifelse(Plaque_morphology %in% c("Tiny", "Cloudy"),
                             "Tiny/Cloudy", Plaque_morphology)),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(Plaque_morphology))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  # scale_x_discrete(labels = rep(c(0, "", 1, "", 2, "", 4, "", 8, ""),
  #                               times = c(1, 5, 1, 5, 1, 5, 1, 5, 1, 5))) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_manual(
    breaks = c("Clear", "Clear?", "Tiny/Cloudy", "None", NA),
    labels = c("Clear", "Unconfirmed", "Tiny/Cloudy", "None", "Phage Extinct"),
    values = c(scales::seq_gradient_pal(low = "#132B43", high = "#56B1F7")(
      seq(0, 1, length.out = 4)), "gray"),
    name = "Plaque Morphology") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Population") +
  guides(fill = guide_legend(nrow = 1)) +
  NULL
dev.off()

png("plots/Plaquemorph_grid_phagetreats_reordered.png", width = 10, height = 3.5, 
    units = "in", res = 150)
ggplot(data = mutate(
  filter(resis_isols, Bact_treat_P == "P"),
  Plaque_morphology = ifelse(Plaque_morphology %in% c("Tiny", "Cloudy"),
                             "Tiny/Cloudy", Plaque_morphology),
  Phage_pop = LETTERS[as.numeric(Phage_pop)-2]),
  aes(x = Bact_isol_forplotting, 
      y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(Plaque_morphology))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  # scale_x_discrete(labels = rep(c(0, "", 1, "", 2, "", 4, "", 8, ""),
  #                               times = c(1, 5, 1, 5, 1, 5, 1, 5, 1, 5))) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        #axis.text.x = element_text(hjust = 0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(.05, "lines")) +
  scale_fill_manual(
    breaks = c("Clear", "Tiny/Cloudy", "None", NA),
    labels = c("Clear", "Tiny/Cloudy", "None", "Phage Extinct"),
    values = c(scales::seq_gradient_pal(low = "#132B43", high = "#56B1F7")(
      seq(0, 1, length.out = 3)), "gray"),
    name = "Plaque Morphology") +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Bact_timepoint * Phage_pop,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint") +
  guides(fill = guide_legend(nrow = 1)) +
  NULL
dev.off()

ggplot(data = filter(resis_isols_bactsum, Bact_treat_P == "P"),
       aes(x = as.factor(Bact_timepoint), y = EOP_avg, group = Bact_pop)) +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.3) +
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ .,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"))) +
  labs(y = "log10(Susceptibility)") +
  scale_y_log10() +
  theme_bw()

resis_isols_phagesum <-
  summarize(group_by(resis_isols, 
                     Phage_treat, Phage_timepoint, Phage_pop, Phage_isol),
            maxpfu = ifelse(any(!is.na(c(pfu_ml, pfu_ml2))),
                            max(pfu_ml, pfu_ml2, na.rm = TRUE),
                            NA))
min(resis_isols_phagesum$maxpfu, na.rm = TRUE)

#Plot phage titers to see how good each phage stock was
ggplot(mutate(group_by(resis_isols,
                       Phage_treat, Phage_timepoint, Phage_pop, Phage_isol,
                       Bact_treat, Bact_timepoint, Bact_pop, Bact_isol),
              pfuplot = 
                ifelse(any(!is.na(c(pfu_ml[bd != TRUE], pfu_ml2[bd2 != TRUE]))),
                       max(pfu_ml[bd != TRUE], pfu_ml2[bd2 != TRUE], na.rm = TRUE),
                       NA)),
       aes(x = Phage_pop, y = pfuplot, 
           color = as.factor(Phage_isol), shape = as.factor(Plaque_morphology))) +
  geom_point() +
  ggh4x::facet_nested(. ~ Phage_treat * Phage_timepoint) +
  scale_y_log10()

ggplot(resis_isols_phagesum,
       aes(x = Phage_timepoint, y = maxpfu)) +
  geom_point() +
  ggh4x::facet_nested(. ~ Phage_treat) +
  scale_y_log10()

# Pigments ----
pig <- read.csv("./clean_data/pig.csv", stringsAsFactors = FALSE,
                 colClasses = "character")
pig <- mutate(
  pig,
  Bact_treat = factor(Bact_treat, levels = c("Anc", "NN", "CNP", "NCP", "CP")),
  Bact_treat_C = factor(Bact_treat_C, levels = c("Anc", "NC", "C")),
  Bact_treat_P = factor(Bact_treat_P, levels = c("Anc", "NP", "P")),
  Measurements = as.numeric(Measurements),
  Meas_norm = as.numeric(Meas_norm))
pig <- group_by(pig, file, Wavelength)

#Do batch correction
pig$Meas_norm_norm <- NA
for(wave in unique(pig$Wavelength)) {
  mydat <- pig[pig$Wavelength == wave & !is.na(pig$Meas_norm), ]
  temp_lmer <- lmer(
    data = mydat,
    Meas_norm ~ (1|file) + Bact_treat + Bact_timepoint +
      (1|Bact_treat:Bact_pop) + (1|Bact_treat:Bact_pop:Bact_isol))
  pig$Meas_norm_norm[pig$Wavelength == wave & !is.na(pig$Meas_norm)] <-
    mydat$Meas_norm -
    predict(temp_lmer,
            newdata = mydat,
            re.form = ~ (1|file),
            random.only = TRUE)
}

#Summarize
pig_sum <-
  summarize(group_by(pig, Wavelength, Bact_pop, Bact_timepoint, 
                     Bact_treat, Bact_treat_C, Bact_treat_P),
            Meas_norm_avg = mean(Meas_norm),
            Meas_norm_sd = sd(Meas_norm),
            Meas_norm_norm_avg = mean(Meas_norm_norm),
            Meas_norm_norm_sd = sd(Meas_norm_norm))

#Plots
ggplot(data = filter(pig, Bact_pop != "Blank"), 
       aes(x = as.numeric(Wavelength), y = Meas_norm, group = Well, color = Bact_treat)) +
  geom_line() +
  facet_grid(Bact_pop ~ Bact_timepoint)

ggplot(data = filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm, color = Bact_treat)) +
  geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.2)) +
  facet_grid(Wavelength ~ Bact_pop, scales = "free")

ggplot(data = filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm_norm, color = Bact_treat)) +
  geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.2)) +
  facet_grid(Wavelength ~ Bact_pop, scales = "free")

png("./plots/pigments.png", width = 8, height = 6, units = "in",
    res = 150)
ggplot(data = filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.6)) +
  facet_grid(Wavelength ~ ., scales = "free",
             labeller = labeller(Wavelength = c("405" = "405nm (pyoverdine)",
                                                "695" = "695nm (pyocyanin)"))) +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Absorbance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

png("./plots/pigments_norm.png", width = 8, height = 6, units = "in",
    res = 150)
ggplot(data = filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.6)) +
  facet_grid(Wavelength ~ ., scales = "free",
             labeller = labeller(Wavelength = c("405" = "405nm (pyoverdine)",
                                                "695" = "695nm (pyocyanin)"))) +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Absorbance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

ggplot(data = filter(pig_sum, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm_avg, color = Bact_treat)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.3)) +
  facet_grid(Wavelength ~ ., scales = "free")

# Biofilm ----

#Data
biofilm <- import_blockmeasures(
  list.files(path = "./biofilm_data/", full.names = T, pattern = ".*Data.csv$"),
  startrow = 26, endrow = 34, endcol = "M")
biofilm <- trans_wide_to_tidy(biofilm, id_cols = "block_name",
                              values_to = "biofilm_abs")
biofilm$block_name <- gsub(" Data$", "", biofilm$block_name)

#Designs
biofilm_designs <- import_blockdesigns(
  list.files(path = "./biofilm_data/", full.names = T, pattern = ".*Design.csv$"),
  join_as_cols = FALSE)
biofilm_designs$block_name <- gsub(" Design$", "", biofilm_designs$block_name)
#Change blanks to NA
biofilm_designs <- mutate(biofilm_designs,
                          Designs = ifelse(Designs == "W", NA, Designs))
biofilm_designs <- separate_tidy(biofilm_designs,
                                 "Designs",
                                 sep = ", ",
                                 into = c("Bact_pop", "Bact_timepoint", "Bact_treat", "Bact_isol"))

#Merge
biofilm <- merge_dfs(biofilm, biofilm_designs, drop = TRUE)

biofilm <- mutate(group_by(biofilm, block_name),
                  biofilm_norm = biofilm_abs - min(biofilm_abs), #almost always the Blank
                  biofilm_norm = ifelse(biofilm_norm == 0, 0.001, biofilm_norm)) 
biofilm <- filter(biofilm, Bact_pop != "Blank")

#Adjust columns for cross-compatibility
biofilm <- mutate(biofilm,
              Bact_isol = gsub("#", "", Bact_isol),
              Bact_pop = gsub("P", "", Bact_pop),
              Bact_timepoint = gsub("T", "", Bact_timepoint),
              Bact_treat_C = factor(
                ifelse(Bact_treat %in% c("NCP", "NN"), "NC",
                       ifelse(Bact_treat %in% c("CNP", "CP"), "C",
                              Bact_treat)),
                levels = c("Anc", "pilA", "NC", "C")),
              Bact_treat_P = factor(
                ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                       ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                              Bact_treat)),
                levels = c("Anc", "pilA", "NP", "P")),
              Bact_treat = factor(Bact_treat,
                                  levels = c("Anc", "pilA", "NN", "CNP", "NCP", "CP")))

#Based on contamination, exclude NCNP Pop 9 T4 and later, CNP Pop 9 T2 and later
biofilm <- filter(biofilm,
              !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 2) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))
#Based on weird blank values (and we already don't have much data anyway)
# exclude the whole 24-07-12 batch
biofilm <- dplyr::filter(ungroup(biofilm), 
                  block_name != "biofilm_data/2024-07-12 Abe Biofilm P9")

#Do stats/batch correction
biofilm_factorized <- mutate(
  biofilm,
  Bact_timepoint <- relevel(factor(Bact_timepoint), ref = "0"),
  Bact_treat <- relevel(factor(Bact_treat), ref = "Anc"))

biofilm$biofilm_norm_log10 <- log10(biofilm$biofilm_norm)
biofilm_lmer <- lmer(
  data = biofilm,
  biofilm_norm_log10 ~ (1|block_name) + Bact_treat + Bact_timepoint +
    (1|Bact_treat:Bact_pop) + (1|Bact_treat:Bact_pop:Bact_isol))
biofilm$biofilm_norm_log10_norm <-
  biofilm$biofilm_norm_log10 -
  predict(biofilm_lmer,
          newdata = biofilm,
          re.form = ~ (1|block_name),
          random.only = TRUE)

ggplot(data = biofilm, 
       aes(x = Bact_timepoint, y = biofilm_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.6)) +
  scale_y_log10() +
  facet_grid(~block_name) +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Absorbance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

ggplot(data = biofilm, 
       aes(x = Bact_timepoint, y = 10**biofilm_norm_log10_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.6)) +
  scale_y_log10() +
  facet_grid(~block_name) +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Absorbance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

pbio <- ggplot(data = biofilm, 
       aes(x = Bact_timepoint, y = biofilm_abs, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Absorbance") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

png("./plots/biofilm_nopops_nolog.png", width = 8, height = 6, units = "in",
    res = 150)
print(pbio)
dev.off()

png("./plots/biofilm_nopops.png", width = 8, height = 6, units = "in",
    res = 150)
print(pbio +
        scale_y_log10())
dev.off()
  

# Twitching ----
twitching <- read.csv("Twitching Motility Data.csv")

#Exclude bad isolates
twitching <-
  filter(twitching,
         !(Date == "2024-02-16" & Bact_pop == 9 &
             Bact_treat == "Anc" & Bact_timepoint == 0 & Bact_isol == 2),
         !(Date == "2024-02-16" & Bact_pop == 8 &
             Bact_treat == "NN" & Bact_timepoint == 1 & Bact_isol == 2),
         !(Date == "2024-02-16" & Bact_pop == 8 &
             Bact_treat == "CP" & Bact_timepoint == 8 & Bact_isol == 2),
         !(Date == "2024-02-03" & Bact_pop == 7 &
             Bact_treat == "NCP" & Bact_timepoint == 1 & Bact_isol == 2))

#Exclude contaminated populations
twitching <- 
  filter(twitching,
         !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 2) |
             (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

twitching <- mutate(
  twitching,
  cm_day = (width/2 + height/2)/2,
  Bact_treat_C = factor(ifelse(Bact_treat %in% c("NCP", "NN"), "NC",
                        ifelse(Bact_treat %in% c("CNP", "CP"), "C",
                               Bact_treat)),
                        levels = c("Anc", "NC", "C", "pilA")),
  Bact_treat_P = factor(ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                        ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                               Bact_treat)),
                        c("Anc", "NP", "P", "pilA")),
  Bact_treat = factor(Bact_treat, 
                      levels = c("Anc", "NN", "CNP", "NCP", "CP", "pilA")))

#Calculate batch-corrected rate
ref_Anc = mean(twitching$cm_day[twitching$Bact_treat == "Anc"])
ref_pilA = mean(twitching$cm_day[twitching$Bact_treat == "pilA"])

twitching <- mutate(
  group_by(twitching, Date),
  cm_day_norm = mynormalize(x = cm_day,
                            b_input = mean(cm_day[Bact_treat == "Anc"]),
                            b_output = ref_Anc,
                            a_input = mean(cm_day[Bact_treat == "pilA"]),
                            a_output = ref_pilA))

#Summarize by isol (which should only affect pilA)
twitching_isols <- summarize(group_by(twitching,
                                Date, Bact_pop, Bact_timepoint,
                                Bact_treat, Bact_treat_C, Bact_treat_P,
                                Bact_isol),
                       cm_day = mean(cm_day),
                       cm_day_norm = mean(cm_day_norm))


ggplot(data = twitching_isols,
       aes(x = Bact_timepoint, y = cm_day, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 5, alpha = 0.7) +
  facet_grid(~Bact_pop) +
  labs(x = "Timepoint", y = "Twitching motility (cm/day)",
       subtitle = "Population") +
  
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 20))

ptwitching_norm <-
  ggplot(data = twitching_isols,
       aes(x = Bact_timepoint, y = cm_day_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 5, alpha = 0.7) +
  facet_grid(~as.numeric(Bact_pop)) +
  labs(x = "Timepoint", y = "Batch-corrected Twitching Motility (cm/day)",
       subtitle = "Population") +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 20))
print(ptwitching_norm)

summarize(group_by(twitching_isols, Bact_treat_P, Bact_timepoint),
          across(cm_day_norm,
                 list(min = min, max = max, mean = mean)))

png("./plots/motility_norm_nopops.png", width = 7, height = 5, units = "in",
    res = 150)               
ptwitching_norm2 <-
  ggplot(data = twitching_isols,
         aes(x = Bact_timepoint, y = cm_day_norm, 
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5,
             position = position_dodge(width = 0.6)) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nTwitching Motility (cm/day)") +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 20))
print(ptwitching_norm2)
dev.off()

png("./plots/motility_norm_nopops_NConly.png", width = 7, height = 5, units = "in",
    res = 150)               
ggplot(data = filter(twitching_isols, Bact_treat_C != "C"),
       aes(x = Bact_timepoint, y = cm_day_norm, 
           color = Bact_treat_P, fill = Bact_treat_P)) +
  geom_point(size = 4, alpha = 0.5,
             position = position_dodge(width = 0.3)) +
  labs(x = "Timepoint", y = "Batch-corrected Twitching Motility\n(cm/day)") +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "No pilus control")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "No pilus control")) +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 20))
dev.off()

if(F) {
  ggplot(data = filter(twitching, Bact_pop == 10),
         aes(x = Bact_timepoint, y = cm_day_norm,
             color = Bact_treat, fill = Bact_treat, shape = Bact_isol)) +
    geom_point(size = 4, alpha = 0.5,
               position = position_dodge(width = 0.5)) +
    guides(shape = guide_legend(override.aes = list(fill = "black")))
}

#Growth curve plots ----
if(!"growth_curve_metrics.csv" %in% list.files()) {
  source("Growthcurve_analysis.R")
  setwd("..") #return wd to parent dir where this script needs it
}

gcdat_sum <- read.csv("growth_curve_metrics.csv", 
                      colClasses = c("Bact_pop" = "character", 
                                     "Bact_timepoint" = "character",
                                     "Bact_isol" = "character"))

png("./plots/max_dens_norm.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = gcdat_sum,
       aes(x = Bact_timepoint, y = max_dens_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ as.numeric(Bact_pop),
             labeller = as_labeller(function(x) {LETTERS[as.numeric(x)-2]})) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nMaximum Density (OD600)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16)) +
  labs(subtitle = "Replicate Population")
dev.off()

png("./plots/max_dens_norm_nopops.png", width = 7, height = 4, units = "in",
    res = 150)
f3d <- 
  ggplot(data = gcdat_sum,
         aes(x = Bact_timepoint, y = max_dens_norm,
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nMaximum Density (OD600)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
f3d
dev.off()

png("./plots/max_percap_norm.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = gcdat_sum,
       aes(x = Bact_timepoint, y = max_percap_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ as.numeric(Bact_pop),
             labeller = as_labeller(function(x) {LETTERS[as.numeric(x)-2]})) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nMaximum Growth Rate (/hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16)) +
  labs(subtitle = "Replicate Population")
dev.off()

png("./plots/max_percap_norm_nopops.png", width = 7, height = 4, units = "in",
    res = 150)
fsmaxpercap <- 
  ggplot(data = gcdat_sum,
         aes(x = Bact_timepoint, y = max_percap_norm,
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nMaximum Growth Rate (/hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
fsmaxpercap
dev.off()

png("./plots/auc_norm.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = gcdat_sum,
       aes(x = Bact_timepoint, y = auc_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ as.numeric(Bact_pop),
             labeller = as_labeller(function(x) {LETTERS[as.numeric(x)-2]})) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nArea Under the Curve\n(OD600 hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16)) +
  labs(subtitle = "Replicate Population")
dev.off()

png("./plots/auc_norm_nopops.png", width = 7, height = 4, units = "in",
    res = 150)
fsauc <- 
  ggplot(data = gcdat_sum,
         aes(x = Bact_timepoint, y = auc_norm,
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nArea Under the Curve\n(OD600 hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
fsauc
dev.off()

png("./plots/lag_time.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = gcdat_sum,
       aes(x = Bact_timepoint, y = lag,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ Bact_pop) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Lag Time (hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
dev.off()

png("./plots/lag_norm.png", width = 10, height = 4, units = "in",
    res = 150)
ggplot(data = gcdat_sum,
       aes(x = Bact_timepoint, y = lag_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(~ as.numeric(Bact_pop),
             labeller = as_labeller(function(x) {LETTERS[as.numeric(x)-2]})) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nLag time (hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16)) +
  labs(subtitle = "Replicate Population")
dev.off()

png("./plots/lag_norm_nopops.png", width = 7, height = 4, units = "in",
    res = 150)
f3c <- 
  ggplot(data = gcdat_sum,
         aes(x = Bact_timepoint, y = lag_norm,
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Batch-corrected\nLag time (hr)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
f3c
dev.off()


png("./plots/figS_gcmetrics.png", width = 13, height = 9, units = "in",
    res = 150)
cowplot::plot_grid(
  cowplot::plot_grid(
    f3c + guides(color = "none", fill = "none", shape = "none"), 
    fsmaxpercap + guides(color = "none", fill = "none", shape = "none"),
    f3d + guides(color = "none", fill = "none", shape = "none"),
    fsauc + guides(color = "none", fill = "none", shape = "none"),
    nrow = 2, labels = "AUTO",
    align = "hv", axis = "tblr"),
  cowplot::get_legend(f3c),
  ncol = 2, rel_widths = c(1, 0.3))
dev.off()

# Multivariate analyses ----
### Merge ----
isol_data <- 
  Reduce(
    full_join,
    list(
      dplyr::rename(mutate(ungroup(twitching_isols), 
                           Bact_pop = as.character(Bact_pop)), 
                    twitch_date = Date),
      mutate(ungroup(resis_isols_bactsum_undup),
             Bact_timepoint = as.character(Bact_timepoint),
             Bact_pop = as.character(Bact_pop),
             Bact_isol = as.character(Bact_isol),
             Bact_treat = plyr::revalue(Bact_treat, c("NCNP" = "NN"), 
                                        warn_missing = FALSE)),
      pivot_wider(
        select(
          filter(ungroup(pig), Bact_pop != "Blank", Wavelength %in% c(405, 695)),
          !c(file, Well, Measurements)),
        id_cols = c("Bact_pop", "Bact_timepoint", "Bact_treat", 
                    "Bact_treat_C", "Bact_treat_P", "Bact_isol"),
        names_from = "Wavelength", values_from = Meas_norm_norm,
        names_prefix = "wavelength_"),
      dplyr::rename(gcdat_sum, gc_date = date),
      select(mutate(filter(ungroup(biofilm), Bact_treat != "pilA"),
                    biofilm_block_name = block_name),
             !block_name & !Well)
    )
  )

isol_data <- mutate(
  isol_data,
  resis = 1-clearing_avg,
  resis_cat = case_when(
    Bact_isol == "pilA" ~ "resis",
    Plaque_morphology %in% c("Cloudy", "Tiny") ~ "part",
    Plaque_morphology == "None" ~ "resis",
    Plaque_morphology %in% c("Clear", "Clear?") ~ "suscep",
    TRUE ~ "nomatch"),
  resis_cat = factor(resis_cat, levels = c("suscep", "part", "resis", "nomatch")),
  twitch_cat = case_when(
    cm_day_norm < 0.5 ~ "none",
    cm_day_norm > 0.7 ~ "full",
    cm_day_norm >= 0.5 & cm_day_norm <= 0.7 ~ "intermediate"))

write.csv(isol_data, "aggregate_isolate_data.csv", row.names = FALSE)

### Plots ----
ggpairs(data = isol_data,
        columns = which(colnames(isol_data) %in% 
                          c("resis", #"cm_day_norm",  
                            "wavelength_405", "wavelength_695",
                            "max_percap_norm", "max_dens_norm", "auc_norm",
                            "lag_norm", "biofilm_norm_log10_norm")),
        upper = list(continuous = "smooth"),
        lower = list(continuous = "points"))


png("./plots/motility_vs_resis.png", width = 6, height = 5, units = "in",
    res = 150)
pmotil_resis <-
  ggplot(
    data = mutate(
      ungroup(isol_data),
      resis_cat = plyr::revalue(resis_cat,
                                c("suscep" = "Susceptible", 
                                  "part" = "Partially Susceptible", 
                                  "resis" = "Resistant"))),
    aes(x = resis_cat, y = cm_day_norm,
        color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 4, alpha = 0.4, position = position_dodge(width = .75)) +
  labs(x = "", y = "Batch-corrected\nTwitching Motility (cm/day)") +
  scale_x_discrete(limits = c("Susceptible", "Partially Susceptible",
                              "Resistant")) +
  scale_color_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6, 3)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P", "pilA"),
                    labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  scale_shape_manual(values = c(25, 21, 23, 24), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C", "pilA"),
                     labels = c("Ancestor", "No", "Yes", "pilA knockout")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 20))
print(pmotil_resis)
dev.off()

png("./plots/fig3_motility_costs_combo.png", width = 11, height = 10, units = "in",
    res = 150)
cowplot::plot_grid(
  cowplot::plot_grid(
    ptwitching_norm2 + guides(color = "none", fill = "none", shape = "none") +
      theme(axis.title.y = element_text(size = 16)),
    pmotil_resis + guides(color = "none", fill = "none", shape = "none"),
    f3c + guides(color = "none", fill = "none", shape = "none") +
      theme(axis.title.y = element_text(size = 16)), 
    f3d + guides(color = "none", fill = "none", shape = "none") +
      theme(axis.title.y = element_text(size = 16)),
    align = "h", axis = "tb", nrow = 2, labels = "AUTO"),
  cowplot::get_legend(ptwitching_norm2),
  rel_widths = c(1, 0.3))
dev.off()

###Motility stats ----
tmp <- mutate(isol_data,
              Bact_treat = relevel(as.factor(Bact_treat), ref = "NN"),
              Bact_treat_P = relevel(as.factor(Bact_treat_P), ref = "NP"),
              resis_cat = relevel(as.factor(resis_cat), ref = "suscep"))

twitch_lmer <- 
  lmer(data = filter(tmp, Bact_timepoint != "pilA" & Bact_timepoint != 0),
     cm_day_norm ~ Bact_treat + Bact_treat:as.numeric(Bact_timepoint) + (1|Bact_pop),
     REML = FALSE)
twitch_lmer_red <- 
  lmer(data = filter(tmp, Bact_timepoint != "pilA" & Bact_timepoint != 0),
       cm_day_norm ~ as.numeric(Bact_timepoint) + (1|Bact_pop),
       REML = FALSE)
anova(twitch_lmer_red, twitch_lmer)
twitch_lmer
summary(twitch_lmer)

twitch_lmer_byres <- 
  lmer(data = filter(tmp, Bact_timepoint != "pilA" & Bact_timepoint != 0),
       cm_day_norm ~ Bact_treat_P*resis_cat + Bact_treat_C + (1|Bact_pop))
twitch_lmer_byres_redresis <- 
  lmer(data = filter(tmp, Bact_timepoint != "pilA" & Bact_timepoint != 0),
       cm_day_norm ~ Bact_treat_P + Bact_treat_C + (1|Bact_pop))
twitch_lmer_byres_redtreatC <- 
  lmer(data = filter(tmp, Bact_timepoint != "pilA" & Bact_timepoint != 0),
       cm_day_norm ~ Bact_treat_P*resis_cat + (1|Bact_pop))
anova(twitch_lmer_byres_redresis, twitch_lmer_byres)
anova(twitch_lmer_byres_redtreatC, twitch_lmer_byres)
twitch_lmer_byres
summary(twitch_lmer_byres)


###Run PCA ----
pca_dat <- isol_data

pca_cols <- which(colnames(pca_dat) %in%
                    c(#"resis", #"cm_day_norm",
                    "wavelength_405", "wavelength_695",
                    "max_percap_norm", "max_dens_norm", "auc_norm", "lag_norm",
                    "biofilm_norm_log10_norm"))

pca_dat <- pca_dat[complete.cases(pca_dat[, pca_cols]), ]

isol_pca <- prcomp(pca_dat[, pca_cols], center = TRUE, scale = TRUE, retx = TRUE)

#Merge pca with data
pca_dat <- cbind(as.data.frame(isol_pca$x),
                 pca_dat)

#Rename rotations for better plotting
replace_names = c("wavelength_405" = "pyoverdin",
                  "wavelength_695" = "pyocyanin",
                  "max_percap_norm" = "r",
                  "max_dens_norm" = "k",
                  "auc_norm" = "auc",
                  "lag_norm" = "lag",
                  "biofilm_norm_log10_norm" = "biofilm")
row.names(isol_pca$rotation) <- plyr::revalue(
  x = row.names(isol_pca$rotation), replace = replace_names)
  
arrow_len <- 3
pca_base <- ggplot(pca_dat, aes(x = PC1, y = PC2)) +
  theme_bw() +
  labs(x = paste("PC1 (", 
                 round((100*((isol_pca$sdev)**2)/
                          sum((isol_pca$sdev)**2))[1], 1),
                 "%)", sep = ""),
       y = paste("PC2 (", 
                 round((100*((isol_pca$sdev)**2)/
                          sum((isol_pca$sdev)**2))[2], 1),
                 "%)", sep = "")) +
  theme(axis.title = element_text(size = 36),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 32),
        plot.title = element_text(size = 36),
        strip.text = element_text(size = 36))

png("./plots/PCA_time.png", width = 10, height = 7, units = "in", res = 150)
pca_base + 
  geom_polygon(stat = "ellipse", aes(fill = Bact_timepoint), alpha = 0.2) +
  geom_point(aes(color = Bact_timepoint), size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = row.names(isol_pca$rotation)),
                           size = 8, alpha = .8, color = "red4", 
                           nudge_x = c(1, 1, 0, 0, 0, 0, 1),
                           nudge_y = c(0, 0, 0, 0, 0, 0.2, 0),
                           min.segment.length = unit(1, "native")) +
  scale_color_manual(name = "Timepoint", values = c("black", viridisLite::viridis(n = 4))) +
  scale_fill_manual(name = "Timepoint", values = c("black", viridisLite::viridis(n = 4)))
dev.off()

png("./plots/PCA_phagetreat&bacttreat_bytime.png", width = 23, height = 8, units = "in", res = 150)
pca_base + 
  #Code for shared polygons across all facets
  # geom_polygon(
  #   data = select(pca_dat, !Bact_timepoint), 
  #   inherit.aes = FALSE, stat = "ellipse", 
  #   aes(x = PC1, y = PC2, fill = Bact_treat), 
  #   alpha = 0.2) +
  geom_polygon(stat = "ellipse", aes(fill = Bact_treat), alpha = 0.2) +
  geom_point(aes(color = Bact_treat, fill = Bact_treat), 
             size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .4, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = cbind(Bact_timepoint = "0",
                                        vars = row.names(isol_pca$rotation),
                                        as.data.frame(isol_pca$rotation)),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = vars),
                           inherit.aes = FALSE,
                           nudge_x = c(1, 1, 0, 0, 0, 0, 0.6),
                           nudge_y = c(0, 0, 0, 0, 0, 0.2, 0),
                           size = 8, alpha = .5, color = "red4", seed = 8,
                           min.segment.length = unit(1, "native")) +
  facet_grid(~ Bact_timepoint,
             labeller = labeller(Bact_timepoint = function(x){paste("Timepoint", x)})) +
  scale_color_manual(name = "Treatments",
                     breaks = c("Anc", "NN", "CNP", "NCP", "CP"),
                     labels = c("Ancestor", "No-phage, Pseudomonas-only",
                                "No-phage, all bacteria",
                                "Phage-added, Pseudomonas-only",
                                "Phage-added, all bacteria"),
                     values = my_cols[c(8, 7, 5, 1, 3)]) +
  scale_fill_manual(name = "Treatments",
                     breaks = c("Anc", "NN", "CNP", "NCP", "CP"),
                     labels = c("Ancestor", "No-phage, Pseudomonas-only",
                                "No-phage, all bacteria",
                                "Phage-added, Pseudomonas-only",
                                "Phage-added, all bacteria"),
                     values = my_cols[c(8, 7, 5, 1, 3)]) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1, byrow = TRUE),
         fill = guide_legend(ncol = 1, byrow = TRUE))
dev.off()

png("./plots/PCA_phagetreat.png", width = 11, height = 8, units = "in", res = 150)
pca_base + 
  geom_point(aes(color = Bact_treat_P), size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = row.names(isol_pca$rotation)),
                           size = 8, alpha = .8, color = "red4", seed = 8,
                           min.segment.length = unit(1, "native")) +
  scale_color_manual(name = "Phage\nTreatment", 
                     breaks = c("Anc", "NP", "P"),
                     values = my_cols[c(8, 1, 2)],
                     labels = c("Ancestor", "No Phage", "Phage"))
dev.off()

png("./plots/PCA_resis.png", width = 11, height = 6, units = "in", res = 150)
pca_base + 
  geom_point(aes(color = as.factor(resis_cat)), size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = row.names(isol_pca$rotation)),
                           size = 8, alpha = .8, color = "red4", seed = 8,
                           min.segment.length = unit(1, "native")) +
  scale_color_manual(name = "Resistance", 
                     breaks = c("suscep", "part", "resis"),
                     values = my_cols[c(8, 1, 2)],
                     labels = c("Susceptible", "Partial Resistance", "Full Resistance"))
dev.off()                        

#phage treat & bact treat             
pca_base + 
  geom_polygon(stat = "ellipse", aes(fill = Bact_treat), alpha = 0.2) +
  geom_point(aes(color = Bact_treat, fill = Bact_treat), 
             size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = row.names(isol_pca$rotation)),
                           size = 8, alpha = .8, color = "red4", seed = 8,
                           min.segment.length = unit(1, "native")) +
  scale_color_manual(name = "Treatments",
                     breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                     values = my_cols[c(8, 7, 1, 5, 3)]) +
  scale_fill_manual(name = "Treatments",
                    breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                    values = my_cols[c(8, 7, 1, 5, 3)])

#time by treat
pca_base + 
  geom_polygon(stat = "ellipse", aes(fill = Bact_timepoint), alpha = 0.2) +
  # geom_polygon(
  #   data = select(pca_dat, !Bact_treat), 
  #   inherit.aes = FALSE, stat = "ellipse", 
  #   aes(x = PC1, y = PC2, fill = Bact_timepoint), 
  #   alpha = 0.2) +
  geom_point(aes(color = Bact_timepoint), size = 5, alpha = 0.7) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  facet_grid(~ Bact_treat) +
  # ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
  #                          aes(x = arrow_len*PC1, y = arrow_len*PC2,
  #                              label = row.names(isol_pca$rotation)),
  #                          size = 8, alpha = .8, color = "red4", seed = 8,
  #                          min.segment.length = unit(1, "native")) +
  scale_color_manual(name = "Timepoint", 
                     breaks = c(0, 1, 2, 4, 8),
                     values = my_cols[c(8, 1, 2, 3, 7)]) +
  scale_fill_manual(name = "Timepoint", 
                    breaks = c(0, 1, 2, 4, 8),
                    values = my_cols[c(8, 1, 2, 3, 7)])


ggplot(data = filter(isol_data, Bact_treat != "pilA"),
       aes(x = Bact_timepoint, y = max_dens_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(as.numeric(Bact_pop) ~ resis_cat) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Timepoint (days)", y = "Maximum Density (OD600)") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16))

ggplot(data = filter(isol_data, Bact_treat != "pilA"),
       aes(x = resis_cat, y = max_dens_norm,
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                     breaks = c("Anc", "NP", "P"),
                     labels = c("Ancestor", "No", "Yes")) +
  scale_fill_manual(values = my_cols[c(8, 2, 6)], name = "Phage Added",
                    breaks = c("Anc", "NP", "P"),
                    labels = c("Ancestor", "No", "Yes")) +
  scale_shape_manual(values = c(25, 21, 23), name = "Other Bacteria Added",
                     breaks = c("Anc", "NC", "C"),
                     labels = c("Ancestor", "No", "Yes")) +
  guides(shape = guide_legend(override.aes = list(fill = "black")))

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
prelimdens_NC <- read.csv("2021-07-28_PA14H6_Coevol.csv", stringsAsFactors = F)
colnames(prelimdens_NC)[5:13] <- paste("-", 0:8, sep = "")

#Pivot longer
prelimdens_NC_lng <- as.data.frame(
  tidyr::pivot_longer(prelimdens_NC, cols = starts_with("-"),
                      names_to = "Dilution", values_to = "count"))
prelimdens_NC_lng <- prelimdens_NC_lng[prelimdens_NC_lng$count != "", ]

#For pop +2 at timepoint 1 bc density was too high to measure estimate it
# as 100 cfu at -8
prelimdens_NC_lng$count[prelimdens_NC_lng$Plating == "PA" & prelimdens_NC_lng$Treat == "+" &
                      prelimdens_NC_lng$Population == 2 & prelimdens_NC_lng$Timepoint == 1 &
                      prelimdens_NC_lng$Dilution == "-8"] <- 100

#Calculate densities & frequencies
prelimdens_NC_lng <- group_by(prelimdens_NC_lng, Timepoint, Plating, Treat, Population)
prelimdens_NC_sum <- summarize(
  prelimdens_NC_lng,
  dens_dil = max(as.numeric(Dilution[!(count %in% c("clumped", "lacy", "mbc", "too")) &
                                       !grepl("~", count)])),
  freq_dil = max(as.numeric(Dilution[!(count %in% c("clumped", "lacy", "mbc", "too"))])),
  dens_val = count[as.numeric(Dilution) == dens_dil],
  freq_val = count[as.numeric(Dilution) == freq_dil],
  dens = NA, freq_L = NA)
for (i in 1:nrow(prelimdens_NC_sum)) {
  #First calculate density
  if (!grepl("[,LS]", prelimdens_NC_sum$dens_val[i])) {
    #If there's not a comma, L, or S, just calculate density (frequency of L is 1)
    prelimdens_NC_sum$dens[i] <- 
      100*as.numeric(prelimdens_NC_sum$dens_val[i])/10**prelimdens_NC_sum$dens_dil[i]
  } else {
    #If there's a comma, split, sum, calculate density
    temp_d <- strsplit(prelimdens_NC_sum$dens_val[i], split = ",")[[1]]
    temp_d <- gsub("[SL]", "", temp_d)
    prelimdens_NC_sum$dens[i] <-
      100*sum(as.numeric(temp_d))/10**prelimdens_NC_sum$dens_dil[i]
  }
  #Now calculate frequency
  if (!grepl("[,LS]", prelimdens_NC_sum$freq_val[i])) {
    #If there's not a comma, L, or S, frequency of L is 1
    prelimdens_NC_sum$freq_L[i] <- 1
  } else {
    #If there's a comma, L, or S, split, calculate frequency
    temp_f <- strsplit(prelimdens_NC_sum$freq_val[i], split = ",")[[1]]
    temp_f <- gsub("~", "", temp_f)
    types <- substr(temp_f, nchar(temp_f), nchar(temp_f))
    temp_f <- substr(temp_f, 1, nchar(temp_f)-1)
    if (length(types) == 1) {
      if (types[1] == "L") {prelimdens_NC_sum$freq_L[i] <- 1
      } else if (types[1] == "S") {prelimdens_NC_sum$freq_L[i] <- 0
      }
    } else {
      prelimdens_NC_sum$freq_L[i] <- 
        as.numeric(temp_f[which(types == "L")])/sum(as.numeric(temp_f))
    }
  }
}

#Add pseudo-data on initial inoculation densities
prelimdens_NC_sum <- 
  full_join(prelimdens_NC_sum,
            data.frame("dens" = c(rep(10**7, 12), rep(10**5, 6)),
                       "Timepoint" = rep(0, 18),
                       "Plating" = c(rep("PA", 12), rep("H6", 6)),
                       "Treat" = c(rep("-", 6), rep("+", 12)),
                       "Population" = rep(1:6, 3),
                       "freq_L" = rep(1, 18)))

prelimdens_NC_sum <- rename(prelimdens_NC_sum,
                            Strain = Plating,
                            Phage_treat = Treat,
                            dens_ml = dens)
prelimdens_NC_sum <- prelimdens_NC_sum[, !names(prelimdens_NC_sum) %in% 
                                         c("dens_dil", "freq_dil", "dens_val", "freq_val")]

prelimdens_NC_sum <- cbind(data.frame(Comm_treat = "-"),
                           prelimdens_NC_sum)

#Set counts below detection to limit of detection
prelimdens_NC_sum$bd <- ifelse(prelimdens_NC_sum$dens_ml == 0, TRUE, FALSE)
for (i in 1:nrow(prelimdens_NC_sum)) {
  if (prelimdens_NC_sum$dens_ml[i] == 0) {
    if(prelimdens_NC_sum$Comm_treat[i] == "+") {
      prelimdens_NC_sum$dens_ml[i] <- 1000/2.5
    } else if (prelimdens_NC_sum$Comm_treat[i] == "-") {
      prelimdens_NC_sum$dens_ml[i] <- 1000/5
    }
  }
}

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
prelimdens_C <- read.csv("2021-12-17_community_experiment.csv",
                      stringsAsFactors = FALSE)

#Split PA data to calculate freq_L and total count
# (note that for now we've summed in when we counted frequency from
# a lower dilution because it only throws off the density by 1 or 2 cfu
# but not accurate)
prelimdens_C$freq_L <- NA
for (i in which(prelimdens_C$Strain == "PA")) {
  #Check if there are small/large freq's
  if(grepl("[,LS]", prelimdens_C$Count[i])) {
    temp_d <- strsplit(prelimdens_C$Count[i], split = ",")[[1]]
    temp_e <- as.numeric(gsub("[SL]", "", temp_d))
    names(temp_e) <- substr(temp_d, nchar(temp_d), nchar(temp_d))
    
    prelimdens_C$Count[i] <- sum(temp_e)
    if (length(temp_e) > 1) {
      prelimdens_C$freq_L[i] <- temp_e["L"]/sum(temp_e)
    } else {
      if(names(temp_e) == "L") {
        prelimdens_C$freq_L[i] <- 1
      } else if (names(temp_e) == "S") {
        prelimdens_C$freq_L[i] <- 0
      }
    }
  } else {
    prelimdens_C$freq_L[i] <- 1
  }
}

#Because I think the small colony morphs that appeared in the
# +phage +competitors on Day 6 were probably mis-counts, set
# those frequencies to 1
prelimdens_C$freq_L[prelimdens_C$Phage_treat == "+" & 
                   prelimdens_C$Timepoint == 6 &
                   prelimdens_C$Strain == "PA"] <- 1

#Calculate density
prelimdens_C$dens_ml <- (1000/2.5) * 
  as.numeric(prelimdens_C$Count)/(10**as.numeric(prelimdens_C$Dilution))

#Fix T0 data by pseudo-replicating for plotting
prelimdens_C <- 
  rbind(
    do.call("rbind", replicate(12, prelimdens_C[prelimdens_C$Timepoint == 0, ],
                               simplify = FALSE)),
    prelimdens_C[prelimdens_C$Timepoint != 0, ])
prelimdens_C$Phage_treat[1:(12*5)] <- rep(c("+", "-"), each = 6*5)
prelimdens_C$Pop[1:(12*5)] <- rep(rep(1:6, each = 5), 2)
prelimdens_C$Count[prelimdens_C$Timepoint == 0 & prelimdens_C$Strain == "H6" &
                  prelimdens_C$Phage_treat == "-"] <- NA
prelimdens_C$dens_ml[prelimdens_C$Timepoint == 0 & prelimdens_C$Strain == "H6" &
                    prelimdens_C$Phage_treat == "-"] <- NA
prelimdens_C <- prelimdens_C[!is.na(prelimdens_C$dens_ml), ]

prelimdens_C <- rename(prelimdens_C,
                       Population = Pop)
prelimdens_C <- prelimdens_C[, !names(prelimdens_C) %in% 
                               c("Count", "Dilution")]

prelimdens_C <- cbind(data.frame(Comm_treat = "+"),
                      prelimdens_C)

#Set counts below detection to limit of detection
prelimdens_C$bd <- ifelse(prelimdens_C$dens_ml == 0, TRUE, FALSE)
for (i in 1:nrow(prelimdens_C)) {
  if (prelimdens_C$dens_ml[i] == 0) {
    if(prelimdens_C$Comm_treat[i] == "+") {
      prelimdens_C$dens_ml[i] <- 1000/2.5
    } else if (prelimdens_C$Comm_treat[i] == "-") {
      prelimdens_C$dens_ml[i] <- 1000/5
    }
  }
}

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

#TODO - reconfirm that all points labeled as 0's are truly 0's
#       and not mis-platings

#Read file
density <- read.csv("2022-09-19_community_experiment_2.csv",
                    stringsAsFactors = FALSE, na.strings = "")
density <- filter(density, !is.na(Timepoint))

#Adjust for consistent style across data
density <- mutate(density,
                  treat_P = ifelse(Phage_treat == "+", "P", "NP"),
                  treat_C = ifelse(Comm_treat == "+", "C", "NC"))

#Calculate total density, large freq
# for total density, note counts are for 2.5 uL
density <- mutate(group_by(density, Timepoint, treat_P, treat_C, Pop, Strain),
              cfu_ml = ifelse(!is.na(Count),
                              400*Count/10**Dilution,
                              400*sum(Small_count/10**Small_dilution, 
                                      Large_count/10**Large_dilution, 
                                      na.rm = TRUE)),
              #Initial densities were estimated for the inoculation
              # source, which was diluted 10:10000 for phage or
              # 100:10000 for bacteria into the actual flask
              cfu_ml = ifelse(Timepoint == 0 & Pop == "H6",
                              10/10000 * Count/10**Dilution,
                              cfu_ml),
              cfu_ml = ifelse(Timepoint == 0 & Pop != "H6",
                              100/10000 * Count/10**Dilution,
                              cfu_ml),
              S_cfu_ml = ifelse(!is.na(Small_count),
                                400*Small_count/10**Small_dilution,
                                0),
              L_cfu_ml = ifelse(!is.na(Large_count),
                                400*Large_count/10**Large_dilution,
                                ifelse(is.na(Small_count), cfu_ml, 0)),
              L_freq = L_cfu_ml/cfu_ml
)
#Fix rows where no data was collected because of mis-plating
density <- mutate(density,
              cfu_ml = ifelse(all(is.na(Count), is.na(Small_count),
                                  is.na(Large_count)),
                              NA, cfu_ml))
#Handle below-detection points
density <- mutate(density,
              bd = ifelse(cfu_ml == 0, TRUE, FALSE),
              cfu_ml = ifelse(cfu_ml == 0,
                              400 * 1/10**Dilution,
                              cfu_ml))
density <- as.data.frame(density)

density <- dup_rows(mydf = density, from_vals = NA, to_vals = c("P", "NP"),
                    mut_col = "treat_P")
density <- dup_rows(mydf = density, from_vals = NA, to_vals = c("C", "NC"),
                    mut_col = "treat_C")

#Change factor orders
density <- mutate(density,
                  treat_C = factor(treat_C, levels = c("NC", "C")),
                  treat_P = factor(treat_P, levels = c("NP", "P")))

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

png("plots/density_legendside.png", width = 7, height = 5, units = "in",
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

png("plots/density_phage.png", width = 7, height = 5, units = "in",
    res = 150)
ggplot(
  data = filter(density, Strain == "H6", Pop %in% 3:10, !is.na(cfu_ml),
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

#Merge with population data
density_mod <- 
  filter(mutate(density, 
                Bact_timepoint = Timepoint,
                Bact_pop = Pop,
                Bact_treat_P = treat_P,
                Bact_treat_C = treat_C),
         Strain == "H6",
         Bact_pop %in% 3:10)
change_before <- function(x) {
  if(length(x) < 2) {return(NA)
  } else {return(c(NA, x[2:length(x)]/x[1:(length(x)-1)]))}
}
change_after <- function(x) {
  if(length(x) < 2) {return(NA)
  } else {return(c(x[2:length(x)]/x[1:(length(x)-1)], NA))}
}
density_mod <- mutate(
  group_by(density_mod, Bact_pop, Bact_treat_P, Bact_treat_C),
  change_before = change_before(cfu_ml[order(as.numeric(Timepoint))]),
  change_after = change_after(cfu_ml[order(as.numeric(Timepoint))]))

if(FALSE) {
  ggplot(data = density_mod, aes(x = Timepoint, y = change_after, group = Bact_pop)) + 
    geom_line() + 
    scale_y_log10() +
    facet_grid(Bact_treat_C ~ Bact_treat_P)
  
  ggplot(data = density_mod, 
         aes(x = Timepoint, y = cfu_ml, group = Bact_pop)) + 
    geom_line(aes(color = change_after < 0.1)) + 
    scale_y_log10() +
    facet_grid(Bact_treat_C ~ Bact_treat_P)
  
  ggplot(data = density_mod, aes(x = change_after)) + 
    geom_histogram() + 
    scale_x_log10() +
    facet_grid(Bact_treat_C ~ Bact_treat_P) +
    geom_vline(xintercept = 1, color = "red")
  
  summary(lm(data = filter(density_mod, Bact_timepoint > 1, Bact_treat_P == "P"),
             cfu_ml ~ Bact_treat_C + Bact_timepoint:Bact_treat_C))
}

# Resistance data ----

### EOP data ----
EOP <- read.csv("2023-05-20 EOP Timeshift data.csv", stringsAsFactors = F)

#ID below-detection counts, and calculate pfu/mL (note 1.5 uL spots)
EOP <- mutate(EOP, 
              bd = ifelse(pfu == 0, TRUE, FALSE),
              pfu = ifelse(pfu == 0, 1, pfu),
              pfu_ml = 1000/1.5 * pfu / 10**dilution)

#Deal with T0 labels
# When it's T0 with a specific treat, the T0 should be matched
#  to that treat name
# When it's T0 with T0, those rows need to be duplicated across all
#  the treats so they plot for all four treats 
EOP <- mutate(EOP,
              Bact_treat = ifelse(Bact_treat == "Anc" & Phage_treat != "Anc",
                                  Phage_treat, Bact_treat),
              Phage_treat = ifelse(Phage_treat == "Anc" & Bact_treat != "Anc",
                                   Bact_treat, Phage_treat))
EOP <- dup_rows(EOP, "Bact_treat", "Anc", c("CP", "NCP", "CNP", "NCNP"))
EOP <- mutate(EOP,
              Phage_treat = ifelse(Phage_treat == "Anc", Bact_treat, Phage_treat))
EOP <- separate_wider_regex(EOP,
                            cols = "Bact_treat",
                            patterns = c("Bact_treat_C" = "^N?C", 
                                         "Bact_treat_P" = "N?P$"),
                            cols_remove = FALSE)
EOP <- separate_wider_regex(EOP,
                            cols = "Phage_treat",
                            patterns = c("Phage_treat_C" = "^N?C", 
                                         "Phage_treat_P" = "N?P$"),
                            cols_remove = FALSE)
EOP$Bact_treat_C <- relevel(as.factor(EOP$Bact_treat_C), ref = "NC")
EOP$Bact_treat_P <- relevel(as.factor(EOP$Bact_treat_P), ref = "NP")

#Based on contamination, exclude NCNP Pop 9 T4 and later, CNP Pop 9 T2 and later
EOP <- filter(EOP, 
              !((Bact_pop == 9 & Bact_treat == "NCNP" & Bact_timepoint >= 4) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

#average technical replicates together (should just be Ancestors)
EOP_isols <- summarize(group_by(EOP, 
                                Date,
                                Bact_pop, Bact_treat, Bact_treat_C,
                                Bact_treat_P, Bact_timepoint, Bact_isol,
                                Phage_pop, Phage_treat, Phage_treat_C,
                                Phage_treat_P, Phage_timepoint, Phage_isol),
                       bd = any(bd == TRUE),
                       pfu_ml = 10**(mean(log10(pfu_ml))))

#Calculate EOP (relative to highest pfu/mL for ea phage stock)
EOP_isols <- mutate(group_by(EOP_isols, 
                             Date, Phage_pop, Phage_treat, 
                             Phage_timepoint, Phage_isol),
              EOP = pfu_ml/max(pfu_ml),
              EOP_plot = EOP)
EOP_isols$EOP_plot[EOP_isols$bd == TRUE] <- 
  max(EOP_isols$EOP[EOP_isols$bd == TRUE], na.rm = TRUE)

#Calculate average across pops
EOP_sum <- summarize(group_by(EOP_isols,
                              Bact_treat, Bact_treat_C,
                              Bact_treat_P, Bact_timepoint,
                              Phage_treat, Phage_treat_C,
                              Phage_treat_P, Phage_timepoint),
                     EOP_avg = 10**mean(log10(EOP), na.rm = TRUE),
                     EOP_plot_avg = 10**mean(log10(EOP_plot), na.rm = TRUE))

#Prepare to merge
EOP_isols <- mutate(EOP_isols, clearing = ifelse(bd == TRUE, 0, 1))

### Spot test data ----
spots <- read.csv("2023-05-24 Spot Timeshift data.csv", stringsAsFactors = F)

if(any(spots$Bact_treat != spots$Phage_treat &
       spots$Bact_treat != "Anc" &
       spots$Phage_treat != "Anc")) {
  stop("treatment columns are not identical")}

#Test plot before manipulation
ggplot(data = spots,
       aes(x = as.factor(paste(Bact_timepoint, Bact_isol)), 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = clearing)) +
  facet_grid(Bact_treat * Phage_treat ~ Bact_pop) +
  # labs(x = "Bacterial timepoint", y = "Phage timepoint",
  #      subtitle = "Population") +
  theme(strip.text = element_text(size = 8))

#Deal with T0 labels
# When it's T0 with a specific treat, the T0 should be matched
#  to that treat name
# When it's T0 with T0, those rows need to be duplicated across all
#  the treats so they plot for all four treats 
spots <- mutate(spots,
              Bact_treat = ifelse(Bact_treat == "Anc" & Phage_treat != "Anc",
                                  Phage_treat, Bact_treat),
              Phage_treat = ifelse(Phage_treat == "Anc" & Bact_treat != "Anc",
                                   Bact_treat, Phage_treat))
spots <- dup_rows(spots, "Bact_treat", "Anc", c("CP", "NCP", "CNP", "NCNP"))
spots <- mutate(spots,
              Phage_treat = ifelse(Phage_treat == "Anc", Bact_treat, Phage_treat))
spots <- separate_wider_regex(spots,
                            cols = "Bact_treat",
                            patterns = c("Bact_treat_C" = "^N?C", 
                                         "Bact_treat_P" = "N?P$"),
                            cols_remove = FALSE)
spots <- separate_wider_regex(spots,
                            cols = "Phage_treat",
                            patterns = c("Phage_treat_C" = "^N?C", 
                                         "Phage_treat_P" = "N?P$"),
                            cols_remove = FALSE)
spots$Bact_treat_C <- relevel(as.factor(spots$Bact_treat_C), ref = "NC")
spots$Bact_treat_P <- relevel(as.factor(spots$Bact_treat_P), ref = "NP")

#Based on contamination, exclude NCNP Pop 9 T4 and later, CNP Pop 9 T2 and later
spots <- filter(spots, 
              !((Bact_pop == 9 & Bact_treat == "NCNP" & Bact_timepoint >= 4) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

#average technical replicates together (should just be Ancestors)
spots_isols <- summarize(group_by(spots, 
                                Date,
                                Bact_pop, Bact_treat, Bact_treat_C,
                                Bact_treat_P, Bact_timepoint, Bact_isol,
                                Phage_pop, Phage_treat, Phage_treat_C,
                                Phage_treat_P, Phage_timepoint, Phage_isol),
                       clearing = mean(clearing))

### Cross streak data ----
cs_isols <- read.csv("William cross-streak data.csv", fileEncoding = "UTF-8")

cs_isols <- separate_wider_regex(cs_isols,
                              cols = "Bact_treat",
                              patterns = c("Bact_treat_C" = "^N?C", 
                                           "Bact_treat_P" = "N?P$"),
                              cols_remove = FALSE)
cs_isols <- separate_wider_regex(cs_isols,
                              cols = "Phage_treat",
                              patterns = c("Phage_treat_C" = "^N?C", 
                                           "Phage_treat_P" = "N?P$"),
                              cols_remove = FALSE)

cs_isols <- pivot_longer(cs_isols, cols = starts_with("Bact_timepoint"),
                         names_to = "Bact_timepoint", 
                         names_prefix = "Bact_timepoint.T",
                         names_transform = as.integer,
                         values_to = "clearing",
                         values_transform = function(x) {as.numeric(x == "S")})

#cs_isols includes no data from CNP or NCNP treats, so there's no
# contamination to exclude

### Merged resistance data ----

#Merge different resistance data
resis_isols <- full_join(cs_isols, full_join(spots_isols, EOP_isols))
resis_isols <- mutate(resis_isols,
                      Bact_treat = factor(Bact_treat),
                      Bact_treat_C = relevel(factor(Bact_treat_C), "NC"),
                      Bact_treat_P = factor(Bact_treat_P))

#Write copy to file
write.csv(resis_isols, "aggregate_resis_data.csv", row.names = FALSE)

#Calculate average for ea bact isol across phages
resis_isols_bactsum <- summarize(
  group_by(resis_isols,
           Bact_treat, Bact_treat_C, Bact_treat_P, Bact_timepoint,
           Bact_pop, Bact_isol),
  clearing_avg = mean(clearing, na.rm = TRUE))
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
         Bact_pop, Bact_isol, clearing_avg))


#Calculate average across pops
resis_sum_treats <- summarize(group_by(resis_isols,
                                Bact_treat, Bact_treat_C,
                                Bact_treat_P, Bact_timepoint,
                                Phage_treat, Phage_treat_C,
                                Phage_treat_P, Phage_timepoint),
                       clearing_avg = mean(clearing, na.rm = TRUE))

png("plots/resis_grid.png", width = 8, height = 5, units = "in", res = 150)
ggplot(data = resis_isols,
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing))) +
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
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Phage_pop * Bact_timepoint,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Population\nBacterial timepoint") +
  NULL
dev.off()

png("plots/resis_grid_avg.png", width = 5.5, height = 4, units = "in", res = 150)
ggplot(data = resis_sum_treats,
       aes(x = as.factor(Bact_timepoint), 
           y = as.factor(Phage_timepoint))) +
  geom_tile(aes(fill = clearing_avg)) +
  facet_grid(Bact_treat_C ~ Bact_treat_P,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only Pseudomonas", C = "All Bacteria"))) +
  labs(x = "Bacterial timepoint", y = "Phage timepoint") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 11.5),
        legend.title = element_text(size = 16)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7",
                         name = "Frequency of\nsusceptibility")
dev.off()

png("plots/resis_grid_phagetreats.png", width = 10, height = 3.5, 
    units = "in", res = 150)
ggplot(data = filter(resis_isols, Bact_treat_P == "P"),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing))) +
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
  ggh4x::facet_nested(Bact_treat_P * Bact_treat_C ~ Phage_pop * Bact_timepoint,
                      labeller = labeller(
                        Bact_treat_P = c(NP = "No Phage", P = "Phage Added"),
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria"))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Population\nBacterial timepoint") +
  NULL
dev.off()

png("plots/resis_grid_phagetreats2.png", width = 10, height = 3.5, 
    units = "in", res = 150)
ggplot(data = filter(resis_isols, Bact_treat_P == "P"),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing))) +
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

png("plots/resis_grid_phagetreats_split1.png", width = 6, height = 3.5, 
    units = "in", res = 150)
ggplot(data = filter(resis_isols, Bact_treat_P == "P", Bact_timepoint %in% c(0, 1, 2)),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing))) +
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
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria\nAdded"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Replicate population") +
  NULL
dev.off()

png("plots/resis_grid_phagetreats_split2.png", width = 4.5, height = 3.5, 
    units = "in", res = 150)
ggplot(data = filter(resis_isols, Bact_treat_P == "P", Bact_timepoint %in% c(4, 8)),
       aes(x = Bact_isol, 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  geom_tile(aes(fill = as.factor(clearing))) +
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
                        Bact_treat_C = c(NC = "Only\nPseudomonas", C = "All Bacteria\nAdded"),
                        Bact_timepoint = function(x) paste("Bacterial timepoint", x))) +
  labs(x = "Bacterial Isolate", y = "Phage timepoint",
       subtitle = "Replicate population") +
  NULL
dev.off()

p1 <- ggplot(data = filter(resis_isols, Bact_treat_P == "P",
                     Bact_treat_C == "NC"),
       aes(x = as.factor(paste(Bact_timepoint, Bact_isol)), 
           y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  facet_wrap(~Phage_pop, ncol = 2,
             labeller = as_labeller(function(x){paste("Pop", x)})) +
  geom_tile(aes(fill = as.factor(clearing))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  scale_x_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        axis.text.x = element_text(hjust = 0.2)) +
  scale_fill_manual(breaks = c(1, 0, NA),
                    labels = c("Susceptible", "Resistant", "Phage Extinct"),
                    values = c("#132B43", "#56B1F7", "gray"),
                    name = "") +
  labs(x = "Bacterial timepoint", y = "Phage timepoint",
       subtitle = "Only Pseudomonas") +
  theme(legend.position = "right") +
  NULL

p2 <- ggplot(data = filter(resis_isols, Bact_treat_P == "P",
                           Bact_treat_C == "C"),
             aes(x = as.factor(paste(Bact_timepoint, Bact_isol)), 
                 y = as.factor(paste(Phage_timepoint, Phage_isol)))) +
  facet_wrap(~Phage_pop, ncol = 2,
             labeller = as_labeller(function(x){paste("Pop", x)})) +
  geom_tile(aes(fill = as.factor(clearing))) +
  scale_y_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  scale_x_discrete(labels = c(0, "", 1, "", 2, "", 4, "", 8, "")) +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text.y = element_text(vjust = -0.1),
        axis.text.x = element_text(hjust = 0.2)) +
  scale_fill_manual(breaks = c(1, 0, NA),
                    labels = c("Susceptible", "Resistant", "Phage Extinct"),
                    values = c("#132B43", "#56B1F7", "gray"),
                    name = "") +
  labs(x = "Bacterial timepoint", y = "Phage timepoint",
       subtitle = "All Bacteria") +
  theme(legend.position = "right") +
  NULL

png("plots/resis_grid_phagetreats_wrapped.png", width = 6.5, height = 6, 
    units = "in", res = 150)
plot_grid(p1 + guides(fill = "none"), 
          p2 + guides(fill = "none"), 
          get_legend(p1),
          nrow = 1, rel_widths = c(1, 1, 0.5))
dev.off()

### Population dynamics vs resistance ----

#Calculate average resistance w/in pops
resis_sum_pops <- summarize(group_by(resis_isols,
                                     Bact_treat, Bact_treat_C,
                                     Bact_treat_P, Bact_timepoint,
                                     Bact_pop,
                                     Phage_treat, Phage_treat_C,
                                     Phage_treat_P, Phage_timepoint,
                                     Phage_pop),
                            clearing_avg = mean(clearing, na.rm = TRUE))

#Merge with population data
pop_resis_merge <- 
  select(
    full_join(
      filter(density_mod, Bact_timepoint %in% c(0, 1, 2, 4, 8)),
      filter(resis_sum_pops, Bact_timepoint == Phage_timepoint),
      by = c("Bact_timepoint", "Bact_pop", "Bact_treat_P", "Bact_treat_C")),
    Strain, cfu_ml, bd, Bact_timepoint, Bact_treat_P, Bact_treat_C, Bact_pop,
    Phage_timepoint, Phage_treat_P, Phage_treat_C, Phage_pop,
    clearing_avg, change_before, change_after)

png("./plots/phagedens_v_resisfreq.png", width = 5, height = 4, units = "in",
    res = 150)
ggplot(data = filter(pop_resis_merge, Phage_treat_P == "P"),
       aes(y = cfu_ml, x = 1-clearing_avg, color = as.factor(Bact_timepoint))) +
  geom_point() +
  scale_y_log10() +
  labs(y = "pfu/mL", x = "Frequency of resistance") +
  scale_color_manual(name = "Timepoint", values = my_cols[c(8, 1, 2, 3, 7)]) +
  theme_bw()
dev.off()
  
png("./plots/changebefore_v_resisfreq.png", width = 5, height = 4, units = "in",
    res = 150)
ggplot(data = filter(pop_resis_merge, Phage_treat_P == "P"),
       aes(y = change_before, x = 1-clearing_avg, color = as.factor(Bact_timepoint))) +
  geom_point() +
  scale_y_log10() +
  labs(y = "Fold-change in pfu/mL from day before", x = "Frequency of resistance") +
  scale_color_manual(name = "Timepoint", values = my_cols[c(8, 1, 2, 3, 7)]) +
  theme_bw()
dev.off()

png("./plots/changeafter_v_resisfreq.png", width = 5, height = 4, units = "in",
    res = 150)
ggplot(data = filter(pop_resis_merge, Phage_treat_P == "P"),
       aes(y = change_after, x = 1-clearing_avg, color = as.factor(Bact_timepoint))) +
  geom_point() +
  scale_y_log10() +
  labs(y = "Fold-change in pfu/mL to day after", x = "Frequency of resistance") +
  scale_color_manual(name = "Timepoint", values = my_cols[c(8, 1, 2, 3, 7)]) +
  theme_bw()
dev.off()

# Pigments ----

pig1 <- import_blockmeasures(
  list.files(path = "./pigment_data/", full.names = T, pattern = ".*Production.csv$"),
  startrow = 42, endrow = 50, endcol = "M",
  metadata = list("Wavelength" = c(41, "A")))
pig2 <- import_blockmeasures(
  list.files(path = "./pigment_data/", full.names = T, pattern = ".*Production.csv$"),
  startrow = 57, endrow = 65, endcol = "M",
  metadata = list("Wavelength" = c(56, "A")))

colnames(pig1)[1] <- "file"
colnames(pig2)[1] <- "file"
pig12 <- trans_wide_to_tidy(list(pig1, pig2), id_cols = c("file", "Wavelength"),
                            values_to_numeric = FALSE)
pig12 <- lapply(
  pig12, 
  FUN = function(x) {x$Wavelength <- gsub(".*:", "", x$Wavelength); return(x)})

pigspectrum <- read_wides(
  list.files(path = "./pigment_data/", full.names = T, pattern = ".*Production.csv$"),
  startrow = 69, endrow = 90)
pigspectrum <- trans_wide_to_tidy(pigspectrum, id_cols = c("file", "Wavelength"),
                                  values_to_numeric = FALSE)

pig <- merge_dfs(pig12, pigspectrum, collapse = TRUE)
pig$file <- regmatches(pig$file, regexpr("\\d{4}-\\d{2}-\\d{2}", pig$file))

#Read designs
designs <-
  read_blocks(
    list.files(path = "./pigment_data/", full.names = T, pattern = ".*Design.csv$"))
designs <- trans_block_to_wide(designs)
designs <- trans_wide_to_tidy(designs, id_cols = "block_name",
                              values_to_numeric = FALSE)
designs <- mutate(designs, 
                  Measurements = ifelse(Measurements == "", NA, Measurements))
designs <- separate_tidy(
  designs, "Measurements", sep = ", ",
  into = c("Bact_pop", "Bact_timepoint", "Bact_treat", "Bact_isol"),
  coerce_NA = FALSE)
colnames(designs)[1] <- "file"
designs$file <- regmatches(designs$file, regexpr("\\d{4}-\\d{2}-\\d{2}", designs$file))

#Merge designs
pig <- merge_dfs(pig, designs, drop = TRUE)

#Adjust columns for cross-compatibility
pig <- mutate(pig,
              Bact_isol = gsub("#", "", Bact_isol),
              Bact_pop = gsub("P", "", Bact_pop),
              Bact_timepoint = gsub("T", "", Bact_timepoint),
              Bact_treat_C = factor(
                ifelse(Bact_treat %in% c("NCP", "NN"), "NC",
                                    ifelse(Bact_treat %in% c("CNP", "CP"), "C",
                                           Bact_treat)),
                levels = c("Anc", "NC", "C")),
              Bact_treat_P = factor(
                ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                                    ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                                           Bact_treat)),
                levels = c("Anc", "NP", "P")),
              Bact_treat = factor(Bact_treat,
                                  levels = c("Anc", "NN", "CNP", "NCP", "CP")))

#Based on contamination, exclude NCNP Pop 9 T4 and later, CNP Pop 9 T2 and later
pig <- filter(pig,
              !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 4) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

#Subtract blank
pig$Measurements[pig$Measurements == "OVRFLW"] <- NA
pig$Measurements <- as.numeric(pig$Measurements)
pig <- mutate(group_by(pig, file, Wavelength),
              Meas_norm = Measurements - Measurements[Bact_pop == "Blank"])
pig <- filter(pig, Bact_isol != "Blank")

#Do batch correction
# pig$Meas_norm2 <- NA
# for(wave in unique(pig$Wavelength)) {
#   #TODO
#   #something like this
#   #see: https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
#   twitching_lmer <- 
#     lmer(data = twitching,
#          cm_day ~ (1|Date) + Bact_treat + Bact_timepoint)
#   
#   twitching$cm_day_norm <-
#     twitching$cm_day -
#     predict(twitching_lmer,
#             newdata = twitching,
#             re.form = ~ (1|Date),
#             random.only = TRUE)
# }

#Summarize
pig_sum <-
  summarize(group_by(pig, Wavelength, Bact_pop, Bact_timepoint, 
                     Bact_treat, Bact_treat_C, Bact_treat_P),
            Meas_norm_avg = mean(Meas_norm),
            Meas_norm_sd = sd(Meas_norm))

#Plots
ggplot(data = filter(pig, Bact_pop != "Blank"), 
       aes(x = as.numeric(Wavelength), y = Meas_norm, group = Well, color = Bact_treat)) +
  geom_line() +
  facet_grid(Bact_pop ~ Bact_timepoint)

ggplot(data = filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm, color = Bact_treat)) +
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

ggplot(data = filter(pig_sum, Bact_pop != "Blank", Wavelength %in% c(405, 695)), 
       aes(x = Bact_timepoint, y = Meas_norm_avg, color = Bact_treat)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.3)) +
  facet_grid(Wavelength ~ ., scales = "free")

# Twitching ----
twitching <- read.csv("Twitching Motility Data.csv")

#Exclude bad isolates
twitching <-
  filter(twitching,
         !(Date == "2024-02-16" & Bact_pop == 9 &
             Bact_treat == "Anc" & Bact_timepoint == 0 & Bact_isol == 2),
         !(Date == "2024-02-16" & Bact_pop == 8 &
             Bact_treat == "NN" & Bact_timepoint == 1 & Bact_isol == 2))

#Exclude contaminated populations
twitching <- 
  filter(twitching,
         !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 4) |
             (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

twitching <- mutate(
  twitching,
  cm_day = (width/2 + height/2)/2,
  Bact_treat_C = factor(ifelse(Bact_treat %in% c("NCP", "NN"), "NC",
                        ifelse(Bact_treat %in% c("CNP", "CP"), "C",
                               Bact_treat))),
  Bact_treat_P = factor(ifelse(Bact_treat %in% c("NN", "CNP"), "NP",
                        ifelse(Bact_treat %in% c("NCP", "CP"), "P",
                               Bact_treat))),
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

png("./plots/motility.png", width = 8, height = 6, units = "in",
    res = 150)               
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
dev.off()

png("./plots/motility_norm.png", width = 8, height = 5, units = "in",
    res = 150)               
ptwitching_norm <-
  ggplot(data = twitching_isols,
       aes(x = Bact_timepoint, y = cm_day_norm, 
           color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 5, alpha = 0.7) +
  facet_grid(~as.numeric(Bact_pop)) +
  labs(x = "Timepoint", y = "Normalized Twitching Motility (cm/day)",
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
dev.off()

png("./plots/motility_norm_nopops.png", width = 7, height = 5, units = "in",
    res = 150)               
ptwitching_norm2 <-
  ggplot(data = twitching_isols,
         aes(x = Bact_timepoint, y = cm_day_norm, 
             color = Bact_treat_P, fill = Bact_treat_P, shape = Bact_treat_C)) +
  geom_point(size = 4, alpha = 0.5,
             position = position_dodge(width = 0.5)) +
  labs(x = "Timepoint", y = "Normalized Twitching Motility\n(cm/day)") +
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

# Multivariate analyses ----
isol_data <- full_join(
  dplyr::select(mutate(twitching, Bact_pop = as.character(Bact_pop)), 
         !c(Date, width, height)),
  mutate(resis_isols_bactsum_undup,
         Bact_timepoint = as.character(Bact_timepoint),
         Bact_pop = as.character(Bact_pop),
         Bact_isol = as.character(Bact_isol),
         Bact_treat = plyr::revalue(Bact_treat, c("NCNP" = "NN"), 
                                    warn_missing = FALSE)))
isol_data <- full_join(
  isol_data,
  pivot_wider(
    select(
      filter(pig, Bact_pop != "Blank", Wavelength %in% c(405, 695)),
      !c(file, Well, Measurements)),
    id_cols = c("Bact_pop", "Bact_timepoint", "Bact_treat", 
                "Bact_treat_C", "Bact_treat_P", "Bact_isol"),
    names_from = "Wavelength", values_from = Meas_norm,
    names_prefix = "wavelength_"))
if(!"growth_curve_metrics.csv" %in% list.files()) {
  source("Growthcurve_analysis.R")
  setwd("..") #return wd to parent dir where this script needs it
}
isol_data <- full_join(
  isol_data,
  read.csv("growth_curve_metrics.csv", 
           colClasses = c("Bact_pop" = "character", 
                          "Bact_timepoint" = "character",
                          "Bact_isol" = "character")))
isol_data <- mutate(isol_data,
                    resis = 1-clearing_avg)


ggpairs(data = isol_data,
        columns = which(colnames(isol_data) %in% 
                          c("resis", #"cm_day_norm",  
                            "wavelength_405", "wavelength_695",
                            "max_percap", "max_dens", "auc")),
        upper = list(continuous = "points"))

png("./plots/maxdens_405nm_treat.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = isol_data,
       aes(x = max_dens, y = wavelength_405, color = Bact_treat_P)) +
  geom_point() +
  scale_color_manual(name = "Phage\nTreatment", 
                     breaks = c("Anc", "NP", "P", "pilA"),
                     values = my_cols[c(8, 1, 2, 7)],
                     labels = c("Ancestor", "No Phage", "Phage", "pilA"))
dev.off()

png("./plots/maxdens_405nm_time.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = isol_data,
       aes(x = max_dens, y = wavelength_405, color = Bact_timepoint)) +
  geom_point() +
  scale_color_manual(name = "Timepoint", 
                     breaks = c(0, 1, 2, 4, 8, "pilA"),
                     values = my_cols[c(8, 1, 2, 3, 7, 4)])
dev.off()

png("./plots/maxdens_695nm_treat.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = isol_data,
       aes(x = max_dens, y = wavelength_695, color = Bact_treat_P)) +
  geom_point() +
  scale_color_manual(name = "Phage\nTreatment", 
                     breaks = c("Anc", "NP", "P", "pilA"),
                     values = my_cols[c(8, 1, 2, 7)],
                     labels = c("Ancestor", "No Phage", "Phage", "pilA"))
dev.off()

png("./plots/maxdens_695nm_time.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = isol_data,
       aes(x = max_dens, y = wavelength_695, color = Bact_timepoint)) +
  geom_point() +
  scale_color_manual(name = "Timepoint", 
                     breaks = c(0, 1, 2, 4, 8, "pilA"),
                     values = my_cols[c(8, 1, 2, 3, 7, 4)])
dev.off()


png("./plots/motility_norm_sensitivesmarked.png", width = 8, height = 6, units = "in",
    res = 150)
ptwitching_norm +
  geom_point(data = filter(isol_data, clearing_avg > 0, !is.na(cm_day_norm)),
             aes(x = as.factor(Bact_timepoint), y = cm_day_norm),
             pch = 4, color = "black", size = 3)
dev.off()

png("./plots/motility_norm_nopops_sensitivesmarked.png", width = 8, height = 6, units = "in",
    res = 150)
ptwitching_norm2 +
  geom_point(data = filter(isol_data, clearing_avg > 0, !is.na(cm_day_norm)),
             aes(x = Bact_timepoint, y = cm_day_norm,
                 group = paste(Bact_treat_P, Bact_treat_C)),
             pch = 4, color = "black", size = 3,
             position = position_dodge(width = 0.5))
dev.off()


##Run PCA
pca_dat <- isol_data

pca_cols <- which(colnames(pca_dat) %in%
                    c(#"resis", #"cm_day_norm",  
                    "wavelength_405", "wavelength_695",
                    "max_percap", "max_dens", "auc"))

pca_dat <- pca_dat[complete.cases(pca_dat[, pca_cols]), ]

isol_pca <- prcomp(pca_dat[, pca_cols], center = TRUE, scale = TRUE, retx = TRUE)

#Merge pca with data
pca_dat <- cbind(as.data.frame(isol_pca$x),
                 pca_dat)

#Rename rotations for better plotting
replace_names = c("wavelength_405" = "w405",
                  "wavelength_695" = "w695",
                  "max_percap" = "r",
                  "max_dens" = "k")
row.names(isol_pca$rotation) <- plyr::revalue(
  x = row.names(isol_pca$rotation), replace = replace_names)
  
arrow_len <- 3
pca_base <- ggplot(pca_dat, aes(x = PC1, y = PC2)) +
  geom_segment(data = as.data.frame(isol_pca$rotation),
               aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               alpha = .5, lwd = 2, color = "red4") +
  ggrepel::geom_text_repel(data = as.data.frame(isol_pca$rotation),
                           aes(x = arrow_len*PC1, y = arrow_len*PC2,
                               label = row.names(isol_pca$rotation)),
                           size = 8, alpha = .8, color = "red4", seed = 8,
                           min.segment.length = unit(1, "native")) +
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
        plot.title = element_text(size = 36))

png("./plots/PCA.png", width = 10, height = 7, units = "in", res = 150)
pca_base + 
  geom_point(aes(color = Bact_timepoint), size = 5, alpha = 0.7) +
  scale_color_manual(name = "Timepoint", 
                     breaks = c(0, 1, 2, 4, 8),
                     values = my_cols[c(8, 1, 2, 3, 7)])
dev.off()

png("./plots/PCA2.png", width = 11, height = 8, units = "in", res = 150)
pca_base + 
  geom_point(aes(color = Bact_treat_P), size = 5, alpha = 0.7) +
  scale_color_manual(name = "Phage\nTreatment", 
                     breaks = c("Anc", "NP", "P"),
                     values = my_cols[c(8, 1, 2)],
                     labels = c("Ancestor", "No Phage", "Phage"))
dev.off()
                        
                        
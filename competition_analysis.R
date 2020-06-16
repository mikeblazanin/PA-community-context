library(ggplot2)
library(tidyr)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

comp_data_0613 <- read.csv("2020-06-13 Competition.csv",
                           stringsAsFactors = FALSE)

comp_data_0613 <- pivot_longer(comp_data_0613,
                               cols = ends_with(".dil"),
                               names_to = "dilution_media",
                               values_to = "dilution")
comp_data_0613 <- pivot_longer(comp_data_0613,
                               cols = matches(".[cp]fu"),
                               names_to = "cfu_media",
                               values_to = "cfu.pfu")
#Drop the .dil and .cfu
comp_data_0613$dilution_media <- 
  sapply(comp_data_0613$dilution_media,
         FUN = function(x) {strsplit(x, split = "\\.")[[1]][1]})
comp_data_0613$cfu_media <- 
  sapply(comp_data_0613$cfu_media,
         FUN = function(x) {strsplit(x, split = "\\.")[[1]][1]})

#Keep only rows with data
comp_data_0613 <- comp_data_0613[
  comp_data_0613$dilution_media == comp_data_0613$cfu_media &
    complete.cases(comp_data_0613),]

#Rename columns
comp_data_0613$Population <- comp_data_0613$cfu_media
comp_data_0613 <- subset(comp_data_0613, 
                         select = -c(dilution_media, cfu_media))

#Calculate density per mL
comp_data_0613$cfu.pfu.mL <- comp_data_0613$cfu.pfu/(10**comp_data_0613$dilution)

#Calculate frequency
comp_data_0613 <- group_by(comp_data_0613,
                           Bacteria.Treat, Phage.Treat, Timepoint)
comp_data_0613 <- mutate(comp_data_0613,
                         freq = cfu.pfu.mL/sum(cfu.pfu.mL[
                           Population %in% c("Pseudo", "Entero", "Staph")]))
comp_data_0613$freq[comp_data_0613$Population == "H6"] <- NA

#Reorder factors for plotting
comp_data_0613$Population <- factor(comp_data_0613$Population,
                                  levels = c("Pseudo", "Entero", "Staph", "H6"))
comp_data_0613$Bacteria.Treat <- factor(comp_data_0613$Bacteria.Treat,
                                        levels = c("PA", "EF", "SA",
                                                   "PA EF", "PA SA", "EF SA",
                                                   "PA SA EF"))

ggplot(data = comp_data_0613, 
       aes(x = as.factor(Timepoint), y = cfu.pfu.mL, 
           color = as.factor(Phage.Treat), group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

ggplot(data = comp_data_0613, 
       aes(x = as.factor(Timepoint), y = freq, 
           color = as.factor(Phage.Treat), group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  #scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])


library(ggplot2)
library(tidyr)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

comp_data <- list()
comp_data[[1]] <- read.csv("2020-06-13 Competition.csv",
                           stringsAsFactors = FALSE)
comp_data[[2]] <- read.csv("2020-07-27 Competition.csv",
                           stringsAsFactors = FALSE)
names(comp_data) <- c("comp_data_0613", "comp_data_0727") 

for (i in 1:length(comp_data)) {
  #Pivot the dilutions & cfu/pfu counts to be tidy
  comp_data[[i]] <- pivot_longer(comp_data[[i]],
                                 cols = ends_with(".dil"),
                                 names_to = "dilution_media",
                                 values_to = "dilution")
  comp_data[[i]] <- pivot_longer(comp_data[[i]],
               cols = matches(".[cp]fu"),
               names_to = "cfu_media",
               values_to = "cfu.pfu")
  
  #Drop the .dil and .cfu from colnames
  comp_data[[i]]$dilution_media <- 
    sapply(comp_data[[i]]$dilution_media,
           FUN = function(x) {strsplit(x, split = "\\.")[[1]][1]})
  comp_data[[i]]$cfu_media <- 
    sapply(comp_data[[i]]$cfu_media,
           FUN = function(x) {strsplit(x, split = "\\.")[[1]][1]})

  #Drop rows where dilution media & cfu media aren't same
  comp_data[[i]] <- comp_data[[i]][
    comp_data[[i]]$dilution_media == comp_data[[i]]$cfu_media, ]
  
  #Make colony_type column
  if (names(comp_data)[i] == "comp_data_0727") {
    comp_data[[i]]$colony_type <- comp_data[[i]]$Pseudo.colony.type
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "Pseudo" &
                                 comp_data[[i]]$colony_type == "g"] <- "S"
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "Pseudo" &
                                 comp_data[[i]]$colony_type == "y"] <- "W"
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "Staph"] <- "S"
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "Entero"] <- "S"
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "H6"] <- "NA"
    
    comp_data[[i]] <- subset(comp_data[[i]], select = -Pseudo.colony.type)
  } else {
    comp_data[[i]]$colony_type <- "S"
    comp_data[[i]]$colony_type[comp_data[[i]]$cfu_media == "H6"] <- "NA"
  }
  
  #Keep only rows with data
  comp_data[[i]] <- comp_data[[i]][complete.cases(comp_data[[i]]),]
  
  #Rename columns
  comp_data[[i]]$Population <- comp_data[[i]]$cfu_media
  comp_data[[i]] <- subset(comp_data[[i]], 
                           select = -c(dilution_media, cfu_media))

  #Add Contam column if none is there
  if (!"Contam" %in% colnames(comp_data[[i]])) {
    comp_data[[i]]$Contam <- 0
  }
  
  #Calculate density per ml
  ##TODO are the data from the first batch effective dilution or dilution plated!?
  comp_data[[i]]$cfu.pfu.mL <- comp_data[[i]]$cfu.pfu/(10**comp_data[[i]]$dilution)
  if (names(comp_data)[i] == "comp_data_0727") {
    comp_data[[i]]$cfu.pfu.mL <- 100 * comp_data[[i]]$cfu.pfu/(10**comp_data[[i]]$dilution)
  }
  
  #Average out replicate platings
  comp_data[[i]] <- group_by(comp_data[[i]],
                             Bacteria.Treat, Phage.Treat, Timepoint, 
                             Population, colony_type)
  comp_data[[i]] <- summarize(comp_data[[i]],
                              cfu.pfu.mL_avg = mean(cfu.pfu.mL),
                              Contam = max(Contam))
  
  #Reorder factors for plotting
  comp_data[[i]]$Population <- factor(comp_data[[i]]$Population,
                                      levels = c("Pseudo", "Entero", "Staph", "H6"))
  comp_data[[i]]$Bacteria.Treat <- factor(comp_data[[i]]$Bacteria.Treat,
                                          levels = c("P", "S", "E", 
                                                     "PS", "PE", "SE",
                                                     "PSE"))
}

ggplot(data = comp_data[[1]], aes(x = as.factor(Timepoint), 
                                  y = cfu.pfu.mL_avg+1, 
                                  color = as.factor(Phage.Treat), 
                                  group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

ggplot(data = comp_data[[2]], aes(x = as.factor(Timepoint), 
                                  y = cfu.pfu.mL_avg+1, 
                                  color = as.factor(Phage.Treat), 
                                  shape = colony_type)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

#Calculate frequency
comp_data_0613 <- group_by(comp_data_0613,
                           Bacteria.Treat, Phage.Treat, Timepoint)
comp_data_0613 <- mutate(comp_data_0613,
                         freq = cfu.pfu.mL/sum(cfu.pfu.mL[
                           Population %in% c("Pseudo", "Entero", "Staph")]))
comp_data_0613$freq[comp_data_0613$Population == "H6"] <- NA

                                                   

ggplot(data = comp_data_0613, 
       aes(x = as.factor(Timepoint), y = cfu.pfu.mL, 
           color = as.factor(Phage.Treat), group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

ggplot(data = comp_data_0613[comp_data_0613$Timepoint == 1, ], 
       aes(x = as.factor(Timepoint), y = cfu.pfu.mL, 
           color = as.factor(Phage.Treat), group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

ggplot(data = comp_data_0613[comp_data_0613$Population != "H6", ], 
       aes(x = as.factor(Timepoint), y = freq, 
           color = as.factor(Phage.Treat), group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  #scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])


#Load ----
library(ggplot2)
library(tidyr)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Get data & analyze ----
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
  
  #Calculate total populations
  if (i == 1) {comp_data_sum <- list()}
  comp_data[[i]] <- group_by(comp_data[[i]],
                             Bacteria.Treat, Phage.Treat, Timepoint, 
                             Population)
  comp_data_sum[[i]] <- summarize(comp_data[[i]],
                                  cfu.pfu.mL_avg_all = sum(cfu.pfu.mL_avg),
                                  Contam = max(Contam))
}

#Incorporate T0 for 07-2020 competition ----
t0_0727 <- read.csv("2020-07-27 Competition T0.csv",
                    stringsAsFactors = F)
t0_0727$cfu.pfu.ml <- 100*t0_0727$cfu/10**t0_0727$Dilution
t0_0727 <- group_by(t0_0727, Bact)
t0_0727_sum <- summarize(t0_0727,
                         cfu.pfu.ml_avg = mean(cfu.pfu.ml))

t0_data <- data.frame(Bacteria.Treat = c("P", "P", "S", "E",
                                         "PS", "PS", "PE", "PE", "SE",
                                         "PSE", "PSE"),
                      Phage.Treat = c(1, 0, 0, 0,
                                      1, 0, 1, 0, 0,
                                      1, 0),
                      Timepoint = 0,
                      Pseudo_dens = NA,
                      Staph_dens = NA,
                      Entero_dens = NA,
                      H6_dens = NA,
                      Contam = 0,
                      stringsAsFactors = F)

for (i in 1:nrow(t0_data)) {
  if ("P" %in% strsplit(t0_data$Bacteria.Treat[i], split = "")[[1]]) {
    t0_data$Pseudo_dens[i] <- 
      0.1 * t0_0727_sum$cfu.pfu.ml_avg[t0_0727_sum$Bact == "Pseudo"]/
      nchar(t0_data$Bacteria.Treat[i])
  }
  if ("S" %in% strsplit(t0_data$Bacteria.Treat[i], split = "")[[1]]) {
    t0_data$Staph_dens[i] <- 
      0.1 * t0_0727_sum$cfu.pfu.ml_avg[t0_0727_sum$Bact == "Staph"]/
      nchar(t0_data$Bacteria.Treat[i])
  }
  if ("E" %in% strsplit(t0_data$Bacteria.Treat[i], split = "")[[1]]) {
    t0_data$Entero_dens[i] <- 
      0.1 * t0_0727_sum$cfu.pfu.ml_avg[t0_0727_sum$Bact == "Entero"]/
      nchar(t0_data$Bacteria.Treat[i])
  }
  if (t0_data$Phage.Treat[i] == 1) {
    t0_data$H6_dens[i] <- 
      0.01 * t0_0727_sum$cfu.pfu.ml_avg[t0_0727_sum$Bact == "H6"]
  }
}

#Melt T0 data
t0_tidy <- pivot_longer(t0_data, cols = ends_with("dens"),
                        names_to = "Population", values_to = "cfu.pfu.mL_avg_all",
                        values_drop_na = TRUE)

#Drop _dens from Pop names
t0_tidy$Population <- gsub(pattern = "_dens", replacement = "", 
                    x = t0_tidy$Population)

#Merge with later data
comp_data_sum[[2]] <- full_join(t0_tidy, comp_data_sum[[2]])

#Calculate frequency ----
for (i in 1:length(comp_data_sum)) {
  comp_data_sum[[i]] <- group_by(comp_data_sum[[i]],
                                 Bacteria.Treat, Phage.Treat, Timepoint)
  comp_data_sum[[i]] <- mutate(comp_data_sum[[i]],
                               freq = cfu.pfu.mL_avg_all/sum(cfu.pfu.mL_avg_all[
                                 Population %in% c("Pseudo", "Entero", "Staph")]))
  comp_data_sum[[i]]$freq[comp_data_sum[[i]]$Population == "H6"] <- NA


  #Factor tweaking for plotting
  comp_data_sum[[i]]$Population <- plyr::revalue(comp_data_sum[[i]]$Population, 
                                   c("Pseudo" = "P. aeruginosa",
                                     "Entero" = "E. faecalis",
                                     "Staph" = "S. aureus",
                                     "H6" = "Phage H6"))
  comp_data_sum[[i]]$Population <- factor(comp_data_sum[[i]]$Population,
                          levels = c("P. aeruginosa", "E. faecalis", 
                                     "S. aureus", "Phage H6"))
  comp_data_sum[[i]]$Phage.Treat <- factor(comp_data_sum[[i]]$Phage.Treat)
  comp_data_sum[[i]]$Phage.Treat <- plyr::revalue(comp_data_sum[[i]]$Phage.Treat, 
                                  c("1" = "Present",
                                    "0" = "Absent"))
  comp_data_sum[[i]]$Bacteria.Treat <- factor(comp_data_sum[[i]]$Bacteria.Treat,
                              levels = c("P", "S", "E", "PS", "PE", "SE", "PSE"))
}


#Plots ----
#Competition #1
ggplot(data = comp_data_sum[[1]], aes(x = as.factor(Timepoint), 
                                  y = cfu.pfu.mL_avg_all+1, 
                                  color = as.factor(Phage.Treat), 
                                  group = Phage.Treat)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)])

#Competition #2, all data
ggplot(data = comp_data_sum[[2]], 
       aes(x = as.factor(Timepoint), 
           y = cfu.pfu.mL_avg_all+1, 
           color = as.factor(Phage.Treat),
           group = Phage.Treat,
           shape = as.factor(Contam))) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)]) +
  labs(x = "Timepoint", y = "Density (cfu or pfu/mL)")

#Competition #2, Non-contaminated data
temp <- comp_data_sum[[2]][comp_data_sum[[2]]$Contam == 0, ]

#Density
ggplot(data = temp, 
       aes(x = as.factor(Timepoint), 
           y = cfu.pfu.mL_avg_all+1, 
           color = as.factor(Phage.Treat),
           group = Phage.Treat)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)]) +
  labs(x = "Timepoint", y = "Density (cfu or pfu/mL)",
       subtitle = "Bacterial Community Combination",
       color = "Phage") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  NULL
  
#Frequency (points)
ggplot(data = temp[temp$Population != "Phage H6", ],
       aes(x = as.factor(Timepoint),
           y = freq,
           color = as.factor(Phage.Treat),
           group = as.factor(Phage.Treat))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line() +
  scale_y_continuous(trans = "log10") +
  facet_grid(Population ~ Bacteria.Treat) +
  theme_bw() +
  scale_color_manual(values = my_cols[c(1, 8)]) +
  labs(x = "Timepoint", y = "Frequency",
       subtitle = "Bacterial Community Combination",
       color = "Phage") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  NULL

#Frequency (stacked bars)
ggplot(data = temp[temp$Population != "Phage H6", ],
       aes(x = as.factor(Timepoint), y = freq)) +
  geom_col(aes(fill = Population)) +
  #scale_y_continuous(trans = "log10") +
  facet_grid(Phage.Treat ~ Bacteria.Treat) +
  theme_bw() +
  scale_fill_manual(values = my_cols[c(3, 1, 2)]) +
  labs(x = "Timepoint", y = "Frequency",
       subtitle = "Bacterial Community Combination",
       fill = "Bacterial species") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  NULL

## Old Plots ----
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


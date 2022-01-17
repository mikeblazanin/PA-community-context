library(ggplot2)
library(tidyr)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#PA14 H6 coevol ----
evol_dens <- read.csv("2021-07-28_PA14H6_Coevol.csv", stringsAsFactors = F)
colnames(evol_dens)[5:13] <- paste("-", 0:8, sep = "")

#Pivot longer
evol_dens_lng <- as.data.frame(
  tidyr::pivot_longer(evol_dens, cols = starts_with("-"),
                      names_to = "Dilution", values_to = "count"))
evol_dens_lng <- evol_dens_lng[evol_dens_lng$count != "", ]

#For pop +2 at timepoint 1 bc density was too high to measure estimate it
# as 100 cfu at -8
evol_dens_lng$count[evol_dens_lng$Plating == "PA" & evol_dens_lng$Treat == "+" &
              evol_dens_lng$Population == 2 & evol_dens_lng$Timepoint == 1 &
                evol_dens_lng$Dilution == "-8"] <- 100

#Calculate densities & frequencies
evol_dens_lng <- group_by(evol_dens_lng, Timepoint, Plating, Treat, Population)
evol_dens_sum <- summarize(
  evol_dens_lng,
  dens_dil = max(as.numeric(Dilution[!(count %in% c("clumped", "lacy", "mbc", "too")) &
                            !grepl("~", count)])),
  freq_dil = max(as.numeric(Dilution[!(count %in% c("clumped", "lacy", "mbc", "too"))])),
  dens_val = count[as.numeric(Dilution) == dens_dil],
  freq_val = count[as.numeric(Dilution) == freq_dil],
  dens = NA, freq_L = NA)
for (i in 1:nrow(evol_dens_sum)) {
  #First calculate density
  if (!grepl("[,LS]", evol_dens_sum$dens_val[i])) {
    #If there's not a comma, L, or S, just calculate density (frequency of L is 1)
    evol_dens_sum$dens[i] <- 
      100*as.numeric(evol_dens_sum$dens_val[i])/10**evol_dens_sum$dens_dil[i]
  } else {
    #If there's a comma, split, sum, calculate density
    temp_d <- strsplit(evol_dens_sum$dens_val[i], split = ",")[[1]]
    temp_d <- gsub("[SL]", "", temp_d)
    evol_dens_sum$dens[i] <-
      100*sum(as.numeric(temp_d))/10**evol_dens_sum$dens_dil[i]
  }
  #Now calculate frequency
  if (!grepl("[,LS]", evol_dens_sum$freq_val[i])) {
    #If there's not a comma, L, or S, frequency of L is 1
    evol_dens_sum$freq_L[i] <- 1
  } else {
    #If there's a comma, L, or S, split, calculate frequency
    temp_f <- strsplit(evol_dens_sum$freq_val[i], split = ",")[[1]]
    temp_f <- gsub("~", "", temp_f)
    types <- substr(temp_f, nchar(temp_f), nchar(temp_f))
    temp_f <- substr(temp_f, 1, nchar(temp_f)-1)
    if (length(types) == 1) {
      if (types[1] == "L") {evol_dens_sum$freq_L[i] <- 1
      } else if (types[1] == "S") {evol_dens_sum$freq_L[i] <- 0
      }
    } else {
      evol_dens_sum$freq_L[i] <- 
        as.numeric(temp_f[which(types == "L")])/sum(as.numeric(temp_f))
    }
  }
}

#Add pseudo-data on initial inoculation densities
evol_dens_sum <- 
  full_join(evol_dens_sum,
            data.frame("dens" = c(rep(10**7, 12), rep(10**5, 6)),
                       "Timepoint" = rep(0, 18),
                       "Plating" = c(rep("PA", 12), rep("H6", 6)),
                       "Treat" = c(rep("-", 6), rep("+", 12)),
                       "Population" = rep(1:6, 3),
                       "freq_L" = rep(1, 18)))
                          

                               




#Community co-passaging ----
comm_data <- read.csv("2021-12-17_community_experiment.csv",
                      stringsAsFactors = FALSE)

#Split PA data to calculate freq_L and total count
# (note that for now we've summed in when we counted frequency from
# a lower dilution because it only throws off the density by 1 or 2 cfu
# but not accurate)
comm_data$freq_L <- NA
for (i in which(comm_data$Strain == "PA")) {
  #Check if there are small/large freq's
  if(grepl("[,LS]", comm_data$Count[i])) {
    temp_d <- strsplit(comm_data$Count[i], split = ",")[[1]]
    temp_e <- as.numeric(gsub("[SL]", "", temp_d))
    names(temp_e) <- substr(temp_d, nchar(temp_d), nchar(temp_d))
    
    comm_data$Count[i] <- sum(temp_e)
    if (length(temp_e) > 1) {
      comm_data$freq_L[i] <- temp_e["L"]/sum(temp_e)
    } else {
      if(names(temp_e) == "L") {
        comm_data$freq_L[i] <- 1
      } else if (names(temp_e) == "S") {
        comm_data$freq_L[i] <- 0
      }
    }
  } else {
    comm_data$freq_L[i] <- 1
  }
}

#Because I think the small colony morphs that appeared in the
# +phage +competitors on Day 6 were probably mis-counts, set
# those frequencies to 1
comm_data$freq_L[comm_data$Phage_treat == "+" & 
                   comm_data$Timepoint == 6 &
                   comm_data$Strain == "PA"] <- 1

#Calculate density
comm_data$dens_ml <- (1000/2.5) * 
  as.numeric(comm_data$Count)/(10**as.numeric(comm_data$Dilution))

#Fix T0 data by pseudo-replicating for plotting
comm_data <- 
  rbind(
    do.call("rbind", replicate(12, comm_data[comm_data$Timepoint == 0, ],
                               simplify = FALSE)),
    comm_data[comm_data$Timepoint != 0, ])
comm_data$Phage_treat[1:(12*5)] <- rep(c("+", "-"), each = 6*5)
comm_data$Pop[1:(12*5)] <- rep(rep(1:6, each = 5), 2)
comm_data$Count[comm_data$Timepoint == 0 & comm_data$Strain == "H6" &
                  comm_data$Phage_treat == "-"] <- NA
comm_data$dens_ml[comm_data$Timepoint == 0 & comm_data$Strain == "H6" &
                    comm_data$Phage_treat == "-"] <- NA
comm_data <- comm_data[!is.na(comm_data$dens_ml), ]

#Merge datasets ----
evol_dens_sum <- rename(evol_dens_sum,
                        Strain = Plating,
                        Phage_treat = Treat,
                        dens_ml = dens)
evol_dens_sum <- evol_dens_sum[, !names(evol_dens_sum) %in% 
                                 c("dens_dil", "freq_dil", "dens_val", "freq_val")]
comm_data <- rename(comm_data,
                    Population = Pop)
comm_data <- comm_data[, !names(comm_data) %in% 
                                 c("Count", "Dilution")]
evol_dens_sum <- cbind(data.frame(Comm_treat = "-"),
                       evol_dens_sum)
comm_data <- cbind(data.frame(Comm_treat = "+"),
                   comm_data)
dens_data <- full_join(evol_dens_sum, comm_data)

#Set counts below detection to limit of detection
dens_data$bd <- ifelse(dens_data$dens_ml == 0, TRUE, FALSE)
for (i in 1:nrow(dens_data)) {
  if (dens_data$dens_ml[i] == 0) {
    if(dens_data$Comm_treat[i] == "+") {
      dens_data$dens_ml[i] <- 1000/2.5
    } else if (dens_data$Comm_treat[i] == "-") {
      dens_data$dens_ml[i] <- 1000/5
    }
  }
}

#Summarize pops for treat-level means
dens_data_sum <- 
  summarize(group_by(dens_data, Comm_treat, Phage_treat, 
                     Timepoint, Strain),
            dens_ml = 10**mean(log10(dens_ml+1)), freq_L = mean(freq_L))


#Plots ----
png("dens_copassaging.png", width = 5, height = 5, units = "in", res = 150)
ggplot(data = dens_data,
       aes(x = Timepoint, y = dens_ml+1, color = Strain)) +
  #geom_point() + 
  geom_line(aes(group = paste(Population, Strain)), alpha = 0.5) +
  geom_line(data = dens_data_sum, lwd = 2, alpha = 0.8) +
  geom_hline(data = data.frame(Comm_treat = c("+", "+", "-", "-"),
                               Phage_treat = c("+", "-", "+", "-"),
                               limit = c(1000/2.5, 1000/2.5, 1000/5, 1000/5)),
             aes(yintercept = limit), lty = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "Competitors Added",
                                                "-" = "No Competitors"))) +
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
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  NULL
dev.off()

png("colony_morph_freq.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = dens_data[dens_data$Strain == "PA", ],
       aes(x = Timepoint, y = freq_L)) +
  geom_line(aes(group = paste(Population, Strain))) +
  geom_line(data = dens_data_sum[dens_data_sum$Strain == "PA", ], lwd = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "Competitors Added",
                                                "-" = "No Competitors"))) +
  xlim(0, NA) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Frequency of Large Colony Morphology") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text =  element_text(size = 12.5)) +
  NULL
dev.off()

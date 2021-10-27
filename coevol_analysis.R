library(ggplot2)
library(tidyr)
library(dplyr)

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
                          
evol_dens_sum_sum <- 
  summarize(group_by(evol_dens_sum, Timepoint, Plating, Treat),
            dens = 10**mean(log10(dens+1)), freq_L = mean(freq_L))
                               

ggplot(data = evol_dens_sum,
       aes(x = Timepoint, y = dens+1, color = Plating)) +
  #geom_point() + 
  geom_line(aes(group = paste(Population, Plating))) +
  geom_line(data = evol_dens_sum_sum, lwd = 2) +
  facet_grid(. ~ Treat,
             labeller = as_labeller(c("+" = "Phage Added",
                                      "-" = "No Phage"))) +
  scale_y_continuous(trans = "log10") +
  xlim(0, NA) +
  labs(y = "Density (cfu/mL or pfu/mL)") +
  theme_bw() +
  NULL

ggplot(data = evol_dens_sum[evol_dens_sum$Plating == "PA", ],
       aes(x = Timepoint, y = freq_L)) +
  geom_line(aes(group = paste(Population, Plating))) +
  geom_line(data = evol_dens_sum_sum[evol_dens_sum_sum$Plating == "PA", ], lwd = 2) +
  facet_grid(. ~ Treat,
             labeller = as_labeller(c("+" = "Phage Added",
                                   "-" = "No Phage"))) +
  xlim(0, NA) +
  labs(y = "Frequency of Large Colony Morphology") +
  theme_bw() +
  NULL

ggplot(data = evol_dens_sum[evol_dens_sum$Plating == "PA", ],
       aes(x = dens, y = freq_L)) +
  geom_line(aes(group = Population, color = as.factor(Population))) +
  scale_x_continuous(trans = "log10") +
  facet_grid(Treat ~ .) +
  labs(y = "Frequency of Large Colony Morphology") +
  NULL

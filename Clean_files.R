## Setup ----
library(dplyr)
library(tidyr)
library(gcplyr)

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
prelimdens_NC <- read.csv("./raw_data/2021-07-28_PA14H6_Coevol.csv", stringsAsFactors = F)
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

write.csv(prelimdens_NC_sum, "./clean_data/prelimdens_NC_sum.csv",
          row.names = FALSE)

## Prelim experiment: population data from run with competitors ----
prelimdens_C <- read.csv("./raw_data/2021-12-17_community_experiment.csv",
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

write.csv(prelimdens_C, "./clean_data/prelimdens_C.csv", row.names = FALSE)

## Population data ----

#Read file
density <- read.csv("./raw_data/2022-09-19_community_experiment_2.csv",
                    stringsAsFactors = FALSE, na.strings = "")
density <- filter(density, !is.na(Timepoint))

#Adjust for consistent style across data
density <- mutate(density,
                  treat_P = ifelse(Phage_treat == "+", "P", "NP"),
                  treat_C = ifelse(Comm_treat == "+", "C", "NC"))

#Cross ancestral rows into each other treatment
density <- dup_rows(mydf = density, from_vals = NA, to_vals = c("P", "NP"),
                    mut_col = "treat_P")
density <- dup_rows(mydf = density, from_vals = NA, to_vals = c("C", "NC"),
                    mut_col = "treat_C")


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
                  cfu_ml = case_when(
                    Timepoint != 0 ~ cfu_ml,
                    Timepoint == 0 & Strain == "H6" & treat_P == "NP" ~ 0,
                    Timepoint == 0 & Strain == "H6" & treat_P == "P" ~
                      10/10000 * Count/10**Dilution,
                    Timepoint == 0 & (Strain == "PA" | treat_C == "C") ~ 
                      100/10000 * Count/10**Dilution,
                    Timepoint == 0 & Strain != "PA" & treat_C != "C" ~ 0),
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
                                  cfu_ml),
                  log_cfuml = log10(cfu_ml))
density <- as.data.frame(density)

write.csv(density, "./clean_data/density.csv", row.names = FALSE)

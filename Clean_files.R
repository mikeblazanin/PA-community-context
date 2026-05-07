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

# Resistance data ----

#Resistance data consists of 4 datasets:
#  A first set of combinations where we used EOP for everything
#  A second set of combinations where we used spot tests for everything
#  A third set of combinations where we used cross-streaks for everything
#  A fourth set where we used EOP again (EOP2) to re-test anything putatively
#   sensitive in the spots/cross-streaks

### EOP data ----
EOP <- read.csv("./raw_data/2023-05-20 EOP Timeshift data.csv", stringsAsFactors = F)

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
              !((Bact_pop == 9 & Bact_treat == "NCNP" & Bact_timepoint >= 2) |
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

### EOP data II ----
EOP2 <- read.csv("./raw_data/2024-04-11_EOP_2_results.csv", stringsAsFactors = F)

#rename columns for eventual merging
EOP2 <- mutate(EOP2, pfu2 = pfu, dilution2 = dilution, 
               Estimated_pfu2 = Estimated_pfu)
EOP2 <- select(EOP2, -pfu, -dilution, -Estimated_pfu, -Bacteria, -Phage)

#ID below-detection counts, and calculate pfu/mL (note 1.5 uL spots)
EOP2 <- mutate(EOP2, 
               bd2 = ifelse(pfu2 == 0, TRUE, FALSE),
               pfu2 = ifelse(pfu2 == 0, 1, pfu2),
               pfu_ml2 = 1000/1.5 * pfu2 / 10**dilution2)

#Deal with T0 labels
# When it's T0 with a specific treat, the T0 should be matched
#  to that treat name
# When it's T0 with T0, those rows need to be duplicated across all
#  the treats so they plot for all four treats 
EOP2 <- mutate(EOP2,
               Bact_treat = ifelse(Bact_treat == "Anc" & Phage_treat != "Anc",
                                   Phage_treat, Bact_treat),
               Phage_treat = ifelse(Phage_treat == "Anc" & Bact_treat != "Anc",
                                    Bact_treat, Phage_treat))
EOP2 <- dup_rows(EOP2, "Bact_treat", "Anc", c("CP", "NCP"))
EOP2 <- mutate(EOP2,
               Phage_treat = ifelse(Phage_treat == "Anc", Bact_treat, Phage_treat))
EOP2 <- separate_wider_regex(EOP2,
                             cols = "Bact_treat",
                             patterns = c("Bact_treat_C" = "^N?C", 
                                          "Bact_treat_P" = "N?P$"),
                             cols_remove = FALSE)
EOP2 <- separate_wider_regex(EOP2,
                             cols = "Phage_treat",
                             patterns = c("Phage_treat_C" = "^N?C", 
                                          "Phage_treat_P" = "N?P$"),
                             cols_remove = FALSE)
EOP2$Bact_treat_C <- relevel(as.factor(EOP2$Bact_treat_C), ref = "NC")

#Calculate EOP (relative to highest pfu/mL for ea phage stock)
EOP2 <- mutate(group_by(EOP2, 
                        Phage_pop, Phage_treat, 
                        Phage_timepoint, Phage_isol),
               EOP2 = pfu_ml2/max(pfu_ml2),
               EOP2_plot = EOP2,
               clearing2 = ifelse(bd2 == TRUE, 0, 1))
EOP2$EOP2_plot[EOP2$EOP2_plot < max(EOP2$EOP2[EOP2$bd2 == TRUE], na.rm = TRUE)] <-
  max(EOP2$EOP2[EOP2$bd2 == TRUE], na.rm = TRUE)

### Spot test data ----
spots <- read.csv("./raw_data/2023-05-24 Spot Timeshift data.csv", stringsAsFactors = F)

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
                !((Bact_pop == 9 & Bact_treat == "NCNP" & Bact_timepoint >= 2) |
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
cs_isols <- read.csv("./raw_data/William cross-streak data.csv", fileEncoding = "UTF-8")

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
resis_isols <- 
  Reduce(
    full_join,
    list(mutate(EOP2, run_EOP2 = TRUE), 
         mutate(cs_isols, run_cs = TRUE),
         mutate(spots_isols, run_spot = TRUE),
         mutate(EOP_isols, run_EOP = TRUE)),
    right = TRUE)
resis_isols <- mutate(ungroup(resis_isols),
                      Bact_treat = factor(Bact_treat),
                      Bact_treat_C = relevel(factor(Bact_treat_C), "NC"),
                      Bact_treat_P = factor(Bact_treat_P))

#Pull earlier results into gaps in later results
resis_isols <- mutate(
  resis_isols,
  clearing2 = ifelse(is.na(clearing2), clearing, clearing2),
  EOP2_plot = ifelse(is.na(EOP2_plot), 
                     ifelse(clearing2 == 0, min(EOP2_plot, na.rm = TRUE), NA),
                     #Note: Pop 3 T0 #3 couldn't be retested so it gets set to NA here
                     EOP2_plot),
  Plaque_morphology = ifelse(is.na(Plaque_morphology), 
                             ifelse(clearing2 == 0, "None", "Clear?"),
                             Plaque_morphology))

#Write copy to file
write.csv(resis_isols, "./clean_data/resis_isols.csv", row.names = FALSE)

# Pigments ----

pig1 <- import_blockmeasures(
  list.files(path = "./raw_data/pigment_data/", full.names = T, pattern = ".*Production.csv$"),
  startrow = 42, endrow = 50, endcol = "M",
  metadata = list("Wavelength" = c(41, "A")))
pig2 <- import_blockmeasures(
  list.files(path = "./raw_data/pigment_data/", full.names = T, pattern = ".*Production.csv$"),
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
  list.files(path = "./raw_data/pigment_data/", full.names = T, pattern = ".*Production.csv$"),
  startrow = 69, endrow = 90)
pigspectrum <- trans_wide_to_tidy(pigspectrum, id_cols = c("file", "Wavelength"),
                                  values_to_numeric = FALSE)

pig <- merge_dfs(pig12, pigspectrum, collapse = TRUE)
pig$file <- regmatches(pig$file, regexpr("\\d{4}-\\d{2}-\\d{2}", pig$file))

#Read designs
designs <-
  read_blocks(
    list.files(path = "./raw_data/pigment_data/", full.names = T, pattern = ".*Design.csv$"))
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
              !((Bact_pop == 9 & Bact_treat == "NN" & Bact_timepoint >= 2) |
                  (Bact_pop == 9 & Bact_treat == "CNP" & Bact_timepoint >= 2)))

#Subtract blank
pig$Measurements[pig$Measurements == "OVRFLW"] <- NA
pig$Measurements <- as.numeric(pig$Measurements)
pig <- mutate(group_by(pig, file, Wavelength),
              Meas_norm = Measurements - Measurements[Bact_pop == "Blank"])
pig <- filter(pig, Bact_isol != "Blank")

#Write copy to file
write.csv(pig, "./clean_data/pig.csv", row.names = FALSE)

# Biofilm ----

#Data
biofilm <- import_blockmeasures(
  list.files(path = "./raw_data/biofilm_data/", full.names = T, pattern = ".*Data.csv$"),
  startrow = 26, endrow = 34, endcol = "M")
biofilm <- trans_wide_to_tidy(biofilm, id_cols = "block_name",
                              values_to = "biofilm_abs")
biofilm$block_name <- gsub(" Data$", "", biofilm$block_name)

#Designs
biofilm_designs <- import_blockdesigns(
  list.files(path = "./raw_data/biofilm_data/", full.names = T, pattern = ".*Design.csv$"),
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

write.csv(biofilm, "./clean_data/biofilm.csv", row.names = FALSE)
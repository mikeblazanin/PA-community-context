library(ggplot2)
library(dplyr)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

## Population data ----

#TODO - reconfirm that all points labeled as 0's are truly 0's
#       and not mis-platings

#Read file
density <- read.csv("2022-09-19_community_experiment_2.csv")
density <- filter(density, !is.na(Timepoint))

#Calculate total density, large freq
# for total density, note counts are for 2.5 uL
density <- mutate(group_by(density, Timepoint, Phage_treat, Comm_treat, Pop, Strain),
              cfu_ml = ifelse(!is.na(Count),
                              400*Count/10**Dilution,
                              400*sum(Small_count/10**Small_dilution, 
                                      Large_count/10**Large_dilution, 
                                      na.rm = TRUE)),
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

density_sum <- summarize(group_by(density,
                              Timepoint, Phage_treat, Comm_treat, Strain),
                     cfu_ml = 10**mean(log10(cfu_ml), na.rm = TRUE),
                     L_freq = mean(L_freq, na.rm = TRUE))

png("William_dens.png", width = 5, height = 5, units = "in",
    res = 150)
ggplot(data = density[density$Pop %in% 3:10 & !is.na(density$cfu_ml), ],
       aes(x = Timepoint, y = cfu_ml, color = Strain)) +
  #geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
  geom_line(data = density_sum, lwd = 2, alpha = 0.8) +
  geom_hline(yintercept = 400, lty = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "Competitors Added",
                                                "-" = "No Competitors"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), limits = c(0, NA)) +
  scale_y_continuous(trans = "log10",
                     breaks = 10**c(2, 6, 10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
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

png("William_morph_freq.png", width = 5, height = 4, units = "in", res = 150)
ggplot(data = density[density$Strain == "PA", ],
       aes(x = Timepoint, y = L_freq)) +
  geom_line(aes(group = paste(Pop, Strain))) +
  geom_line(data = density_sum[density_sum$Strain == "PA", ], lwd = 2) +
  facet_grid(Comm_treat ~ Phage_treat,
             labeller = labeller(Phage_treat = c("+" = "Phage Added",
                                                 "-" = "No Phage"),
                                 Comm_treat = c("+" = "Competitors Added",
                                                "-" = "No Competitors"))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9), limits = c(0, NA)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "Frequency of Large Colony Morphology") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text =  element_text(size = 12.5)) +
  NULL
dev.off()


## EOP data ----

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
EOP <- dup_rows(EOP, "Bact_treat", "Anc", c("CP", "NCP", "CNP", "NCNP"))
EOP <- mutate(EOP,
              Phage_treat = ifelse(Phage_treat == "Anc", Bact_treat, Phage_treat))

#average technical replicates together (should just be Ancestors)
EOP_isols <- summarize(group_by(EOP, Bact_pop, Bact_treat, Bact_timepoint,
                                Phage_pop, Phage_treat, Phage_timepoint),
                       bd = any(bd == TRUE),
                       pfu_ml = 10**(mean(log10(pfu_ml))))

#Calculate EOP (relative to highest pfu/mL for ea phage stock)
EOP_isols <- mutate(group_by(EOP_isols, Phage_pop, Phage_treat, Phage_timepoint),
              EOP = pfu_ml/max(pfu_ml))

#Plot
ggplot(data = filter(EOP_isols, Phage_pop == 5),
       aes(x = Bact_timepoint, y = EOP, color = as.factor(Phage_timepoint))) +
  geom_point(aes(shape = bd)) + 
  geom_line(aes(group = paste(Phage_timepoint, Phage_pop))) +
  facet_wrap(~Bact_treat) +
  scale_y_log10() +
  scale_shape_manual(values = c(16, 4), breaks = c(FALSE, TRUE)) +
  guides(shape = "none")

ggplot(data = EOP_isols,
       aes(x = as.factor(Bact_timepoint), 
           y = as.factor(Phage_timepoint))) +
  geom_tile(aes(fill = log10(EOP))) +
  facet_grid(Phage_treat ~ Phage_pop) +
  labs(x = "Bacterial timepoint", y = "Phage timepoint")


       
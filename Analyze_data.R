library(ggplot2)
library(dplyr)

EOP <- read.csv("2023-05-20 EOP Timeshift data.csv", stringsAsFactors = F)

#ID below-detection counts, and calculate pfu/mL
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


       
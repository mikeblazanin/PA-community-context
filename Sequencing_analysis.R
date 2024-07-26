library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

mutations <- read.table("mutations_compare_all.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE, quote = "")

mutations_wide <- pivot_wider(mutations,
                              names_from = title, values_from = frequency)

mutations_wide$mut_id <- 1:nrow(mutations_wide)

mutations_wide <- mutate(group_by(mutations_wide, mut_id),
                         num_anc = sum(!is.na(c(`01_B_T0_3`, `02_B_T0_4`,
                                              `03_B_T0_5`, `04_B_T0_6`,
                                              `05_B_T0_7`, `06_B_T0_8`,
                                              `07_B_T0_9`, `08_B_T0_10`))))

mutations_long <- pivot_longer(mutations_wide,
                               names_to = "title",
                               values_to = "frequency",
                               cols = `01_B_T0_3`:`58_P_T8_CP10`)
mutations_long <- filter(mutations_long, !is.na(frequency))
mutations_long <- mutate(mutations_long,
                         title = ifelse(title == "42_P_T0", "42_P_T0_Anc", title))

mutations_long <- separate_wider_delim(
  mutations_long, title, "_",
  names = c("Samp_id", "Samp_type", "Timepoint", "TreatPop"))
mutations_long <- mutate(
  mutations_long,
  Timepoint = gsub("T", "", Timepoint),
  Treat = ifelse(Timepoint == 0, "Anc", str_extract(TreatPop, "CP|NCP|CNP|NN")),
  Pop = str_extract(TreatPop, "[:digit:]{1,2}"))

mut_sum <- summarize(group_by(mutations_long, Samp_type, Timepoint, Treat, Pop, type),
                     n_type = n(),
                     fixed = sum(frequency == 1),
                     polymorphic = sum(frequency > 0 & frequency < 1))

png("./plots/muts_per_pop.png", width = 7, height = 4, units = "in",
    res = 150)              
ggplot(filter(mut_sum, Samp_type == "B"),
       aes(x = as.numeric(Pop), y = n_type, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(~ Treat,
             labeller = labeller(
               Treat = c(Anc = "Ancestor", 
                         "NN" = "Only Pseudomonas\nNo Phage",
                         "CNP" = "All Bacteria\nNo Phage",
                         "NCP" = "Only Pseudomonas\nPhage Added",
                         "CP" = "All Bacteria\nPhage Added"))) +
  labs(x = "Population", y = "Number of mutations")
dev.off()

png("./plots/muts_in_ancestors.png", width = 6, height = 6, units = "in",
    res = 150)
ggplot(filter(mutations_long, Samp_type == "B"), 
       aes(x = as.numeric(Pop), y = frequency)) +
  geom_point(position = position_jitter(width = 0.1),
             aes(color = num_anc > 1),
             alpha = 0.5) +
  facet_grid(Treat ~ .,
             labeller = labeller(
               Treat = c(Anc = "Ancestor", 
                         "NN" = "Only Pseudomonas\nNo Phage",
                         "CNP" = "All Bacteria\nNo Phage",
                         "NCP" = "Only Pseudomonas\nPhage Added",
                         "CP" = "All Bacteria\nPhage Added"))) +
  scale_color_manual("Appears in\nMultiple\nAncestors",
                     breaks = c(TRUE, FALSE),
                     values = my_cols[c(2, 8)]) +
  labs(x = "Population")
dev.off()

ggplot(filter(mutations_long, Samp_type == "B", num_anc <= 1), 
       aes(x = as.numeric(Pop), y = frequency)) +
  geom_point(position = position_jitter(width = 0.1),
             aes(color = num_anc > 0),
             alpha = 0.5) +
  facet_grid(Treat ~ .)
  #scale_y_log10() +
  # scale_color_manual("Appears in\nMultiple\nAncestors",
  #                    breaks = c(TRUE, FALSE),
  #                    values = my_cols)

mut_sum2 <- summarize(
  group_by(filter(mutations_long, Samp_type == "B"),
          across(c(-Samp_id, -Samp_type, -Timepoint, -TreatPop,
                                    -frequency, -Treat, -Pop))),
  total_freq = sum(frequency),
  max_freq = max(frequency),
  nhits = n(),
  total_freq_NN = sum(frequency[Treat == "NN"]),
  total_freq_NCP = sum(frequency[Treat == "NCP"]),
  total_freq_CNP = sum(frequency[Treat == "CNP"]),
  total_freq_CP = sum(frequency[Treat == "CP"]),
  totalfreq_diff = total_freq_CP + total_freq_NCP - total_freq_CNP - total_freq_NN,
  totalfreq_frac = (total_freq_CP + total_freq_NCP)/total_freq,
  max_freq_NN = max(c(0, frequency[Treat == "NN"])),
  max_freq_NCP = max(c(0, frequency[Treat == "NCP"])),
  max_freq_CNP = max(c(0, frequency[Treat == "CNP"])),
  max_freq_CP = max(c(0, frequency[Treat == "CP"])),
  max_freq_P = max(max_freq_NCP, max_freq_CP),
  max_freq_NP = max(max_freq_CNP, max_freq_NN),
  max_freq_C = max(max_freq_CP, max_freq_CNP),
  max_freq_NC = max(max_freq_NCP, max_freq_NN))

mutlong_B <- mutate(
  group_by(filter(mutations_long, Samp_type == "B"),
           across(c(-Samp_id, -Samp_type, -Timepoint, -TreatPop,
                    -frequency, -Treat, -Pop))),
  total_freq = sum(frequency),
  max_freq = max(frequency),
  nhits = n(),
  total_freq_NN = sum(frequency[Treat == "NN"]),
  total_freq_NCP = sum(frequency[Treat == "NCP"]),
  total_freq_CNP = sum(frequency[Treat == "CNP"]),
  total_freq_CP = sum(frequency[Treat == "CP"]),
  totalfreq_diff = total_freq_CP + total_freq_NCP - total_freq_CNP - total_freq_NN,
  totalfreq_frac = (total_freq_CP + total_freq_NCP)/total_freq,
  max_freq_NN = max(c(0, frequency[Treat == "NN"])),
  max_freq_NCP = max(c(0, frequency[Treat == "NCP"])),
  max_freq_CNP = max(c(0, frequency[Treat == "CNP"])),
  max_freq_CP = max(c(0, frequency[Treat == "CP"])),
  max_freq_P = max(max_freq_NCP, max_freq_CP),
  max_freq_NP = max(max_freq_CNP, max_freq_NN),
  max_freq_C = max(max_freq_CP, max_freq_CNP),
  max_freq_NC = max(max_freq_NCP, max_freq_NN))

write.table(
  filter(mut_sum2, num_anc == 0, max_freq > 0.3,
         seq_id %in% c("NZ_KE136335.1", "NZ_KE136336.1", "NZ_KE136338.1")),
  "mutations_of_interest.tsv", sep = "\t", row.names = FALSE)

write.table(
  filter(mutlong_B, num_anc == 0, max_freq > 0.3,
         seq_id %in% c("NZ_KE136335.1", "NZ_KE136336.1", "NZ_KE136338.1")),
  "mutations_of_interest_long.tsv", sep = "\t", row.names = FALSE)

mut_of_interest <- filter(mutlong_B, max_freq_NP == 0)

ggplot(data = mut_of_interest, aes(x = Treat, y = frequency)) + 
  geom_point(aes(color = mutation_category),
             position = position_jitter(width = 0.2),
             alpha = 0.5) +
  scale_y_log10()

mut_of_interest2 <- filter(mutlong_B, max_freq_NP == 0,
                           frequency >= 0.15, Samp_type == "B")

ggplot(data = mut_of_interest2,
       aes(x = Treat, y = frequency, fill = mutation_category)) +
  geom_col(position = "stack")

mut_of_interest3 <- filter(mut_of_interest2, TreatPop != "NCP7",
                           grepl("pil|Pil", gene_product) == TRUE)

mut_of_interest3 <- mutate(
  mut_of_interest3,
  mutation_category2 = ifelse(
    mutation_category == "small_indel",
    ifelse((nchar(ref_seq) - nchar(new_seq)) %% 3 == 0, 
           "in-frame indel", 
           "frameshift indel"),
    mutation_category))
                                                       

ggplot(data = mut_of_interest3,
       aes(x = Pop, y = 1, fill = mutation_category2)) +
  geom_col(position = "stack") +
  facet_grid(~Treat)


mut_of_interest3_sum_mutcag <- 
  summarize(
    group_by(mut_of_interest3,
             Samp_id, Samp_type, Timepoint, TreatPop, Treat, Pop,
             mutation_category),
    nmuts = n())

mut_of_interest3_sum_mutcag <-
  left_join(expand.grid(
    Treat = unique(mut_of_interest3_sum_mutcag$Treat),
    Pop = unique(mut_of_interest3_sum_mutcag$Pop),
    mutation_category = unique(mut_of_interest3_sum_mutcag$mutation_category)),
    mut_of_interest3_sum_mutcag)
mut_of_interest3_sum_mutcag <- mutate(mut_of_interest3_sum_mutcag,
                                      nmuts = ifelse(is.na(nmuts), 0, nmuts))

ggplot(data = mut_of_interest3_sum_mutcag,
       aes(x = Treat, y = nmuts, color = mutation_category)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  facet_wrap(~Pop)

ggplot(data = mut_of_interest3_sum_mutcag,
       aes(x = Treat, y = nmuts, color = mutation_category)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  facet_wrap(~Pop)



pop_sum_ofinterest <- summarize(
  group_by(filter(mutlong_B, Samp_type == "B",
                  max_freq_NP == 0),
           Samp_id, Samp_type, Timepoint, TreatPop, Treat, Pop),
  nmuts = n(),
  nfixed = sum(frequency == 1))

ggplot(data = pop_sum_ofinterest, aes(x = Treat, y = nfixed)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  coord_cartesian(ylim = c(NA, 5))


ggplot(data = mutate(mut_sum2,
                     pilus_related = grepl("pil|Pil", gene_product)),
       aes(x = (total_freq_NN + total_freq_CNP)/16,
           y = (total_freq_NCP + total_freq_CP)/16)) +
  geom_point(aes(size = nhits, color = pilus_related), alpha = 0.5) +
  scale_size_continuous(name = "Number of\noccurrences", range = c(0, 4)) +
  labs(x = "Average frequency in\nno-phage populations",
       y = "Average frequency in\nphage-added populations") +
  scale_y_log10()

png("./plots/mutation_frequency.png", width = 6, height = 6, units = "in",
    res = 150)
ggplot(data = mutate(mut_sum2,
                     pilus_related = grepl("pil|Pil", gene_product)),
       aes(y = max_freq_P, x = max_freq_NP)) +
  geom_point(aes(size = total_freq/32, color = pilus_related), alpha = 0.5) +
  scale_size_continuous(name = "Average frequency", range = c(1, 4)) +
  labs(x = "Max frequency in\nno-phage populations",
       y = "Max frequency in\nphage-added populations")
dev.off()

png("./plots/mutation_frequency_simple.png", width = 6, height = 3.2, units = "in",
    res = 150)
ggplot(data = mutate(mut_sum2,
                     pilus_related = grepl("pil|Pil", gene_product)),
       aes(y = max_freq_P, x = max_freq_NP)) +
  geom_point(aes(color = pilus_related), alpha = 0.3, size = 2.5) +
  labs(x = "Max frequency in\nno-phage populations",
       y = "Max frequency in\nphage-added populations") +
  scale_color_manual(name = "Location of\nbacterial mutation",
                     breaks = c(TRUE, FALSE),
                     values = my_cols[c(6, 5)],
                     labels = c("Mutation in\nphage receptor",
                                "Mutation in\nother gene")) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme_bw() +
  guides(color = guide_legend(byrow = TRUE)) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(0.1, "in"))
dev.off()

png("./plots/mutation_frequency_NCvsC_amongP.png", width = 6, height = 5, units = "in",
    res = 150)
ggplot(data = mutate(filter(mut_sum2, max_freq_NP == 0),
                     pilus_related = grepl("pil|Pil", gene_product)),
       aes(y = max_freq_C, x = max_freq_NC)) +
  geom_point(aes(color = pilus_related), alpha = 0.5, size = 2.5) +
  labs(x = "Max frequency in\nonly-Pseudomonas populations",
       y = "Max frequency in\nall bacteria populations",
       subtitle = "Among phage-added-only mutations") +
  scale_color_manual(name = "Pilus-related\nmutation",
                     breaks = c(FALSE, TRUE),
                     values = my_cols[c(5, 6)]) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.subtitle = element_text(size = 16))
dev.off()

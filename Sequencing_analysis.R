##TODO ----
# somehow translate 100% frequency deletions into 0's for all the other mutations?
# plot mutation trajectories excluding those that appear in control treatments
# Plot only pilus & color pilus mutants in trajectories by gene or mutation
# Add other gene categories in trajectories plot
# Investigate other top freq muts and color in categories
# ML to predict treatment? phage density?
# ML on mutisol to id treatment? sensitivity to phage? twitching?
# PCA to look at mutisol
# Notes for PCA
#   perhaps cast mutation_type into binary columns so it can be used too?
#   rather than mutation_type, just divide into larger categories:
#     breaking - deletion, frameshift, nonsense
#     modifying - inframe, nonsyn
#     intergenic - intergenic
#     silent - synonymous
#   convert deletions into mutations at each gene covered by the deletion
#   by referencing the annotations file
# Calculate diversity metrics of each pop-timepoint
# Compare mutations from isolates to mutations in mutpop
#   (i.e. where are all the background/ancestral mutations in isolates?)
# Nucelotide diversity, fst, parallel selection, correlated coevol mutation, 
#   outlier analysis
# Discriminant analysis to look for diff genes



##Setup ----
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Global settings
make_slowplots <- FALSE

##Coverage analysis ----
covdat <- read.table("all_depths.txt", sep = "\t", comment.char = "",
                     stringsAsFactors = FALSE, header = TRUE)
colnames(covdat)[1:2] <- c("title", "seq_id")
covdat <- separate_wider_delim(covdat, 
                               cols = title,
                               delim = "_",
                               names = c(NA, "seqrun", "SeqCoastid", "seqtype",
                                         "Timepoint", "Pop", "Treat", "Isol"),
                               cols_remove = FALSE)
covdat <- mutate(covdat,
                 seq_id = gsub("e82cc082e40344a9", "EF", seq_id),
                 seq_id = gsub("6a49b35b19134f7e", "SA", seq_id),
                 seq_id = gsub("PP203295", "H6_1", seq_id),
                 seq_id = gsub("contig_", "", seq_id),
                 Timepoint = as.numeric(gsub("T", "", Timepoint)),
                 Pop = gsub("P", "", Pop),
                 Isol = gsub("#", "", Isol),
                 Treat_C = factor(
                   ifelse(Treat %in% c("NCP", "NN"), "NC",
                          ifelse(Treat %in% c("CNP", "CP"), "C",
                                 Treat)),
                   levels = c("Anc", "NC", "C")),
                 Treat_P = factor(
                   ifelse(Treat %in% c("NN", "CNP"), "NP",
                          ifelse(Treat %in% c("NCP", "CP"), "P",
                                 Treat)),
                   levels = c("Anc", "NP", "P")),
                 Treat = factor(Treat,
                                levels = c("Anc", "NN", "CNP", "NCP", "CP")))


ggplot(data = filter(covdat, Isol == "Pop", seqtype == "PA"),
       aes(x = seq_id, y = meandepth)) +
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0)) +
  scale_y_log10()

ggplot(data = filter(covdat, seq_id == "PA14_1", Isol == "Pop"),
       aes(x = Timepoint, y = meandepth, color = Treat)) +
  geom_point(alpha = 0.5, position = position_dodge(width = 0.3)) +
  scale_y_log10()

head(sort(filter(covdat, seq_id == "PA14_1", Isol == "Pop")$meandepth))

ggplot(data = filter(covdat, Isol != "Pop"),
       aes(x = paste0("T", Timepoint, " ", Treat, " P", Pop, " ", Isol),
           y = meandepth)) +
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0))
filter(covdat, seq_id == "PA14_1", Isol != "Pop",
       meandepth == min(filter(covdat, seq_id == "PA14_1", Isol != "Pop")$meandepth))
  
##Initial data reading & organization ----
mutations <- read.table("mutations_compare_all.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE,
                        quote = '', #don't allow quotes
                        comment.char = "") #don't allow comments since we use # for isolates

#Missing samples:
# 4777_25 - AX T0 P0 Anc #0 - expected bc reference
# 5174_134 - PA T0 P3 Anc #1 - expected bc reference

#Drop columns that contain no info
mutations <- select(mutations,
                    !clone & !codon_position_is_indeterminate &
                    !mutator_status & !population & !time & !treatment)

#Rename seq_id values
mutations <- mutate(mutations,
                    seq_id = gsub("e82cc082e40344a9", "EF", seq_id),
                    seq_id = gsub("6a49b35b19134f7e", "SA", seq_id),
                    seq_id = gsub("PP203295", "H6_1", seq_id),
                    seq_id = gsub("contig_", "", seq_id))

#Rename frequency
colnames(mutations)[colnames(mutations) == "frequency"] <- "freq"

#Number each unique mutation (shared across samples)
mutations <- 
  mutate(group_by(mutations,
                  across(-c(freq, multiple_polymorphic_SNPs_in_same_codon,
                            new_read_count, new_read_count_basis,
                            ref_read_count, ref_read_count_basis,
                            title))),
         mut_id = cur_group_id())

mutations <- separate_wider_delim(
  mutations,
  cols = title,
  delim = "_",
  names = c("seqrun", "SeqCoastid", "seqtype",
            "Timepoint", "Pop", "Treat", "Isol"),
  cols_remove = FALSE)
mutations <- mutate(mutations,
                    Timepoint = as.numeric(gsub("T", "", Timepoint)),
                    Pop = gsub("P", "", Pop),
                    Isol = gsub("#", "", Isol),
                    Treat_C = factor(
                      ifelse(Treat %in% c("NCP", "NN"), "NC",
                             ifelse(Treat %in% c("CNP", "CP"), "C",
                                    Treat)),
                      levels = c("Anc", "NC", "C")),
                    Treat_P = factor(
                      ifelse(Treat %in% c("NN", "CNP"), "NP",
                             ifelse(Treat %in% c("NCP", "CP"), "P",
                                    Treat)),
                      levels = c("Anc", "NP", "P")),
                    Treat = factor(Treat,
                                   levels = c("Anc", "NN", "CNP", "NCP", "CP")))

#Split samples by type
mutpop <- filter(mutations, Isol == "Pop")
mutisol <- filter(mutations, Isol != "Pop")

##Population sample analyses ----

#Filter out fixed-nonvariate mutations
#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
mutpop <- filter(mutpop, !mut_id %in% c("2958", "3198", "7885", "11212"))

mutpop <- mutate(
  mutpop,
  gene_category = case_when(
    grepl("pil", gene_name, ignore.case = TRUE) ~ "pilus",
    grepl("CupB3", gene_name) ~ "pilus",
    grepl("pilus", gene_name) ~ "pilus",
    grepl("Wsp", gene_name) ~ "biofilm",
    grepl("MorA", gene_name) ~ "biofilm",
    grepl("RhlR", gene_name) ~ "QS",
    grepl("LasR", gene_name) ~ "QS",
    grepl("LasI", gene_name) ~ "QS",
    grepl("quorum", gene_name) ~ "QS",
    grepl("LasA", gene_name) ~ "LasA/LasB",
    grepl("LasB", gene_name) ~ "LasA/LasB",
    grepl("Mut", gene_name) ~ "mutator",
    grepl("Tps", gene_name) ~ "TPS system",
    grepl("type VI", gene_name) ~ "T6SS",
    grepl("Type VI", gene_name) ~ "T6SS",
    grepl("Tse", gene_name) ~ "T6SS",
    grepl("flagella", gene_name) ~ "flagella",
    grepl("Thioredoxin family protein/Helix-turn-helix domain-containing protein", gene_name) ~ "mystery_1",
    TRUE ~ "other"))

#For population samples, fill in all the missing mutation-sample combinations
# with NA
mutpop <- complete(
  ungroup(mutpop),
  #Create the nesting defining each unique mutation
  nesting(aa_new_seq, aa_position, aa_ref_seq, codon_new_seq, codon_number, 
          codon_position, codon_ref_seq, gene_name, gene_position, gene_product, 
          gene_strand, genes_inactivated, genes_overlapping, genes_promoter, 
          insert_position, locus_tag, locus_tags_inactivated, 
          locus_tags_overlapping, locus_tags_promoter, mutation_category,
          new_seq, position, position_end, position_start,
          ref_seq, repeat_length, repeat_new_copies, repeat_ref_copies, 
          repeat_seq, seq_id, size, snp_type, transl_table, type, mut_id, 
          gene_category),
  #Create the nesting defining each sample
  nesting(seqrun, SeqCoastid, seqtype, Timepoint, Pop, Treat, Isol, title,
          Treat_C, Treat_P)
  #Columns left out (and therefore filled in by NA by default):
  # multiple_polymorphic_SNPs_in_same_codon,
  # new_read_count, new_read_count_basis, 
  # ref_read_count, ref_read_count_basis,
  # freq
)

##Code for checking unique samples
# temp <- summarize(group_by(mutpop,
#                         seqrun, SeqCoastid, seqtype, Timepoint, Pop, Treat, Isol, title,
#                         Treat_C, Treat_P),
#                num_muts = n())
# View(temp[order(temp$seqtype, temp$Timepoint, temp$Treat, as.numeric(temp$Pop)), ])
                        

#Duplicating ancestral H6 rows into every pop & treatment
mutpop <-
  bind_rows(
    filter(mutpop, !(seqtype == "H6" & Pop == 0)),
    crossing(
      Pop = as.character(unique(filter(mutpop, seqtype == "H6", Pop != 0)$Pop)),
      Treat = as.character(unique(filter(mutpop, seqtype == "H6", Treat != "Anc")$Treat)),
      select(filter(mutpop, seqtype == "H6", Pop == 0),
             -Pop, -Treat)))
#Duplicating ancestral PA rows into every treatment
mutpop <-
  bind_rows(
    filter(mutpop, !(seqtype == "PA" & Treat == "Anc")),
    crossing(
      Treat = as.character(unique(filter(mutpop, seqtype == "PA", Treat != "Anc")$Treat)),
      select(filter(mutpop, seqtype == "PA", Treat == "Anc"),
             -Treat)))

#Fix treatment columns post-duplication
mutpop <- mutate(mutpop,
                 freq0 = ifelse(is.na(freq), 0, freq),
                 freq01 = ifelse(freq0 == 0, 0.01, freq0),
                 Treat_C = factor(
                   ifelse(Treat %in% c("NCP", "NN"), "NC",
                          ifelse(Treat %in% c("CNP", "CP"), "C",
                                 Treat)),
                   levels = c("Anc", "NC", "C")),
                 Treat_P = factor(
                   ifelse(Treat %in% c("NN", "CNP"), "NP",
                          ifelse(Treat %in% c("NCP", "CP"), "P",
                                 Treat)),
                   levels = c("Anc", "NP", "P")),
                 Treat = factor(Treat,
                                levels = c("Anc", "NN", "CNP", "NCP", "CP")))

mutpop <- mutate(group_by(mutpop, seqtype, mut_id, Treat, Pop),
                 #everabove_025 = any(freq0 > 0.25, na.rm = TRUE),
                 everabove_010 = any(freq0 > 0.1, na.rm = TRUE),
                 everabove_005 = any(freq0 > 0.05, na.rm = TRUE),
                 everabove_010_noanc = any(freq0[Timepoint != 0] > 0.1, na.rm = TRUE),
                 everabove_005_noanc = any(freq0[Timepoint != 0] > 0.05, na.rm = TRUE),
                 fixed = all(freq0 == 1))

# test <- filter(mutpop, Pop == 5, Treat == "NCP", everabove_010 == TRUE,
#                seqtype == "PA", seq_id == "PA14_1")
# test <- mutate(
#   group_by(test, mut_id, seqtype, Treat, Pop),
#   freq_change = (freq01 - lag(freq01, order_by = Timepoint)),
#   abs_freq_rate =
#     (freq01 - lag(freq01, order_by = Timepoint))/
#     (Timepoint - lag(Timepoint, order_by = Timepoint)),
#   rel_freq_rate = abs_freq_rate/lag(freq01, order_by = Timepoint))
# 
# View(select(arrange(test, mut_id, Timepoint),
#             mut_id, Timepoint, freq, freq0, freq01, freq_change, abs_freq_rate,
#             rel_freq_rate))

mutpop <- mutate(
  group_by(mutpop, mut_id, seqtype, Treat, Pop),
  freq_change = (freq01 - lag(freq01, order_by = Timepoint)),
  abs_freq_rate =
    (freq01 - lag(freq01, order_by = Timepoint))/
    (Timepoint - lag(Timepoint, order_by = Timepoint)),
  rel_freq_rate = abs_freq_rate/lag(freq01, order_by = Timepoint))

# top_muts <- slice_max(
#   group_by(
#     filter(mutpop, 
#            everabove_010_noanc == TRUE, seqtype == "PA", seq_id == "PA14_1"),
#     Pop, Treat, Timepoint),
#   order_by = freq0, n = 3)
# View(arrange(select(top_muts, Treat, Pop, Timepoint, mut_id, freq0, 
#                     gene_name, mutation_category, type, snp_type, everything()),
#              Treat, Pop, Timepoint))

# top_muts2 <- slice_max(
#   group_by(
#     filter(mutpop,
#            !is.na(rel_freq_rate), everabove_010_noanc == TRUE, seqtype == "PA", seq_id == "PA14_1"),
#     Pop, Treat, Timepoint),
#   order_by = rel_freq_rate, n = 10)
# View(arrange(select(top_muts2, Treat, Pop, Timepoint, mut_id, rel_freq_rate, freq0,
#                     gene_name, position, mutation_category, type, snp_type, everything()),
#              Treat, Pop, Timepoint))
# sort(table(top_muts2$gene_name), decreasing = TRUE)

# temp <- filter(mutpop,
#             seqtype == "PA", seq_id == "PA14_1",
#             gene_name == "Thioredoxin family protein/Helix-turn-helix domain-containing protein")
# View(temp)
# table(temp$position)
#Thioredoxin family protein runs from 4629026	to 4629529
#Helix-turn-helix domain-containing protein runs from 4629721	to 4630503
#Gap is 4629530 - 4629720

ggplot(
  data = filter(mutpop, seqtype == "PA", gene_category == "mystery_1",
                seq_id %in% c("PA14_1")),
  aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
      color = as.factor(mut_id))) +
  geom_line() +
  facet_grid(Treat ~ Pop) +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
  theme_bw() +
  NULL

#Top hits that I haven't yet understood the function of:
#Thioredoxin family protein/Helix-turn-helix domain-containing protein 
# intergenic SNPs located between 4629614 and 4629639 are the top hit
#PA2778 family cysteine peptidase located around 6494573
#Amino-acid carrier protein AlsT/L-asparaginase I
#DUF4340 domain-containing protein located around 1387230
#Exo-alpha-sialidase

mut_sum <- summarize(
  group_by(filter(mutpop, seqtype == "PA", seq_id == "PA14_1"),
           aa_new_seq, aa_position, aa_ref_seq, codon_new_seq, codon_number, 
           codon_position, codon_ref_seq, gene_name, gene_position, gene_product, 
           gene_strand, genes_inactivated, genes_overlapping, genes_promoter, 
           insert_position, locus_tag, locus_tags_inactivated, 
           locus_tags_overlapping, locus_tags_promoter, mutation_category,
           new_seq, position, position_end, position_start,
           ref_seq, repeat_length, repeat_new_copies, repeat_ref_copies, 
           repeat_seq, seq_id, size, snp_type, transl_table, type, mut_id, 
           gene_category, Treat, Treat_C, Treat_P),
  avg_freq = mean(freq0))

#Pilus mutations
mutpop_pil <- filter(mutpop, seqtype == "PA", seq_id == "PA14_1",
                  everabove_010_noanc == TRUE, gene_category == "pilus",
                  gene_name != "type 4b pilus Flp biogenesis protein TadZ")
mutpop_pil <- mutate(mutpop_pil, 
                     lfreq01 = log10(freq01),
                     freq_cat = as.numeric(cut(freq0,
                                    breaks = c(0, 0.03, 0.1, 0.32, 1),
                                    include.lowest = TRUE)))

mutpop_pil_w <- pivot_wider(
  select(mutpop_pil,
         -seqrun, -SeqCoastid, -title,
         -multiple_polymorphic_SNPs_in_same_codon,
         -new_read_count, -new_read_count_basis, 
         -ref_read_count, -ref_read_count_basis,
         -freq, -freq0, -freq01, -freq_cat,
         -fixed, -freq_change, -abs_freq_rate, -rel_freq_rate),
  names_from = Timepoint, values_from = lfreq01, names_prefix = "T")

View(arrange(select(mutpop_pil_w, T0, T1, T2, T4, T8, Treat, Pop, mut_id, gene_name,
            everything()),
            Pop, gene_name, mut_id, Treat))

mutpop_pil_comp <-
  summarize(group_by(mutpop_pil,
                     across(-c(Treat, Treat_C, Treat_P, title,
                               seqrun, SeqCoastid, title,
                               multiple_polymorphic_SNPs_in_same_codon,
                               new_read_count, new_read_count_basis, 
                               ref_read_count, ref_read_count_basis,
                               freq, freq0, freq01, fixed, freq_change, 
                               abs_freq_rate, rel_freq_rate, 
                               lfreq01, freq_cat))),
            diff_f = ifelse(any(Treat == "CP"), freq_cat[Treat == "CP"], 1) - 
              ifelse(any(Treat == "NCP"), freq_cat[Treat == "NCP"], 1))

print(
  ggplot(
    data = mutpop_pil_comp,
    aes(x = as.numeric(Timepoint), y = diff_f, group = mut_id, 
        color = gene_name)) +
    #geom_point(aes(alpha = gene_category)) +
    geom_line(position = position_jitter(width = 0, height = 0.2)) +
    facet_grid( ~ Pop) +
    #coord_cartesian(ylim = c(0.1, 1)) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    # scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
    #                    values = my_cols[c(6, 5, 3, 8)]) +
    # scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
    #                    values = c(0.8, 0.8, 0.8, 0.15)) +
    theme_bw() +
    NULL
)



#Pilus-only PCA
pilus_pcadat <- select(
  filter(ungroup(mutpop), seqtype == "PA", seq_id == "PA14_1", 
         gene_category == "pilus"),
  Timepoint, Pop, Treat, freq0, gene_name)
pilus_pcadat <-
  filter(group_by(pilus_pcadat, gene_name), any(freq0 >= 0.1))
pilus_pcadat$Treat[pilus_pcadat$Timepoint == 0] <- "Anc"
pilus_pcadat <-
  summarize(group_by(pilus_pcadat, 
                     Timepoint, Pop, Treat, gene_name),
            pooledfreq = max(freq0),
            pooledfreq = ifelse(pooledfreq < 0.1, 0, pooledfreq))
pilus_pcadat <- pivot_wider(ungroup(pilus_pcadat),
    names_from = gene_name, values_from = pooledfreq)
pilus_pca <- prcomp(select(pilus_pcadat, !c(Timepoint, Pop, Treat)),
                    retx = TRUE)
pilus_pcadat <- cbind(as.data.frame(pilus_pca$x), pilus_pcadat)

#All-data PCA
pcadat <- select(
  filter(ungroup(mutpop), seqtype == "PA", seq_id == "PA14_1"),
  Timepoint, Pop, Treat, freq0, gene_name)
pcadat <- filter(group_by(pcadat, gene_name), any(freq0 >= 0.1))
pcadat$Treat[pcadat$Timepoint == 0] <- "Anc"
pcadat <-
  summarize(group_by(pcadat, Timepoint, Pop, Treat, gene_name),
            pooledfreq = max(freq0),
            pooledfreq = ifelse(pooledfreq < 0.1, 0, pooledfreq))
pcadat <- pivot_wider(ungroup(pcadat),
                      names_from = gene_name, values_from = pooledfreq)
mypca <- prcomp(select(pcadat, !c(Timepoint, Pop, Treat)),
                retx = TRUE)
pcadat <- cbind(as.data.frame(mypca$x), pcadat)


png("./plots/pca_pilusmuts.png", width = 10, height = 3, units = "in",
    res = 150)
ggplot(pilus_pcadat, aes(x = PC1, y = PC2, color = as.factor(Treat))) +
  geom_polygon(data = filter(pilus_pcadat, Treat %in% c("Anc", "NCP", "CP")),
               stat = "ellipse", aes(fill = as.factor(Treat)), alpha = 0.2) +
  geom_point() +
  labs(x = paste("PC1 (", 
                 round((100*((pilus_pca$sdev)**2)/
                          sum((pilus_pca$sdev)**2))[1], 1),
                 "%)", sep = ""),
       y = paste("PC2 (", 
                 round((100*((pilus_pca$sdev)**2)/
                          sum((pilus_pca$sdev)**2))[2], 1),
                 "%)", sep = "")) +
  scale_color_manual(name = "Treatments",
                     breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                     values = my_cols[c(8, 7, 1, 5, 3)]) +
  scale_fill_manual(name = "Treatments",
                    breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                    values = my_cols[c(8, 7, 1, 5, 3)]) +
  facet_grid(~ Timepoint) +
  theme_bw() +
  NULL
dev.off()

png("./plots/pca_allmuts.png", width = 10, height = 3, units = "in",
    res = 150)
ggplot(pcadat, aes(x = PC1, y = PC2, color = as.factor(Treat))) +
  geom_polygon(stat = "ellipse", aes(fill = as.factor(Treat)), alpha = 0.2) +
  geom_point() +
  labs(x = paste("PC1 (", 
                 round((100*((mypca$sdev)**2)/
                          sum((mypca$sdev)**2))[1], 1),
                 "%)", sep = ""),
       y = paste("PC2 (", 
                 round((100*((mypca$sdev)**2)/
                          sum((mypca$sdev)**2))[2], 1),
                 "%)", sep = "")) +
  scale_color_manual(name = "Treatments",
                     breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                     values = my_cols[c(8, 7, 1, 5, 3)]) +
  scale_fill_manual(name = "Treatments",
                    breaks = c("Anc", "NN", "NCP", "CNP", "CP"),
                    values = my_cols[c(8, 7, 1, 5, 3)]) +
  facet_grid(~ Timepoint) +
  theme_bw() +
  NULL
dev.off()
         
if (make_slowplots) {
  png("./plots/manhattan.png", width = 8, height = 6, units = "in",
      res = 150)
  print(
    ggplot(data = mut_sum,
         aes(x = position, y = avg_freq, color = gene_category)) +
    geom_point() +
    facet_grid(Treat ~ .)
  )
  dev.off()

  ggplot(data = filter(mutpop, freq0 > 0, seqtype == "PA", seq_id == "PA14_1"), 
         aes(x = freq0)) +
    geom_histogram() +
    facet_grid(Treat * Timepoint ~ Pop, scales = "free_y") +
    scale_x_log10()


  # png("./plots/mut_trajectories_above025.png", width = 12, height = 7, units = "in",
  #     res = 150)
  # ggplot(
  #   data = filter(mutpop, seqtype == "PA", everabove_025 == TRUE, seq_id %in% c("PA14_1")),
  #   aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
  #       color = gene_category)) +
  #   geom_point(aes(alpha = gene_category)) +
  #   geom_line(aes(alpha = gene_category)) +
  #   facet_grid(Treat ~ Pop) +
  #   #coord_cartesian(ylim = c(0.1, 1)) +
  #   scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
  #   scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
  #                      values = my_cols[c(6, 5, 3, 8)]) +
  #   scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
  #                      values = c(0.7, 0.7, 0.7, 0.3)) +
  #   theme_bw() +
  #   NULL
  # dev.off()

  png("./plots/mut_trajectories_above010.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", everabove_010 == TRUE, 
                    seq_id %in% c("PA14_1")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_category)) +
      #geom_point(aes(alpha = gene_category)) +
      geom_line(aes(alpha = gene_category)) +
      facet_grid(Treat ~ Pop) +
      #coord_cartesian(ylim = c(0.1, 1)) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = my_cols[c(6, 5, 3, 8)]) +
      scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = c(0.8, 0.8, 0.8, 0.15)) +
      theme_bw() +
      NULL
  )
  dev.off()

  png("./plots/mut_trajectories_above005.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
  ggplot(
    data = filter(mutpop, seqtype == "PA", everabove_005 == TRUE, 
                  seq_id %in% c("PA14_1")),
    aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
        color = gene_category)) +
    #geom_point(aes(alpha = gene_category)) +
    geom_line(aes(alpha = gene_category)) +
    facet_grid(Treat ~ Pop) +
    #coord_cartesian(ylim = c(0.1, 1)) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                       values = my_cols[c(6, 5, 3, 8)]) +
    scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                       values = c(0.8, 0.8, 0.8, 0.15)) +
    theme_bw() +
    NULL
  )
  dev.off()

  everabove_010_noanc

  png("./plots/mut_trajectories_above010_noanc.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
  ggplot(
    data = filter(mutpop, seqtype == "PA", everabove_010_noanc == TRUE, 
                  seq_id %in% c("PA14_1")),
    aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
        color = gene_category)) +
    #geom_point(aes(alpha = gene_category)) +
    geom_line(lwd = 1.1, aes(alpha = gene_category)) +
    facet_grid(Treat ~ Pop) +
    #coord_cartesian(ylim = c(0.1, 1)) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    scale_color_manual(breaks = c("biofilm", "flagella", "LasA/LasB", "QS",
                                  "T6SS", "pilus", "TPS system", "other"),
                       values = my_cols) +
    scale_alpha_manual(breaks = c("biofilm", "flagella", "LasA/LasB", "QS",
                                  "T6SS", "pilus", "TPS system", "other"),
                       values = c(rep(0.6, 7), 0.1)) +
    theme_bw() +
    NULL
  )
  dev.off()

  png("./plots/mut_trajectories_above005_noanc.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", everabove_005_noanc == TRUE, 
                    seq_id %in% c("PA14_1")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_category)) +
      #geom_point(aes(alpha = gene_category)) +
      geom_line(aes(alpha = gene_category)) +
      facet_grid(Treat ~ Pop) +
      #coord_cartesian(ylim = c(0.1, 1)) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = my_cols[c(6, 5, 3, 8)]) +
      scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = c(0.8, 0.8, 0.8, 0.15)) +
      theme_bw() +
      NULL
  )
  dev.off()

  png("./plots/mut_trajectories_all_noother.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, gene_category != "other",
                    seqtype == "PA", seq_id %in% c("PA14_1")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_category)) +
      #geom_point(aes(alpha = gene_category)) +
      geom_line(aes(alpha = gene_category)) +
      facet_grid(Treat ~ Pop) +
      #coord_cartesian(ylim = c(0.1, 1)) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = my_cols[c(6, 5, 3, 8)]) +
      scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = c(0.8, 0.8, 0.8, 0.15)) +
      theme_bw() +
      NULL
  )
  dev.off()

  png("./plots/mut_trajectories_all.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", seq_id %in% c("PA14_1")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_category)) +
      #geom_point(aes(alpha = gene_category)) +
      geom_line(aes(alpha = gene_category)) +
      facet_grid(Treat ~ Pop) +
      #coord_cartesian(ylim = c(0.1, 1)) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = my_cols[c(6, 5, 3, 8)]) +
      scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                         values = c(0.8, 0.8, 0.8, 0.15)) +
      theme_bw() +
      NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_pilus.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", seq_id == "PA14_1",
                    everabove_010_noanc == TRUE, gene_category == "pilus",
                    gene_name != "type 4b pilus Flp biogenesis protein TadZ"),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_name)) +
      #geom_point(aes(alpha = gene_category)) +
      geom_line() +
      facet_grid(Treat ~ Pop) +
      #coord_cartesian(ylim = c(0.1, 1)) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      # scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
      #                    values = my_cols[c(6, 5, 3, 8)]) +
      # scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
      #                    values = c(0.8, 0.8, 0.8, 0.15)) +
      theme_bw() +
      NULL
  )
  dev.off()

  png("./plots/mut_changes_above010_noanc.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
  ggplot(
    data = filter(mutpop, seqtype == "PA", everabove_010_noanc == TRUE, 
                  seq_id %in% c("PA14_1")),
    aes(x = as.numeric(Timepoint), y = rel_freq_rate, group = mut_id, 
        color = gene_category)) +
    geom_point(aes(alpha = gene_category)) +
    geom_line(aes(alpha = gene_category)) +
    facet_grid(Treat ~ Pop) +
    #coord_cartesian(ylim = c(0.1, 1)) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    scale_color_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                       values = my_cols[c(6, 5, 3, 8)]) +
    scale_alpha_manual(breaks = c("pilus", "biofilm", "mutator", "other"),
                       values = c(0.8, 0.8, 0.8, 0.15)) +
    theme_bw() +
    NULL
  )
  dev.off()

  png("./plots/mut_freqs.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(data = filter(mutpop, #gene_category != "other", 
                         seq_id %in% c("PA14_1"), seqtype == "PA"),
           aes(x = as.factor(Timepoint), y = freq0, fill = gene_category)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(Treat ~ Pop) +
      # scale_fill_manual(breaks = c("biofilm", "flagella", "LasA/LasB", "QS",
      #                               "T6SS", "pilus", "TPS system", "other"),
      #                    values = my_cols) +
      NULL
  )
  dev.off()

  ggplot(
    data = filter(mutpop, seqtype == "PA", gene_category == "pilus",
                  seq_id %in% c("PA14_1")),
    aes(x = as.numeric(Timepoint), y = freq0, group = gene_name, 
        color = as.factor(gene_name))) +
    geom_line() +
    facet_grid(Treat ~ Pop) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    theme_bw() +
    NULL
}


# temp <- filter(mutpop, fixed == TRUE)
# View(temp[order(temp$mut_id, temp$Timepoint, temp$Treat, as.numeric(temp$Pop)),
#           c("freq", "mut_id", "Timepoint", "Treat", "Pop")])
# 
# temp <- filter(mutpop, seqtype == "PA", Isol == "Pop", Treat == "NCP",
#                Pop == 3, everabove_025 == TRUE, seq_id == "PA14_1")
# View(temp[order(temp$mut_id, temp$Timepoint), 
#           c("freq", "mut_id", "Timepoint")])

##Isolate analyses ----

##Note: five isolates have only the four shared fixed mutations 
##      (i.e. they're ancestral)
# 5174_58 - PA T1 P4 CP #3 - susceptible
# 5174_64 - PA T4 P4 NCP #3 - susceptible
# 5174_69 - PA T1 P5 NCP #1 - susceptible
# 5174_70 - PA T1 P5 NCP #2 - resistant
# 5174_129 - PA T8 P10 CP #1 - resistant

#Filter out sample that had 0 coverage (so mutation is a full-genome deletion)
mutisol <- filter(mutisol, title != "5174_67_PA_T8_P4_NCP_#5")

mutisol <- left_join(
  mutisol,
  select(
    mutate(read.csv("aggregate_isolate_data.csv", stringsAsFactors = FALSE,
                    colClasses = "character"),
           Treat = Bact_treat, Treat_C = Bact_treat_C, Treat_P = Bact_treat_P,
           Timepoint = Bact_timepoint, Pop = Bact_pop, Isol = Bact_isol),
    Treat, Treat_C, Treat_P, Timepoint, Pop, Isol, 
    clearing_avg, EOP_avg, Plaque_morphology, resis, resis_cat,
    cm_day, cm_day_norm, twitch_cat,
    wavelength_405, wavelength_695,
    max_percap, max_dens, lag, auc, max_percap_time, max_percap_dens, min_dens,
    max_percap_norm, max_dens_norm, auc_norm, lag_norm,
    biofilm_abs, biofilm_norm, biofilm_norm_log10, 
    biofilm_norm_log10_norm, biofilm_block_name))

mutisol <- mutate(
  mutisol,
  gene_category = case_when(
    grepl("pil", gene_name, ignore.case = TRUE) ~ "pilus",
    grepl("CupB3", gene_name) ~ "pilus",
    grepl("Wsp", gene_name) ~ "biofilm",
    grepl("MorA", gene_name) ~ "biofilm",
    grepl("LasR", gene_name) ~ "biofilm",
    grepl("Mut", gene_name) ~ "mutator",
    TRUE ~ "other"),
  mutation_category2 = case_when(
    mutation_category == "small_indel" & grepl("intergenic", gene_position) ~
      "small_indel_intergenic",
    mutation_category == "small_indel" & grepl("coding", gene_position) &
      (nchar(ref_seq) - nchar(new_seq)) %% 3 == 0 ~ "small_indel_inframe",
    mutation_category == "small_indel" & grepl("coding", gene_position) &
      (nchar(ref_seq) - nchar(new_seq)) %% 3 != 0 ~ "small_indel_frameshift",
    TRUE ~ mutation_category),
  mutation_category3 = case_when(
    mutation_category2 %in% 
      c("large_deletion", "small_indel_frameshift", "snp_nonsense") ~ "loss of function",
    mutation_category2 %in% c("small_indel_inframe") ~ mutation_category2,
    mutation_category2 %in% c("small_indel_intergenic", "snp_intergenic") ~ "intergenic",
    mutation_category2 %in% c("snp_nonsense", "snp_nonsynonymous") ~ "other change",
    mutation_category2 %in% c("snp_synonymous") ~ "silent",
    TRUE ~ "NA"),
  gene_short = coalesce(
    # 1) Capitalized short symbol at end, e.g. "… PilT" or "… PilY1"
    str_extract(gene_name, "([A-Z][A-Za-z0-9]+)$"),
    # 2) Gene token immediately following "protein"
    str_extract(gene_name, "(?i)\\bPil[A-Za-z0-9]*(?=$|/)"),
    # 3) As a final fallback use the provided gene_name column
    gene_name),
  concise_mut = case_when(
    mut_id %in% c("2958", "3198", "7885", "11212") ~ "", #(ancestral mutations)
    grepl("snp", mutation_category) ~
      ifelse(grepl("intergenic", gene_position),
             paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
             paste0(gene_short, " (", aa_ref_seq, "->", aa_new_seq, ")")),
    mutation_category == "small_indel" & type == "DEL" ~
      ifelse(grepl("intergenic", gene_position),
             paste0(gene_short, " (intergenic ", size, "bp DEL)"),
             paste0(gene_short, " (", size, "bp DEL)")),
    mutation_category == "small_indel" & type == "INS" ~
      ifelse(grepl("intergenic", gene_position),
             paste0(gene_short, " (intergenic ", length(new_seq), "bp INS)"),
             paste0(gene_short, " (", length(new_seq), "bp INS)")),
    mutation_category == "small_indel" & type == "SUB" ~
      ifelse(grepl("intergenic", gene_position),
             paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
             paste0(gene_short, " (nt ", ref_seq, "->", new_seq, ")")),
    mutation_category == "large_deletion" ~
      paste0(gene_short, " (", size, "bp DEL)"),
    TRUE ~ "nomatch")
  )

nice_isoltable <- 
  pivot_wider(
    select(
      mutate(ungroup(mutisol),
             mut_combin = ifelse(
               mut_id %in% c("2958", "3198", "7885", "11212"),
               NA,
               paste0(mutation_category2, " [", concise_mut, "]"))),
      Timepoint, Treat, Pop, Isol, Treat_C, Treat_P,
      clearing_avg, EOP_avg, Plaque_morphology, resis_cat, 
      twitch_cat,
      gene_category, mut_combin),
    names_from = gene_category,
    names_prefix = "mut_in_",
    values_from = mut_combin,
    values_fn = function(x) {
      if(!all(is.na(x))) {paste0(x[!is.na(x)], collapse = "; ")} else {(NA)}})

write.csv(x = nice_isoltable[order(nice_isoltable$resis_cat,
                                   nice_isoltable$Timepoint,
                                   nice_isoltable$Treat_P, 
                                   nice_isoltable$Treat_C,
                                   nice_isoltable$Pop,
                                   nice_isoltable$Isol), ],
          file = "nice_isolmuts.csv", row.names = FALSE)

#Filter out fixed-nonvariate mutations
#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
mutisol <- filter(mutisol, !mut_id %in% c("2958", "3198", "7885", "11212"))

#Note: we sequenced more NCP isols than CP
                    
png("./plots/isol_mut_types.png", width = 8, height = 8, units = "in",
    res = 150)
ggplot(data = mutisol,
       aes(x = Treat, fill = mutation_category2)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(~ gene_category)
dev.off()

ggplot(data = filter(mutisol),
       aes(x = resis_cat, fill = mutation_category2)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(Treat ~ gene_category)

ggplot(data = filter(mutisol),
       aes(x = resis_cat, fill = mutation_category2)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(Treat ~ gene_category)

ggplot(data = filter(mutisol),
       aes(x = resis_cat, fill = mutation_category3)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(gene_category ~ Treat)

ggplot(data = filter(mutisol),
       aes(x = resis_cat, fill = mutation_category3)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(gene_category ~ Treat)

ggplot(data = filter(mutisol, gene_category == "pilus"),
       aes(x = twitch_cat, fill = mutation_category3)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(gene_category ~ Treat)

ggplot(data = filter(mutisol, gene_category == "pilus"),
       aes(x = twitch_cat, fill = mutation_category3)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(gene_category ~ Treat)


mutisol_sum <- summarize(
  group_by(mutisol, Timepoint, Treat, Pop, Isol, Treat_C, Treat_P),
  nmuts = n(),
  npil = sum(gene_category == "pilus"))

ggplot(data = mutisol_sum,
       aes(x = Timepoint, y = nmuts)) +
  geom_point(alpha = 0.5,
             position = position_jitter(width = 0.2, height = 0)) +
  scale_y_log10() +
  facet_grid(~Treat)

ggplot(data = mutisol_sum,
       aes(x = Timepoint, y = npil)) +
  geom_point(alpha = 0.5,
             position = position_jitter(width = 0.2, height = 0)) +
  facet_grid(~Treat)

View(filter(mutisol[order(mutisol$gene_name, mutisol$mut_id, 
                          mutisol$Timepoint, mutisol$Treat, 
                          as.numeric(mutisol$Pop), mutisol$Isol), ],
            !mut_id %in% c("2958", "3198", "7885", "11212")))

View(filter(mutisol, 
            !mut_id %in% c("2958", "3198", "7885", "11212"),
            new_seq == ref_seq))


##Old code ----
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
  filter(mutlong_B,
         seq_id %in% c("NZ_KE136335.1", "NZ_KE136336.1", "NZ_KE136338.1")),
  "all_PA_mutations_long.tsv", sep = "\t", row.names = FALSE)

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

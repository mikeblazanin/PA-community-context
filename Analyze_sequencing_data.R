##TODO ----
#Top priority:
# nothing left!

#Maybe:
# create gene_category2 which ID's type IV pilus genes more-specifically
# adjust coverage analysis for repetitive regions
# Get total # reads (not just aligned reads) & use as denom for coverage analy
# why are cov-based densities of other bact higher in NC treats?
# Look for dift fates of phages depending on initial inocula 3-10 & pilus muts
# Is there anything in common among the populations where phage went ~extinct?
#  Some way those 5 pops cluster?
# Diversity metrics
# Make mutation trajectory plots for non-PA data

#Other ideas:
# check Ian's ideas about KEGG/BRITE
# compare pilus gene size to # of muts
# Plot mutations in N samples vs freq plot
#  (n reflects mut accessibility, freq reflects selection)
# PCA mutisol
# PCA but casting paste(gene_name, mutation_cat3 a la mutisol)
# Diversity metrics of each pop-timepoint
# fst
# parallel selection - Toby
# outlier analysis - Chelsea/Toby
# ML on mutpop to predict treatment or phage density
# ML on mutisol to predict treatment or resistance or twitching
# plot mutation trajectories excluding those that appear in control treatments
# Discriminant analysis to look for diff genes
# ID mutations that are moving together
# correlated coevol mutation
# Use CVTK or lofreq to call mutations
#  Brennan et al. PNAS 2022
#  Buffalo & Coop. PNAS 2020

##Setup ----
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(lme4)
library(ggh4x)
library(emmeans)
library(pbkrtest)
library(lmerTest)
library(stringdist)
library(BiocManager)
#BiocManager::install("clusterProfiler")
#BiocManager::install("GO.db")
library(clusterProfiler)
library(ggrepel)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Extended version of Carto "Safe" palette
my_cols2 <- c("#888888", "#117733", "#332288", "#661100","#882255",
              "#88CCEE", "#CC6677", "#DDCC77",  
            "#AA4499", "#44AA99","#999933",  "#6699CC",
             "#E69F00")

#Global settings
make_slowplots <- TRUE

##Coverage analysis ----
covdat <- read.csv("./clean_data/covdat.csv", stringsAsFactors = FALSE)
covdat <- 
  as_tibble(
    mutate(covdat,
           seqrun = as.character(seqrun),
           SeqCoastid = as.character(SeqCoastid),
           Pop = as.character(Pop),
           Treat = factor(Treat, levels = c("Anc", "NN", "CNP", "NCP", "CP")),
           Treat_C = factor(Treat_C, levels = c("Anc", "NC", "C")),
           Treat_P = factor(Treat_P, levels = c("Anc", "NP", "P"))))

covbycontig <- 
  summarize(
    group_by(
      filter(covdat, Isol == "Pop", seqtype == "PA"),
      Strain, seq_id, Timepoint, Pop, Treat, Treat_C, Treat_P, endpos,
      included_in_treat),
    totbases = meandepth * endpos)
covbycontig <- mutate(group_by(covbycontig, Timepoint, Pop, 
                               Treat, Treat_C, Treat_P, included_in_treat),
                      rel_freq = totbases/sum(totbases))
covbysp <- 
  summarize(group_by(covbycontig, Strain, Timepoint, Pop, Treat, Treat_C, Treat_P,
                     included_in_treat),
            totbases = weighted.mean(x = totbases, w = endpos))
covbysp <-
  mutate(group_by(covbysp, Timepoint, Pop, Treat, Treat_C, Treat_P,
                  included_in_treat),
         rel_freq = totbases/sum(totbases))


#Incorporate cfu/pfu data
cfudensity <- read.csv("clean_data/density.csv", stringsAsFactors = FALSE)

#Change factor orders
cfudensity <- mutate(cfudensity,
                     Treat_C = factor(treat_C, levels = c("NC", "C")),
                     Treat_P = factor(treat_P, levels = c("NP", "P")),
                     Treat = paste0(Treat_C, Treat_P),
                     Treat = ifelse(Treat == "NCNP", "NN", Treat),
                     Pop = as.character(Pop))
cfudensity <- select(cfudensity,
                     Timepoint, Pop, Strain,
                     Treat, Treat_C, Treat_P,
                     cfu_ml)
#We're missing T2 CP3 CP4 CP5 AX
cfudensity <- mutate(group_by(cfudensity,
                              Timepoint, Pop, Treat, Treat_C, Treat_P),
                     rel_freq_cfu = cfu_ml/sum(cfu_ml, na.rm = TRUE))
covbysp <- left_join(covbysp, cfudensity)
covbysp <- mutate(covbysp,
                  rel_freq_cfu = ifelse(is.na(rel_freq_cfu), 0, rel_freq_cfu))

if(make_slowplots) {
  png("./plots/coverage_by_contig.png", width = 12, height = 4, units = "in",
      res = 150)
  print(
    ggplot(data = filter(covdat, Isol == "Pop", seqtype == "PA"),
         aes(x = seq_id, y = log10(meandepth))) +
    geom_point(alpha = 0.2) +
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_x_discrete(limits = c("PA14_1", "H6_1",
                              paste0("SA_", 1:2),
                              paste0("EF_", 1:4),
                              paste0("AX1_", 1:71))) +
    labs(y = "Mean coverage depth", x = "Contig from reference genomes",
         subtitle = "Metagenomic sequencing") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4))
  )
  dev.off()
  
  png("./plots/coverage_by_contig_phageseq.png", width = 12, height = 4, units = "in",
      res = 150)
  print(
    ggplot(data = filter(covdat, Isol == "Pop", seqtype == "H6"),
         aes(x = seq_id, y = log10(meandepth))) +
    geom_point(alpha = 0.2) +
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_x_discrete(limits = c("PA14_1", "H6_1",
                                paste0("SA_", 1:2),
                                paste0("EF_", 1:4),
                                paste0("AX1_", 1:71))) +
    labs(y = "Mean coverage depth", x = "Contig from reference genomes",
         subtitle = "Viral filtrate sequencing") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4))
  )
  dev.off()
  
  #Coverage PA by timepoint and treat
  ggplot(data = filter(covdat, seqtype == "PA", seq_id == "PA14_1", Isol == "Pop"),
         aes(x = Timepoint, y = meandepth, color = Treat)) +
    geom_point(alpha = 0.5, position = position_dodge(width = 0.5),
               size = 2) +
    scale_y_log10(limits = c(100, NA)) +
    scale_color_manual(name = "Treatments",
                       breaks = c("Anc", "NN", "CNP", "NCP", "CP"),
                       labels = c("Ancestor", "No-phage, Pseudomonas-only",
                                  "No-phage, all bacteria",
                                  "Phage-added, Pseudomonas-only",
                                  "Phage-added, all bacteria"),
                       values = my_cols[c(8, 7, 5, 1, 3)]) +
    labs(x = "Timepoint (days)", y = "Mean coverage depth aligned to\nPseudomonas genome",
         subtitle = "Metagenomic sequencing") +
    theme_bw()

  #Coverage by contig by treatment
  ggplot(data = mutate(filter(covdat, Isol == "Pop", seqtype == "PA"),
                       seq_id = factor(seq_id,
                                       levels = c("PA14_1", "H6_1",
                                                  paste0("SA_", 1:2),
                                                  paste0("EF_", 1:4),
                                                  paste0("AX1_", 1:71)))),
         aes(x = Treat, y = log10(meandepth))) +
    geom_point(alpha = 0.2) +
    facet_wrap(~seq_id, nrow = 8) +
    scale_y_continuous(labels = math_format(10^.x)) +
    #scale_color_discrete(palette = "Set2") +
    labs(y = "Mean coverage depth", x = "Treatment",
         subtitle = "Metagenomic sequencing") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4))

  png("./plots/density_fromseq.png", width = 10, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(covbysp, Treat != "Anc"),
      aes(x = Timepoint, y = log10(rel_freq), color = Strain)) +
      geom_point(position = position_jitter(width = 0.05, height = 0),
                 alpha = 0.3) + 
      geom_line(aes(group = paste(Pop, Strain)), alpha = 0.5) +
      facet_grid(Treat_C ~ Treat_P,
                 labeller = labeller(Treat_P = c("P" = "Phage Added",
                                                 "NP" = "No Phage"),
                                     Treat_C = c("C" = "All Bacteria",
                                                 "NC" = "Only Pseudomonas"))) +
      scale_x_continuous(breaks = c(0, 3, 6, 9), minor_breaks = 0:9, 
                         limits = c(0, NA)) +
      scale_y_continuous(labels = math_format(10^.x)) +
      labs(y = "Frequency", x = "Timepoint (days)") +
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
            legend.background = element_rect(color = "black")) +
      NULL)
  dev.off()
  
  p1 <- ggplot(data = arrange(
    filter(covbysp, Timepoint != 0,
           !(Timepoint == 2 & Treat == "CP" & Pop %in% c(3, 4, 5) & Strain == "AX"),
           included_in_treat == TRUE),
    Timepoint),
    aes(x = log10(rel_freq), y = log10(rel_freq_cfu), 
        color = Strain, fill = Strain, shape = included_in_treat)) +
    geom_point(alpha = 0.7, show.legend = TRUE, size = 3) +
    scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                       labels = c("P. aeruginosa", "S. aureus",
                                  "E. faecalis", "A. xylosoxidans",
                                  "Phage H6"),
                       values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_fill_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                      labels = c("P. aeruginosa", "S. aureus",
                                 "E. faecalis", "A. xylosoxidans",
                                 "Phage H6"),
                      values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_shape_manual(name = "Included in treatment",
                       limits = c(TRUE, FALSE),
                       values = c(18, 16),
                       labels = c("Yes", "No")) +
    scale_x_continuous(labels = math_format(10^.x)) +
    scale_y_continuous(labels = math_format(10^.x)) +
    geom_abline(slope = 1, intercept = 0) +
    # guides(color = guide_legend(override.aes = list(size = 2)),
    #        shape = guide_legend(override.aes = list(size = 2))) +
    labs(x = "", y = "Frequency (CFU or PFU)") +
    theme_bw()
  p2 <- ggplot(data = arrange(
    filter(covbysp, Timepoint != 0,
           included_in_treat == FALSE,
           rel_freq != 0),
    Timepoint),
    aes(y = log10(rel_freq), x = 0, 
        color = Strain, fill = Strain, shape = included_in_treat)) +
    geom_point(position = position_dodge(width = 1), alpha = 0.7, size = 3) +
    coord_flip() +
    scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                       labels = c("P. aeruginosa", "S. aureus",
                                  "E. faecalis", "A. xylosoxidans",
                                  "Phage H6"),
                       values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_fill_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                      labels = c("P. aeruginosa", "S. aureus",
                                 "E. faecalis", "A. xylosoxidans",
                                 "Phage H6"),
                      values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_shape_manual(name = "Included in treatment",
                       limits = c(TRUE, FALSE),
                       values = c(18, 16),
                       labels = c("Yes", "No")) +
    scale_y_continuous(labels = math_format(10^.x)) +
    scale_x_continuous(breaks = 0, limits = c(-1, 1)) +
    labs(y = "Frequency (metagenomics)", x = "") +
    theme_bw() +
    guides(color = "none", fill = "none", shape = "none")
  p3 <- ggplot(data = arrange(
    filter(covbysp, Timepoint != 0,
           included_in_treat == FALSE,
           rel_freq == 0),
    Timepoint),
    aes(x = 0, y = 0, color = Strain, fill = Strain, shape = included_in_treat)) +
    coord_flip() +
    geom_dotplot(aes(group = Strain), position = position_dodge(1),
                 alpha = 0.7, stackdir = "center", dotsize = 5) +
    # geom_point(position = 
    #              position_jitterdodge(dodge.width = 0.5, jitter.height = 0.5,
    #                                   jitter.width = 0), 
    #            alpha = 0.7) +
    scale_color_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                       labels = c("P. aeruginosa", "S. aureus",
                                  "E. faecalis", "A. xylosoxidans",
                                  "Phage H6"),
                       values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_fill_manual(limits = c("PA", "SA", "EF", "AX", "H6"),
                      labels = c("P. aeruginosa", "S. aureus",
                                 "E. faecalis", "A. xylosoxidans",
                                 "Phage H6"),
                      values = my_cols[c(3, 1, 2, 7, 8)]) +
    scale_shape_manual(name = "Included in treatment",
                       limits = c(TRUE, FALSE),
                       values = c(18, 16),
                       labels = c("Yes", "No")) +
    scale_x_continuous(breaks = 0, limits = c(-1, 1)) +
    scale_y_continuous(breaks = 0, limits = c(-1, 1)) +
    labs(y = "", x = "") +
    theme_bw() +
    guides(color = "none", fill = "none", shape = "none")
  png("./plots/density_comparison.png", width = 10.5, height = 7, units = "in",
      res = 150)
  print(
    cowplot::plot_grid(NULL, 
                     p1 + guides(color = "none", fill = "none", shape = "none"), 
                     cowplot::get_legend(p1),
                     p3, 
                     p2, 
                     NULL,
                     nrow = 2,
                     rel_heights = c(1, 0.2),
                     rel_widths = c(0.2, 1, 0.25),
                     align = "hv", axis = "tblr")
  )
  dev.off()
}

##Read in clean data ----
source("Define_colgroups_sequencing.R")

term2gene <- read.csv("./clean_data/term2gene.csv", stringsAsFactors = FALSE)
term2name <- read.csv("./clean_data/term2name.csv", stringsAsFactors = FALSE,
                       row.names = 1)

mutations <- read.csv("./clean_data/mutations_PAseq_PA14seqid.csv", 
                      stringsAsFactors = FALSE, colClasses = "character")
mutationsd <- read.csv("./clean_data/mutationsd_PAseq_PA14seqid.csv", 
                       stringsAsFactors = FALSE)

mutations <- filter(mutations, seqtype == "PA", seq_id == "PA14_1")
mutationsd <- filter(mutationsd, seqtype == "PA", seq_id == "PA14_1")

mutations <-
  as_tibble(
    mutate(mutations,
           Treat = factor(Treat, levels = c("Anc", "NN", "CNP", "NCP", "CP")),
           Treat_C = factor(Treat_C, levels = c("Anc", "NC", "C")),
           Treat_P = as.factor(Treat_P),
           mut_id = as.integer(mut_id),
           Timepoint = as.numeric(Timepoint),
           rpt_bit_score = as.numeric(rpt_bit_score)))
mutations <- dplyr::rename(mutations, `e-value_score` = `e.value_score`)

mutationsd <-
  as_tibble(
    mutate(mutationsd,
           Treat = factor(Treat, levels = c("Anc", "NN", "CNP", "NCP", "CP")),
           Pop = as.character(Pop),
           seqrun = as.character(seqrun),
           SeqCoastid = as.character(SeqCoastid),
           Treat_C = factor(Treat_C, levels = c("Anc", "NC", "C")),
           Treat_P = as.factor(Treat_P)))
mutationsd <- dplyr::rename(mutationsd, `e-value_score` = `e.value_score`)


#Split samples by type ----
mutpop <- filter(mutations, Isol == "Pop")
mutisol <- filter(mutations, Isol != "Pop")
mutpopd <- filter(mutationsd, Isol == "Pop")
mutisold <- filter(mutationsd, Isol != "Pop")

##Population sample analyses ----

### Setup ----
#For population samples, fill in all the missing mutation-sample combinations
# with NA
mutpop <- complete(
  ungroup(mutpop),
  nesting(!!!syms(intersect(mutcols, colnames(mutpop)))),
  nesting(!!!syms(intersect(sampcols, colnames(mutpop))))
  #Columns left out (and therefore filled in by NA by default):
  # mut_and_samp_cols
)
mutpopd <- complete(
  ungroup(mutpopd),
  nesting(!!!syms(intersect(mutcols, colnames(mutpopd)))),
  nesting(!!!syms(intersect(sampcols, colnames(mutpopd))))
  #Columns left out (and therefore filled in by NA by default):
  # mut_and_samp_cols
)

#Duplicating ancestral PA rows into every treatment
mutpop <-
  bind_rows(
    filter(mutpop, !Treat == "Anc"),
    crossing(
      Treat = as.character(unique(filter(mutpop, Treat != "Anc")$Treat)),
      select(filter(mutpop, Treat == "Anc"),
             -Treat)))
#Drop ancestral from mutpopd
mutpopd <- filter(mutpopd, Treat != "Anc")

#Fix treatment columns post-duplication
mutpop <- mutate(mutpop,
                 freq = as.numeric(freq),
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
mutpopd <- mutate(mutpopd,
                  freq = as.numeric(freq),
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

mutpop <- mutate(group_by(mutpop, 
                          across(any_of(c(mutcols, "seqtype", "Treat", "Pop")))),
                 everabove_010 = any(freq0 > 0.1, na.rm = TRUE),
                 fixed = all(freq0 == 1))
mutpopd <- mutate(group_by(mutpopd, 
                          across(any_of(c(mutcols, "seqtype", "Treat", "Pop")))),
                 everabove_010 = any(freq0 > 0.1, na.rm = TRUE),
                 fixed = all(freq0 == 1))

### Filter out fixed-nonvariate mutations ----
mutpop <- mutate(group_by(mutpop,
                          across(any_of(mutcols)),
                          Pop),
                 freq_popavg = mean(freq0),
                 freq_popsd = sd(freq0))
mutpopd <- mutate(group_by(mutpopd,
                          across(any_of(mutcols)),
                          Pop),
                 freq_popavg = mean(freq0),
                 freq_popsd = sd(freq0))

#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
ids_out <- unique(filter(mutpop, freq_popavg == 1, freq_popsd == 0)$mut_id)
ids_outd <- unique(filter(mutpopd, freq_popavg == 1, freq_popsd == 0)$mut_id)
ids_out
ids_outd
mutpop <- filter(mutpop, !mut_id %in% ids_out)
mutpopd <- filter(mutpopd, !mut_id %in% ids_outd)

### Analysis of mutation trajectories ----
if(make_slowplots) {
  png("./plots/mut_freq_distr.png", width = 10, height = 8,
      units = "in", res = 150)
  print(
    ggplot(data = 
             mutate(
               rbind(
                 filter(mutpop, 
                        freq0 > 0, Timepoint != 0),
                 mutate(filter(mutpop, 
                               freq0 > 0, Timepoint == 0, Treat == "NN"),
                        Treat = "Anc")),
               Pop = as.numeric(Pop),
               Treat = factor(Treat, levels = c("Anc", "NN", "CNP", "NCP", "CP"))),
           aes(x = freq0, fill = Treat)) +
      geom_histogram() +
      facet_nested(Treat * Timepoint ~ Pop, scales = "free_y",
                   labeller = labeller(Timepoint = function(x) {paste0("T", x)},
                                       Pop = function(x) {paste0("Rep", x)})) +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75)) +
      scale_fill_manual(name = "Treatments",
                        breaks = c("Anc", "NN", "CNP", "NCP", "CP"),
                        labels = c("Ancestor", "No-phage, Pseudomonas-only",
                                   "No-phage, all bacteria",
                                   "Phage-added, Pseudomonas-only",
                                   "Phage-added, all bacteria"),
                        values = my_cols[c(8, 7, 5, 1, 3)]) +
      labs(x = "Frequency", y = "Mutation count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      NULL)
  dev.off()
  
  png("./plots/mut_trajectories_above010.png", 
      width = 12, height = 7, units = "in",
      res = 150)
  f5a <- ggplot(
      data = mutate(filter(mutpop, everabove_010 == TRUE),
                    Pop = as.numeric(Pop)),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_catpilnot)) +
      geom_line(lwd = 1.2, alpha = 0.6) +
      facet_nested(Treat_P * Treat_C ~ Pop,
                 labeller = labeller(
                   Pop = function(x) {paste("Rep", LETTERS[as.numeric(x)-2])},
                   Treat_C = c("NC" = "Pseudo.-only",
                             "C" = "All bacteria"),
                   Treat_P = c("NP" = "No-phage",
                               "P" = "Phage-added"))) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
      scale_color_manual(breaks = c("pilus", "other"),
                         values = my_cols[c(6, 8)],
                         name = "Gene function") +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 18),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 16),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15)) +
      NULL
  print(f5a)
  dev.off()
  
  png("./plots/mut_trajectories.png", 
      width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = mutate(mutpop, Pop = as.numeric(Pop)),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_catpilnot)) +
      geom_line(lwd = 1.2, alpha = 0.5) +
      facet_nested(Treat_P * Treat_C ~ Pop,
                   labeller = labeller(
                     Pop = function(x) {paste("Rep", LETTERS[as.numeric(x)-2])},
                     Treat_C = c("NC" = "Pseudo.-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "other"),
                         values = my_cols[c(6, 8)],
                         name = "Gene function") +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 18),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 16),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15)) +
      NULL
  )
  dev.off()
  
  pildat <- mutate(
    group_by(filter(mutpop,
                    everabove_010 == TRUE, gene_category == "pilus"),
             gene_name, gene_short),
    Pop = as.numeric(Pop))
  pildat <- mutate(
    ungroup(pildat),
    gene_name = 
      factor(gene_name,
             levels = 
               arrange(
                 summarize(group_by(pildat, gene_name), freqsum = sum(freq0)), 
                 desc(freqsum))$gene_name),
    gene_short = 
      factor(gene_short,
             levels = 
               arrange(
                 summarize(group_by(pildat, gene_short), freqsum = sum(freq0)), 
                 desc(freqsum))$gene_short)
  )
  pildat <- mutate(
    group_by(pildat, gene_short),
    gene_num = cur_group_id(),
    gene_name_num = paste0(cur_group_id(), ": ", gene_short))
  pildat <- mutate(
    ungroup(pildat),
    gene_name_num = factor(
      gene_name_num,
      levels = arrange(summarize(group_by(pildat, gene_num, gene_name_num)),
                       gene_num)$gene_name_num))
  pildat_sum <- summarize(
    group_by(pildat, gene_num, gene_name_num, Pop, Treat, Treat_C, Treat_P, 
             mut_id, mutation_category3),
    max_freq = max(freq0),
    max_timepoint = Timepoint[which.max(freq0)])
  png("./plots/mut_trajectories_pilus.png", width = 16, height = 6, units = "in",
      res = 150)
  f5b <- 
    ggplot(
      data = pildat,
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_name_num, lty = mutation_category3)) +
    geom_line(lwd = 1.05, alpha = 0.6) +
    geom_text_repel(data = pildat_sum,
                    aes(x = max_timepoint, y = max_freq, label = gene_num),
                    nudge_y = 0.05, size = 4, show.legend = FALSE) +
    # geom_text(data = pildat_sum,
    #           aes(x = max_timepoint, y = max_freq+0.05, label = gene_num),
    #           size = 4, show.legend = FALSE) +
    facet_nested(Treat_P * Treat_C ~ Pop,
                 labeller = labeller(
                   Pop = function(x) {paste("Rep", LETTERS[as.numeric(x)-2])},
                   Treat_C = c("NC" = "Pseudo.-only",
                               "C" = "All bacteria"),
                   Treat_P = c("NP" = "No-phage",
                               "P" = "Phage-added"))) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8), limits = c(0, 9)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
    scale_color_discrete(name = "Pilus gene") +
    scale_linetype_manual(
      name = "Putative mutation effect",
      breaks = c("loss of function", "nonsyn change", "intergenic", "silent"),
      labels = c("loss of function", "nonsyn. change", "intergenic", "silent"),
      values = c(1, 6, 2, 3)) +
    labs(x = "Timepoint (days)", y = "Frequency") +
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 18),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 16),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15)) +
    NULL
  print(f5b)
  dev.off()
  
  png("./plots/fig4_mut_trajectories_combo.png", width = 14, height = 11, units = "in",
      res = 150)
  print(
    cowplot::plot_grid(f5a, f5b, rel_heights = c(1, 0.6), 
                     ncol = 1, labels = "AUTO",
                     align = "hv", axis = "lr")
  )
  dev.off()
  
  png("./plots/mut_freqs_pilus.png", width = 10, height = 10, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(pildat, freq0 > 0),
      aes(x = gene_name, y = freq0, color = as.factor(Pop), group = Pop,
          shape = mutation_category3)) +
      geom_point(alpha = 0.5, position = position_dodge(width = 0.5)) +
      facet_nested(Treat_P * Treat_C ~ .,
                   labeller = labeller(
                     Treat_C = c("NC" = "Pseudomonas-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_color_manual(
        name = "Population", values = my_cols,
        breaks = 3:10, labels =  LETTERS[(3:10)-2]) +
      scale_fill_manual(
        name = "Population", values = my_cols,
        breaks = 3:10, labels =  LETTERS[(3:10)-2]) +
      scale_shape_manual(
        name = "Putative mutation effect",
        breaks = c("loss of function", "nonsyn change", "intergenic", "silent"),
        labels = c("loss of function", "nonsyn. change", "intergenic", "silent"),
        values = c(15, 16, 17, 18)) +
      labs(x = "Gene name", y = "Frequency") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
            plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "npc")) +
      NULL
  )
  dev.off()
}

###Analyze pilus category frequencies ----
gcat_freq <- summarize(
  group_by(filter(mutpop, mutation_category != "snp_synonymous"), 
           gene_catpilnot, across(any_of(sampcols))),
  total_freq = sum(freq0),
  indep_freq = 1-prod(1-freq0))

gcatfreq_sum <-
  summarize(group_by(filter(gcat_freq, Timepoint != 0),
                     gene_catpilnot, Pop, Treat, Treat_C, Treat_P),
            avgindep_freq = mean(indep_freq),
            avgtotal_freq = mean(total_freq))

####Plots ----
ggplot(data = filter(gcat_freq, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = total_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(Pop, Timepoint))) +
  facet_grid(Treat_P ~ gene_catpilnot) +
  #scale_y_log10() +
  NULL

ggplot(data = filter(gcat_freq, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(Pop, Timepoint))) +
  facet_grid(Treat_P ~ gene_catpilnot) +
  #scale_y_log10() +
  NULL

ggplot(filter(gcatfreq_sum, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = avgindep_freq)) +
  geom_line(aes(group = Pop)) +
  facet_grid(~Treat_P)
ggplot(filter(gcatfreq_sum, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = avgtotal_freq)) +
  geom_line(aes(group = Pop)) +
  facet_grid(~Treat_P)

png("./plots/gene_freq_pop_diff_indep_PNP.png", width = 6, height = 5,
    units = "in", res = 150)
fs17b <-
  ggplot(data = filter(gcat_freq, 
                     gene_catpilnot == "pilus", Timepoint == 8),
       aes(x = Treat_P, y = indep_freq)) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop))) +
  facet_grid(~ Treat_C,
             labeller = labeller(
               Treat_C = c("NC" = "Only\nPseudo.",
                           "C" = "All\nBacteria"))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "Cumulative frequency\nTimepoint 8\n(assuming indep.)") +
  scale_x_discrete(breaks = c("NP", "P"),
                   labels = c("No Phage", "Phage Added")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13)) +
  NULL
fs17b
dev.off()

png("./plots/gene_freq_pop_diff_indep_CNC.png", width = 6, height = 5,
    units = "in", res = 150)
pindepfreq_bytreat <- 
  ggplot(data = filter(gcat_freq, 
                       gene_catpilnot == "pilus", Timepoint != 0),
         aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop, Timepoint))) +
  facet_grid(~ Treat_P,
             labeller = labeller(
               Treat_P = c("NP" = "No Phage",
                           "P" = "Phage\nAdded"))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Timepoint\n(days)", values = my_cols[c(1:3, 7)]) +
  labs(x = "", y = "Cumulative pilus mut. freq.\n(assuming indep.)") +
  scale_x_discrete(breaks = c("NC", "C"),
                   labels = c("Only Pseudo.", "All Bacteria")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  NULL
pindepfreq_bytreat
dev.off()

png("./plots/gene_freq_pop_diff_total_PNP.png", width = 6, height = 5,
    units = "in", res = 150)
fs17a <-
  ggplot(data = filter(gcat_freq, 
                     gene_catpilnot == "pilus", Timepoint == 8),
       aes(x = Treat_P, y = total_freq)) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop))) +
  facet_grid(~ Treat_C,
               labeller = labeller(
                 Treat_C = c("NC" = "Only\nPseudo.",
                             "C" = "All\nBacteria"))) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "Cumulative frequency\nTimepoint 8\n(assuming exclus.)") +
  scale_x_discrete(breaks = c("NP", "P"),
                   labels = c("No Phage", "Phage Added")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13)) +
  NULL
fs17a
dev.off()

png("./plots/gene_freq_pop_diff_total_CNC.png", width = 6, height = 5,
    units = "in", res = 150)
f6b <- 
  ggplot(data = filter(gcat_freq, 
                     gene_catpilnot == "pilus", Timepoint != 0),
       aes(x = Treat_C, y = total_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop, Timepoint))) +
  facet_grid(~ Treat_P,
               labeller = labeller(
                 Treat_P = c("NP" = "No Phage",
                             "P" = "Phage\nAdded"))) +
  #scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Timepoint\n(days)", values = my_cols[c(1:3, 7)]) +
  labs(x = "", y = "Cumulative pilus mut. freq.\n(assuming exclus.)") +
  scale_x_discrete(breaks = c("NC", "C"),
                   labels = c("Only Pseudo.", "All Bacteria")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  NULL
f6b
dev.off()

png("./plots/figS17_gene_freq_pop_diff_PNP.png", width = 12, height = 5,
    units = "in", res = 150)
cowplot::plot_grid(fs17a, fs17b,
                   nrow = 1, align = "hv", axis = "tblr", labels = "AUTO")
dev.off()



ggplot(data = filter(gcat_freq, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = total_freq, color = gene_catpilnot)) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop))) +
  facet_grid(Treat_P ~ Timepoint) +
  #scale_y_log10() +
  NULL

ggplot(data = filter(gcat_freq, gene_catpilnot == "pilus"),
       aes(x = Treat_C, y = indep_freq, color = gene_catpilnot)) +
  geom_point() +
  geom_line(aes(group = paste(gene_catpilnot, Pop))) +
  facet_grid(Treat_P ~ Timepoint) +
  #scale_y_log10() +
  NULL


####Stats ----

#Effect of C
totalfreq_full_lme <- 
  lmer(total_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
             data = filter(gcat_freq, gene_catpilnot == "pilus",
                           Treat_P == "P", Timepoint != 0),
       REML = FALSE)
totalfreq_reduced_lme <-
  lmer(total_freq ~ as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(totalfreq_reduced_lme, totalfreq_full_lme)
totalfreq_full_lme

indepfreq_full_lme <- 
  lmer(indep_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
indepfreq_reduced_lme <- 
  lmer(indep_freq ~ as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(indepfreq_reduced_lme, indepfreq_full_lme)
indepfreq_full_lme

#Effect of P
totalfreq_full_lmeP <- 
  lmer(total_freq ~ Treat_P + Treat_C + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Timepoint == 8),
       REML = FALSE)
totalfreq_reduced_lmeP <-
  lmer(total_freq ~ Treat_C + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Timepoint == 8),
       REML = FALSE)
anova(totalfreq_reduced_lmeP, totalfreq_full_lmeP)
totalfreq_full_lmeP

indepfreq_full_lmeP <- 
  lmer(indep_freq ~ Treat_P + Treat_C + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Timepoint == 8),
       REML = FALSE)
indepfreq_reduced_lmeP <- 
  lmer(indep_freq ~ Treat_C + (1|Pop),
       data = filter(gcat_freq, gene_catpilnot == "pilus",
                     Timepoint == 8),
       REML = FALSE)
anova(indepfreq_reduced_lmeP, indepfreq_full_lmeP)
indepfreq_full_lme


### GO enrichment analysis ----

#GO enrichment analysis using cluterProfiler in R
#GSEA analysis looks for terms that are biased up/down across all genes w/ term
#enrichment looks for overrepresentation of terms in genes at the top
# -- this one is what we use

#Calculate gene frequencies
gnamefreq <- summarize(
  group_by(
    filter(mutpopd, !is.na(Locus.Tag), mutation_category != "snp_synonymous"), 
    Locus.Tag, Product,
    across(any_of(sampcols))),
  indep_freq = 1-prod(1-freq0))

#There are no missing combos to fill in (filled in previously)

non0_mean <- function(x, cond) {
  x <- x[cond]
  if(length(x) > 0) {return(mean(x))} else {return(0)}
}
gnamefreq_sum <- 
  summarize(
    group_by(filter(gnamefreq, indep_freq >= 0.1), Locus.Tag, Product),
    PvsNP_non0mean = 
      non0_mean(indep_freq, Treat_P == "P" & Timepoint == 8) -
      non0_mean(indep_freq, Treat_P == "NP" & Timepoint == 8),
    CvsNC_inNP_non0mean =
      non0_mean(indep_freq, Treat_C == "C" & Treat_P == "NP" & Timepoint == 8) -
      non0_mean(indep_freq, Treat_C == "NC" & Treat_P == "NP" & Timepoint == 8),
    CvsNC_inP_non0mean =
      non0_mean(indep_freq, Treat_C == "C" & Treat_P == "P") -
      non0_mean(indep_freq, Treat_C == "NC" & Treat_P == "P")
  )

gochangecutoff <- 0.15
if(!any(gnamefreq_sum$PvsNP_non0mean > gochangecutoff)) {
  goPvsNP <- NULL
} else {
  goPvsNP <-
    cbind(Treat = "P", Context = "All",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, PvsNP_non0mean > gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
if(!any(gnamefreq_sum$CvsNC_inNP_non0mean > gochangecutoff)) {
  goCvsNC_inNP <- NULL
} else {
  goCvsNC_inNP <-
    cbind(Treat = "C", Context = "NP",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, CvsNC_inNP_non0mean > gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
if(!any(gnamefreq_sum$CvsNC_inP_non0mean > gochangecutoff)) {
  goCvsNC_inP <- NULL
} else {
  goCvsNC_inP <-
    cbind(Treat = "C", Context = "P",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, CvsNC_inP_non0mean > gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
if(!any(gnamefreq_sum$PvsNP_non0mean < -gochangecutoff)) {
  goNPvsP <- NULL
} else {
  goNPvsP <-
    cbind(Treat = "NP", Context = "All",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, PvsNP_non0mean < -gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
if(!any(gnamefreq_sum$CvsNC_inNP_non0mean < -gochangecutoff)) {
  goNCvsC_inNP <- NULL
} else {
  goNCvsC_inNP <-
    cbind(Treat = "NC", Context = "NP",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, CvsNC_inNP_non0mean < -gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
if(!any(gnamefreq_sum$CvsNC_inP_non0mean < -gochangecutoff)) {
  goNCvsC_inP <- NULL
} else {
  goNCvsC_inP <-
    cbind(Treat = "NC", Context = "P",
          as.data.frame(
              enricher(
                gene = filter(gnamefreq_sum, CvsNC_inP_non0mean < -gochangecutoff)$Locus.Tag,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 2)))
}
gores <- do.call(rbind,
                 list(goCvsNC_inNP, goCvsNC_inP, 
                      goNCvsC_inNP, goNCvsC_inP,
                      goPvsNP, goNPvsP))
gores <- filter(gores, pvalue <= 0.05)
gores <- mutate(gores,
                FoldPlot = ifelse(Treat %in% c("NC", "NP"),
                                  -FoldEnrichment,
                                  FoldEnrichment),
                TreatPlot = ifelse(grepl("P", Treat), "P", "C"))

write.csv(filter(gores, !grepl("GO:", Description), Count > 1),
          "./clean_data/GOresults.csv", row.names = FALSE)


f7a <- ggplot(
  mutate(
    arrange(
      filter(gores, !grepl("GO:", Description), pvalue < 0.01, 
             TreatPlot == "P", Count > 1),
      FoldPlot),
    Description = factor(Description, levels = Description)),
  aes(x = FoldPlot, y = Description, size = -log10(p.adjust))) +
  geom_point() +
  scale_size_continuous(breaks = c(4, 3, 2, 1), 
                        labels = c("0.0001", "0.001", "0.01", "0.1"),
                        limits = c(4, 0.5), range = c(0, 5),
                        name = "Adjusted p-value") +
  # scale_color_manual(breaks = c(TRUE, FALSE),
  #                    values = my_cols[c(6, 5)]) +
  geom_vline(xintercept = 0, lty = 2) +
  #guides(color = "none") +
  labs(x = "Fold change in phage-added treatment vs no-phage", y = "") +
  theme_bw() +
  NULL

f7b <- ggplot(
  mutate(
    arrange(filter(gores, !grepl("GO:", Description), pvalue < 0.01, 
                   TreatPlot == "C", Count > 1),
            FoldPlot),
    Description = factor(Description, levels = Description)), 
  aes(x = FoldPlot, y = Description, size = -log10(p.adjust))) +
  geom_point() +
  scale_size_continuous(breaks = c(3, 2, 1), labels = c("0.001", "0.01", "0.1"),
                        limits = c(3, 0.5), range = c(0, 5),
                        name = "Adjusted p-value") +
  facet_grid(Context ~ ., scales = "free_y",
             labeller = labeller(
               Context = c("NP" = "In no-phage\ntreatment",
                           "P" = "In phage-added\ntreatment"))) +
  geom_vline(xintercept = 0, lty = 2) +
  #guides(color = "none") +
  labs(x = "Fold change in all bacteria treatment vs Pseudomonas-only", y = "") +
  theme_bw() +
  NULL
  
png("./plots/GOhits.png", width = 12, height = 12,
    units = "in", res = 150)
cowplot::plot_grid(f7a, f7b, nrow = 2, align = "v", axis = "lr",
                     labels = c("A", "B"))
dev.off()

### Mutation frequencies vs population densities ----
densH6 <- filter(cfudensity, Strain == "H6", !Pop %in% c(1, 2))
densH6 <- dplyr::rename(densH6,
                        H6_pfuml = cfu_ml)

gcat_freqw <- pivot_wider(select(ungroup(gcat_freq),
                                 gene_catpilnot, Timepoint, Pop,
                                 Treat, Treat_C, Treat_P,
                                 total_freq, indep_freq),
                          names_from = gene_catpilnot,
                          values_from = c(total_freq, indep_freq))
gcat_freqw_sum <- summarize(group_by(filter(gcat_freqw, Timepoint != 0),
                                     Pop, Treat, Treat_C, Treat_P),
                            across(contains("freq"), mean))
                            
dens_freq <- left_join(gcat_freqw, 
                       mutate(densH6, Pop = as.character(Pop)))
dens_freq <- mutate(dens_freq,
                    H6_pfuml = ifelse(is.na(H6_pfuml) & Treat_P == "NP",
                                      400,
                                      H6_pfuml))

pindepfreq_H6dens <- 
  ggplot(data = filter(dens_freq, Timepoint != 0,
                     Treat %in% c("CP", "NCP")),
       aes(x = indep_freq_pilus, y = log10(H6_pfuml), 
           color = Treat_C, fill = Treat_C, shape = Treat_C)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = math_format(10^.x),
                     breaks = c(2, 4, 6, 8, 10)) +
  geom_smooth(method = "lm", show.legend = FALSE) +
  scale_color_manual(name = "Community Treatment", 
                     breaks = c("NC", "C"),
                     labels = c("Only Pseudo.", "All Bacteria"),
                     values = my_cols[c(3, 7)]) +
  scale_fill_manual(name = "Community Treatment", 
                    breaks = c("NC", "C"),
                    labels = c("Only Pseudo.", "All Bacteria"),
                    values = my_cols[c(3, 7)]) +
  scale_shape_manual(name = "Community Treatment",
                     breaks = c("NC", "C"),
                     labels = c("Only Pseudo.", "All Bacteria"),
                     values = c(21, 23)) +
  labs(x = "Cumulative pilus mut. freq.\n(assuming indep.)",
       y = "Phage population density (pfu/mL)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

f6a <- 
  ggplot(data = filter(dens_freq, Timepoint != 0,
                     Treat %in% c("CP", "NCP")),
       aes(x = total_freq_pilus, y = log10(H6_pfuml), 
           color = Treat_C, fill = Treat_C, shape = Treat_C)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = math_format(10^.x),
                     breaks = c(2, 4, 6, 8, 10)) +
  geom_smooth(method = "lm", show.legend = FALSE) +
  scale_color_manual(name = "Community Treatment", 
                     breaks = c("NC", "C"),
                     labels = c("Only Pseudo.", "All Bacteria"),
                     values = my_cols[c(3, 7)]) +
  scale_fill_manual(name = "Community Treatment", 
                    breaks = c("NC", "C"),
                    labels = c("Only Pseudo.", "All Bacteria"),
                    values = my_cols[c(3, 7)]) +
  scale_shape_manual(name = "Community Treatment",
                     breaks = c("NC", "C"),
                     labels = c("Only Pseudo.", "All Bacteria"),
                     values = c(21, 23)) +
  labs(x = "Cumulative pilus mut. freq.\n(assuming exclus.)",
       y = "Phage population density (pfu/mL)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

png("./plots/fig5_density_vs_muts.png", width = 12, height = 5,
    units = "in", res = 150)
cowplot::plot_grid(
  f6a + labs(x = "Total pilus mut. freq", y = "Phage population density (pfu/mL)"),
  f6b + labs(y = "Total pilus mut. freq"),
  nrow = 1, rel_widths = c(1, 0.8),
  align = "hv", axis = "tb", labels = "AUTO"
)
dev.off()

png("./plots/figS_density_vs_muts_all.png", width = 12, height = 10,
    units = "in", res = 150)
cowplot::plot_grid(
  cowplot::plot_grid(
    f6a + guides(color = "none", fill = "none", shape = "none"),
    pindepfreq_H6dens + guides(color = "none", fill = "none", shape = "none"), 
    cowplot::get_legend(f6a),
    rel_widths = c(1, 1, 0.5),
    nrow = 1, labels = c("A", "B", ""),
    align = "hv", axis = "tblr"),
  cowplot::plot_grid(f6b + guides(color = "none"),
                     pindepfreq_bytreat + guides(color = "none"),
                     cowplot::get_legend(f6b),
                     rel_widths = c(1, 1, 0.5),
                     nrow = 1, align = "h", axis = "tb",
                     labels = c("C", "D")),
  nrow = 2
)
dev.off()

densfreq_indeplm <- 
  lm(data = filter(dens_freq, Timepoint != 0, Treat %in% c("CP", "NCP")),
     log10(H6_pfuml) ~ Treat + indep_freq_pilus:Treat)
densfreq_totallm <- 
  lm(data = filter(dens_freq, Timepoint != 0, Treat %in% c("CP", "NCP")),
     log10(H6_pfuml) ~ Treat + total_freq_pilus:Treat)

anova(densfreq_indeplm)
anova(densfreq_totallm)

summary(densfreq_indeplm)
summary(densfreq_totallm)



##Isolate analyses ----

### Isolate samples merge with phenotype data ----
#Filter out sample that had 0 coverage (so mutation is a full-genome deletion)
mutisol <- filter(mutisol, title != "5174_67_PA_T8_P4_NCP_#5")

mutisol <- left_join(
  mutisol,
  select(
    mutate(
      filter(read.csv("./clean_data/aggregate_isolate_data.csv", stringsAsFactors = FALSE,
                      colClasses = "character"),
             Bact_timepoint != "pilA"),
      Treat = Bact_treat, 
      Treat_C = Bact_treat_C, 
      Treat_P = Bact_treat_P,
      Timepoint = as.numeric(Bact_timepoint), 
      Pop = Bact_pop, Isol = Bact_isol),
    any_of(sampcols))
)

### Make nice table ----

##Note: five isolates have only the four shared fixed mutations 
##      (i.e. they're ancestral)
# 5174_58 - PA T1 P4 CP #3 - susceptible
# 5174_64 - PA T4 P4 NCP #3 - susceptible
# 5174_69 - PA T1 P5 NCP #1 - susceptible
# 5174_70 - PA T1 P5 NCP #2 - resistant
# 5174_129 - PA T8 P10 CP #1 - resistant

nice_isoltable <- 
  pivot_wider(
    select(
      mutate(ungroup(mutisol),
             mut_combin = ifelse(
               mut_id %in% ids_out,
               NA,
               paste0(mutation_category2, " [", concise_mut, "]")),
             genecat_tmp = ifelse(
               gene_category %in% c("pilus", "mutator", "biofilm"),
               gene_category,
               "other")),
      Timepoint, Treat, Pop, Isol, Treat_C, Treat_P,
      clearing_avg, EOP_avg, Plaque_morphology, resis_cat, 
      twitch_cat,
      genecat_tmp, mut_combin),
    names_from = genecat_tmp,
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
          file = "./clean_data/nice_isolmuts.csv", row.names = FALSE)

### Create plottable frequency of mutations ----
isol_freq_sum <- 
  summarize(
    group_by(
      mutate(mutisol,
             gene_category = case_when(
               gene_category == "pilus" ~ "pilus",
               gene_category == "biofilm" ~ "biofilm",
               TRUE ~ "other")),
      gene_category, across(any_of(sampcols)), mutation_category3),
    mut_genes = {
      mygenes <- gene_short[!mut_id %in% ids_out]
      case_when(
        any(mutation_category_manual == "intergenic") ~ "intergenic",
        length(mygenes) == 1 ~ paste(mygenes, collapse = " "),
        length(mygenes) == 0 ~ "none",
        length(mygenes) > 1 ~ "multiple",
        TRUE ~ "tbd")
    })

#Fill in missing combos
isol_freq_sum <-
  complete(
    ungroup(isol_freq_sum),
    nesting(!!!syms(intersect(sampcols, colnames(isol_freq_sum)))),
    gene_category,
    #Columns left out (and therefore filled in by "none"):
    # mut_genes
    fill = list("mut_genes" = "none", "mutation_category3" = "none")
  )

#Manually recode large deletion to "multiple"
isol_freq_sum <- 
  mutate(isol_freq_sum,
         mut_genes = ifelse(mut_genes == "pilU|pilT",
                            "multiple",
                            mut_genes))

pil_order <- unique(isol_freq_sum[isol_freq_sum$gene_category == "pilus", ]$mut_genes)
pil_order <-
  c("none", "multiple", #"intergenic",
    pil_order[!pil_order %in% c("none", "multiple", "intergenic")])


### Plots of mutation types ----
if(make_slowplots) {
  png("./plots/isolmut_pilus_freq.png",
      width = 9, height = 5, units = "in", res = 150)
  print(
    ggplot(
      mutate(
        filter(isol_freq_sum, gene_category == "pilus"),
        mut_genes = factor(mut_genes, levels = pil_order)),
      aes(x = resis_cat)) +
      geom_bar(aes(fill = mut_genes),position = "fill") +
      facet_nested(. ~ Treat_P * Treat_C,
                   labeller = labeller(
                     Treat_C = c("NC" = "Pseudomonas-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_x_discrete(limits = c("suscep", "part", "resis"),
                       labels = c("Susceptible", "Partial", "Resistant")) +
      scale_fill_manual(
        name = "Pilus Genes\nMutated",
        limits = pil_order,
        values = c("gray70", "black", #"gray30",
                   viridis_pal()(length(pil_order)-2))) +
      labs(x = "Isolates' Resistance to Phage",
           y = "Frequency") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_blank(),
            axis.title = element_text(size = 18),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 18),
            strip.text = element_text(size = 18)) +
      geom_text(data = summarize(group_by(isol_freq_sum, resis_cat, Treat_C, Treat_P),
                                 n = length(unique(title))),
                aes(y = 0.1, label = paste0("n = ", n)), size = 5) +
      NULL)
  dev.off()
  
  png("./plots/isolmut_pilus_muttypes.png",
      width = 9, height = 5, units = "in", res = 150)
  print(
    ggplot(
      mutate(
        filter(isol_freq_sum, gene_category == "pilus"),
        mutation_category3 = factor(
          mutation_category3,
          levels = c("none", "silent", "nonsyn change", "loss of function"))),
      aes(x = resis_cat)) +
      geom_bar(aes(fill = mutation_category3),position = "fill") +
      facet_nested(. ~ Treat_P * Treat_C,
                   labeller = labeller(
                     Treat_C = c("NC" = "Pseudomonas-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_x_discrete(limits = c("suscep", "part", "resis"),
                       labels = c("Susceptible", "Partial", "Resistant")) +
      scale_fill_manual(
        name = "Pilus Gene(s)\nMutation Effect",
        limits = c("none", "silent", "nonsyn change", "loss of function"),
        values = brewer_pal(palette = "Blues")(5)[-1]) +
      labs(x = "Isolates' Resistance to Phage",
           y = "Frequency") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_blank(),
            axis.title = element_text(size = 18),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 18),
            strip.text = element_text(size = 18)) +
      geom_text(data = summarize(group_by(isol_freq_sum, resis_cat, Treat_C, Treat_P),
                                 n = length(unique(title))),
                aes(y = 0.1, label = paste0("n = ", n)), size = 5) +
      NULL)
  dev.off()
}


#Filter out fixed-nonvariate mutations
#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
mutisol <- filter(mutisol, !mut_id %in% ids_out)

#Note: we sequenced more NCP isols than CP

png("./plots/isol_pilmuts_byresis.png", width = 8, height = 10, units = "in",
    res = 150)
print(
  ggplot(
    data = filter(mutisol, gene_category == "pilus"),
    aes(x = gene_name, y = resis_cat, 
        color = mutation_category3, fill = mutation_category3,
        group = Pop)) +
    geom_point(alpha = 0.5, 
               position = position_jitterdodge(jitter.width = 0,
                                               jitter.height = 0.05,
                                               dodge.width = 0.6)) +
    facet_grid(Treat ~ .) +
    scale_shape_manual(values = 21:24) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    NULL
)
dev.off()

png("./plots/isol_pilmuts_bytwitch.png", width = 8, height = 10, units = "in",
    res = 150)
print(
  ggplot(
    data = filter(mutisol, gene_category == "pilus"),
    aes(x = gene_name, y = twitch_cat, 
        color = mutation_category3, fill = mutation_category3,
        group = Pop)) +
    geom_point(alpha = 0.5, 
               position = position_jitterdodge(jitter.width = 0,
                                               jitter.height = 0.05,
                                               dodge.width = 0.6)) +
    facet_grid(Treat ~ .) +
    scale_shape_manual(values = 21:24) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    NULL
)
dev.off()

print(
  ggplot(
    data = filter(mutisol, gene_category == "pilus"),
    aes(x = gene_name, y = as.numeric(auc_norm), color = Pop, fill = Pop,
        group = Pop, shape = mutation_category3)) +
    geom_point(alpha = 0.5) +
    facet_grid(~ Treat) +
    scale_shape_manual(values = 21:24) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    NULL
)

print(
  ggplot(
    data = filter(mutisol, gene_category == "pilus"),
    aes(x = gene_name, y = as.numeric(lag_norm), color = Pop, fill = Pop,
        group = Pop, shape = mutation_category3)) +
    geom_point(alpha = 0.5) +
    facet_grid(~ Treat) +
    scale_shape_manual(values = 21:24) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    NULL
)

png("./plots/isol_pilstruc_byresis.png", width = 10, height = 7, units = "in",
    res = 150)
print(
  ggplot(
    data = filter(mutisol, gene_category == "pilus"),
    aes(x = gene_name, y = resis_cat, color = Pop, fill = Pop,
        group = Pop, shape = mutation_category3)) +
    geom_point(alpha = 0.5, 
               position = position_jitterdodge(jitter.width = 0,
                                               jitter.height = 0.05,
                                               dodge.width = 0.6)) +
    facet_grid(Treat ~ pilmut_struc, scales = "free_x") +
    scale_shape_manual(values = 21:24) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    NULL
)
dev.off()

## Comparing isolate vs population sequencing ----

#For isolate samples, fill in all the missing mutation-sample combinations
# with NA
mutcomp <- 
  complete(
    ungroup(mutate(mutisol, freq = 1)),
    #Create the nesting defining each unique mutation
    nesting(!!!syms(intersect(mutcols, colnames(mutisol)))),
    #Create the nesting defining each population
    nesting(!!!syms(intersect(sampcols, colnames(mutisol))))
    #Columns left out (and therefore filled in by NA by default):
    # mut_and_samp_cols
  )

mutcomp <-
  select(mutcomp,
         -any_of(c(setdiff(mut_and_samp_cols, "freq"),
                   "title", "seqrun", "SeqCoastid")))

mutcomp <- 
  mutate(ungroup(mutcomp),
         freq0 = ifelse(is.na(freq), 0, freq),
         resis_cat_num = recode(resis_cat,
                                resis = 1, part = 0.5, suscep = 0))

mutcomp <-
  summarize(
    group_by(mutcomp,
             across(any_of(c(mutcols,
                             "seqtype", "Timepoint", "Pop",
                             "Treat", "Treat_C", "Treat_P")))),
    isolfreq = mean(freq0),
    avgresis = mean(resis_cat_num))

mutcomp <-
  full_join(
    mutcomp,
    select(filter(mutpop),
           -c(multiple_polymorphic_SNPs_in_same_codon,
              new_read_count, new_read_count_basis,
              ref_read_count, ref_read_count_basis,
              seqrun, SeqCoastid,
              Isol, title)))

#Re-run crossing to fill in all the missing mutation-sample combinations
# with NA
mutcomp <- 
  complete(
    ungroup(mutcomp),
    #Create the nesting defining each unique mutation
    nesting(!!!syms(intersect(mutcols, colnames(mutcomp)))),
    #Create the nesting defining each sample
    nesting(!!!syms(intersect(sampcols, colnames(mutcomp)))),
    #Columns left out (and therefore filled in by NA by default):
    # mut_and_samp_cols
  )

mutcomp <- 
  mutate(ungroup(mutcomp),
         freq0 = ifelse(is.na(freq), 0, freq))

mutcomp <- mutate(
  group_by(mutcomp, seqtype, mut_id, Treat, Pop),
  everabove_010_orisolmut = 
    any((!is.na(everabove_010) & everabove_010 == TRUE) | 
          (!is.na(isolfreq) & isolfreq > 0)))

tmp <- filter(mutcomp,
              everabove_010_orisolmut == TRUE, 
              gene_category == "pilus", Treat %in% c("NCP", "CP"))

png("./plots/isolmut_vs_popmut_pilus.png", width = 16, height = 8, units = "in",
    res = 150)
ggplot(tmp,
       aes(x = Timepoint, y = freq0, color = gene_name, fill = gene_name,
           shape = as.factor(avgresis), group = mut_id)) +
  geom_line() +
  geom_point(alpha = 0.5, aes(size = isolfreq)) +
  facet_grid(Treat ~ Pop) +
  scale_size_continuous(breaks = c(0, 0.5, 1), range = c(0, 4))
dev.off()


tmp <- mutate(filter(mutcomp,
              everabove_010_orisolmut == TRUE, 
              Treat %in% c("NCP", "CP")),
              Pop = as.numeric(Pop))

png("./plots/isolmut_vs_popmut.png", width = 16, height = 8, units = "in",
    res = 150)
ggplot(tmp,
       aes(x = Timepoint, y = freq0, 
           color = gene_catpilnot,
           group = mut_id)) +
  geom_line() +
  geom_point(alpha = 0.5, aes(size = isolfreq)) +
  facet_nested(Treat_P * Treat_C ~ Pop,
               labeller = labeller(
                 Pop = function(x) {paste("Rep", LETTERS[as.numeric(x)-2])},
                 Treat_C = c("NC" = "Pseudomonas-only",
                             "C" = "All bacteria"),
                 Treat_P = c("NP" = "No-phage",
                             "P" = "Phage-added"))) +
  scale_size_continuous(name = "Frequency among\nisolate sequencing",
                        breaks = c(0, 0.5, 1),
                        range = c(0, 4)) +
  scale_color_manual(breaks = c("pilus", "other"),
                     values = my_cols[c(6, 8)],
                     name = "Gene function") +
  labs(y = "Frequency in population sequencing")
dev.off()


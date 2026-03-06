##TODO ----
#Best ideas
# Get total # reads (not just aligned reads)
# run isols in polym mode, try to figure out contamination (where are bg muts
#  (in isols!?!)
# compare pilus gene size to # of muts
# Plot mutations in N samples vs freq plot
#  (n reflects mut accessibility, freq reflects selection)
# why are cov-based densities of other bact higher in NC treats?
# PCA mutisol
# Convert big deletions into deletions at ea gene via annotations file
# PCA but casting paste(gene_name, mutation_cat3 a la mutisol)
# Diversity metrics of each pop-timepoint
# Look at pilus structure & where genes act in pilus construction
# ID genes mutated more highly in one treat vs others
# Look for dift fates of phages depending on initial inocula 3-10 & pilus muts
# Is there anything in common among the populations where phage went ~extinct?
#  Some way those 5 pops cluster?
# How often are the same gene targeted with different mutations?
# fst
# parallel selection - Toby
# outlier analysis - Chelsea/Toby

#Good ideas
# ML on mutpop to predict treatment or phage density
# ML on mutisol to predict treatment or resistance or twitching
# Mutations associated with phage density

#Other ideas
# plot mutation trajectories excluding those that appear in control treatments
# Discriminant analysis to look for diff genes
# ID mutations that are moving together
# correlated coevol mutation


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
#library(BiocManager)
#BiocManager::install("clusterProfiler")
#library(clusterProfiler)
#BiocManager::install("AnnotationHub")
#library(AnnotationHub)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

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
                 Strain = case_when(
                   str_starts(seq_id, "EF_") ~ "EF",
                   str_starts(seq_id, "SA_") ~ "SA",
                   str_starts(seq_id, "AX1_") ~ "AX",
                   seq_id == "H6_1" ~ "H6",
                   seq_id == "PA14_1" ~ "PA",
                   TRUE ~ "unknown"),
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
                                levels = c("Anc", "NN", "CNP", "NCP", "CP")),
                 included_in_treat = 
                   Strain == "PA" |
                   (Strain %in% c("SA", "EF", "AX") & Treat_C == "C") |
                   (Strain == "H6" & Treat_P == "P"))

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

cfudensity <- readRDS("density.rds")
cfudensity <- mutate(cfudensity,
                     Treat_C = treat_C, Treat_P = treat_P,
                     Treat = paste0(Treat_C, Treat_P),
                     Treat = ifelse(Treat == "NCNP", "NN", Treat),
                     Pop = as.character(Pop))
#We're missing T2 CP3 CP4 CP5 AX
cfudensity <- mutate(group_by(cfudensity,
                              Timepoint, Pop, Treat, Treat_C, Treat_P),
                     rel_freq_cfu = cfu_ml/sum(cfu_ml, na.rm = TRUE))
covbysp <- left_join(covbysp,
                    select(cfudensity,
                           Timepoint, Pop, Treat, Treat_C, Treat_P,
                           Strain, cfu_ml, rel_freq_cfu))
covbysp <- mutate(covbysp,
                  rel_freq_cfu = ifelse(is.na(rel_freq_cfu), 0, rel_freq_cfu))

if(make_slowplots) {
  png("./plots/coverage_by_contig.png", width = 12, height = 4, units = "in",
      res = 150)
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
  dev.off()
  
  png("./plots/coverage_by_contig_phageseq.png", width = 12, height = 4, units = "in",
      res = 150)
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
  dev.off()
  
  png("./plots/coverage_PAby_TimepointTreat.png", width = 6, height = 3, 
      units = "in", res = 150)
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
  dev.off()
  
  png("./plots/coverage_bycontig_bytreat.png", width = 15, height = 12, units = "in",
      res = 150)
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
  dev.off()
    
  png("./plots/density_fromseq.png", width = 10, height = 7, units = "in",
      res = 150)
  ##TODO: use total reads as denom, not total aligned reads
  print(
    ggplot(
      data = filter(covbysp, Treat != "Anc"),
      aes(x = Timepoint, y = log10(rel_freq_noAX), color = Strain)) +
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
  png("./plots/density_fromseq.png", width = 10.5, height = 7, units = "in",
      res = 150)
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
  dev.off()
}


##Initial breseq data reading & organization ----
mutcols <- 
  c("aa_new_seq", "aa_position", "aa_ref_seq", "codon_new_seq", "codon_number", 
    "codon_position", "codon_ref_seq", "gene_name", "gene_position", 
    "gene_product", "gene_strand", "genes_inactivated", "genes_overlapping", 
    "genes_promoter", "insert_position", "locus_tag", "locus_tags_inactivated",
    "locus_tags_overlapping", "locus_tags_promoter", "mutation_category",
    "new_seq", "position", "position_end", "position_start",
    "ref_seq", "repeat_length", "repeat_new_copies", "repeat_ref_copies",
    "repeat_seq", "seq_id", "size", "snp_type", "transl_table", "type", "mut_id",
    "gene_category", "mutation_category2", "mutation_category3",
    "pilmut_struc", "gene_short", "concise_mut", "mut_combin")
sampcols <- c("seqrun", "SeqCoastid", "seqtype", "Timepoint", "Pop", "Treat", 
              "Isol", "title", "Treat_C", "Treat_P",
              "clearing_avg", "EOP_avg", "Plaque_morphology", "resis", 
              "resis_cat", "cm_day", "cm_day_norm", "twitch_cat",
              "wavelength_405", "wavelength_695",
              "max_percap", "max_dens", "lag", "auc", "max_percap_time", 
              "max_percap_dens", "min_dens",
              "max_percap_norm", "max_dens_norm", "auc_norm", "lag_norm",
              "biofilm_abs", "biofilm_norm", "biofilm_norm_log10",
              "biofilm_norm_log10_norm", "biofilm_block_name",
              "resis_cat_num")
mut_and_samp_cols <- 
  c("multiple_polymorphic_SNPs_in_same_codon",
    "new_read_count", "new_read_count_basis",
    "ref_read_count", "ref_read_count_basis",
    "freq", "freq0", "freq01",
    "everabove_025", "everabove_010", "everabove_005", 
    "everabove_010_noanc", "everabove_005_noanc", "fixed",
    "freq_change", "abs_freq_rate", "rel_freq_rate", "lfreq01", "freq_cat",
    "isolfreq", "everabove_010_orisolmut")

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
mutations <- mutate(group_by(mutations, across(any_of(mutcols))),
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


##Read in annotations from PA annotations file for later use ----
anno <- read.csv("PA14_annotated.csv")

##Shared analysis labeling ----
mutations <- mutate(
  ungroup(mutations),
  gene_category = case_when(
    grepl("CupB3", gene_name) ~ "pilus",
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
    grepl("T2SS", gene_name) ~ "T2SS",
    grepl("type II secretion system", gene_name) ~ "T2SS",
    grepl("flagella", gene_name) ~ "flagella",
    grepl("Thioredoxin family protein/Helix-turn-helix domain-containing protein", gene_name) ~ "mystery_1",
    grepl("FimX", gene_name) ~ "pilus",
    grepl("pil", gene_name, ignore.case = TRUE) ~ "pilus",
    grepl("pilus", gene_name) ~ "pilus",
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
    mutation_category2 %in% c("small_indel_intergenic", "snp_intergenic") ~ "intergenic",
    mutation_category2 %in% c("small_indel_inframe", "snp_nonsynonymous") ~ "nonsyn change",
    mutation_category2 %in% c("snp_synonymous") ~ "silent",
    TRUE ~ "NA"),
  pilmut_struc = ifelse(
    gene_category == "pilus",
    case_when(
      grepl("4b", gene_name) ~ "4b",
      grepl("IVB", gene_name) ~ "4b",
      grepl("TadZ", gene_name) ~ "unknown",
      grepl("PilT", gene_name) ~ "motor",
      grepl("PilB", gene_name) ~ "motor",
      grepl("PilQ", gene_name) ~ "OMP",
      grepl("PilS", gene_name) ~ "regulatory",
      grepl("PilM", gene_name) ~ "alignment_sc",
      grepl("PilC", gene_name) ~ "motor",
      grepl("pilC", gene_name) ~ "motor",
      grepl("PilN", gene_name) ~ "alignment_sc",
      grepl("PilZ", gene_name) ~ "motor",
      grepl("TadG", gene_name) ~ "unknown",
      grepl("PilU", gene_name) ~ "motor",
      grepl("PilR", gene_name) ~ "regulatory",
      grepl("PilO", gene_name) ~ "alignment_sc",
      grepl("FimV", gene_name) ~ "OMP",
      grepl("PilJ", gene_name) ~ "regulatory",
      grepl("PilF", gene_name) ~ "OMP",
      grepl("PilY1", gene_name) ~ "fiber",
      grepl("CupB3", gene_name) ~ "usher_chaperone_assembly",
      grepl("PilV", gene_name) ~ "fiber",
      grepl("PilG", gene_name) ~ "regulatory",
      grepl("PilK", gene_name) ~ "tbd",
      grepl("PilP", gene_name) ~ "alignment_sc",
      grepl("PilD", gene_name) ~ "usher_chaperone_assembly",
      grepl("PilE", gene_name) ~ "fiber",
      grepl("PilX", gene_name) ~ "fiber",
      grepl("FleP", gene_name) ~ "regulatory",
      grepl("PilM", gene_name) ~ "alignment",
      grepl("PilR", gene_name) ~ "regulatory",
      grepl("PilP", gene_name) ~ "alignment_sc",
      grepl("FimL", gene_name) ~ "regulatory",
      TRUE ~ "tbd"),
    NA),
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
             paste0(gene_short, " (intergenic ", nchar(new_seq), "bp INS)"),
             paste0(gene_short, " (", nchar(new_seq), "bp INS)")),
    mutation_category == "small_indel" & type == "SUB" ~
      ifelse(grepl("intergenic", gene_position),
             paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
             paste0(gene_short, " (nt ", ref_seq, "->", new_seq, ")")),
    mutation_category == "large_deletion" ~
      paste0(gene_short, " (", size, "bp DEL)"),
    TRUE ~ "nomatch")
)

#Split samples by type ----
mutpop <- filter(mutations, Isol == "Pop", seqtype == "PA")
mutisol <- filter(mutations, Isol != "Pop", seqtype == "PA")

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

mutpop <- mutate(group_by(mutpop, 
                          across(any_of(c(mutcols, "seqtype", "Treat", "Pop")))),
                 everabove_010 = any(freq0 > 0.1, na.rm = TRUE),
                 everabove_010_noanc = any(freq0[Timepoint != 0] > 0.1, na.rm = TRUE),
                 fixed = all(freq0 == 1))

### Filter out fixed-nonvariate mutations ----
mutpop <- mutate(group_by(mutpop,
                          across(any_of(mutcols)),
                          Pop),
                 freq_popavg = mean(freq0),
                 freq_popsd = sd(freq0))

unique(filter(mutpop, freq_popavg == 1, freq_popsd == 0)$mut_id)
#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
mutpop <- filter(
  mutpop,
  !mut_id %in% unique(filter(mutpop, freq_popavg == 1, freq_popsd == 0)$mut_id))

###Analyze gene category frequencies ----

##TODO - remove synonymous mutations

gcat_freq <- summarize(
  group_by(
    filter(
      mutpop, seq_id == "PA14_1", seqtype == "PA",
      gene_name != "type 4b pilus Flp biogenesis protein TadZ"), 
    gene_category,
    across(any_of(sampcols))),
  total_freq = sum(freq0),
  max_freq = max(freq0),
  indep_freq = 1-prod(1-freq0),
  avg_freq = mean(freq0))

View(summarize(group_by(gcat_freq, Treat, gene_category),
          avg_total = mean(total_freq)))

ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1")),
       aes(x = Treat_C, y = total_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(Pop, Timepoint))) +
  facet_grid(Treat_P ~ gene_category) +
  #scale_y_log10() +
  NULL

ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1")),
       aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(Pop, Timepoint))) +
  facet_grid(Treat_P ~ gene_category) +
  #scale_y_log10() +
  NULL

png("./plots/gene_freq_pop_diff_total_PNP.png", width = 20, height = 4,
    units = "in", res = 150)
ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1"),
                     Timepoint == 8),
       aes(x = Treat_P, y = total_freq)) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop))) +
  facet_nested(~ gene_category + Treat_C,
               labeller = labeller(
                 Treat_C = c("NC" = "Only\nPseudo.",
                             "C" = "All\nBacteria"))) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "Total frequency\nTimepoint 8\n(assuming exclus.)") +
  scale_x_discrete(breaks = c("NP", "P"),
                   labels = c("No Phage", "Phage Added")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13)) +
  NULL
dev.off()

png("./plots/gene_freq_pop_diff_total_CNC.png", width = 20, height = 4,
    units = "in", res = 150)
ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1"),
                     Timepoint != 0),
       aes(x = Treat_C, y = total_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop, Timepoint))) +
  facet_nested(~ gene_category + Treat_P,
               labeller = labeller(
                 Treat_P = c("NP" = "No Phage",
                             "P" = "Phage\nAdded"))) +
  #scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Timepoint\n(days)", values = my_cols[c(1:3, 7)]) +
  labs(x = "", y = "Total frequency\n(assuming exclus.)") +
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
dev.off()

png("./plots/gene_freq_pop_diff_indep_PNP.png", width = 20, height = 4,
    units = "in", res = 150)
ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1"),
                     Timepoint == 8),
       aes(x = Treat_P, y = indep_freq)) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop))) +
  facet_nested(~ gene_category + Treat_C,
             labeller = labeller(
               Treat_C = c("NC" = "Only\nPseudo.",
                           "C" = "All\nBacteria"))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "Total frequency\nTimepoint 8\n(assuming indep.)") +
  scale_x_discrete(breaks = c("NP", "P"),
                   labels = c("No Phage", "Phage Added")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 13)) +
  NULL
dev.off()

png("./plots/gene_freq_pop_diff_indep_CNC.png", width = 20, height = 4,
    units = "in", res = 150)
ggplot(data = filter(gcat_freq, !gene_category %in% c("other", "mystery_1"),
                     Timepoint != 0),
       aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop, Timepoint))) +
  facet_nested(~ gene_category + Treat_P,
               labeller = labeller(
                 Treat_P = c("NP" = "No Phage",
                             "P" = "Phage\nAdded"))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Timepoint\n(days)", values = my_cols[c(1:3, 7)]) +
  labs(x = "", y = "Total frequency\n(assuming indep.)") +
  scale_x_discrete(breaks = c("NC", "C"),
                   labels = c("Only Pseudo.", "All Bacteria")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  NULL
dev.off()

ggplot(data = filter(gcat_freq, gene_category == "pilus"),
       aes(x = Treat_C, y = total_freq, color = gene_category)) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop))) +
  facet_grid(Treat_P ~ Timepoint) +
  #scale_y_log10() +
  NULL

ggplot(data = filter(gcat_freq, gene_category == "pilus"),
       aes(x = Treat_C, y = indep_freq, color = Pop)) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop))) +
  facet_grid(Treat_P ~ Timepoint) +
  #scale_y_log10() +
  NULL

ggplot(data = filter(gcat_freq, gene_category == "pilus"),
       aes(x = Treat_C, y = max_freq, color = Pop)) +
  geom_point() +
  geom_line(aes(group = paste(gene_category, Pop))) +
  facet_grid(Treat_P ~ Timepoint) +
  #scale_y_log10() +
  NULL


totalfreq_full_lme <- 
  lmer(total_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
             data = filter(gcat_freq, gene_category == "pilus",
                           Treat_P == "P", Timepoint != 0),
       REML = FALSE)
totalfreq_reduced_lme <-
  lmer(total_freq ~ as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_category == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(totalfreq_reduced_lme, totalfreq_full_lme)
totalfreq_full_lme

indepfreq_full_lme <- 
  lmer(indep_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_category == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
indepfreq_reduced_lme <- 
  lmer(indep_freq ~ as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_category == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(indepfreq_reduced_lme, indepfreq_full_lme)
indepfreq_full_lme

maxfreq_full_lme <- 
  lmer(max_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_category == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
maxfreq_reduced_lme <- 
  lmer(max_freq ~ as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_freq, gene_category == "pilus",
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(maxfreq_reduced_lme, maxfreq_full_lme)
maxfreq_full_lme

gcat_type_freq_pil <- summarize(
  group_by(
    filter(
      mutpop, seq_id == "PA14_1", seqtype == "PA",
      gene_category == "pilus",
      gene_name != "type 4b pilus Flp biogenesis protein TadZ"), 
    gene_category, mutation_category3,
    across(any_of(sampcols))),
  total_freq = sum(freq0),
  indep_freq = 1-prod(1-freq0),
  avg_freq = mean(freq0))

ggplot(data = gcat_type_freq_pil,
       aes(x = Treat_C, y = indep_freq, color = mutation_category3)) +
  geom_point() +
  geom_line(aes(group = paste(mutation_category3, Pop))) +
  facet_grid(Treat_P * mutation_category3 ~ Timepoint) +
  #scale_y_log10() +
  NULL

totalfreq_pil_full_lme <- 
  lmer(total_freq ~ 
         Treat_C + mutation_category3 + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_type_freq_pil, 
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
totalfreq_pil_reduced_lme <-
  lmer(total_freq ~ mutation_category3 + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_type_freq_pil, 
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(totalfreq_pil_reduced_lme, totalfreq_pil_full_lme)
totalfreq_pil_full_lme

indepfreq_pil_full_lme <- 
  lmer(indep_freq ~ Treat_C + mutation_category3 + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_type_freq_pil, 
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
indepfreq_pil_reduced_lme <- 
  lmer(indep_freq ~ mutation_category3 + as.factor(Timepoint) + (1|Pop),
       data = filter(gcat_type_freq_pil, 
                     Treat_P == "P", Timepoint != 0),
       REML = FALSE)
anova(indepfreq_pil_reduced_lme, indepfreq_pil_full_lme)
indepfreq_pil_full_lme

###Analyze unique gene frequencies ----
gnamefreq <- summarize(
  group_by(
    filter(
      mutpop, seq_id == "PA14_1", seqtype == "PA"), 
    gene_name, gene_product, gene_category,
    across(any_of(sampcols))),
  total_freq = sum(freq0),
  max_freq = max(freq0),
  indep_freq = 1-prod(1-freq0),
  avg_freq = mean(freq0))

gnamefreq_subset <- 
  filter(group_by(gnamefreq, gene_name),
         any(total_freq > 0.1, indep_freq > 0.1))

#fill in missing combinations (doesn't actually fill in any combos bc
# they were all already filled in previously)
gnamefreq_subset <- 
  complete(ungroup(gnamefreq_subset),
           nesting(!!!syms(intersect(sampcols, colnames(gnamefreq_subset)))),
           gene_name,
           fill = list("total_freq" = 0, "max_freq" = 0,
                       "indep_freq" = 0, "avg_freq" = 0))

# #Demo examples
# testdat <- filter(gnamefreq_subset, Timepoint != 0, Treat_P == "P",
#                   gene_name == "EamA domain-containing protein")
# ggplot(data = testdat,
#        aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
#   geom_point() +
#   geom_line(aes(group = paste(Timepoint, Pop))) +
#   facet_grid(~ Treat_P)
# testlm_full <- lmer(indep_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
#                     data = testdat, REML = FALSE)
# testlm_red <- lmer(indep_freq ~ as.factor(Timepoint) + (1|Pop),
#                    data = testdat, REML = FALSE)
# anova(testlm_full, testlm_red)
# coef(summary(testlm_full))
# 
# testdat <- filter(gnamefreq_subset, Timepoint == 8,
#                   gene_name == "EamA domain-containing protein")
# testdat <- filter(gnamefreq_subset, Timepoint == 8,
#                   gene_name == "3-oxoadipate CoA-transferase")
# ggplot(data = testdat,
#        aes(x = Treat_C, y = indep_freq, color = as.factor(Timepoint))) +
#   geom_point() +
#   geom_line(aes(group = paste(Timepoint, Pop))) +
#   facet_grid(~ Treat_P)
# testlm <- lmer(
#   indep_freq ~ Treat_C + Treat_P + (1|Pop),
#   data = testdat,
#   REML = FALSE)
# testem <- emmeans(testlm, ~ Treat_C * Treat_P)
# emm_within_P <- as.data.frame(contrast(testem, method = "pairwise", by = "Treat_P"))
# emm_within_C <- as.data.frame(contrast(testem, method = "pairwise", by = "Treat_C"))

gnamefreq_subset_sum <-
  group_modify(
    group_by(filter(gnamefreq_subset, Timepoint != 0),
             gene_name, gene_category),
    ~ {
      mtc <- myTryCatch({
        out <- 
          data.frame(
            NCvC_inP_pval = NA,
            NCvC_inP_est = NA,
            NCvC_inNP_pval = NA,
            NCvC_inNP_est = NA,
            NPvP_inC_pval = NA,
            NPvP_inC_est = NA,
            NPvP_inNC_pval = NA,
            NPvP_inNC_est = NA)
        
        if(any(filter(.x, Treat_P == "P")$indep_freq > 0)) {
          #First test for significant effect of C vs NC within P treats
          # across all timepoints
          full_lme_P <- lmer(
            indep_freq ~ Treat_C + as.factor(Timepoint) + (1|Pop),
            data = filter(.x, Treat_P == "P"),
            REML = FALSE)
          red_lme_P <- lmer(
            indep_freq ~ as.factor(Timepoint) + (1|Pop),
            data = filter(.x, Treat_P == "P"),
            REML = FALSE)
          
          out$NCvC_inP_pval <- anova(red_lme_P, full_lme_P)$`Pr(>Chisq)`[2]
          out$NCvC_inP_est <- coef(summary(full_lme_P))["Treat_CC", "Estimate"]
        }
        
        if(any(filter(.x, Timepoint == 8)$indep_freq > 0)) {
          #Then do NC vs C, P vs NP contrasts for Timepoint 8
          full_lme_T8 <- lmer(
            indep_freq ~ Treat_C*Treat_P + (1|Pop),
            data = filter(.x, Timepoint == 8),
            REML = FALSE)
          full_emm_T8 <- emmeans(full_lme_T8, ~ Treat_C * Treat_P)
          emm_within_P <- 
            as.data.frame(contrast(full_emm_T8, method = "pairwise", by = "Treat_P"))
          emm_within_C <- 
            as.data.frame(contrast(full_emm_T8, method = "pairwise", by = "Treat_C"))
          
          out$NCvC_inNP_pval <- emm_within_P[emm_within_P$Treat_P == "NP", "p.value"]
          out$NCvC_inNP_est <- emm_within_P[emm_within_P$Treat_P == "NP", "estimate"]
          out$NPvP_inC_pval <- emm_within_C[emm_within_C$Treat_C == "C", "p.value"]
          out$NPvP_inC_est <- emm_within_C[emm_within_C$Treat_C == "C", "estimate"]
          out$NPvP_inNC_pval <- emm_within_C[emm_within_C$Treat_C == "NC", "p.value"]
          out$NPvP_inNC_est <- emm_within_C[emm_within_C$Treat_C == "NC", "estimate"]  
        }
        
        out
      })

      if(is.null(mtc$error)) {
        return(
          cbind(mtc$value,
                "warn" = ifelse(is.null(mtc$warning), NA, as.character(mtc$warning)),
                "error" = NA))
      } else {
        return(
          cbind(data.frame(
            NCvC_inP_pval = NA,
            NCvC_inP_est = NA,
            NCvC_inNP_pval = NA,
            NCvC_inNP_est = NA,
            NPvP_inC_pval = NA,
            NPvP_inC_est = NA,
            NPvP_inNC_pval = NA,
            NPvP_inNC_est = NA),
            "warn" = ifelse(is.null(mtc$warning), NA, as.character(mtc$warning)),
            "error" = as.character(mtc$warning)))
      }
    })

gnfs_sum_lng <- 
  pivot_longer(
    filter(gnamefreq_subset_sum,
           if_any(ends_with("pval"), ~ !is.na(.) & . < 0.01)),
    cols = ends_with(c("pval", "est")),
    names_to = c("comparison", "in_treat", ".value"),
    names_pattern = "(.*)_in(.*)_(.*)")

png("./plots/genefreq_byname.png", width = 15, height = 20, units = "in",
    res = 150)
ggplot(data = filter(gnfs_sum_lng, pval < 0.01),
       aes(y = gene_name, x = -est, size = -log10(pval),
           color = comparison, shape = in_treat)) +
  geom_point() +
  facet_grid(gene_category ~ ., space = "free_y", drop = TRUE, scales = "free_y") +
  scale_size_continuous(limits = c(2, NA), range = c(0.5, 4)) +
  scale_shape_manual(values = 15:18) +
  theme(strip.text.y = element_text(angle = 0))
dev.off()

### GO enrichment analysis ----

#GO enrichment analysis using cluterProfiler in R
#org.Ce.eg.db but need something like that for PA14
#Can do a GSEA analysis where you put in a ranked list of
# genes and it does kind of a rolling window approach and
# stops at some point for sig GO terms
#Or can just do a plain GO analysis that looks for enriched
# terms in a hard set of genes
#remove synon
#I think first plan of attack is take my post-stats gene names
# and put them into a simple GO analysis
#Second plan of attack (if needed) is to rank genes by frequency
# change per replicate, then feed that into GSEA
# (because by binning by gene name before GO I'll lose any
# pathways where there are multiple hits but to different
# genes in the pathway)

hub <- AnnotationHub()
q <- query(hub, "Pseudomonas aeruginosa")
id <- q$ah_id[length(q)]

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

# mutpop <- mutate(
#   group_by(mutpop, mut_id, seqtype, Treat, Pop),
#   freq_change = (freq01 - lag(freq01, order_by = Timepoint)),
#   abs_freq_rate =
#     (freq01 - lag(freq01, order_by = Timepoint))/
#     (Timepoint - lag(Timepoint, order_by = Timepoint)),
#   rel_freq_rate = abs_freq_rate/lag(freq01, order_by = Timepoint))

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

if(make_slowplots) {
  png("./plots/mut_trajectories_mystery1.png", width = 10, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", gene_category == "mystery_1",
                    seq_id %in% c("PA14_1")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = as.factor(mut_id))) +
      geom_line() +
      facet_grid(Treat ~ Pop) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      theme_bw() +
      guides(color = "none") +
      NULL)
  dev.off()
}

#Top hits that I haven't yet understood the function of:
#Thioredoxin family protein/Helix-turn-helix domain-containing protein 
# intergenic SNPs located between 4629614 and 4629639 are the top hit
#PA2778 family cysteine peptidase located around 6494573
#Amino-acid carrier protein AlsT/L-asparaginase I
#DUF4340 domain-containing protein located around 1387230
#Exo-alpha-sialidase

mut_sum <- summarize(
  group_by(filter(mutpop, seqtype == "PA", seq_id == "PA14_1"),
           across(any_of(c(mutcols, "Treat", "Treat_C", "Treat_P")))),
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
         any_of(c(
           mutcols, 
           setdiff(sampcols, c("seqrun", "SeqCoastid", "title")),
           "lfreq01"))),
  names_from = Timepoint, values_from = lfreq01, names_prefix = "T")

View(arrange(select(mutpop_pil_w, T0, T1, T2, T4, T8, Treat, Pop, mut_id, gene_name,
            everything()),
            Pop, gene_name, mut_id, Treat))

mutpop_pil_comp <-
  summarize(group_by(mutpop_pil,
                     across(!any_of(
                       c("Treat", "Treat_C", "Treat_P",
                         "seqrun", "SeqCoastid", "title",
                         mut_and_samp_cols)))),
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

if(make_slowplots) {
  png("./plots/pca_pilusmuts.png", width = 10, height = 3, units = "in",
      res = 150)
  print(
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
      NULL)
  dev.off()
  
  png("./plots/pca_allmuts.png", width = 10, height = 3, units = "in",
      res = 150)
  print(
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
      NULL)
  dev.off()
  
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
    geom_line(lwd = 1.1) + #, aes(alpha = gene_category)) +
    facet_grid(Treat ~ Pop) +
    #coord_cartesian(ylim = c(0.1, 1)) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
    # scale_color_manual(breaks = c("biofilm", "flagella", "LasA/LasB", "QS",
    #                               "T6SS", "pilus", "TPS system", "other"),
    #                    values = my_cols) +
    # scale_alpha_manual(breaks = c("biofilm", "flagella", "LasA/LasB", "QS",
    #                               "T6SS", "pilus", "TPS system", "other"),
    #                   values = c(rep(0.6, 7), 0.1)) +
    theme_bw() +
    NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_above010_noanc_base.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = mutate(
        filter(mutpop, 
               seqtype == "PA", everabove_010_noanc == TRUE, 
               seq_id %in% c("PA14_1")),
        Pop = as.numeric(Pop)),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id)) +
      geom_line(lwd = 0.7, alpha = 0.2) +
      facet_nested(Treat_P + Treat_C ~ Pop,
                   labeller = labeller(
                     Treat_P = c("NP" = "No Phage",
                                 "P" = "Phage\nAdded"),
                     Treat_C = c("NC" = "Only\nPseudomonas",
                                 "C" = "All Bacteria"),
                     Pop = function(x) {paste("Pop.", LETTERS[as.numeric(x)-2])})) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      theme(axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_blank(),
            strip.text.y = element_text(size = 16),
            strip.text.x = element_text(size = 14),
            plot.subtitle = element_text(size = 16)) +
      NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_above010_noanc_pilus.png", width = 13.5, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = mutate(
        filter(mutpop, 
               seqtype == "PA", everabove_010_noanc == TRUE, 
               seq_id %in% c("PA14_1")),
        Pop = as.numeric(Pop),
        pilgeneTF = gene_category == "pilus"),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id,
          color = pilgeneTF, alpha = pilgeneTF, lwd = pilgeneTF)) +
      geom_line() +
      facet_nested(Treat_P + Treat_C ~ Pop,
                   labeller = labeller(
                     Treat_P = c("NP" = "No Phage",
                                 "P" = "Phage\nAdded"),
                     Treat_C = c("NC" = "Only\nPseudomonas",
                                 "C" = "All Bacteria"),
                     Pop = function(x) {paste("Pop.", LETTERS[as.numeric(x)-2])})) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(name = "Gene for",
                         limits = c(TRUE, FALSE),
                         labels = c("Pilus", "Other"),
                         values = my_cols[c(1, 8)]) +
      scale_alpha_manual(limits = c(TRUE, FALSE),
                         values = c(0.8, 0.2)) +
      scale_linewidth_manual(breaks = c(TRUE, FALSE),
                             values = c(1.3, 0.7)) +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      theme(axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_blank(),
            strip.text.y = element_text(size = 16),
            strip.text.x = element_text(size = 14),
            plot.subtitle = element_text(size = 16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16)) +
      guides(alpha = "none", linewidth = "none") +
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
  
  png("./plots/mut_freqs_pilus.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seqtype == "PA", seq_id == "PA14_1",
                    everabove_010_noanc == TRUE, gene_category == "pilus",
                    gene_name != "type 4b pilus Flp biogenesis protein TadZ"),
      aes(x = gene_name, y = freq0, color = Pop, group = Pop,
          shape = mutation_category3)) +
      geom_point(alpha = 0.5, position = position_dodge(width = 0.5)) +
      facet_grid(Treat ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      NULL
  )
  dev.off()
  
  tmp <- filter(mutpop, seqtype == "PA", seq_id == "PA14_1",
         everabove_010_noanc == TRUE, gene_category == "pilus",
         gene_name != "type 4b pilus Flp biogenesis protein TadZ")
  View(filter(tmp, is.na(mutation_category3)))

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
  
  png("./plots/mut_trajectories_H6.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seq_id == "H6_1", 
                    !(mutation_category == "large_deletion" & size == "57534")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id)) +
      #geom_point() +
      geom_line() +
      facet_grid(seqtype*Treat ~ Pop) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      theme_bw() +
      NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_SA.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, str_starts(seq_id, "SA_") &
                    !(mutation_category == "large_deletion" & size == "57534")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id)) +
      #geom_point() +
      geom_line() +
      facet_grid(seqtype*Treat ~ Pop) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      theme_bw() +
      NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_AX.png", width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = filter(mutpop, seq_id == "H6_1", 
                    !(mutation_category == "large_deletion" & size == "57534")),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id)) +
      #geom_point() +
      geom_line() +
      facet_grid(seqtype*Treat ~ Pop) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      theme_bw() +
      NULL
  )
  dev.off()
}


# temp <- filter(mutpop, fixed == TRUE)
# View(temp[order(temp$mut_id, temp$Timepoint, temp$Treat, as.numeric(temp$Pop)),
#           c("freq", "mut_id", "Timepoint", "Treat", "Pop")])
# 
# temp <- filter(mutpop, seqtype == "PA", Isol == "Pop", Treat == "NCP",
#                Pop == 3, everabove_025 == TRUE, seq_id == "PA14_1")
# View(temp[order(temp$mut_id, temp$Timepoint), 
#           c("freq", "mut_id", "Timepoint")])

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
    select(filter(mutpop,
                  seqtype == "PA", seq_id == "PA14_1"),
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

tmp <- filter(mutcomp, seqtype == "PA", seq_id == "PA14_1",
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


tmp <- filter(mutcomp, seqtype == "PA", seq_id == "PA14_1",
              everabove_010_orisolmut == TRUE, 
              Treat %in% c("NCP", "CP"))

png("./plots/isolmut_vs_popmut.png", width = 16, height = 8, units = "in",
    res = 150)
ggplot(tmp,
       aes(x = Timepoint, y = freq0, color = gene_category,
           group = mut_id)) +
  geom_line() +
  geom_point(alpha = 0.5, aes(size = isolfreq)) +
  facet_grid(Treat ~ Pop) +
  scale_size_continuous(breaks = c(0, 0.5, 1),
                        range = c(0, 4))
dev.off()

# View(
#   select(
#     filter(mutisol,
#            !mut_id %in% c("2958", "3198", "7885", "11212")),
#     Timepoint, Treat, Pop, Isol, Treat_C, Treat_P,
#     clearing_avg, EOP_avg, Plaque_morphology, resis_cat, 
#     twitch_cat,
#     gene_category, concise_mut, freq0))

# ggplot(data = filter(mutisol,
#                      !mut_id %in% c("2958", "3198", "7885", "11212")),
#        aes(x = freq0)) +
#   geom_histogram()


# ggplot(data = mutisol2_sum,
#        aes(x = isolavg_freq, y = freq0, color = gene_name)) +
#   geom_point(alpha = 0.5) +
#   facet_grid(Treat ~ Pop) +
#   theme() +
#   guides(color = "none")



##Isolate analyses ----

### Isolate samples merge with phenotype data ----
#Filter out sample that had 0 coverage (so mutation is a full-genome deletion)
mutisol <- filter(mutisol, title != "5174_67_PA_T8_P4_NCP_#5")

mutisol <- left_join(
  mutisol,
  select(
    mutate(
      filter(read.csv("aggregate_isolate_data.csv", stringsAsFactors = FALSE,
                      colClasses = "character"),
             Bact_timepoint != "pilA"),
      Treat = Bact_treat, 
      Treat_C = Bact_treat_C, 
      Treat_P = Bact_treat_P,
      Timepoint = as.numeric(Bact_timepoint), 
      Pop = Bact_pop, Isol = Bact_isol),
    any_of(sampcols))
)


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

#Create plottable frequency of mutations
isol_freq_sum <- 
  summarize(
    group_by(
      mutate(mutisol,
             gene_category = case_when(
               gene_category == "pilus" ~ "pilus",
               gene_category == "biofilm" ~ "biofilm",
               TRUE ~ "other")),
      gene_category, across(any_of(sampcols))),
    mut_genes = {
      mygenes <- gene_short[!mut_id %in% c("2958", "3198", "7885", "11212")]
      case_when(
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
    fill = list("mut_genes" = "none")
  )

#Manually recode large deletion to "multiple"
isol_freq_sum <- 
  mutate(isol_freq_sum,
         mut_genes = ifelse(
           mut_genes == "[twitching motility protein PilU]–[type IV pilus twitching motility protein PilT]",
         "multiple",
         mut_genes))

pil_order <- unique(isol_freq_sum[isol_freq_sum$gene_category == "pilus", ]$mut_genes)
pil_order <-
  c("none", "multiple", pil_order[!pil_order %in% c("none", "multiple")])

png("./plots/pilus_mut_freq_isols.png",
    width = 9, height = 5, units = "in", res = 150)
ggplot(
  mutate(
    filter(isol_freq_sum, gene_category == "pilus"),
    mut_genes = factor(mut_genes, levels = pil_order)),
  aes(x = resis_cat)) +
  geom_bar(aes(fill = mut_genes), position = "fill") +
  scale_x_discrete(limits = c("suscep", "part", "resis"),
                   labels = c("Susceptible", "Partial", "Resistant")) +
  scale_fill_manual(
    name = "Pilus Genes\nMutated",
    limits = pil_order,
    values = c("gray70", "black",
               viridis_pal()(length(pil_order)-2))) +
  labs(x = "Isolates' Resistance to Phage",
       y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18)) +
  geom_text(data = summarize(group_by(isol_freq_sum, resis_cat), 
                             n = length(unique(title))),
            aes(y = 0.1, label = paste0("n = ", n)), size = 5) +
  NULL
dev.off()


png("./plots/all_mut_freq_isols.png",
    width = 14, height = 5, units = "in", res = 150)
ggplot(isol_freq_sum,
  aes(x = resis_cat)) +
  geom_bar(aes(fill = mut_genes), position = "fill") +
  facet_grid(~ gene_category) +
  scale_x_discrete(limits = c("suscep", "part", "resis"),
                   labels = c("Susceptible", "Partial", "Resistant")) +
  # scale_fill_manual(
  #   name = "Pilus Genes\nMutated",
  #   limits = pil_order,
  #   values = c("gray70", "black",
  #              viridis_pal()(length(pil_order)-2))) +
  labs(x = "Isolates' Resistance to Phage",
       y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18)) +
  geom_text(data = summarize(group_by(isol_freq_sum, resis_cat), 
                             n = length(unique(title))),
            aes(y = 0.1, label = paste0("n = ", n)), size = 5) +
  guides(fill = "none") +
  NULL
dev.off()

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

png("./plots/isol_mut_types3.png", width = 8, height = 6, units = "in",
    res = 150)
ggplot(data = filter(mutisol),
       aes(x = resis_cat, fill = mutation_category3)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(gene_category ~ Treat)
dev.off()

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

png("./plots/isol_pilmuts_byresis.png", width = 8, height = 10, units = "in",
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
    aes(x = gene_name, y = twitch_cat, color = Pop, fill = Pop,
        group = Pop, shape = mutation_category3)) +
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


##Old code ----
if(FALSE) {
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
}
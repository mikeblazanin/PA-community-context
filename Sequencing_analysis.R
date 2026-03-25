##TODO ----
#Top priority:
# fix gene cat analy to remove synonymous mutations
# Update gene cat stats using old & new pilus category
# Correlate phage freq vs freq of old & new pilus category
# fix gene_short for intergenic/multihit cases

#Maybe:
# create gene_category2 which ID's type IV pilus genes more-specifically
# adjust coverage analysis for repetitive regions
# Get total # reads (not just aligned reads) & use as denom for coverage analy
# run isols in polym mode, try to figure out contamination (where are bg muts
#  (in isols!?!)
# why are cov-based densities of other bact higher in NC treats?
# Look for dift fates of phages depending on initial inocula 3-10 & pilus muts
# Is there anything in common among the populations where phage went ~extinct?
#  Some way those 5 pops cluster?

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

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Extended version of Carto "Safe" palette
my_cols2 <- c("#888888", "#117733", "#332288", "#661100","#882255",
              "#88CCEE", "#CC6677", "#DDCC77",  
            "#AA4499", "#44AA99","#999933",  "#6699CC",
             "#E69F00")

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
  c(
    #columns from breseq:
    "aa_new_seq", "aa_position", "aa_ref_seq", "codon_new_seq", "codon_number", 
    "codon_position", "codon_ref_seq", "gene_name", "gene_position", 
    "gene_product", "gene_strand", "genes_inactivated", "genes_overlapping", 
    "genes_promoter", "insert_position", "locus_tag", "locus_tags_inactivated",
    "locus_tags_overlapping", "locus_tags_promoter", "mutation_category",
    "new_seq", "position", "position_end", "position_start",
    "ref_seq", "repeat_length", "repeat_new_copies", "repeat_ref_copies",
    "repeat_seq", "seq_id", "size", "snp_type", "transl_table", "type", 
    #columns from my analyses:
    "mut_id", "gene_category", "mutation_category2", "mutation_category3",
    "pilmut_struc", "gene_short", "concise_mut", "mut_combin",
    #Columns from eggnog and bakta:
    "Sequence.Id", "Type", "Start", "Stop", "Strand", "Locus.Tag", "Gene", 
    "Product", "DbXrefs", "seed_ortholog", "evalue", "score", "eggNOG_OGs", 
    "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", 
    "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", 
    "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
    #Columns from InterPro
    "seqMD5", "seq_len", "analysis", "signat_accession", "signat_description",
    "start", "stop", "e-value_score", "match_status", "query_date",
    "InterPro_accession", "InterPro_annotations", "InterPro_GOs", 
    "Pathways", "InterPro_GOs_clean",
    #Columns from anno merged df
    "GOsmerge",
    #Columns from BLAST for repeat regions
    "rpt_hit_id", "rpt_bit_score", "rpt_match_seq")
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

#Drop non-PA results
nonPAmuts <- filter(mutations, seqtype != "PA" | seq_id != "PA14_1")
mutations <- filter(mutations, seqtype == "PA", seq_id == "PA14_1")

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

##Incorporate annotations ----
bakta <- read.csv("PA14_annotated.csv")
bakta$Locus.Tag[bakta$Locus.Tag == ""] <-
  paste0("PEANK_MB", 
        formatC(1:sum(bakta$Locus.Tag == ""),
                width = floor(log10(sum(bakta$Locus.Tag == "")))+1,
                flag = "0"))

eggnog <- read.delim("PA14_eggnog_out3.emapper.annotations",
                     comment = "#",
                     header = FALSE)
colnames(eggnog) <-
  gsub("#", "", 
       read.delim("PA14_eggnog_out2.emapper.annotations",
                  comment = "", header = FALSE,
                  skip = 4, nrows = 1))
colnames(eggnog)[colnames(eggnog) == "query"] <- "Locus.Tag"

interpro <- read.table("PA14_interpro_out.tsv", sep = "\t",
                       stringsAsFactors = FALSE, header = FALSE,
                       quote = "", #don't allow quotes
                       comment.char = "") #don't allow comments
colnames(interpro) <- 
  c("Locus.Tag", "seqMD5", "seq_len", "analysis", "signat_accession",
    "signat_description", "start", "stop", "e-value_score",
    "match_status", "query_date", "InterPro_accession",
    "InterPro_annotations", "InterPro_GOs", "Pathways")
interpro_sum <- 
  summarize(group_by(interpro, Locus.Tag),
            across(everything(),
                   ~ paste(., collapse = ";")))

anno = full_join(bakta, eggnog)
anno = full_join(anno, interpro_sum)

anno <- mutate(ungroup(anno),
               GOsmerge = paste0(GOs, "|", InterPro_GOs),
               GOsmerge = gsub("NA", "-", GOsmerge),
               GOsmerge = gsub(";|,", "|", GOsmerge),
               GOsmerge = gsub("(-\\|)*|(\\|-)*$", "", GOsmerge),
               GOsmerge = gsub("\\([^\\|]*\\)", "", GOsmerge))
anno$GOsmerge <-
  sapply(X = anno$GOsmerge,
         FUN = function(x) {paste0(unique(strsplit(x, split = "\\|")[[1]]),
                                   collapse = "|")})

#Make term to gene mapping from background annotation
term2gene <- 
  select(
    separate_longer_delim(filter(ungroup(anno), !is.na(GOsmerge), GOsmerge != "-"),
                          cols = GOsmerge, delim = "|"),
    GOsmerge, Locus.Tag)
#Get indirect annotations
term2gene <- buildGOmap(term2gene)

#Add indirect annotations
anno <- full_join(anno,
                  (summarize(group_by(term2gene, Locus.Tag),
                             GOsmerge_ind = paste(GOsmerge, collapse = "|"))))

#Pull GO term names for interpretation
term2name <- data.frame(
  GOsmerge_ind = AnnotationDbi::keys(GO.db::GO.db), 
  name = AnnotationDbi::Term(AnnotationDbi::keys(GO.db::GO.db)))


#Manually determine Locus.Tag values because locus_tag from breseq
# is not a match
mutations$Locus.Tag <- "unmatched"
mutations$mutation_category_manual <- NA
mutations$Locus.Tag[mutations$size == 6541479] <- NA #no coverage sample
for (i in which(!is.na(mutations$Locus.Tag))) {
  matches <- anno[(anno$Start <= mutations$position_start[i] & 
                     anno$Stop >= mutations$position_start[i]) |
                    (anno$Start <= mutations$position_end[i] &
                       anno$Stop >= mutations$position_end[i]), ]
  if(nrow(matches) == 0) {
    mutations$mutation_category_manual[i] <- "intergenic"
    promgenes <- c()
    #preceding gene (in promoter if on the + strand)
    before <- max(which(anno$Stop < mutations$position_start[i]))
    if(anno$Strand[before] == "-" & 
       (mutations$position_start[i] - anno$Stop[before]) <= 150) {
      promgenes <- before
    } else {promgenes <- c(NA)}
    #following gene (in promoter if on the - strand)
    after <- min(which(anno$Start > mutations$position_end[i]))
    if(anno$Strand[after] == "+" &
       (anno$Start[after] - mutations$position_end[i]) <= 150) {
      promgenes <- c(promgenes, after)
    } else {promgenes <- c(promgenes, NA)}
    if(any(!is.na(promgenes))) {promgenes <- anno$Locus.Tag[promgenes]}
    mutations$Locus.Tag[i] <- paste(promgenes, collapse = "|")
    #mutations$Locus.Tag[i] <- paste(promgenes[promgenes != ""], collapse = "|")
  } else if(nrow(matches) == 1) {
    mutations$mutation_category_manual[i] <- "genic"
    mutations$Locus.Tag[i] <- matches$Locus.Tag[1]
  } else { #nrow(matches) > 1
    #myrows <- which(!is.na(matches$Locus.Tag) & matches$Locus.Tag != "")
    #If only one match has a locus tag, use that one
    # if(length(myrows) == 1) {
    #   mutations$Locus.Tag[i] <- matches$Locus.Tag
    # } else {
    mutations$mutation_category_manual[i] <- "multigene"
    mutations$Locus.Tag[i] <- paste(matches$Locus.Tag, collapse = "|")
    #}
  } 
}

#Checks
# View(mutations[grepl("\\|", mutations$mutation_category) & 
#                  !grepl("\\|", mutations$Locus.Tag), ])
# View(filter(mutations, Locus.Tag == "unmatched"))
# View(select(filter(mutations, grepl("\\|", Locus.Tag)), 
#             gene_name, gene_position, gene_product, gene_strand, genes_promoter, 
#             mutation_category, position, position_end, position_start, snp_type, 
#             type, mut_id, mutation_category2, mutation_category3, Locus.Tag))
# View(select(filter(mutations, mutation_category3 == "intergenic"), 
#             gene_name, gene_position, gene_product, gene_strand, genes_promoter, 
#             mutation_category, position, position_end, position_start, snp_type, 
#             type, mut_id, mutation_category2, mutation_category3, Locus.Tag))
# View(select(filter(mutations, mutation_category3 == "intergenic", 
#                    (Locus.Tag != "" & genes_promoter == "")), 
#             gene_name, gene_position, gene_product, gene_strand, genes_promoter, 
#             mutation_category, position, position_end, position_start, snp_type, 
#             type, mut_id, mutation_category2, mutation_category3, Locus.Tag))

# tmp <- summarize(group_by(mutations, across(any_of(mutcols))))
# tmp <- mutate(ungroup(tmp),
#               nsep = rowSums(across(everything(),
#                                     ~ str_detect(.x, "[|/]")),
#                              na.rm = TRUE),
#               ncomma = rowSums(across(everything(),
#                                       ~ str_detect(.x, ",")),
#                                na.rm = TRUE))
# tmp <- tmp[order(tmp$nsep, tmp$ncomma, decreasing = TRUE), ]
# write.csv(tmp, "tmp.csv", row.names = FALSE)

### Split multihit rows ----
###Notes:

#Across all mutation types:
# gene_name, genes_inactivated, genes_overlapping, genes_promoter
#   use _ in gene names where originally contained / or ,
# gene_product
#   uses the original gene names with / or , in them
mutationsd <- list(
  filter(mutations, mutation_category_manual == "multigene"),
  filter(mutations, mutation_category_manual == "intergenic"),
  filter(mutations, 
         #grepl("\\|", Locus.Tag) & 
           !mutation_category_manual %in% c("multigene", "intergenic"))
  # filter(mutations, 
  #        !grepl("\\|", Locus.Tag) & 
  #          !mutation_category_manual %in% c("multigene", "intergenic"))
)
if(sum(sapply(FUN = nrow, X = mutationsd)) != nrow(mutations)) {
  stop("missing some cases in mutations splitting")}

##multigene hits (bc of overlapping genes) have | in 16 columns
# they use , in genes_overlapping
# also will have to re-make concise_mut
# they may also have / or , in names in gene_product (as all types do)
if(any(str_count(mutationsd[[1]]$genes_overlapping, fixed(",")) > 1)) {
  stop("Too many commas in some multigene$genes_overlapping")}
mutationsd[[1]]$genes_overlapping <- gsub(",", "|", mutationsd[[1]]$genes_overlapping)
mutationsd[[1]] <- separate_longer_delim(
  mutationsd[[1]],
  cols = all_of(c("aa_new_seq", "aa_position", "aa_ref_seq", "codon_new_seq",
                  "codon_number", "codon_position", "codon_ref_seq", "gene_name",
                  "gene_position", "gene_product", "gene_strand", "locus_tag",
                  "mutation_category", "snp_type", "transl_table", #"mutation_category2",
                  "Locus.Tag",
                  "genes_overlapping")),
  delim = "|")


##intergenic hits have / in 4 columns, but gene_product also uses it
# within gene names so can't split that one easily
# of course use | in 1 column: "Locus.Tag"
# use , in the genes_promoter column
# they also have gene_position with a 'intergenic (pos/pos)' setup
# Note also that gene_name and gene_product for intergenic hits contain
#  both neighboring genes, even if it's only in the promoter of one/none
# they may also have / or , in names in gene_product (as all types do)
mutationsd[[2]]$genes_promoter <- gsub(",", "|", mutationsd[[2]]$genes_promoter)
mutationsd[[2]]$gene_position <- gsub("/", ")/intergenic (", mutationsd[[2]]$gene_position)
mutationsd[[2]] <- separate_longer_delim(
  mutationsd[[2]],
  cols = all_of(c("gene_name", "gene_strand", "locus_tag",
                  "Locus.Tag",
                  "genes_promoter", "gene_position")),
  delim = stringr::regex("[\\|/]"))

##Additional cases we basically only want to split Locus.Tag
# this includes coding indels, but they use / in gene_position denoting start/end
# rather than gene1/gene2
# use , in genes_overlapping (if span multiple genes)
# they may also have / or , in names in gene_product (as all types do)
mutationsd[[3]] <- separate_longer_delim(
  mutationsd[[3]], cols = all_of(c("Locus.Tag")), delim = "|")

##Merge
mutationsd <- do.call(rbind, mutationsd)

#Merge
mutationsd <- 
  left_join(
    mutate(mutationsd,
           Locus.Tag = ifelse(Locus.Tag == "NA" | Locus.Tag == "unmatched",
                              NA, Locus.Tag)), 
    anno,
    na_matches = "never")

#Check for name mismatches
mutationsd_sum <- summarize(group_by(mutationsd, across(any_of(mutcols))))
mutationsd_sum <- mutate(group_by(mutationsd_sum, across(any_of(mutcols))),
                         sd_gn = stringdist(gene_name, Product),
                         sd_gpd = stringdist(gene_product, Product),
                         sd_go = stringdist(genes_overlapping, Product),
                         sd_gpm = stringdist(genes_promoter, Product),
                         nas = sum(is.na(sd_gn), is.na(sd_gpd), is.na(sd_go), is.na(sd_gpm)),
                         sdmin = min(sd_gn, sd_gpd, sd_go, sd_gpm))
# any(!mutationsd_sum$nas %in% c(0, 4))
# View(mutationsd_sum[order(mutationsd_sum$sdmin, decreasing = TRUE),
#                     c("gene_name", "gene_product", "genes_overlapping",
#                       "genes_promoter", "Gene", "Product",
#                       "sd_gn", "sd_gpd", "sd_go", "sd_gpm", "sdmin")])

#Overwrite original mutations frame with re-joined version
mutations <-
  summarize(
    group_by(mutationsd, across(any_of(c("mut_id", sampcols)))),
    across(.cols = -any_of(sampcols),
           .fns = function(x) {
             if(length(unique(x)) == 1) {as.character(x[1])
             } else {paste(x, collapse = "|")
             }}))

##Shared analysis labeling ----
catmutate <- function(df) {
  return(
    mutate(
      ungroup(df),
      mutation_category2 = case_when(
        mutation_category == "small_indel" & grepl("intergenic", gene_position) ~
          "small_indel_intergenic",
        mutation_category == "small_indel" & grepl("coding", gene_position) &
          (nchar(ref_seq) - nchar(new_seq)) %% 3 == 0 ~ "small_indel_inframe",
        mutation_category == "small_indel" & grepl("coding", gene_position) &
          (nchar(ref_seq) - nchar(new_seq)) %% 3 != 0 ~ "small_indel_frameshift",
        grepl("intergenic", gene_position) ~ paste0(mutation_category, "_intergenic"),
        TRUE ~ mutation_category),
      mutation_category3 = case_when(
        mutation_category2 %in% 
          c("large_deletion", "small_indel_frameshift", "snp_nonsense") ~ "loss of function",
        grepl("intergenic", mutation_category2) ~ "intergenic",
        mutation_category2 %in% c("small_indel_inframe", "snp_nonsynonymous") ~ "nonsyn change",
        mutation_category2 %in% c("snp_synonymous") ~ "silent",
        grepl("\\|", mutation_category2) ~ "multigene",
        TRUE ~ "other"),
      gene_category = case_when(
        grepl("CupB3", Product) ~ "pilus",
        grepl("Wsp", Product) ~ "biofilm",
        grepl("MorA", Product) ~ "biofilm",
        grepl("RhlR", Product) ~ "QS",
        grepl("LasR", Product) ~ "QS",
        grepl("LasI", Product) ~ "QS",
        grepl("quorum", Product) ~ "QS",
        grepl("LasA", Product) ~ "LasA_LasB",
        grepl("LasB", Product) ~ "LasA_LasB",
        grepl("Mut", Product) ~ "mutator",
        grepl("Tps", Product) ~ "TPS system",
        grepl("type VI", Product) ~ "T6SS",
        grepl("Type VI", Product) ~ "T6SS",
        grepl("Tse", Product) ~ "T6SS",
        grepl("T2SS", Product) ~ "T2SS",
        grepl("type II secretion system", Product) ~ "T2SS",
        grepl("flagella", Product) ~ "flagella",
        grepl("Thioredoxin family protein/Helix-turn-helix domain-containing protein", Product) ~ "mystery_1",
        grepl("FimX", Product) ~ "pilus",
        grepl("pil", Product, ignore.case = TRUE) ~ "pilus",
        grepl("pilus", Product) ~ "pilus",
        grepl("CRISPR", Product) ~ "CRISPR",
        grepl(paste0(term2name$GOsmerge_ind[grep("pilus|Pilus|pili", term2name$name)],
                     collapse = "|"),
              GOsmerge_ind) ~ "pilus",
        TRUE ~ "other"),
      gene_catpilnot = ifelse(gene_category == "pilus", "pilus", "other"),
      pilmut_struc = ifelse(
        gene_category == "pilus",
        case_when(
          grepl("4b", Product) ~ "4b",
          grepl("IVB", Product) ~ "4b",
          grepl("TadZ", Product) ~ "unknown",
          grepl("PilT", Product) ~ "motor",
          grepl("PilB", Product) ~ "motor",
          grepl("PilQ", Product) ~ "OMP",
          grepl("PilS", Product) ~ "regulatory",
          grepl("PilM", Product) ~ "alignment_sc",
          grepl("PilC", Product) ~ "motor",
          grepl("pilC", Product) ~ "motor",
          grepl("PilN", Product) ~ "alignment_sc",
          grepl("PilZ", Product) ~ "motor",
          grepl("TadG", Product) ~ "unknown",
          grepl("PilU", Product) ~ "motor",
          grepl("PilR", Product) ~ "regulatory",
          grepl("PilO", Product) ~ "alignment_sc",
          grepl("FimV", Product) ~ "OMP",
          grepl("PilJ", Product) ~ "regulatory",
          grepl("PilF", Product) ~ "OMP",
          grepl("PilY1", Product) ~ "fiber",
          grepl("CupB3", Product) ~ "usher_chaperone_assembly",
          grepl("PilV", Product) ~ "fiber",
          grepl("PilG", Product) ~ "regulatory",
          grepl("PilK", Product) ~ "tbd",
          grepl("PilP", Product) ~ "alignment_sc",
          grepl("PilD", Product) ~ "usher_chaperone_assembly",
          grepl("PilE", Product) ~ "fiber",
          grepl("PilX", Product) ~ "fiber",
          grepl("FleP", Product) ~ "regulatory",
          grepl("PilM", Product) ~ "alignment",
          grepl("PilR", Product) ~ "regulatory",
          grepl("PilP", Product) ~ "alignment_sc",
          grepl("FimL", Product) ~ "regulatory",
          TRUE ~ "tbd"),
        NA))
  )
}
#TODO
# need to fix gene_short for intergenic/multihit cases bc it's just
#  pulling the last bit instead of each last bit
#   gene_short = coalesce(
#     # 1) Capitalized short symbol at end, e.g. "… PilT" or "… PilY1"
#     str_extract(gene_name, "([A-Z][A-Za-z0-9]+)$"),
#     # 2) Gene token immediately following "protein"
#     str_extract(gene_name, "(?i)\\bPil[A-Za-z0-9]*(?=$|/)"),
#     # 3) As a final fallback use the provided gene_name column
#     gene_name),
#   concise_mut = case_when(
#     mut_id %in% c("2958", "3198", "7885", "11212") ~ "", #(ancestral mutations)
#     grepl("snp", mutation_category) ~
#       ifelse(grepl("intergenic", gene_position),
#              paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
#              paste0(gene_short, " (", aa_ref_seq, "->", aa_new_seq, ")")),
#     mutation_category == "small_indel" & type == "DEL" ~
#       ifelse(grepl("intergenic", gene_position),
#              paste0(gene_short, " (intergenic ", size, "bp DEL)"),
#              paste0(gene_short, " (", size, "bp DEL)")),
#     mutation_category == "small_indel" & type == "INS" ~
#       ifelse(grepl("intergenic", gene_position),
#              paste0(gene_short, " (intergenic ", nchar(new_seq), "bp INS)"),
#              paste0(gene_short, " (", nchar(new_seq), "bp INS)")),
#     mutation_category == "small_indel" & type == "SUB" ~
#       ifelse(grepl("intergenic", gene_position),
#              paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
#              paste0(gene_short, " (nt ", ref_seq, "->", new_seq, ")")),
#     mutation_category == "large_deletion" ~
#       paste0(gene_short, " (", size, "bp DEL)"),
#     TRUE ~ "nomatch")
# )

mutations <- catmutate(mutations)
mutationsd <- catmutate(mutationsd)

#Check for annotations among 'known' genes
tmp <- 
  summarize(group_by(
    left_join(
      separate_longer_delim(
        filter(
          summarize(group_by(mutationsd,
                             Locus.Tag, gene_category, Product, GOsmerge_ind)),
          !is.na(GOsmerge_ind), GOsmerge_ind != "-"),
        cols = GOsmerge_ind, delim = "|"),
      term2name),
    Locus.Tag, gene_category, Product, GOsmerge_ind, name))
#View(filter(tmp, gene_category == "pilus"))
#View(filter(tmp, grepl("pil", name)))

#Incorporate repetitive regions ----
rptdat <- rbind(read.csv("PA14_PA14_BLAST.csv"),
                #there were no hits for PA14 vs H6
                read.csv("PA14_AX1_BLAST.csv"),
                read.csv("PA14_EF_BLAST.csv"),
                read.csv("PA14_SA_BLAST.csv"))

#drop trivial full self-hit
rptdat <- filter(rptdat, !(alignment.length == 6541479 & pct.identity == 100))

rptdat$rpt_hit_id <- as.character(1:nrow(rptdat))

#Rename query (seq_id) and subject values
rptdat <- mutate(rptdat,
                 query = gsub("contig_", "", query),
                 seq_id = query,
                 subject = gsub("contig_", "", subject),
                 subject = gsub("e82cc082e40344a9", "EF", subject),
                 subject = gsub("6a49b35b19134f7e", "SA", subject),
                 subject = gsub("PP203295", "H6_1", subject))

if(FALSE) {
  rptcov <- rep(0, n = max(rptdat$query.end))
  for(i in 1:nrow(rptdat)) {
    rptcov[rptdat$query.start[i]:rptdat$query.end[i]] <-
      rptcov[rptdat$query.start[i]:rptdat$query.end[i]] + 1
  }
  ggplot(data = rptdat,
         aes(xmin = query.start, xmax = query.end, 
             ymin = 0, ymax = 1,
             fill = subject)) +
    geom_rect(alpha = 0.5) +
    guides(fill = "none")
  ggplot(data = rptcov, aes(x = 1:length(rptcov), y = cov)) + geom_line()
}

id_rpthit <- function(mutdf, rptdf) {
  mutdf$rpt_hit_id <- NA
  for(i in 1:nrow(mutdf)) {
    matches <- filter(rptdf,
                      seq_id == mutdf$seq_id[i] &
                        ((query.start <= mutdf$position_start[i] &
                            query.end >= mutdf$position_start[i]) |
                           (query.start <= mutdf$position_end[i] &
                              query.end >= mutdf$position_end[i])))
    
    if(nrow(matches) > 0) {
      mutdf$rpt_hit_id[i] <-
        as.character(matches$rpt_hit_id[which.max(matches$bit.score)])
    } else {
      mutdf$rpt_hit_id[i] <- "none"
    }
  }
  return(mutdf)
}

mutations <- 
  left_join(id_rpthit(mutations, rptdat),
            select(mutate(rptdat, 
                          rpt_bit_score = bit.score,
                          rpt_match_seq = subject),
                   seq_id, rpt_hit_id, rpt_bit_score, rpt_match_seq))
mutationsd <- 
  left_join(id_rpthit(mutationsd, rptdat),
            select(mutate(rptdat, 
                          rpt_bit_score = bit.score,
                          rpt_match_seq = subject),
                   seq_id, rpt_hit_id, rpt_bit_score, rpt_match_seq))

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
mutpopd <- mutate(group_by(mutpopd, 
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
                                     Pop = function(x) {paste0("Pop", x)})) +
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
    NULL
  dev.off()
  
  png("./plots/mut_trajectories_above010.png", 
      width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = mutate(filter(mutpop, everabove_010 == TRUE),
                    Pop = as.numeric(Pop)),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_catpilnot)) +
      geom_line(aes(alpha = gene_catpilnot)) +
      facet_nested(Treat_P * Treat_C ~ Pop,
                 labeller = labeller(
                   Pop = function(x) {paste0("Pop", LETTERS[as.numeric(x)-2])},
                   Treat_C = c("NC" = "Pseudomonas-only",
                             "C" = "All bacteria"),
                   Treat_P = c("NP" = "No-phage",
                               "P" = "Phage-added"))) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "other"),
                         values = my_cols[c(6, 8)],
                         name = "Gene function") +
      scale_alpha_manual(breaks = c("pilus", "other"),
                         values = c(0.5, 0.2),
                         name = "Gene function") +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      NULL
  )
  dev.off()
  
  png("./plots/mut_trajectories_above010_noanc.png", 
      width = 12, height = 7, units = "in",
      res = 150)
  print(
    ggplot(
      data = mutate(filter(mutpop, everabove_010_noanc == TRUE),
                    Pop = as.numeric(Pop)),
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_catpilnot)) +
      geom_line(aes(alpha = gene_catpilnot)) +
      facet_nested(Treat_P * Treat_C ~ Pop,
                   labeller = labeller(
                     Pop = function(x) {paste0("Pop", LETTERS[as.numeric(x)-2])},
                     Treat_C = c("NC" = "Pseudomonas-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8)) +
      scale_color_manual(breaks = c("pilus", "other"),
                         values = my_cols[c(6, 8)],
                         name = "Gene function") +
      scale_alpha_manual(breaks = c("pilus", "other"),
                         values = c(0.5, 0.2),
                         name = "Gene function") +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      NULL
  )
  dev.off()
  
  pildat <- mutate(
    group_by(filter(mutpop,
                    everabove_010 == TRUE, gene_category == "pilus"),
             gene_name),
    Pop = as.numeric(Pop))
  pildat <- mutate(
    ungroup(pildat),
    gene_name = 
      factor(gene_name,
             levels = 
               arrange(
                 summarize(group_by(pildat, gene_name), freqsum = sum(freq0)), 
                 desc(freqsum))$gene_name)
  )
  pildat <- mutate(
    group_by(pildat, gene_name),
    gene_num = cur_group_id(),
    gene_name_num = paste0(cur_group_id(), ": ", gene_name))
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
  print(
    ggplot(
      data = pildat,
      aes(x = as.numeric(Timepoint), y = freq0, group = mut_id, 
          color = gene_name_num, lty = mutation_category3)) +
      geom_line(lwd = 1.05, alpha = 0.6) +
      geom_text(data = pildat_sum,
                aes(x = max_timepoint, y = max_freq+0.05, label = gene_num),
                size = 4, show.legend = FALSE) +
      facet_nested(Treat_P * Treat_C ~ Pop,
                   labeller = labeller(
                     Pop = function(x) {paste0("Pop", LETTERS[as.numeric(x)-2])},
                     Treat_C = c("NC" = "Pseudomonas-only",
                                 "C" = "All bacteria"),
                     Treat_P = c("NP" = "No-phage",
                                 "P" = "Phage-added"))) +
      scale_x_continuous(breaks = c(0, 1, 2, 4, 8), limits = c(0, 9)) +
      scale_color_manual(name = "Gene name", values = my_cols2) +
      scale_linetype_manual(
        name = "Putative mutation effect",
        breaks = c("loss of function", "nonsyn change", "intergenic", "silent"),
        labels = c("loss of function", "nonsyn. change", "intergenic", "silent"),
        values = c(1, 6, 2, 3)) +
      labs(x = "Timepoint (days)", y = "Frequency") +
      theme_bw() +
      NULL
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

###Analyze gene category frequencies CHANGE TO ONLY PILUS ----
gcat_freq <- summarize(
  group_by(
    filter(
      mutpop, gene_name != "type 4b pilus Flp biogenesis protein TadZ"), 
    gene_category,
    across(any_of(sampcols))),
  total_freq = sum(freq0),
  max_freq = max(freq0),
  indep_freq = 1-prod(1-freq0),
  avg_freq = mean(freq0))

# View(summarize(group_by(gcat_freq, Treat, gene_category),
#           avg_total = mean(total_freq)))

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
      mutpop, gene_category == "pilus",
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

### GO enrichment analysis ----

#GO enrichment analysis using cluterProfiler in R
#GSEA analysis looks for terms that are biased up/down across all genes w/ term
#enrichment looks for overrepresentation of terms in genes at the top
# -- this one is what we use

#Calculate gene frequencies
gnamefreq <- summarize(
  group_by(
    filter(mutpopd, !is.na(Locus.Tag)), 
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

goPvsNP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, PvsNP_non0mean > 0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
      pvalue <= 0.05)
goCvsNC_inNP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, CvsNC_inNP_non0mean > 0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
      pvalue <= 0.05)
goCvsNC_inP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, CvsNC_inP_non0mean > 0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
      pvalue <= 0.05)
goNPvsP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, PvsNP_non0mean < -0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
    pvalue <= 0.05)
goNCvsC_inNP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, CvsNC_inNP_non0mean < -0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
    pvalue <= 0.05)
goNCvsC_inP <-
  filter(
    as.data.frame(
      enricher(
        gene = filter(gnamefreq_sum, CvsNC_inP_non0mean < -0.15)$Locus.Tag,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.5,
        qvalueCutoff = 0.5,
        minGSSize = 2)),
    pvalue <= 0.05)


goPvsNP$Description <- factor(goPvsNP$Description,
                              levels = goPvsNP$Description)
ggplot(filter(goPvsNP, pvalue < 0.01),
       aes(x = FoldEnrichment, y = Description, size = pvalue)) +
  geom_point() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_size_continuous(limits = c(NA, -2), transform = "log10")
  scale_y_discrete(
    breaks = goPvsNP$Description[order(goPvsNP$FoldEnrichment, decreasing = TRUE)])


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

### Create plottable frequency of mutations ----
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
        any(mutation_category3 == "intergenic") ~ "intergenic",
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
  c("none", "multiple", "intergenic",
    pil_order[!pil_order %in% c("none", "multiple", "intergenic")])

if(make_slowplots) {
  png("./plots/isolmut_pilus_freq.png",
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
      values = c("gray70", "black", "gray30",
                 viridis_pal()(length(pil_order)-3))) +
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
}

### Plots of mutation types ----

#Filter out fixed-nonvariate mutations
#mut_id 2958, 3198, and 7885 are fixed across all pops & treatments & timepoints
#mut_id 11212 is fixed in pop 9 across all treatments & timepoints
mutisol <- filter(mutisol, !mut_id %in% c("2958", "3198", "7885", "11212"))

#Note: we sequenced more NCP isols than CP
                    
if(make_slowplots) {
  png("./plots/isolmut_types3.png", width = 7, height = 5, units = "in",
      res = 150)
  ggplot(
    data = mutate(mutisol,
                  mutation_category2 = factor(
                    mutation_category2,
                    levels = c("large_deletion", "small_indel_frameshift",
                               "snp_nonsense", "small_indel_inframe", 
                               "snp_nonsynonymous", "small_indel_intergenic",
                               "snp_intergenic", "snp_synonymous"))),
    aes(x = resis_cat)) +
    geom_bar(aes(fill = mutation_category2), position = "fill") +
    geom_text(data = 
                summarize(group_by(mutisol, resis_cat, Treat_C, Treat_P, Treat), 
                          n = n()),
              aes(y = 0.1, label = paste0("n = ", n)), size = 5) +
    facet_nested(Treat_P * Treat_C ~ gene_catpilnot,
                 labeller = labeller(
                   Pop = function(x) {paste0("Pop", LETTERS[as.numeric(x)-2])},
                   Treat_C = c("NC" = "Pseudomonas-only",
                               "C" = "All bacteria"),
                   Treat_P = c("NP" = "No-phage",
                               "P" = "Phage-added"))) +
    scale_x_discrete(limits = c("part", "resis"),
                     labels = c("Partial", "Resistant")) +
    scale_fill_manual(values = my_cols[c(8, 1:7)], 
                      name = "Mutation type",
                      labels = c("Large deletion", "Frameshift indel",
                                 "Nonsense SNP", "In-frame indel",
                                 "Nonsynonymous SNP", "Intergenic indel",
                                 "Intergenic SNP", "Synonymous SNP")) +
    labs(x = "Isolates' Resistance to Phage",
         y = "Frequency") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 14, 
                                     hjust = 1, vjust = 0.5),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 18))
  dev.off()
}


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
                 Pop = function(x) {paste0("Pop", LETTERS[as.numeric(x)-2])},
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

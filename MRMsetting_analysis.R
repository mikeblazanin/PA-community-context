##Setup ----
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

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
    "pilmut_struc", "gene_short", "concise_mut", "mut_combin", "gene_catpilnot",
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

muttest1 <- read.table("compare_all_test1.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE,
                        quote = '', #don't allow quotes
                        comment.char = "") #don't allow comments since we use # for isolates
muttest2 <- read.table("compare_all_test2.tsv", sep = "\t",
                       stringsAsFactors = FALSE, header = TRUE,
                       quote = '', #don't allow quotes
                       comment.char = "") #don't allow comments since we use # for isolates
muttest1$run <- "Unfiltered_reads"
muttest2$run <- "fastp_filtered_reads"
muttest <- rbind(muttest1, muttest2)

muttest <- filter(muttest, seq_id == "PA14_contig_1")


#Drop columns that contain no info
muttest <- select(muttest,
                    !clone & !codon_position_is_indeterminate &
                      !mutator_status & !population & !time & !treatment)

#Rename seq_id values
muttest <- mutate(muttest,
                    seq_id = gsub("e82cc082e40344a9", "EF", seq_id),
                    seq_id = gsub("6a49b35b19134f7e", "SA", seq_id),
                    seq_id = gsub("PP203295", "H6_1", seq_id),
                    seq_id = gsub("contig_", "", seq_id))

#Rename frequency
colnames(muttest)[colnames(muttest) == "frequency"] <- "freq"

#Number each unique mutation (shared across samples)
muttest <- mutate(group_by(muttest, across(any_of(mutcols))),
                    mut_id = cur_group_id())

muttest <- separate_wider_delim(
  muttest,
  cols = title,
  delim = "_",
  names = c("seqrun", "SeqCoastid", "seqtype",
            "Timepoint", "Pop", "Treat", "Isol", "MRM"),
  cols_remove = FALSE)

ggplot(muttest,
       aes(x = MRM, y = freq)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = mut_id), alpha = 0.5) +
  scale_y_log10() +
  facet_grid(~run) +
  labs(x = "")

ggplot(muttest,
       aes(x = run, y = freq)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = mut_id), alpha = 0.5) +
  scale_y_log10() +
  facet_grid(~MRM) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "")

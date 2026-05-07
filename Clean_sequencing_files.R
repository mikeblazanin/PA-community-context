##Setup ----
library(tidyr)
library(dplyr)
library(stringr)

##Coverage analysis ----
covdat <- read.table("./raw_data/all_depths_MRM3.txt", sep = "\t", comment.char = "",
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

write.csv(covdat, "./clean_data/covdat.csv", row.names = FALSE)
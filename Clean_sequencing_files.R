##Setup ----
library(tidyr)
library(dplyr)
library(stringr)

source("Define_colgroups_sequencing.R")

##Coverage ----
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

##Initial breseq data reading & organization ----
mutations <- read.table("./raw_data/mutations_compare_all_MRM3.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE,
                        quote = '', #don't allow quotes
                        comment.char = "") #don't allow comments since we use # for isolates

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

runcombos <- unique(paste(mutations$seqrun, mutations$SeqCoastid))
allcombos <- c(paste("4777", 1:58), paste("5174", 1:134))
runcombos[!runcombos %in% allcombos]
allcombos[!allcombos %in% runcombos]

#Missing samples:
# 4777_25 - AX T0 P0 Anc #0 - expected bc reference
# 5174_133 - doesn't exist, SeqCoast skipped a tube ID for some reason
# 5174_134 - PA T0 P3 Anc #1 - expected bc reference


#Note: we used to drop non-PA samples right here

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

#Deal with interpro
interpro <- read.table("./raw_data/PA14_interpro_out.tsv", sep = "\t",
                       stringsAsFactors = FALSE, header = FALSE,
                       quote = "", #don't allow quotes
                       comment.char = "") #don't allow comments
colnames(interpro) <- 
  c("Locus.Tag", "seqMD5", "seq_len", "analysis", "signat_accession",
    "signat_description", "start", "stop", "e-value_score",
    "match_status", "query_date", "InterPro_accession",
    "InterPro_annotations", "InterPro_GOs", "Pathways")

interpro_sum <- 
  summarize(group_by(filter(interpro, InterPro_GOs != "-"), 
                     Locus.Tag),
            across(everything(),
                   ~ paste(., collapse = ";")))

#Do bakta and eggnog
bakta <- read.csv("./raw_data/PA14_bakta_annotated.csv")
bakta$Locus.Tag[bakta$Locus.Tag == ""] <-
  paste0("PEANK_MB", 
         formatC(1:sum(bakta$Locus.Tag == ""),
                 width = floor(log10(sum(bakta$Locus.Tag == "")))+1,
                 flag = "0"))

eggnog <- read.delim("./raw_data/PA14_eggnog_out3.emapper.annotations",
                     comment = "#",
                     header = FALSE)
colnames(eggnog) <-
  gsub("#", "", 
       read.delim("./raw_data/PA14_eggnog_out3.emapper.annotations",
                  comment = "", header = FALSE,
                  skip = 4, nrows = 1))
colnames(eggnog)[colnames(eggnog) == "query"] <- "Locus.Tag"

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
        NA),
      gene_short = case_when(
        !Gene %in% c("", "|NA", "NA|", NA) ~ Gene,
        !is.na(Product) ~ Product,
        TRUE ~ ""),
      concise_mut = case_when(
        mut_id %in% c("1071", "1167", "2787", "5034") ~ "", #(ancestral mutations)
        grepl("snp", mutation_category) ~
          ifelse(mutation_category_manual == "intergenic",
                 paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
                 paste0(gene_short, " (", aa_ref_seq, "->", aa_new_seq, ")")),
        mutation_category == "small_indel" & type == "DEL" ~
          ifelse(mutation_category_manual == "intergenic",
                 paste0(gene_short, " (intergenic ", size, "bp DEL)"),
                 paste0(gene_short, " (", size, "bp DEL)")),
        mutation_category == "small_indel" & type == "INS" ~
          ifelse(mutation_category_manual == "intergenic",
                 paste0(gene_short, " (intergenic ", nchar(new_seq), "bp INS)"),
                 paste0(gene_short, " (", nchar(new_seq), "bp INS)")),
        mutation_category == "small_indel" & type == "SUB" ~
          ifelse(mutation_category_manual == "intergenic",
                 paste0(gene_short, " (intergenic nt ", ref_seq, "->", new_seq, ")"),
                 paste0(gene_short, " (nt ", ref_seq, "->", new_seq, ")")),
        mutation_category == "large_deletion" ~
          paste0(gene_short, " (", size, "bp DEL)"),
        TRUE ~ "nomatch"),
      concise_mut = trimws(concise_mut, which = "both"),
    )
  )
}

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

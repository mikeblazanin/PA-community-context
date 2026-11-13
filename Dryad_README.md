#### File: aggregate_isolate_data.csv

**Description:** A file containing the aggregated results of all phenotypes measured for bacterial isolates.

##### Variables

* twitch_date: the date twitching motility data was collected
* Bact_pop: Bacterial replicate population
* Bact_timepoint: Bacterial timepoint
* Bact_treat: Bacterial treatment (NC for no community, C for community; NP for no phage, P for phage)
* Bact_treat_C: Bacterial community treatment (NC for no community, C for community)
* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage)
* Bact_isol: Bacterial isolate number
* cm_day: centimeters per day measured from twitching motility assay
* cm_day_norm: centimeters per day measured from twitching motility assay, batch corrected
* clearing_avg: whether bacteria was cleared by phages, averaged across all phages tested
* EOP_avg: efficiency of plaquing from phages, averaged across all phages tested
* Plaque_morphology: morphology of phage plaques
* wavelength_405: absorbance of overnight bacterial cultures at 405 nm
* wavelength_695: absorbance of overnight bacterial cultures at 695 nm
* gc_date: date when growth curve was run
* Well: well of plate where growth curve was run
* max_percap: maximum cellular (i.e. per-capita) growth rate from growth curve
* max_dens: maximum density from growth curve (OD600)
* lag: lag time from growth curve
* auc: area under the curve from growth curve
* max_percap_time: time when maximum cellular growth rate was achieved from growth curve
* max_percap_dens: density when maximum cellular growth rate was achieved from growth curve
* min_dens: minimum density from growth curve
* max_percap_norm: maximum cellular (i.e. per-capita) growth rate from growth curve, batch corrected
* max_dens_norm: maximum density from growth curve (OD600), batch corrected
* auc_norm: area under the curve from growth curve, batch corrected
* lag_norm: lag time from growth curve, batch corrected
* biofilm_abs: absorbance corresponding to biofilm formation
* biofilm_norm: absorbance corresponding to biofilm formation, batch corrected
* biofilm_norm_log10: absorbance corresponding to biofilm formation, batch corrected, log10
* biofilm_norm_log10_norm: absorbance corresponding to biofilm formation, batch corrected, log10, batch corrected
* biofilm_block_name: unique name of biofilm batch
* resis: 1 - clearing_avg, the level of resistance
* resis_cat: categorical level of resistance (resistant, partial, or susceptible)
* twitch_cat: categorical level of twitching (full if above 0.7, none if below 0.5, intermediate in-between)

#### File: biofilm.csv

**Description:** 
Raw and batch-corrected absorbance measurements from the biofilm formation assay, with bacterial isolate and treatment metadata.

##### Variables

* block_name: unique name of biofilm batch
* Well: well of plate where isolate was run
* biofilm_abs: absorbance corresponding to biofilm formation
* Bact_pop: Bacterial replicate population
* Bact_timepoint: Bacterial timepoint
* Bact_treat: Bacterial treatment (NC for no community, C for community; NP for no phage, P for phage)
* Bact_treat_C: Bacterial community treatment (NC for no community, C for community)
* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage)
* Bact_isol: Bacterial isolate number
* biofilm_norm: absorbance corresponding to biofilm formation, batch corrected

#### File: density.csv

**Description:** 
Bacterial colony-count data and derived density estimates, including small- and large-colony morphotypes, across sampling timepoints and treatments.

##### Variables

* Timepoint: Sampling timepoint.
* Phage_treat: Phage treatment (NP for no phage, P for phage).
* Comm_treat: Community treatment (NC for no community, C for community).
* Pop: Replicate population.
* Strain: Species counted.
* Count: Number of colonies counted on the plate.
* Dilution: Dilution of the sample plated (log10).
* Small_count: Number of small-colony morphotype colonies counted.
* Small_dilution: Dilution of the sample used to count small-colony morphotype colonies.
* Large_count: Number of large-colony morphotype colonies counted.
* Large_dilution: Dilution of the sample used to count large-colony morphotype colonies.
* treat_P: Phage treatment (NP for no phage, P for phage).
* treat_C: Community treatment (NC for no community, C for community).
* cfu_ml: Estimated colony-forming units per milliliter (CFU/mL) in the original sample.
* S_cfu_ml: Estimated small-colony morphotype CFU per milliliter.
* L_cfu_ml: Estimated large-colony morphotype CFU per milliliter.
* L_freq: Frequency of the large-colony morphotype.
* bd: whether any counts were below detection (i.e. 0 CFUs observed)
* log_cfuml: Log10-transformed CFU per milliliter.

#### File: growth_curve_metrics.csv

**Description:** 
Growth-curve metrics calculated for bacterial isolates, including raw and batch-corrected growth characteristics.

##### Variables

* date: Date the assay was performed.
* Well: Well of the plate in which the sample was assayed.
* Bact_treat: Combined bacterial treatment (community and phage treatment; NC/C and NP/P).
* Bact_treat_C: Bacterial community treatment (NC for no community, C for community).
* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage).
* Bact_timepoint: Bacterial sampling timepoint.
* Bact_pop: Bacterial replicate population.
* Bact_isol: Bacterial isolate number.
* max_percap: Maximum cellular (i.e. per-capita) growth rate from growth curve.
* max_dens: Maximum density from growth curve (OD600).
* lag: Lag time from growth curve.
* auc: Area under the curve from growth curve.
* max_percap_time: Time when maximum cellular growth rate was achieved from growth curve.
* max_percap_dens: Density when maximum cellular growth rate was achieved from growth curve.
* min_dens: Minimum density from growth curve.
* max_percap_norm: Maximum cellular (i.e. per-capita) growth rate from growth curve, batch corrected.
* max_dens_norm: Maximum density from growth curve (OD600), batch corrected.
* auc_norm: Area under the curve from growth curve, batch corrected.
* lag_norm: Lag time from growth curve, batch corrected.

#### File: GOresults.csv

**Description:** 
Results of Gene Ontology (GO) enrichment analyses, including enrichment statistics and the genes associated with each enriched GO term, as output by clusterProfiler.

##### Variables

* Treat: Combined treatment (community and phage treatment; NC/C and NP/P).
* Context: Context or comparison grouping used for the GO enrichment analysis.
* ID: Gene Ontology (GO) term identifier.
* Description: Description/name of the GO term.
* GeneRatio: Ratio of input genes associated with the GO term to input genes represented in the enrichment analysis.
* BgRatio: Ratio of background/universe genes associated with the GO term to background genes represented in the enrichment analysis.
* RichFactor: Number of input genes associated with the GO term divided by the total number of background genes associated with that term.
* FoldEnrichment: Enrichment of the GO term in the input gene set relative to the background; calculated as GeneRatio divided by BgRatio.
* zScore: Standardized enrichment score comparing the observed number of input genes in the term with the expected number under the background distribution.
* pvalue: Unadjusted p-value for enrichment of the GO term.
* p.adjust: P-value adjusted for multiple hypothesis testing.
* qvalue: Q-value estimate for false-discovery-rate control.
* geneID: Gene identifiers from the input gene set that are associated with the GO term.
* Count: Number of input genes associated with the GO term.
* FoldPlot: Fold-enrichment value used for plotting.
* TreatPlot: Treatment/grouping label used for plotting.

#### File: covdat.csv

**Description:** 
Sequencing coverage statistics for reference sequences/contigs, with sequencing sample and treatment metadata.

##### Variables

* seqrun: Sequencing run identifier.
* SeqCoastid: Sequencing sample identifier.
* seqtype: Type of sequencing sample (PA for bacterial, H6 for viral).
* Timepoint: Sampling timepoint.
* Pop: Replicate population.
* Treat: Combined treatment (community and phage treatment; NC/C and NP/P).
* Isol: Isolate number.
* title: Sample or sequencing-library title.
* seq_id: Identifier of the reference sequence/contig.
* startpos: Start coordinate of the reference sequence region summarized.
* endpos: End coordinate of the reference sequence region summarized.
* numreads: Number of reads mapped to the reference sequence region.
* covbases: Number of bases in the reference sequence region covered by at least one read.
* coverage: Percentage of bases in the reference sequence region covered by at least one read.
* meandepth: Mean read depth across the reference sequence region.
* meanbaseq: Mean base quality of bases covering the reference sequence region.
* meanmapq: Mean mapping quality of reads covering the reference sequence region.
* Strain: Microbial strain.
* Treat_C: Community treatment (NC for no community, C for community).
* Treat_P: Phage treatment (NP for no phage, P for phage).
* included_in_treat: Boolean of whether the strain was included in the treatment (and therefore expected to have coverage).

#### File: nice_isolmuts.csv

**Description:** 
Summary of bacterial isolate phenotypes and mutation-category counts or indicators, with treatment and sampling metadata.

##### Variables

* Timepoint: Sampling timepoint.
* Treat: Combined treatment (community and phage treatment; NC/C and NP/P).
* Pop: Replicate population.
* Isol: Isolate number.
* Treat_C: Community treatment (NC for no community, C for community).
* Treat_P: Phage treatment (NP for no phage, P for phage).
* clearing_avg: Whether bacteria was cleared by phages, averaged across all phages tested.
* EOP_avg: Efficiency of plaquing from phages, averaged across all phages tested.
* Plaque_morphology: Morphology of phage plaques.
* resis_cat: Categorical level of resistance (resistant, partial, or susceptible).
* twitch_cat: Categorical level of twitching (full, intermediate, or none).
* mut_in_other: List of mutations assigned to the other gene category.
* mut_in_mutator: List of mutations in genes associated with mutator phenotypes.
* mut_in_pilus: List of mutations in pilus-associated genes.
* mut_in_biofilm: List of mutations in biofilm-associated genes.

#### File: pig.csv

**Description:** 
Plate-reader absorbance measurements collected for the pigment assay for pyocyanin and pyoverdin, with bacterial isolate and treatment metadata.

##### Variables

* file: Name of the source data file.
* Wavelength: Wavelength used for the plate-reader absorbance measurement.
* Well: Well of the plate in which the sample was assayed.
* Measurements: Raw plate-reader measurement.
* Bact_pop: Bacterial replicate population.
* Bact_timepoint: Bacterial sampling timepoint.
* Bact_treat: Combined bacterial treatment (community and phage treatment; NC/C and NP/P).
* Bact_isol: Bacterial isolate number.
* Bact_treat_C: Bacterial community treatment (NC for no community, C for community).
* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage).
* Meas_norm: Plate-reader measurement after normalization.

#### File: prelimdens_C.csv

**Description:** 
Bacterial density and large-colony morphotype frequency estimates for preliminary experiment in community-present condition.

##### Variables

* Comm_treat: Community treatment (NC for no community, C for community).
* Timepoint: Sampling timepoint.
* Phage_treat: Phage treatment (NP for no phage, P for phage).
* Population: Replicate population.
* Strain: Bacterial strain.
* freq_L: Frequency of the large-colony morphotype.
* dens_ml: Bacterial density, expressed as cells or CFU per milliliter.
* bd: whether any counts were below detection (i.e. 0 CFUs observed)

#### File: prelimdens_NC_sum.csv

**Description:** 
Bacterial density and large-colony morphotype frequency estimates for preliminary experiment in no-community (Pseudomonas-only) condition.

##### Variables

* Comm_treat: Community treatment (NC for no community, C for community).
* Timepoint: Sampling timepoint.
* Strain: Bacterial strain.
* Phage_treat: Phage treatment (NP for no phage, P for phage).
* Population: Replicate population.
* dens_ml: Bacterial density, expressed as cells or CFU per milliliter.
* freq_L: Frequency of the large-colony morphotype.
* bd: whether any counts were below detection (i.e. 0 CFUs observed)

#### File: gcdat.csv

**Description:** 

##### Variables

* file:
* Time:
* Well:
* Measurements:
* date:
* Time_hr:
* block_name:
* Bact_pop:
* Bact_timepoint:
* Bact_treat:
* Bact_isol:
* Bact_treat_C:
* Bact_treat_P:
* meas_blanked:
* sm_meas:
* deriv:
* deriv_percap:

#### File: resis_isols.csv

**Description:** 

##### Variables

* Bact_pop:
* Phage_pop:
* Bact_timepoint:
* Bact_treat_C:
* Bact_treat_P:
* Bact_treat:
* Bact_isol:
* Phage_timepoint:
* Phage_treat_C:
* Phage_treat_P:
* Phage_treat:
* Phage_isol:
* Plaque_morphology:
* pfu2:
* dilution2:
* Estimated_pfu2:
* bd2:
* pfu_ml2:
* EOP2:
* EOP2_plot:
* clearing2:
* run_EOP2:
* Date:
* clearing:
* run_cs:
* run_spot:
* bd:
* pfu_ml:
* EOP:
* EOP_plot:
* run_EOP:

#### File: twitching.csv

**Description:** 

##### Variables

* Date:
* Bact_pop:
* Bact_treat:
* Bact_timepoint:
* Bact_isol:
* width:
* height:
* cm_day:
* Bact_treat_C:
* Bact_treat_P:

#### File: term2name.csv

**Description:** 

##### Variables

* :
* GOsmerge_ind:
* name:

#### File: term2gene.csv

**Description:** 

##### Variables

* GOsmerge:
* Locus.Tag:

#### File: mutations_PAseq_PA14seqid.csv

**Description:** 

##### Variables

* mut_id:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* Treat_C:
* Treat_P:
* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* transl_table:
* type:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: mutationsd_PAseq_PA14seqid.csv

**Description:** 

##### Variables

* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* transl_table:
* type:
* mut_id:
* Treat_C:
* Treat_P:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: PA14_interpro_out_sum.tsv

**Description:** 

##### Variables

* Locus.Tag:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:

#### File: mutations_PAseq_otherseqid.csv

**Description:** 

##### Variables

* mut_id:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* Treat_C:
* Treat_P:
* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* transl_table:
* type:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: mutationsd_PAseq_otherseqid.csv

**Description:** 

##### Variables

* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* transl_table:
* type:
* mut_id:
* Treat_C:
* Treat_P:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: mutations_H6seq.csv

**Description:** 

##### Variables

* mut_id:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* Treat_C:
* Treat_P:
* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* transl_table:
* type:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: mutationsd_H6seq.csv

**Description:** 

##### Variables

* aa_new_seq:
* aa_position:
* aa_ref_seq:
* codon_new_seq:
* codon_number:
* codon_position:
* codon_ref_seq:
* freq:
* gene_name:
* gene_position:
* gene_product:
* gene_strand:
* genes_inactivated:
* genes_overlapping:
* genes_promoter:
* insert_position:
* locus_tag:
* locus_tags_inactivated:
* locus_tags_overlapping:
* locus_tags_promoter:
* multiple_polymorphic_SNPs_in_same_codon:
* mutation_category:
* new_read_count:
* new_read_count_basis:
* new_seq:
* position:
* position_end:
* position_start:
* ref_read_count:
* ref_read_count_basis:
* ref_seq:
* repeat_length:
* repeat_new_copies:
* repeat_ref_copies:
* repeat_seq:
* seq_id:
* size:
* snp_type:
* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* transl_table:
* type:
* mut_id:
* Treat_C:
* Treat_P:
* Locus.Tag:
* mutation_category_manual:
* Sequence.Id:
* Type:
* Start:
* Stop:
* Strand:
* Gene:
* Product:
* DbXrefs:
* seed_ortholog:
* evalue:
* score:
* eggNOG_OGs:
* max_annot_lvl:
* COG_category:
* Description:
* Preferred_name:
* GOs:
* EC:
* KEGG_ko:
* KEGG_Pathway:
* KEGG_Module:
* KEGG_Reaction:
* KEGG_rclass:
* BRITE:
* KEGG_TC:
* CAZy:
* BiGG_Reaction:
* PFAMs:
* seqMD5:
* seq_len:
* analysis:
* signat_accession:
* signat_description:
* start:
* stop:
* e-value_score:
* match_status:
* query_date:
* InterPro_accession:
* InterPro_annotations:
* InterPro_GOs:
* Pathways:
* GOsmerge:
* GOsmerge_ind:
* mutation_category2:
* mutation_category3:
* gene_category:
* gene_catpilnot:
* pilmut_struc:
* gene_short:
* concise_mut:
* rpt_hit_id:
* rpt_bit_score:
* rpt_match_seq:

#### File: Analyze_data.R

**Description:** 

#### File: Analyze_sequencing_data.R

**Description:** 

#### File: Define_colgroups_sequencing.R

**Description:** 

#### File: renv.lock

**Description:** 

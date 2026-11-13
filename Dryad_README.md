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

##### Variables

* Timepoint:
* Phage_treat:
* Comm_treat:
* Pop:
* Strain:
* Count:
* Dilution:
* Small_count:
* Small_dilution:
* Large_count:
* Large_dilution:
* treat_P:
* treat_C:
* cfu_ml:
* S_cfu_ml:
* L_cfu_ml:
* L_freq:
* bd:
* log_cfuml:

#### File: growth_curve_metrics.csv

**Description:** 

##### Variables

* date:
* Well:
* Bact_treat:
* Bact_treat_C:
* Bact_treat_P:
* Bact_timepoint:
* Bact_pop:
* Bact_isol:
* max_percap:
* max_dens:
* lag:
* auc:
* max_percap_time:
* max_percap_dens:
* min_dens:
* max_percap_norm:
* max_dens_norm:
* auc_norm:
* lag_norm:

#### File: GOresults.csv

**Description:** 

##### Variables

* Treat:
* Context:
* ID:
* Description:
* GeneRatio:
* BgRatio:
* RichFactor:
* FoldEnrichment:
* zScore:
* pvalue:
* p.adjust:
* qvalue:
* geneID:
* Count:
* FoldPlot:
* TreatPlot:

#### File: covdat.csv

**Description:** 

##### Variables

* seqrun:
* SeqCoastid:
* seqtype:
* Timepoint:
* Pop:
* Treat:
* Isol:
* title:
* seq_id:
* startpos:
* endpos:
* numreads:
* covbases:
* coverage:
* meandepth:
* meanbaseq:
* meanmapq:
* Strain:
* Treat_C:
* Treat_P:
* included_in_treat:

#### File: nice_isolmuts.csv

**Description:** 

##### Variables

* Timepoint:
* Treat:
* Pop:
* Isol:
* Treat_C:
* Treat_P:
* clearing_avg:
* EOP_avg:
* Plaque_morphology:
* resis_cat:
* twitch_cat:
* mut_in_other:
* mut_in_mutator:
* mut_in_pilus:
* mut_in_biofilm:

#### File: pig.csv

**Description:** 

##### Variables

* file:
* Wavelength:
* Well:
* Measurements:
* Bact_pop:
* Bact_timepoint:
* Bact_treat:
* Bact_isol:
* Bact_treat_C:
* Bact_treat_P:
* Meas_norm:

#### File: prelimdens_C.csv

**Description:** 

##### Variables

* Comm_treat:
* Timepoint:
* Phage_treat:
* Population:
* Strain:
* freq_L:
* dens_ml:
* bd:

#### File: prelimdens_NC_sum.csv

**Description:** 

##### Variables

* Comm_treat:
* Timepoint:
* Strain:
* Phage_treat:
* Population:
* dens_ml:
* freq_L:
* bd:

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

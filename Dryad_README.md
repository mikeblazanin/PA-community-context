# Community diversity favors coexistence between bacteria and parasitic bacteriophages

Dataset DOI: [10.5061/dryad.zcrjdfnw7](https://doi.org/10.5061/dryad.zcrjdfnw7)

## Description of the data and file structure

We experimentally evolved *Pseudomonas aeruginosa* bacteria and a *Pseudomonas*-targeting phage with and without a community representing a cystic fibrosis disease pulmonary infection. We isolated evolved Pseudomonas and phages, and quantified infectivity/resistance between bacteria-phage pairs, twitching motility of bacteria, growth curves of bacteria, biofilm formation of bacteria, and pigment production of bacteria. We whole-genome sequenced bacterial isolates, and did metagenomics of our mixed experimental communities of bacteria (with or without phage).

### Files and variables

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

Raw and processed plate-reader measurements collected during bacterial growth-curve assays, with isolate and treatment metadata.

##### Variables

* file: Name of the source data file.

* Time: Time at which the plate-reader measurement was recorded.

* Well: Well of the plate in which the sample was assayed.

* Measurements: Raw plate-reader measurement (OD600).

* date: Date the assay was performed.

* Time_hr: Time since beginning the growth curve at which the plate-reader measurement was recorded, in hours.

* block_name: Unique name of the assay batch.

* Bact_pop: Bacterial replicate population.

* Bact_timepoint: Bacterial sampling timepoint.

* Bact_treat: Combined bacterial treatment (community and phage treatment; NC/C and NP/P).

* Bact_isol: Bacterial isolate number.

* Bact_treat_C: Bacterial community treatment (NC for no community, C for community).

* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage).

* meas_blanked: Plate-reader measurement after blank subtraction.

* sm_meas: Smoothed plate-reader measurement.

* deriv: Derivative of the smoothed measurement with respect to time.

* deriv_percap: Per-capita derivative of the smoothed measurement (derivative divided by measurement).

#### File: resis_isols.csv

**Description:** 

Phage resistance assay measurements for bacterial and phage isolates, including plaque counts, efficiency of plaquing, and bacterial clearing.

##### Variables

* Bact_pop: Bacterial replicate population.

* Phage_pop: Phage replicate population.

* Bact_timepoint: Bacterial sampling timepoint.

* Bact_treat_C: Bacterial community treatment (NC for no community, C for community).

* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage).

* Bact_treat: Combined bacterial treatment (community and phage treatment; NC/C and NP/P).

* Bact_isol: Bacterial isolate number.

* Phage_timepoint: Phage sampling timepoint.

* Phage_treat_C: Phage community treatment (NC for no community, C for community).

* Phage_treat_P: Phage phage treatment (NP for no phage, P for phage).

* Phage_treat: Combined phage treatment (community and phage treatment; NC/C and NP/P).

* Phage_isol: Phage isolate number.

* Plaque_morphology: Morphology of phage plaques.

* pfu2: Number of plaques counted.

* dilution2: Dilution where plaques were counted.

* Estimated_pfu2: Estimated plaque-forming units.

* bd2: Whether PFUs were below detection (i.e. 0 counted)

* pfu_ml2: Estimated plaque-forming units per milliliter (PFU/mL).

* EOP2: Efficiency of plaquing (EOP), relative to the most susceptible host.

* EOP2_plot: EOP2 value formatted for plotting.

* clearing2: Whether the bacterial isolate was cleared at all.

* run_EOP2: Boolean indicating whether data was collected in the second run of EOP assay batches.

* Date: Date the assay was performed.

* clearing: Whether the bacterial isolate was cleared by the tested phage (first run of EOP assay batches)

* run_cs: Boolean indicating whether cross-streak data was collected.

* run_spot: Boolean indicating whether spot test data was collected

* bd: Whether PFUs were below detection (first run of EOP assay batches)

* pfu_ml: Estimated plaque-forming units per milliliter (PFU/mL) (first run of EOP assay batches)

* EOP: Efficiency of plaquing (EOP) of the phage on the bacterial isolate, relative to the most susceptible host (first run of EOP assay batches).

* EOP_plot: EOP value formatted for plotting.

* run_EOP: Boolean indicating whether data was collected in the first run of EOP assay batches.

#### File: twitching.csv

**Description:** 

Measurements from the bacterial twitching motility assay, with bacterial isolate and treatment metadata.

##### Variables

* Date: Date the assay was performed.

* Bact_pop: Bacterial replicate population.

* Bact_treat: Combined bacterial treatment (community and phage treatment; NC/C and NP/P).

* Bact_timepoint: Bacterial sampling timepoint.

* Bact_isol: Bacterial isolate number.

* width: Width of the twitching motility zone measured in the assay.

* height: Height of the twitching motility zone measured in the assay.

* cm_day: Distance of twitching motility centimeters per day, calculated from the assay measurements.

* Bact_treat_C: Bacterial community treatment (NC for no community, C for community).

* Bact_treat_P: Bacterial phage treatment (NP for no phage, P for phage).

#### File: term2name.csv

**Description:** 

Mapping table linking GO term identifiers used in the merged GO annotation analysis to GO term names.

##### Variables

* : rownames for table, formatted as required for clusterProfiler

* GOsmerge_ind: GO term identifier/index used to distinguish entries in the merged GO annotation table.

* name: Name/description associated with the GO term.

#### File: term2gene.csv

**Description:** 

Mapping table linking GO term identifiers to bacterial gene locus tags for GO enrichment analysis.

##### Variables

* GOsmerge: GO term identifier used in the merged GO annotation table.

* Locus.Tag: Locus tag identifying the annotated gene.

#### Files: mutations_PAseq_PA14seqid.csv, mutationsd_PAseq_PA14seqid.csv, mutations_PAseq_otherseqid.csv, mutationsd_PAseq_otherseqid.csv, mutations_H6seq.csv, mutationsd_H6seq.csv

**Description:** 

Mutation calls from breseq for sequencing samples, mapped to the PA14 reference contigs ("PA14seqid") or other contigs ("otherseqid"), from sequencing of DNA extracted for bacteria ("PAseq") or DNA viruses ("H6seq"). Mutations are listed individually with all genes they affect ("mutations" files), or mutations hitting multiple genes have been split into one row per gene ("mutationsd" files).

##### Variables

* mut_id: Unique mutation identifier.

* seqrun: Sequencing run identifier.

* SeqCoastid: Sequencing sample identifier.

* seqtype: Type of sequencing sample (PA for bacterial, H6 for viral).

* Timepoint: Sampling timepoint.

* Pop: Replicate population.

* Treat: Treatment group.

* Isol: Isolate number.

* title: Full name of sequencing sample (seqrun, SeqCoastid, seqtype, Timepoint, Pop, Treat, Isol)

* Treat_C: Community treatment (NC for no community, C for community).

* Treat_P: Phage treatment (NP for no phage, P for phage).

* aa_new_seq: Amino acid encoded by the mutated sequence at the affected position.

* aa_position: Position of the affected amino acid in the protein sequence.

* aa_ref_seq: Reference amino acid at the affected position.

* codon_new_seq: Mutant codon sequence.

* codon_number: Number of the affected codon in the coding sequence.

* codon_position: Position within the affected codon (1, 2, or 3).

* codon_ref_seq: Reference codon sequence.

* freq: Estimated frequency of the mutation in the sample.

* gene_name: Name of the gene affected by the mutation.

* gene_position: Position of the mutation relative to the affected gene.

* gene_product: Annotated product/function of the affected gene.

* gene_strand: Strand of the affected gene.

* genes_inactivated: Gene(s) predicted to be inactivated by the mutation.

* genes_overlapping: Gene(s) overlapping the mutation position.

* genes_promoter: Gene(s) whose promoter region overlaps the mutation.

* insert_position: Position at which an insertion occurs.

* locus_tag: Locus tag of the affected gene.

* locus_tags_inactivated: Locus tag(s) of gene(s) predicted to be inactivated by the mutation.

* locus_tags_overlapping: Locus tag(s) of genes overlapping the mutation position.

* locus_tags_promoter: Locus tag(s) of genes whose promoter region overlaps the mutation.

* multiple_polymorphic_SNPs_in_same_codon: Indicator of whether multiple polymorphic single-nucleotide polymorphisms occur in the same codon.

* mutation_category: Mutation class/category assigned by breseq.

* new_read_count: Number of sequencing reads supporting the alternate (mutant) allele.

* new_read_count_basis: Basis or denominator used to calculate the alternate-allele read count/frequency.

* new_seq: Alternate (mutant) nucleotide sequence.

* position: Reference-genome coordinate at which the mutation begins.

* position_end: End coordinate of the affected region in the reference sequence.

* position_start: Start coordinate of the affected region in the reference sequence.

* ref_read_count: Number of sequencing reads supporting the reference allele.

* ref_read_count_basis: Basis or denominator used to calculate the reference-allele read count/frequency.

* ref_seq: Reference nucleotide sequence.

* repeat_length: Length of the repeat region affected by the mutation.

* repeat_new_copies: Number of repeat copies in the mutated sequence.

* repeat_ref_copies: Number of repeat copies in the reference sequence.

* repeat_seq: Nucleotide sequence of the repeated element.

* seq_id: Identifier of the reference sequence/contig.

* size: Size of the mutation or affected sequence, in base pairs.

* snp_type: Classification of the single-nucleotide polymorphism.

* transl_table: Genetic code/translation table used to translate the coding sequence.

* type: Type of mutation as classified by breseq.

* Locus.Tag: Locus tag identifying the annotated gene.

* mutation_category_manual: Manually assigned mutation category.

* Sequence.Id: Identifier of the sequence/contig containing the annotated feature.

* Type: Type of annotated genomic feature.

* Start: Start coordinate of the annotated feature on the sequence.

* Stop: End coordinate of the annotated feature on the sequence.

* Strand: Strand orientation of the annotated feature.

* Gene: Gene name or identifier.

* Product: Annotated product or function of the gene.

* DbXrefs: Cross-references to external biological databases.

* seed_ortholog: Best seed ortholog used by eggNOG-mapper for functional annotation.

* evalue: E-value for the sequence similarity/orthology assignment.

* score: Score associated with the sequence similarity/orthology assignment.

* eggNOG_OGs: eggNOG orthologous group(s) assigned to the gene.

* max_annot_lvl: Most specific taxonomic level at which the eggNOG annotation was assigned.

* COG_category: Cluster of Orthologous Groups functional category.

* Description: Functional description of the gene or protein assigned by the eggNOG functional annotation.

* Preferred_name: Preferred gene name assigned by the functional annotation.

* GOs: Gene Ontology terms associated with the gene.

* EC: Enzyme Commission number(s) associated with the gene product.

* KEGG_ko: KEGG Orthology identifier(s) assigned to the gene.

* KEGG_Pathway: KEGG pathway identifier(s) associated with the gene.

* KEGG_Module: KEGG module identifier(s) associated with the gene.

* KEGG_Reaction: KEGG reaction identifier(s) associated with the gene.

* KEGG_rclass: KEGG reaction-class identifier(s) associated with the gene.

* BRITE: KEGG BRITE hierarchy annotation(s).

* KEGG_TC: KEGG transporter classification identifier(s).

* CAZy: Carbohydrate-Active enZYmes (CAZy) family annotation(s).

* BiGG_Reaction: BiGG reaction identifier(s) associated with the gene.

* PFAMs: Pfam protein family/domain annotation(s).

* seqMD5: MD5 checksum of the protein sequence.

* seq_len: Length of the protein sequence, in amino acids.

* analysis: InterProScan analysis/member database that detected the protein signature (e.g., Pfam).

* signat_accession: Accession identifier of the detected signature/model in the member database.

* signat_description: Description of the detected signature/model.

* start: Start position of the signature match on the protein sequence (1-based, inclusive).

* stop: End position of the signature match on the protein sequence (1-based, inclusive).

* e-value_score: E-value or score reported for the signature match by the member database analysis.

* match_status: Status of the signature match (T indicates a true match).

* query_date: Date on which the InterProScan query was run.

* InterPro_accession: InterPro entry accession associated with the signature match, when available.

* InterPro_annotations: Description/name of the InterPro entry associated with the signature match, when available.

* InterPro_GOs: Gene Ontology annotations associated with the InterPro match, when available.

* Pathways: Pathway cross-references associated with the InterPro match, when available.

* GOsmerge: GO term identifier used in the merged GO annotation table.

* GOsmerge_ind: GO term identifier/index used to distinguish entries in the merged GO annotation table.

* mutation_category2: Additional mutation category assigned during downstream annotation.

* mutation_category3: Additional mutation category assigned during downstream annotation.

* gene_category: Functional/category grouping assigned to the affected gene.

* gene_catpilnot: Gene category indicating pilus-related versus non-pilus-related genes.

* pilmut_struc: Structural classification of the pilus-associated mutation.

* gene_short: Shortened gene name/label.

* concise_mut: Concise textual description of the mutation.

* rpt_hit_id: Identifier of the repeat/sequence hit reported by the repeat annotation step.

* rpt_bit_score: Bit score for the repeat/sequence hit.

* rpt_match_seq: Sequence matching the repeat/sequence hit.

#### File: PA14_interpro_out_sum.tsv

**Description:** 

InterProScan protein-signature and functional annotation results for the PA14 reference genome.

##### Variables

* Locus.Tag: Locus tag identifying the annotated gene.

* seqMD5: MD5 checksum of the protein sequence.

* seq_len: Length of the protein sequence, in amino acids.

* analysis: InterProScan analysis/member database that detected the protein signature (e.g., Pfam).

* signat_accession: Accession identifier of the detected signature/model in the member database.

* signat_description: Description of the detected signature/model.

* start: Start position of the signature match on the protein sequence (1-based, inclusive).

* stop: End position of the signature match on the protein sequence (1-based, inclusive).

* e-value_score: E-value or score reported for the signature match by the member database analysis.

* match_status: Status of the signature match (T indicates a true match).

* query_date: Date on which the InterProScan query was run.

* InterPro_accession: InterPro entry accession associated with the signature match, when available.

* InterPro_annotations: Description/name of the InterPro entry associated with the signature match, when available.

* InterPro_GOs: Gene Ontology annotations associated with the InterPro match, when available.

* Pathways: Pathway cross-references associated with the InterPro match, when available.

## Code/software

R is needed to re-run analysis files included (described below). All exact version and package dependencies are described in the renv.lock file (described below). In order to re-create the environment where analyses were run, download all files into a single directory. Install R 4.5.1, set the working directory to the directory containing files, and download renv 1.1.8 (as in the renv.lock file). Then run renv::restore()

#### File: renv.lock

**Description:** renv lockfile recording the R package dependencies and package versions used for the analysis.

Once the software environment is restored, run Analyze_data.R to regenerate all phenotypic analyses and plots, or Analyze_sequencing_data.R to regenerate all genotypic analyses and plots.

#### File: Analyze_data.R

**Description:** 

R script used to analyze the experimental phenotype and isolate data.

#### File: Analyze_sequencing_data.R

**Description:** 

R script used to process and analyze sequencing data, including mutation and annotation data.

#### File: Define_colgroups_sequencing.R

**Description:** 

R script defining column groups or annotations used in the sequencing-data analysis/visualization.

## Access information

Other publicly accessible locations of the data:

* [https://github.com/mikeblazanin/PA-community-context](https://github.com/mikeblazanin/PA-community-context)

Data was derived from the following sources:

* None

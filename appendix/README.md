# Appendix

This folder contains all the notebooks and code used to generate the figures for the thesis.

## Contents:

### In House Packages

1. g4predict: A regular expression matching package for identifying putative G-Quadruplexes by the Quadparser method, plus methods for filtering them, including using interval scheduling algorithms to yield the maximum number of non-overlapping PG4s: [https://www.github.com/mparker2/g4predict](https://www.github.com/mparker2/g4predict)
* g4seeqer: DNA and RNA G-Quadruplex prediction using an artificial neural network, plus Cython implementation of g4hunter and tools for interpretation of neural network scores using mutation mapping: [https://www.github.com/mparker2/g4seeqer](https://www.github.com/mparker2/g4seeqer)
* g4netx: A package which identifies all possible overlapping registers of Quadparser PG4s using network analysis: [https://www.github.com/mparker2/g4predict](https://www.github.com/mparker2/g4predict)
* matplotlib_logo: Plotting package for making sequence logos: [https://www.github.com/mparker2/matplotlib_logo](https://www.github.com/mparker2/matplotlib_logo)

### Chapter 2: A Recurrent Neural Network to Predict G Quadruplex Structure

1. Training of the DNA G4Seeqer neural network on g4seq data, plus testing on hold out data and BG4 data. Comparisons to Quadron, G4Hunter and neural net trained on trinucleotide content: [01_g4seeqer_model_training_and_validation.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_2/01_g4seeqer_model_training_and_validation.ipynb)
* Training of RNA G4Seeqer (rG4Seeqer) model on RT pausing data, plus validation on G4RNA database: [02_rg4seeqer_model_training_and_validation.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_2/02_rg4seeqer_model_training_and_validation.ipynb)
* Mutation mapping analysis and interpretation of G4Seeqer results, including triplex/hairpin analyses, loop length analysis, and comparison of G4seeqer scores to UV melting temperatures: [03_mutation_mapping_and_triplex_hairpin.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_2/03_mutation_mapping_and_triplex_hairpin.ipynb)
* Analysis suggesting the applicability of the G4seeqer model trained on human genome to other genomes with different G4 content: [04_human_mouse_pg4_overlap.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_2/04_human_mouse_pg4_overlap.ipynb)

### Chapter 3: Global analysis of Predicted G-Quadruplexes in the Arabidopsis thaliana genome

1. Analysis of PG4 densities in genomes from ENSEMBL plants, plus: [01_plant_g4_density_genome_size.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_3/01_plant_g4_density_genome_size.ipynb)
* Plotting metagene profiles of two and three tetrad PG4 density across Arabidopsis genes: [02_g4_density_metagene_profiles.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_3/02_g4_density_metagene_profiles.ipynb)
* Reverse translation analysis of PG4 enrichment through codon selection in Arabidopsis CDS regions: [03_reverse_translation_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_3/03_reverse_translation_analysis.ipynb)
* Analysis of codons and protein motifs which contribute to the formation of PG4s in Arabidopsis CDS regions, plus analysis of tRNA adaptation index: [04_cds_codon_usage_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_3/04_cds_codon_usage_analysis.ipynb)

### Chapter 4: Global effect of G Quadruplex stabilisation on gene expression

1. Differential Expression analysis of NMM vs DMSO microarray dataset: [01_nmm_v_dmso_microarray_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/01_nmm_v_dmso_microarray_analysis.ipynb)
* Differential Expression analysis of drought stress dataset: [02_drought_stress_microarray_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/02_drought_stress_microarray_analysis.ipynb)
* Differential Expression analysis of DNA damage response dataset: [03_dna_damage_response_microarray_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/03_dna_damage_response_microarray_analysis.ipynb)
* Correlation of NMM and Berberine Differential Expression with PG4 densities: [04_pg4_prediction_and_nmm_correlation.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/04_pg4_prediction_and_nmm_correlation.ipynb)
* Self Organising Map analysis to identify PG4 types enriched in downregulated genes: [05_self_organising_map_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/05_self_organising_map_analysis.ipynb)
* Normalisation of PolII ChIP-chip dataset: [06_polii_chip_microarray_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/06_polii_chip_microarray_analysis.ipynb)
* Pol II occupancy bootstrap metagene plots: [07_polii_occupancy_metagene_profiles.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/07_polii_occupancy_metagene_profiles.ipynb)
* Normalisation of GRO/RNAseq counts and correlation with PG4 density: [08_groseq_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/08_groseq_analysis.ipynb)
* Differential Expression analysis of Maize NMM vs DMSO dataset and overlap with NMM regulated Arabidopsis orthologs: [09_zea_mays_nmm_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_4/09_zea_mays_nmm_analysis.ipynb)

### Chapter 5: Effect of G-Quadruplexes on expression and splicing of the Extensin gene family

1. Differential Expression analysis of NMM vs DMSO RNAseq experiment, plus GOseq Gene Ontology analysis: [01_nmm_dmso_rnaseq_differential_expression.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/01_nmm_dmso_rnaseq_differential_expression.ipynb)
* PG4 enrichment in Gene Ontology groups and comparison with GOseq results from NMM vs DMSO: [02_gene_ontology_pg4_enrichment.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/02_gene_ontology_pg4_enrichment.ipynb)
* Predictions of PG4s in Extensin gene family: [03_extensin_gene_pg4_estimate_table.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/03_extensin_gene_pg4_estimate_table.ipynb)
* Analysis of LRX1 and EXT3 qPCR data: [04_extensin_qpcr.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/04_extensin_qpcr.ipynb)
* Analysis of EXT repeat CD spec data: [05_ext_repeat_cd_spectroscopy.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/05_ext_repeat_cd_spectroscopy.ipynb)
* Visualisation of splice splice consensus sequences using weblogos: [06_exitron_splice_junction_logos.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/06_exitron_splice_junction_logos.ipynb)
* Simulation of expected number of mismapped splice sites using polyester simulated RNAseq data: [07a_splice_simulation.py](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/07a_splice_simulation.ipynb); [07b_splicing_simulation.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/07b_splicing_simulation.ipynb)
* Dotplot for EXT9 and mappability vs spliced read analyses for Extensins: [08_ext9_dotplot_and_mappability.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/08_ext9_dotplot_and_mappability.ipynb)
* Mapping of Sanger sequenced products to TAIR10 genome: [09_sanger_sequencing_splice_variants.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/09_sanger_sequencing_splice_variants.ipynb)
* Limma voom analysis of differential junction usage:[10_differential_splicing_analysis.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/10_differential_splicing_analysis.ipynb)
* Analysis of percent spliced reads in Extensins in NMM vs Control:[11_extensin_percent_spliced_nmm_vs_dmso.ipynb](https://www.github.com/mparker2/mparker_phd_thesis/tree/master/appendix/chapter_5/11_extensin_percent_spliced_nmm_vs_dmso.ipynb)

# Global effect of G Quadruplex stabilisation on gene expression
\label{chap:global_nmm}

## Introduction

As was shown in Chapter 4, the distribution of PG4s around and within genic regions is not uniform. In the *Arabidopsis thaliana* genome, template stranded PG4s are enriched in the 5' UTRs and promoter proximal regions, whilst coding strand PG4s are enriched in 3' UTRs and promoter distal regions. These enrichments do not appear to be explained simply by the GC content of these sequences, nor by the requirement to code for particular amino acids. These features may be deliberately conserved, indicating a biological function for PG4s within gene bodies. As G4s are formed from single stranded DNA, it has been previously hypothesised that coding strand G4s might function to promote transcription by competing with double stranded DNA to produce regions of open chromatin which could easily become transcription bubbles [@Rhodes2015]. G4s in the coding strand also have an opportunity to form in mRNA and regulate stability, splicing or translation. Template stranded G4s which occur downstream of translocating RNA Polymerase II (Pol II), on the other hand, might cause blockages which prevent elongation, causing downregulation of the gene [@Rhodes2015].

Arguably the best *in vivo* evidence for G4 formation was conducted using an antibody raised against G4 structures, referred to as BG4. BG4 was used to visualise G4s in human cancer cells by immunofluorescence [@Biffi2013]. G4 foci were identified at both the telomeres (which are highly GC rich and have been shown to form G4s *in vitro*), and in interstitial regions which contain actively transcribed euchromatin. The G4 density of the cells was also seen to fluctuate throughout the cell cycle, with the greatest number of foci appearing during S phase, when the DNA is decondensed to allow replication to occur.

The BG4 antibody has been further used to conduct ChIPseq experiments in human cell lines: the DNA fragments to which the antibody binds were enriched and subsequently sequenced [@Hansel2016]. BG4 ChIPseq peaks were found to overlap with regions of DNAse sensitvitiy, wich are thought ot correspond to regions of open chromatin. These regions are commonly found around the transcriptional start sites of transcripts and are associated with active transcription and promoter proximal pausing. Interestingly, genes in which BG4 peaks strengthened after treatment with HDAC inhibitors (which cause relaxation of heterochromatin) also saw a corresponding increase in gene expression, suggesting that promoter G4 formation may have a positive impact upon gene expression [@Hansel2016].

Another method by which the effect of G4s on transcription has been studied is through stabilisation with G4 binding ligands. These ligands generally bind to G4s through external hydrophobic pi-stacking above a planar tetrad, or through intercalation between the inner faces of tetrads.  Treatment of yeast species *Saccharomyces cerevisiae* with the G4 stacking ligand N-methyl-mesoporphyrin was shown to upregulate the expression of genes containing coding strand G4s in their promoters [@Hershman2008]. Rodriguez et al. showed that treatment of human cells with Pyridostatin resulted in DNA damage at PG4 containing regions in gene bodies, caused by both replication dependent and transcription dependent damage [@Rodriguez2012]. These genes were also downregulated in expression, suggesting that stabilised G4s arrest Pol II. G4s have also been shown to cause pausing of polymerases *in vitro*, in the presence and absence of G4 binding ligands [@Han1999; @Siddiqui-Jain2002; @Dexheimer2006; @Cogoi2006; @Chambers2015; @Kwok2016].

Here we conduct a global study of G4 stabilisation in the model plant *Arabidopsis thaliana*, using the G4 binding ligand NMM. Previous studies have shown that treatment of Arabidopsis seedlings with NMM cause developmental defects [@Nakagawa2012], suggesting an effect on gene expression. We also identify alterations in Pol II occupancy potentially caused by G4 dense transcribed regions of genes. Interestingly, Mullen et al. showed that genes involved in response to drought stress tended to be more likely to contain PG4s in their gene bodies [@Mullen2012]. Since G4 formation is dependent upon potassium cation concentration, the intracellular concentration of which increases during drought stress, this could constitute a regulatory mechanism of G4s. We also investigate the overlap of NMM regulated genes with drought stress responsive genes.

\newpage

## Methods

### Plant Growth Conditions

Plants for microarray analysis were grown and harvested by Manoj Valluru. For N-Methyl Mesoporphyrin (NMM) (Frontier Scientific, NMM580) treatment microarray experiments, the *Arabidopsis thaliana* Columbia (Col-0) ecotype was used. Seeds were surface sterilised, stratified for 2-3 days at 4°C and sown on vertical plates containing Murashige & Skoog (MS) agar with 0.5% sucrose and 0.8% agar. Plants were then transferred to growth cabinets with 16 hours light at 23°C.

### RNA Extraction & Microarray data generation

RNA extraction was performed by Manoj Valluru. Total nucleic acid isolation protocol was carried out by phenol-chloroform extraction as described by White and Kaper [@White1989]. The resulting pellets were resuspended in sterile water and stored at -80C. The RNA concentration and quality was checked using the NanoDrop 1000 Spectrophotometer (ThermoScientific).

cDNA library preparation and microarray analysis were performed by the Genomics Core Facility in the Sheffield Institute for Translational Neuroscience. An Affymetrix Arabidopsis Gene 1.0 ST array was used. RNA integrity and abundance was measured using an Aligent Bioanalyser 2100. Hybridization and scanning procedures were conducted according to the manufacturer using the Affymetrix Gene Chip hybridisation system.

### Microarray Analysis

NMM Microarray analysis was conducted in R version 3.3.2 using the packages `oligo` version 1.36.1 and `puma` version 3.14.0 [@Carvalho2010; @Pearson2009; @Liu2013]. `oligo/puma` was chosen over `oligo/limma` [@Ritchie2015] for this analysis as NMM treatment appeared to cause large consistent changes in gene expression which violated the assumptions used in Robust Multi-chip Averaging (RMA), namely that most genes do not change their expression and that there is no correlation between average expression and log fold change. As a result of this violation no genes were found to be differentially expressed at FDR threshold of 0.1 when using `limma` after RMA normalisation. Instead, CEL files were read into R using `oligo` but were not normalised using RMA. The `puma` Bayesian probabilistic method was used to normalise data and conduct differential expression analysis. `puma` Probability of Positive Log Ratio (PPLR) values were calculated for each contrast. PPLR tends to zero for genes which have a negative change in expression and to one for genes which have a positive change. Strongly differentially expressed genes were produced using an absolute log2 fold change threshold of 1 and a PPLR of 0.05 (or 0.95 for positively differentially expressed genes), and moderately differentially expressed genes were produced using an absolute log2 fold change threshold of 0.5 and a PPLR of 0.05 (or 0.95 for positively differentially expressed genes). Annotation of microarray data was conducted using the `oligo getNetAffx` function and TAIR10 AGI ids were extracted. Code used for normalisation and differential expression analysis is available in Appendix \ref{nmm_v_dmso_microarray_analysis}.

### Analysis of previously published microarray data

Processed Berberine expression data was downloaded from the supplementary material of Nakagawa et al. 2012. Data was generated using an Affymetrix ATH1 GeneChip array from Col-0 plants grown on MS media containing 12.5μM Berberine for 14 days.

Drought stress microarrays were downloaded from GSE65414 [@Linster2015]. Drought stressed plants were grown on soil for 6 weeks under with 8 hour light period, with normal watering, followed by 10 days drought stress. Data was generated using the Affymetrix Arabidopsis Gene 1.1 ST Array.

Raw drought stress microarray data was processed in R version 3.4.1 using `oligo` version 1.42.0 and `limma` version 3.34.6 [@Carvalho2010; @Ritchie2015]. CEL files were read into R using `oligo` and quantile normalised & median polished using Robust Multi-chip averaging. Linear modelling was then performed using `limma`, and p values were adjusted for multiple testing using Benjamini Hochberg correction. Strongly differentially expressed genes were produced using an absolute log2 fold change threshold of 1 and a FDR of 0.05, and moderately differentially expressed genes were produced using an absolute log2 fold change threshold of 0.5 and a FDR of 0.05. Annotation of microarray data was conducted using the `oligo` `getNetAffx` function and Ensembl annotations were extracted. Code is available in Appendix \ref{drought_stress_microarray_analysis}.

DNA damage stress microarrays were downloaded from GSE23892 [@Bohmdorfer2011]. Plants were grown on half MS media in 24 hour light for 5 days. Gamma irradiated plants were treated with 100 Gy of ionising radiation at a dose rate of 21-34 Gy per minute. Tissue was harvested 1.5 hours after irradiation. Data was generated using an Affymetrix ATH1 Gene Array.

DNA damage stress data was processed in R version 3.4.1 using `affy` version 1.56.0 and `limma` version 3.34.6 [@Carvalho2010; @Ritchie2015]. CEL files were read into R using `affy` and quantile normalised & median polished using Robust Multi-chip averaging. Linear modelling and differentially expressed gene identification was performed using `limma` as described above for drought stress data. Code is available in Appendix \ref{dna_damage_response_microarray_analysis}.

### Genome and Annotations used

All further analyses were performed using the TAIR10 genome [@Initiative2000], downloaded in fasta format from arabidopsis.org, and the Araport11 genome annotation [@Cheng2017], downloaded in GTF format from araport.org. Annotations were filtered to obtain protein coding genes only using the `CGAT` version 0.2 `gtf2gtf` script [@Sims2014]. To obtain sets of genic features such as exons, introns, CDS and UTRs, CGAT `gtf2gtf` was also used. Overlapping exons from different isoforms of the same gene were flattened to produce non-overlapping exons, and bed files of exons, CDS, 5' UTRs and 3' UTRs were generated from these flattened exons using `awk`. Bed files of introns were created using CGAT `gtf2gtf` to generate exon "complementation". Bed files of whole gene bodies were generated using CGAT `gtf2gtf` to merge all intervals into a single interval spanning the entire gene.

### PG4 prediction

G Quadruplex predictions in the TAIR10 genome were carried out using an in-house script (g4predict) which utilises the Quadparser method [@Huppert2005]. Results were filtered using a dynamic programming approach, commonly used in interval scheduling [@Ray], to produce the greatest number of non-overlapping PG4s. Scripts can be found in Appendix \ref{g4predict}. To count G4s per gene, the bed files containing PG4s were overlapped with bed files generated from Araport11 for exon, intron, CDS, 5' UTR, 3' UTR and full gene bodies. This was done using `bedtools intersect` in count mode [@Quinlan2010]. `bedtools` version 2.27.1 was used. PG4s on the template and coding strands of gene features were counted separately. For multi-exon genes, counts for different exons were summed using `awk` scripts, and counts were normalised by length to get PG4 densities per kb. Barplots of average PG4 density for various gene features and gene sets were produced in python using `pandas` version 0.20.3 and `seaborn` version 0.8.1 [@Mckinney2011; @Waskom2014]. Errorbars for these plots are estimated 68% confidence intervals generated using 1000 bootstrapped samples. Statistical hypothesis testing was done using the Mann-Whitney U test. Code is available in Appendix \ref{pg4_prediction_and_nmm_correlation}.

Maximal PG4 densities were calculated using a sliding window of 200bp generated using `bedtools makewindows`, with a step size of 5bp [@Quinlan2010]. `bedtools` version 2.27.1 was used. `bedtools map` was used to count the number of PG4s overlapping each window. The score of the maximum scoring 200bp window overlapping the transcribed body (exons and introns) of a gene was assigned as the maximal PG4 density of the gene, using `bedtools map`. Coding and template strand densities were calculated separately. Pointplots of average expression change during NMM treatment for genesets with different maximal PG4 densities were produced in python using `pandas` and `seaborn` [@Mckinney2011; @Waskom2014]. Errorbars for these plots are estimated 68% confidence intervals generated using 1000 bootstrapped samples. Statistical hypothesis testing was done using the Mann-Whitney U test. Code is available in Appendix \ref{pg4_prediction_and_nmm_correlation}.

For analyses where G4seeqer was used, PG4 predictions were conducted on the TAIR10 genome using the G4seeqer command line tool. A step size of 5bp and G4Hunter threshold of 0.75 were used. All intervals tested using G4seeqer were output regardless of neural network score (i.e. a threshold of 0 was used). Maximum G4seeqer score overlapping each gene body (exons and introns) was calculated using `bedtools map` [@Quinlan2010]. Coding and template strand scores were calculated separately.

### Self Organising Map Analysis

Loop lengths and total loop lengths for each Quadparser 2 tetrad PG4 in the TAIR10 genome were extracted from bed files output by g4predict. PG4s which did not overlap with gene bodies were discarded. Self Organising Maps were trained in R version 3.4.1 using the package `kohonen` version 3.0.6 on loop length data [@Wehrens2007]. 36 clusters were used. To identify enrichment of specific clusters of PG4s in NMM downregulated genes, the total number of PG4s from each cluster overlapping the geneset was calculated, and compared to an expected number of overlaps computed by permuting PG4s amongst all genes. Genes were weighted by length such that a 2kb gene was twice as likely to be assigned PG4s as a 1kb gene. Coding and template strand PG4s were permuted separately. The log fold change between observed and expected overlap was then calculated for each cluster. SOM plots were made in Python  using `matplotlib` version 2.1.0 [@Hunter2007]. Code is available in Appendix \ref{self_organising_map_analysis}.

### Venn diagrams

Venn diagrams of geneset overlaps were produced in Python 3.5 using the package `matplotlib_venn` version 0.11.5. Statistical hypothesis testing of overlaps was conducted using hypergeometric tests with `scipy.stats` version 0.19.1 [@Jones2001].

### Pol II ChIP-tiling array analysis

RNA Polymerase II ChIP-tiling array data was downloaded from GEO Accession GSE21673 [@Chodavarapu2010]. This data was from plants grown under 24hr light on soil for 10-14 days before being harvested. The authors conducted the ChIP using Abcam ab817 Pol II antibody and the GEO entry lists the tiling arrays used as Affymetrix Arabidopsis Tiling 1.0R Array. Pol II occupancy tracks were generated from CEL files in R version 3.3.2 using `STARR` version 1.28.0 [@Zacher2010]. Cyclic lowess method was used for probe intensity normalisation. The enrichment ratio of PolII signal intensity over control was calculated and saved in BigWig format using `rtracklayer` version 1.32.2 [@Lawrence2009]. Code used to normalize ChIP-chip data is available in Appendix \ref{polii_chip_microarray_analysis}. Metagene profiles for all genes were produced using CGAT `bam2geneprofile` version 0.2 [@Sims2014]. Gene profiles of merged exons (without introns) were produced using 100 bins across the gene body, with an upstream and downstream extension of 500bp at 10bp resolution (i.e. binned in 10bp intervals).

To compared Pol II occupancy of G4 containing genes with non-G4 containing genes, genesets with max G4seeqer scores greater than 0.95 and less than 0.05, or maximal PG4 density greater than 2 or equal to zero were used. Metagene profile matrices were read into Python using `pandas` and averaged profiles for each geneset were generated using `seaborn` bootstrapping to estimate central tendency and confidence intervals [@Mckinney2011; @Waskom2014]. Bootstrapped profiles were smoothed using a moving average of 20. 1000 iterations were used for all bootstraps. Profiles were normalised so that the absolute area under the curve was equal to one. Code used to generate metagene profiles is available in Appendix \ref{polii_occupancy_metagene_profiles}.

### GRO/RNA seq analysis

Global Run On (GRO) and RNA sequencing data for 6-d-old Arabidopsis seedlings, grown at 22 °C on half Linsmaier and
Skoog (LS) medium in constant light, was downloaded from the European Nucleotide Archive (ENA), GEO accession GSE83108 [@Hetzel2016]. Quality control analyses were performed using `FastQC` version 0.11.3 and `fastq-screen` version 0.5.1 [@Andrews2010; @Wingett2011]. Mapping to the TAIR10 genome with splice junction annotations from Araport11 [@Initiative2000; @Cheng2017] was conducted using `STAR` version 2.4.2a [@Dobin2013] with default parameters, and generated BAM files were sorted using `samtools` version 1.4.1 [@Li2009]. Exonic read counts per gene were then counted using `featureCounts` version 1.6.1 [@Liao2013]. Read counts were normalised for library depth in R version 3.4.1 [@RCoreTeam2013] using `DESeq2` version 1.18.1 [@Love2014] and log2 transformed to get log counts per million (logCPM). The ratio of GROseq to RNAseq reads was calculated by subtracting the average RNAseq logCPM from the average GROseq logCPM. Scatter plots of RNAseq logCPM vs GROseq logCPM were generated in Python using `seaborn` [@Waskom2014].

To contrast GRO/RNA seq ratios of G4 containing genes with non-G4 containing genes, genesets with max G4seeqer scores greater than 0.95 and less than 0.05, or maximal PG4 density greater than 2 or equal to zero were used. Overlayed histograms and kernel density estimate plots of GRO/RNA seq ratio for these genesets were generated using `seaborn` version 0.8.1 [@Waskom2014]. Statistical hypothesis testing was conducted using Welch's unpaired T-test. Code used for all GROseq analysis is available in Appendix \ref{groseq_analysis}.

### *Zea mays* RNAseq data analysis.

*Z. mays* data was downloaded from the Europeaen Nucleotide Archive accession PRJEB23390 [@Tokan2018]. Maize seedlings were germinated for 5 days at room temperature then grown in 1/4 Reid-York solution for. After 1 day plants were treated with 16 μM NMM, and grown for a further 3 days. Root tissue was collected for RNAseq using an Illumina TruSeq stranded mRNA kit, on an Illumina NextSeq 500 with 80bp paired end reads [@Tokan2018].

Expression counts were generated at the transcript level using `salmon` version 0.8.2 [@Patro2017]. Library type was automatically determined as ISR (Inward facing pairs, stranded, read one on reverse strand). Mean insert size for each sample was assessed to be in the 120-140bp range. Gene level abundance was then aggregated from this using `tximport` version 1.6.0 [@Soneson2015] and differential expression testing was conducted using `edgeR` linear modelling (version 3.20.7) [@Robinson2010; @McCarthy2012]. An absolute log2 fold change threshold of 1 and an FDR of 0.05 was used to identify differentially expressed genes.

Arabidopsis  orthologs of *Z. mays* genes were downloaded from Ensembl plants using `biomaRt` version 2.34.2 [@Durinck2005; @Durinck2009]. Overlaps of orthologous NMM regulated genes from the *Z. mays* RNAseq dataset and the Arabidopsis microarray were produced in python using `pandas` and `matplotlib_venn`. Code used to generate figures is available in Appendix \ref{zea_mays_nmm_analysis}.

\newpage

## Results

#### NMM causes global change in gene expression

In order to test the effect of G4 stabilisation on gene expression in Arabidopsis, we conducted a microarray analysis using RNA from 7 day old seedlings treated with NMM for 6 hours. Control samples were treated with DMSO for 6 hours. Two biological replicates per condition were used for control samples, and three for treated samples. We found 1098 and 858 genes were differentially downregulated and upregulated respectively using a log fold change threshold of 1 (PPLR < 0.05). When a less stringent log fold change threshold of 0.5 was used (PPLR < 0.05), downregulated genes outnumbered upregulated by a ratio of 2:1 (3882 downregulated, 1930 upregulated), suggesting that NMM has a global effect on gene expression unlikely to be caused by a single transcription factor. An MA plot of average gene expression vs log fold change, with lowess curve fitting, showed that there appeared to be a slight skew towards downregulation in genes with higher average expression (Fig \ref{nmm}a). To more clearly visualise this skew we binned genes by average expression quartile (Fig \ref{nmm}b). This showed there was indeed a relationship between average expression and expression change upon NMM treatment, with highly expressed genes tending to be more downregulated by NMM treatment. This global pattern, which violates some of the assumptions that are usually used in microarray normalisation and analysis, potentially suggests a widespread effect of NMM directly upon either transcription or mRNA stability.

\newpage

![**Global effect of NMM on Gene Expression.** **a)** MA plot showing relationship between average gene expression and Log2 fold change in expression upon treatment with NMM. Orange line is lowess curve fit showing slight negative correlation between expression and fold change for genes in the expression range 4-8. **b)** Boxplot of Log2 fold change in expression upon treatment with NMM, cut on quartiles by average expression. Lowest expressed 25% is leftmost, and highest expressed is rightmost. Lower quartile, median, and upper quartile are at 4.6, 5.4, and 6.5, respectively. **c)** Histogram and kernel density estimate showing distribution of Probability of Positive Log Ratio (PPLR) values. PPLRs which tend towards zero represent negatively differentially expressed genes, whilst PPLRs which tend towards one represent positively differentially expressed genes. \label{nmm}](figures/nmm_vs_dmso.svg)

\newpage

To support the hypothesis that NMM alters gene expression through G4 stabilisation, we correlated our results with processed data from a Berberine treatment array (Fig \ref{berberine}a) [@Nakagawa2012]. Berberine is another G4 stabilising drug, but with a very different structure and method of action (intercalation with G4s rather than hydrophobic stacking). Despite the differences in structure of the two drugs, and the very different conditions (plants were grown on Berberine for 14 days, compared with 6 hour treatment of 7 day old seedlings with NMM), the log fold changes from our data correlated well with the Berberine dataset (Pearsons R: 0.43, Spearman’s ρ: 0.44). There was a strong overlap between the genes downregulated by NMM and those downregulated by Berberine (Fig \ref{berberine}b, p=1.1e-36). Both genesets show an greater exonic G4 density on the template strand than genes not regulated by either drug, however genes which are regulated by both drugs had the greatest average exonic PG4. These results suggest that the main effects on gene expression were through G4 interaction, with any off target effects being less significant contributors.

\newpage

![**Comparison of gene expression during NMM treatment with expression during Berberine treatment.** **a)** Scatter plot with regression line showing the correlation in expression change for NMM vs DMSO and Berberine vs Control. Processed Berberine data was taken from supplementary information of Nakagawa et al. 2012, however only differentially regulated genes were reported. **b)** Venn diagram reporting the overlap of genes downregulated by NMM with those downregulated by Berberine. **c)** Bar plot showing the average exonic PG4 densities of NMM and Berberine downregulated genesets, on the coding and template strands, respectively. Bar colours match set colours from Fig 3b. Errorbars are 68% confidence intervals for mean generated using 1000 bootstrapped samples. \label{berberine}](figures/nmm_berberine.svg)

\newpage

### Genes downregulated by NMM are enriched in two tetrad PG4s

We next investigated whether NMM regulated genes were enriched for PG4s. We first measured the density of PG4s in different regions of each gene (i.e. promoters, exons, introns, CDS and UTRs) and looked at the differences between the differentially expressed gene sets at 6 hours NMM treatment. We did not see a strong difference in three tetrad PG4 density in our gene sets, for any gene feature (data not shown). This is likely due to the low density of these PG4s in Arabidopsis. We discovered an enrichment of the 2 tetrad PG4s on the template strand of genes which were down-regulated by NMM, with approximately 10% more genes containing template PG4s than expected for a gene set of that size (p=2e-48) (Fig \ref{nmm_g4_sq}). This enrichment occurred most specifically in the CDS and 5' UTR regions of genes (Fig \ref{nmm_g4}). Genes which were very strongly downregulated by NMM (logFC < -1) tended to contain large numbers of PG4s throughout their exonic bodies, particularly in coding regions and in the 5' UTR, whilst moderately downregulated genes (logFC < -0.5) tended to have greater concentration of PG4s in 5' UTRs.

\newpage

![**NMM downregulated genes are enriched in PG4s** Heatmaps showing fractions of genes containing at least one predicted G4 in their gene body for upregulated genes vs all genes (top row) and downregulated genes vs all genes (bottom row) (For down and upregulated genes, FDR < 0.05 and absolute logFC > 0.5). PG4 predictions for the coding strand are in the left hand column whilst PG4 predictions for template strand are on the right. P values for each heatmap are calculated using Chi-squared tests. Genes downregulated by NMM show a particularly strong enrichment of PG4s, and particularly on the template strand. \label{nmm_g4_sq}](figures/nmm_g4_presence_absence_chisquared.svg)

\newpage

![**Distribution of PG4s in genes differentially regulated by NMM.** Bar plots showing the average PG4 densities of genes up or downregulated by NMM, for **a)** full gene body (exons and introns), **b)** coding regions, **c)** 5' UTR, **d)** 3' UTR, and **e)** introns, respectively. In each figure, left and right panels represent coding and template strand, respectively. Genesets are separated into three categories by strength of regulation: green: not differentially expressed, blue: moderately differentially expressed (PPLR < 0.05, logFC > 0.5), orange: strongly differentially expressed (PPLR < 0.05, logFC > 1). Errorbars are 68% confidence intervals for mean generated using 1000 bootstrapped samples. Genes which are strongly downregulated by NMM tend to have higher PG4 densities on the template strand of coding regions and 5' UTRs, whilst moderately downregulated genes tend to have greater PG4 density on the template strand of their 5' UTRs. \label{nmm_g4}](figures/nmm_regulated_g4_distribution.svg)

\newpage

The PG4 density of genes downregulated by both NMM and Berberine was also calculated (Fig \ref{berberine}c). We found that whilst the gene sets downregulated by either treatment were enriched in PG4s, those which were in the intersection of the two sets had the greatest average exonic PG4 density. This is further evidence that these drugs are regulating gene expression through G4 stabilisation.

Previous studies have shown that NMM is highly selective towards parallel G4s [@Nicoludis2012a; @Nicoludis2012; @Tippana2014; @Kreig2015], and can induce anti-parallel G4s to switch structure [@Nicoludis2012a]. G4s with short loop lengths are more likely to form parallel structures [@Tippana2014]. In order to test if NMM was selective towards G4s with particular loop lengths, we used Self Organising Maps to cluster all predicted Arabidopsis 2 tetrad PG4s into 36 groups, based on the length of each loop 1-3, and the total loop length (Fig \ref{som}b). Each cluster contained between approximately 1000 and 8000 PG4s (Fig \ref{som}a). We then analysed the relative enrichment of each PG4 cluster on each strand, within genes which were downregulated by NMM, compared to the expected number of PG4s from that cluster if PG4s were uniformly distributed across all genes. No particular PG4 cluster was strongly enriched on the coding strand of downregulated genes (Fig \ref{som}c). One cluster, however, was strongly enriched on the template strand (Fig \ref{som}d). This cluster contained PG4s with very short loop lengths of 1-2bp and a total loop length of 5-6bp. This conformed well with our prior expectations, as G4s with short loop lengths are known to form propeller-like parallel G4s of the kind favoured by NMM (Fig \ref{g4_struct}c) [@Nicoludis2012a].

\newpage

![**Self Organising Maps demonstrate NMM downregulated genes contain specific PG4 types.** Self Organising Map (SOM) plots for clustering of Quadparser predicted two tetrad PG4s by loop length. In each figure, each circle represents a cluster of similar PG4s. **a)** SOM plot coloured by cluster size. Each cluster contained between 1000 and 8000 PG4s. **b)** SOM plot coloured by total length in bases of all loops. **c)** and **d)** Log2 fold enrichment of each cluster in the gene bodies of genes downregulated by NMM, on the coding and template strands, respectively. Log2 fold enrichments were generated by comparing actual overlap of PG4s with downregulated genes, with expected overlap when PG4s were permuted amongst all genes. \label{som}](figures/som_enrichment.svg)

\newpage

### NMM downregulation is correlated with maximal G4 density

Previous studies have suggested that G4s cause the largest effect on gene expression when grouped in clusters [@Yoshida2018]. This may be due to an increase in the likelihood of a single G4 being formed at any one time, or through increased polymorphism of G4 formation. To identify whether NMM regulated genes tend to contain G4 clusters, we used a sliding window of 200bp to count two tetrad G4 density across the whole transcribed region of each gene, including introns. Each gene was then assigned the maximum density score for the gene. Genes were then binned by their maximal density, and expression change after NMM treatment was calculated (Fig. \ref{max_g4}). We found that genes with higher maximal G4 density tended to have more negative log fold changes during NMM treatment. This suggests that clusters of G4s do have a stronger effect on gene expression. It is possible that a single region of high G4 density may be sufficient to cause downregulation of an otherwise G4 free gene during NMM treatment.

\newpage

![**NMM regulates genes with large maximal PG4 density** Point plot showing mean log fold change in gene expression during NMM treatment for genes binned by "maximal PG4 density", defined as the greatest concentration of PG4 motifs within a 200bp sliding window anywhere in the body of the gene (i.e. exon or intron). Left and right panels depict coding and template strands, respectively. Errorbars are 68% confidence intervals for mean generated using 1000 bootstrapped samples. \label{max_g4}](figures/nmm_exprs_g4_max_density.svg)

\newpage

### G Quadruplexes may cause downregulation by Pol II stalling

Since transcriptional downregulation by NMM stabilised G4s appears to occur most strongly for genes which have many G4s in the gene body, and because this effect is specific to the template strand, we hypothesised that this could be a result of RNA Polymerase (Pol II) stalling. Pol II is the RNA polymerase which is responsible for transcribing all protein coding genes. The Pol II complex scans along the template strand and uses complementary base pairing to produce an mRNA copy which corresponds to the sequence of the coding strand. Since only the template strand is read directly, this might explain why coding strand G4s do not cause downregulation, since they do not form blockages which prevent the elongation of Pol II. To test whether G4s cause blockages which slow or stall the elongation of Pol II in the abscence of stabilisation by NMM, we reanalysed Pol II ChIP tiling array data [@Chodavarapu2010]. Metagene profiles of the transcriptional start and end sites showed an accumulation of Pol II at the transcriptional termination site (TTS) (Fig \ref{pol_ii}, orange profiles). This was surprising as it is in disagreement with Pol II occupancy profiles in human cell lines, where there is generally a much larger peak of paused Pol II at the start of the gene. We contrasted this result with occupancy profiles for G4 dense genes (Fig \ref{pol_ii}, blue profiles). For genes which contained at least one G4 dense region, measured either by G4Seeqer score greater than 0.95, or by maximum Quadparser G4 density per 200bp of greater than 3, we found that there was greater Pol II occupancy at the TSS and in the TSS proximal part of the gene body. Greater Pol II occupancy can be explained by either of two factors: increased initiation and transcription in G4 dense genes, or slower Pol II elongation. Since G4 dense genes do not have higher expression than non-G4 containing genes, we suggest that template strand G4s cause a reduction in Pol II speed. This may result in abortive transcription or alteration of co-transcriptional processes such as splicing.

\newpage

![**PG4 dense genes have altered RNA Polymerase II occupancy** Metagene profiles showing RNA Polymerase II (Pol II) occupancy across binned exonic gene bodies. Profiles are made up of 500bp upstream region (in 10bp bins), 100 gene body bins, and 500bp downstream region (in 10bp bins). **a)** Metagene profile for genes containing PG4 predicted by G4Seeqer (max G4seeqer score > 0.9), vs. genes containing no G4s (max G4seeqer score < 0.1). **b)** Metagene profile for genes containing two tetrad maximal PG4 density per 200bp of 3 or greater, vs. genes with maximal PG4 density of zero (contain no PG4s). Averages were generated using 1000 bootstrapped samples from each geneset. Bootstrapped samples are normalised such that the absolute area under the curve was equal to one. Errorbars are 68% confidence intervals for bootstrapped means. \label{pol_ii}](figures/polii_occupancy.svg)

\newpage

### G Quadruplex stalling may result in abortive transcription

To test whether template G4 containing genes might produce abortive transcripts that are degraded by the exosome rather than maturing to mRNAs, we analysed publicly available GRO-seq data [@Hetzel2016]. Since GRO-seq captures nascent RNA irrespective of its stability, and RNAseq captures stable RNAs, the ratio between the normalised read counts in GROseq vs RNAseq represents an estimate of the amount of unstable products produced at each gene locus (Fig. \ref{gro}a). We found that the largest difference in ratio was between non-coding and protein coding RNAs, with non-coding RNAs having much greater GRO/RNA ratios (data not shown). This is likely explained by the higher rate of modification of many ncRNAs, e.g. tRNAs, which prevent reverse transcription and sequencing by RNAseq. Other ncRNAs such as natural antisense RNAs may also be unstable and degraded quickly.

G4 predictions were calculated using G4Seeqer and the score of the maximum score region within the transcribed body of each gene was assigned to the gene. A G4 containing set was produced using genes which contained a maximum template strand G4seeqer score of more than 0.95, and a G4 negative set was produced using genes with a maximum score of only 0.05 or less. We found a small but significant positive increase in GRO/RNA ratio for G4 dense genes (t-test p = 0.009), suggesting that some abortive transcripts are produced from these genes (Fig. \ref{gro}b, right panel). In contrast, genes with high scoring G4 regions on the coding strand did not have greater GRO/RNA ratios (p=0.4) (Fig. \ref{gro}b, left panel).

We conducted the same analysis using genes containing three or more two tetrad PG4s in 200bp maximal density clusters, as described above (Fig. \ref{gro}c). For these genes, we did not find any significant increase in GRO/RNA ratio, suggesting that these two tetrad G4s are not sufficiently stable to cause abortive transcription. Since these genes do have higher promoter proximal Pol II occupancy, we suggest that elongation occurs more slowly across the genes, but does not cause premature termination.

\newpage

![**Log Ratio of GRO / RNA seq counts per million detects abortive transcription of PG4 dense genes.** **a)** Scatter plot showing measured expression in log2 counts per million (logCPM) for each gene from GRO-seq vs. RNAseq datasets. Genes which fall below and to the right of the orange line have positive log2 GRO/RNA ratios. **b)**  Histogram and kernel density estimates of GRO/RNA ratio for genes containing PG4 predicted by G4Seeqer (max G4seeqer score > 0.95) in orange, vs. genes containing no G4s (max G4seeqer score < 0.1) in blue. Left and right panels represent coding and template strand G4 rich genes, respectively. **c)** Histogram and kernel density estimates of GRO/RNA ratio for genes containing two tetrad maximal PG4 density per 200bp of 3 or greater in orange, vs. genes with maximal PG4 density of 0 (contain no PG4s) in blue. Left and right panels represent coding and template strand G4 rich genes, respectively. \label{gro}](figures/gro_g4s.svg)

\newpage

### NMM treatment does not activate DNA damage response pathways

Rodriguez et al. showed that treatment of human cells with Pyridostatin (PDS), a G4 binding agent, caused an increase in the histone marker γH2AX, which is associated with double stranded DNA breaks [@Rodriguez2012]. Cells in gap phases 1 & 2, during which DNA is not being replicated, showed a reduction in double strand breaks when also treated with 5,6-dichloro-1-β-D-ribofuranosylbenzimidazole, an inhibitor of transcription. This suggests that PDS causes transcription and replication dependent DNA damage, presumably resulting from Pol II and replication fork stalling at G4 loci. To see if NMM causes similar DNA damage in *Arabidopsis thaliana*, we looked for transcriptional signatures of damage response. A set of 7 genes: AGO2, PARP1, RPA1E, BRCA1, GRG, RAD51, and RAD17, which have been reported as strong transcriptional markers of DNA damage [@Ryu2018] were tested. Only AGO2 showed upregulation (logFC=1.23, PPLR=1) upon NMM treatment (Fig. \ref{damage}a). We next compared data from a microarray in which seedlings were treated with Gamma Irradiation so see if there was any overlap with genes regulated by NMM. We found a small but significant overlap in upregulated genes (P = 2.9e-29). Gene ontology analysis of the overlapping gene set yielded no clear results, however. In comparison, the total gamma irradiation upregulated set was strongly enriched for genes involved in double strand break repair (FDR = 9.68e-4). Taken together, this data suggests that NMM treatment does not cause major DNA instability or DNA damage response.

\newpage

![**NMM does not activate DNA damage response** **a)** Barplot showing normalised log2 expression of seven DNA damage responsive genes: AGO2, PARP1, RPA1E, BRCA1, GRG, RAD51, and RAD17, in the presence and absence of NMM. Errorbars are standard errors calculated by `puma`. Asterisks denote significant results. **b)** Venn diagram showing the overlap between NMM upregulated genes and genes upregulated in response to gamma irradiation. \label{damage}](figures/dna_damage_response.svg){height=750px}

\newpage

### G4 dense genes are modulated by environmental stress and correlate with NMM treatment

Because G4s require potassium cations or other divalent cations for formation, it has been noted by previous studies that they might form more readily under stress conditions in plants such as drought stress, when the intracellular concentrations of such ions is increased [@Mullen2010; @Mullen2012]. Mullen et al. showed that genes involved in response to drought stress tended to be more likely to contain G-Quadruplexes in their gene bodies. To more closely examine this hypothesis we reanalysed microarrays conducted using RNA from drought stressed plants [@Linster2015] to determine whether differentially expressed genes contained enrichments or depletion of PG4 structures within various gene regions. As with experiments conducted on the NMM microarray, genes were considered to be moderately differentially expressed if they underwent a log change in expression of greater than 0.5 fold (FDR < 0.05). Genes which changed in expression by more than 1 log fold (FDR < 0.05) were considered strongly differentially expressed. 2947 and 491 genes showed moderate or strong upregulation, respectively, and 2572 and 984 genes showed moderate or strong downregulation. These gene sets were then analysed for two tetrad G4 density using the Quadparser method. As with NMM downregulated genes, we found that genes downregulated during drought stress tended to have greater numbers of template strand PG4s in exonic regions. Around 8% more genes than expected contained at least one template PG4 in the gene body (p=1e-11) (Fig \ref{drought_g4_sq}). This enrichment appeared specifically in CDS and 5' UTR regions (Fig. \ref{drought_g4}b-c). Interestingly, we also found that genes which were upregulated during drought stress were more likely to have PG4s in the coding strand of their 3' UTRs (Fig. \ref{drought_g4}d).

\newpage

![**Drought downregulated genes are enriched in PG4s** Heatmaps showing fractions of genes containing at least one predicted G4 in their gene body for drought upregulated genes vs all genes (top row) and drought downregulated genes vs all genes (bottom row) (For down and upregulated genes, FDR < 0.05 and absolute logFC > 0.5). PG4 predictions for the coding strand are in the left hand column whilst PG4 predictions for template strand are on the right. P values for each heatmap are calculated using Chi-squared tests. Genes downregulated by drought stress show a particularly strong enrichment of PG4s, and particularly on the template strand. \label{drought_g4_sq}](figures/drought_g4_presence_absence_chisquared.svg){height=750px}

\newpage

![**Distribution of PG4s in genes differentially regulated by Drought stress.** Bar plots showing the average PG4 densities of genes up or downregulated by Drought stress, for **a)** full gene body (exons and introns), **b)** coding regions, **c)** 5' UTR, **d)** 3' UTR, and **e)** introns, respectively. In each figure, left and right panels represent coding and template strand, respectively. Genesets are separated into three categories by strength of regulation: green: not differentially expressed, blue: moderately differentially expressed (PPLR < 0.05, logFC > 0.5), orange: strongly differentially expressed (PPLR < 0.05, logFC > 1). Errorbars are 68% confidence intervals for mean generated using 1000 bootstrapped samples. Genes which are downregulated by drought stress tend to have higher PG4 densities on the template strand of coding regions and 5' UTRs. Genes which are upregulated by drought stress are more likely to contain PG4s in their 3' UTRs. \label{drought_g4}](figures/drought_regulated_g4_distribution.svg){height=750px}

\newpage

Next we investigated whether there were any similarities in the expression profiles of NMM treated seedlings and drought stressed plants. We found a strong overlap between genes moderately downregulated by NMM and those moderately downregulated by drought stress (p = 1.7e-196) (Fig \ref{nmm_drought}a). Analysis of the PG4 density of these gene sets and the overlap between them showed that the genes which were downregulated in both experiments tended to be the most PG4 rich, particularly in the 5' UTR of the gene (Fig \ref{nmm_drought}b). This suggests that these genes could indeed be regulated through the same mechanism of G4 stabilisation.

\newpage

![**Overlap of genes downregulated by NMM with those downregulated by Drought stress.** **a)** Venn diagram reporting the overlap of genes downregulated by NMM with those downregulated by drought stress. **c)** Bar plot showing the average PG4 densities in 5' UTRs of NMM and drought stress downregulated genesets. Left and right panels show densities on the coding and template strands, respectively. Both  genesets show an greater exonic G4 density on the template strand than genes not regulated by either drug, however genes which are regulated by both drugs had the greatest average PG4 density in 5' UTRs. Bar colours match set colours from Fig 3b. Errorbars are 68% confidence intervals for mean generated using 1000 bootstrapped samples. \label{nmm_drought}](figures/nmm_drought_venn.svg){height=750px}

\newpage

### NMM affects gene expression in *Zea mays*

A recently published study by Tokan et al. provides the first analysis of differential expression during NMM treatment in a crop plant, *Zea mays* [@Tokan2018]. Their paper focussed on the effect of NMM on the expression of transposable elements, however, we reanalysed their data to identify whether protein-coding genes in maize respond to NMM similarly to Arabidopsis. We found a strong overlap between genes differentially regulated in Arabidopsis and their orthologs in maize, in both the up- and downregulated genesets (Fig. \ref{maize_venn}a and b respectively). Furthermore, we found that maize genes which were differentially regulated by NMM tended to have greater two tetrad PG4 densities in their CDSs. Unlike Arabidopsis, however, this increase in PG4 density was not restricted to the template strand, and not only in downregulated genes. Genes which are upregulated by NMM in maize are also more likely to contain CDS PG4s on both the coding and template strands (Fig \ref{maize_pg4}a). This suggests some other mechanism of G4 regulation beyond Pol II stalling. Whilst 5' UTRs of maize genes tend to have high densities of PG4s on the template strand (Fig \ref{maize_pg4}b, right panel), we did not see a greater level of PG4s in either up- or downregulated sets. In fact, downregulated genes appear to have fewer template strand PG4s in the 5'UTR (Fig \ref{maize_pg4}b, right panel). Interestingly, upregulated genes also tended to have greater levels of coding and template PG4s in introns (Fig \ref{maize_pg4}c).

\newpage

![**NMM expression change of Arabidopsis and *Zea mays* orthologs** Venn diagrams showing the overlap between the set of genes **a)** upregulated and **b)** downregulated by NMM in Arabidopsis and the unique set of Arabidopsis orthologs for *Z. mays* genes regulated by NMM. P values were calculated by hypergeometric test compared to the unique set of *Arabidopsis* genes with at least one *Z. mays* ortholog. \label{maize_venn}](figures/zmays_athaliana_nmm_expression_venn.svg){height=750px}

\newpage

![**Distribution of Two tetrad PG4s in *Z. mays* genes regulated by NMM** Bar plots showing the average PG4 densities of *Z. mays* genes up or downregulated by NMM, for **a)** CDSs, **b)** 5' UTRs and **c)** introns, respectively. In each figure, left and right panels represent coding and template strand, respectively. Errorbars are 68% confidence intervals for bootstrapped means. \label{maize_pg4}](figures/zmays_nmm_g4_distribution.svg){height=750px}

\newpage

## Discussion

We have conducted the first in detail analysis of the effects of G4 stabilisation in the higher plant, *Arabidopsis thaliana*. To determine how gene expression is altered by G4 stabilisation, we carried out a microarray analysis of plants treated with the G4 binding drug, NMM. NMM treatment had very strong effects on gene expression, particularly in genes which contained large numbers of parallel two tetrad G4s in the transcribed gene body. Moderately downregulated genes had large numbers of G4s in their 5' UTRs, whilst very strongly downregulated genes contained large numbers of G4s throughout the exonic regions of the gene. This is contrary to evidence from human systems which has suggested G4s influence transcription most when located in the promoter region. Furthermore, we find that the effect on gene expression is entirely strand dependent: only template stranded G4s strongly affect expression. Many previous studies have suggested that two tetrad G4s are not biologically relevant due to their relative instability compared to three tetrad G4s. The majority of these studies have been focussed upon their role in human biology. However, several other studies have shown that two tetrad G4s are stable *in vitro* and at physiological temperatures. In organisms that exist at lower temperatures it is credible that two tetrad G4s may play a role. Indeed, here we find that the strongest effect of NMM treatment is on genes predicted to form only two tetrad G4s.

Since the effect of NMM on gene expression mainly appears to be confined to genes with template stranded G4s, G4s will not be present in the mRNA of downregulated genes. Any direct binding of NMM to G4s must therefore occur in the DNA. The most likely explanation for a template stranded effect is that stabilised G4s interact with the elongating Pol II, which uses the non-coding strand of genes as a template for RNA polymerisation. Since G4s have previously been shown to cause polymerase stalling in vitro, we investigated the Pol II occupancy profiles of G4 containing genes. This data showed that G4 containing genes had higher Pol II density at the TSS proximal end of the gene, and lower density at the TTS. An increase in Pol II occupancy could be explained by one of two factors: more Pol II molecules binding and initiating transcription; or a reduction in Pol II speed. We suggest that G4s in the template strand block the elongation process and cause Pol II to slow down, resulting in a higher Pol II occupancy upstream of the G4 dense region. The Pol II data analysed is captured in the absence of any G4 stabilising drugs, indicating that G4 dependent Pol II slowing is a commonly occurring phenomenon. We hypothesise that two tetrad G4s form naturally in genes, causing Pol II slowing, but still creating full length products. Changes in Pol II speed may have consequential effects for co-transcriptional processes such as mRNA splicing. In some cases blockages may cause premature termination, resulting in truncated products. Analysis of the ratio of GRO/RNAseq read counts suggests that G4 dense genes do indeed produce slightly more unstable products than genes containing no G4s. In the presence of NMM, however, stabilised G4s are likely to become too difficult for the transcription complex to unwind, and cause greater levels of premature termination, the products of which are presumably degraded. The result of this is the dramatic downregulation seen in the microarray.

Previous studies of PG4 localisation in Arabidopsis have highlighted a greater number of two tetrad PG4s in genes annotated as responsive to drought than expected given the distribution across all genes [@Mullen2010; @Mullen2012]. Since intracellular potassium levels are increased during drought stress, it is possible that the stability of G4s could be increased. To investigate this potential G4 dependent regulatory mechanism, we analysed microarray data from drought stressed plants. We found that genes which are downregulated by drought stress contained more PG4s in the template strand, particularly in 5' UTRs and to a lesser extent in coding regions. This result matched closely to the enrichment of G4s in genes downregulated by NMM. Indeed, when we studied the overlap between drought stress downregulated genes and NMM downregulated genes, we found a strong overlap. Furthermore, the genes which were downregulated in both conditions were those which had the greatest PG4 density in their 5' UTRs. We suggest that during drought stress G4s in these genes form more strongly, causing blockages that pause Pol II, downregulating the expression of the gene. Finally, we found that genes upregulated by drought stress tended to contain higher levels of G4s in their 3' UTRs. This effect was not replicated by NMM treatment, suggesting an alternative mechanism of action. Since the 3' UTR is known to be an important regulator of mRNA stability and translation, we speculate these G4s form more strongly in the mRNA during drought stress and recruit some G4 binding factor which could enhance the stability of the mRNA.

\newpage
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJxekZrOU85eWhFTGYzMDR5Ijp7In
RleHQiOiJ3ZSIsInN0YXJ0IjoyMDIxMCwiZW5kIjoyMDIxMH0s
IjVzWEhYN3V3MzZPTmU5aTciOnsidGV4dCI6IlBQTFIiLCJzdG
FydCI6MjA1OTIsImVuZCI6MjA1OTZ9LCJpb2hLM0FaU1hneEVa
aFlhIjp7InRleHQiOiJuIE1BIHBsb3Qgb2YgYXZlcmFnZSBnZW
5lIGV4cHJlc3Npb24gdnMgbG9nIGZvbGQgY2hhbmdlLCB3aXRo
IGxvd2VzcyBjdXJ2ZSBmaXR04oCmIiwic3RhcnQiOjIwOTAzLC
JlbmQiOjIxMTE0fSwiSWt5c2VaWmw1NlFnWHRPdiI6eyJ0ZXh0
IjoiVGhpcyBnbG9iYWwgcGF0dGVybiwgd2hpY2ggdmlvbGF0ZX
Mgc29tZSBvZiB0aGUgYXNzdW1wdGlvbnMgdGhhdCBhcmUgdXN1
YWxseSB1c+KApiIsInN0YXJ0IjoyMTQwOCwiZW5kIjoyMTYzM3
0sInNyMlJ4ZHphTDFRa1Y2Sm0iOnsidGV4dCI6Ik1BIHBsb3Qg
c2hvd2luZyByZWxhdGlvbnNoaXAgYmV0d2VlbiBhdmVyYWdlIG
dlbmUgZXhwcmVzc2lvbiBhbmQgTG9nMiBmb2xkIGNoYW7igKYi
LCJzdGFydCI6MjE2OTksImVuZCI6MjE5NjF9LCIxTGJOcENEdn
ZueGZ2TTNQIjp7InRleHQiOiJCb3RoICBnZW5lc2V0cyBzaG93
IGFuIGdyZWF0ZXIgZXhvbmljIEc0IGRlbnNpdHkgb24gdGhlIH
RlbXBsYXRlIHN0cmFuZCB0aGFuIGdl4oCmIiwic3RhcnQiOjI0
MzkxLCJlbmQiOjI0MzkxfSwiekNXSVdXM3RGaGNrb0QxdyI6ey
J0ZXh0IjoiY29kaW5nIGFuZCB0ZW1wbGF0ZSBzdHJhbmRzIiwi
c3RhcnQiOjI0MzQ4LCJlbmQiOjI0Mzc1fSwieGpNbFFISkdsMl
dHUUtRSCI6eyJ0ZXh0Ijoic3RyaWtpbmciLCJzdGFydCI6MjUx
NDYsImVuZCI6MjUxNDZ9LCJnYVJGVU5TZmJvTVllTTg1Ijp7In
RleHQiOiIqRGlzdHJpYnV0aW9uIG9mIFBHNHMgaW4gZ2VuZXMg
ZGlmZmVyZW50aWFsbHkgcmVndWxhdGVkIGJ5IE5NTS4qKiBCYX
IgcGxvdHMgc2hv4oCmIiwic3RhcnQiOjI2NDg2LCJlbmQiOjI3
NDk3fSwiSXV2Qno0ZkF4Vk13MmU4dyI6eyJ0ZXh0IjoiY29tcG
FyZWQgdG8gcGVybXV0ZWQgcHJvZmlsZXMgYWNyb3NzIGFsbCBn
ZW5lcy4iLCJzdGFydCI6Mjg2NzIsImVuZCI6Mjg3ODJ9LCI2ck
FQeFJhRTFjQmNuVExvIjp7InRleHQiOiJwcm9wZWxsZXItbGlr
ZSBwYXJhbGxlbCBHNHMiLCJzdGFydCI6MjkxNzYsImVuZCI6Mj
kyMDN9LCJROWtHa3l3Z2h2em13TUhEIjp7InRleHQiOiJtYXkg
YmUgc3VmZmljaWVudCB0byBjYXVzZSBkb3ducmVndWxhdGlvbi
BvZiBhbiBvdGhlcndpc2UgRzQgZnJlZSBnZW5lIGR1cmluZyBO
4oCmIiwic3RhcnQiOjMxMDU5LCJlbmQiOjMxMTUxfSwiOFF2bU
RQeDJMWkxFakdCTiI6eyJ0ZXh0IjoiVGhpcyBzdWdnZXN0cyB0
aGF0IGNsdXN0ZXJzIG9mIEc0cyBkbyBoYXZlIGEgc3Ryb25nZX
IgZWZmZWN0IG9uIGdlbmUgZXhwcmVzc2lvIiwic3RhcnQiOjMw
OTIzLCJlbmQiOjMxMDAxfSwicGxESDFuSWFsTHdvejI4cSI6ey
J0ZXh0IjoibWVhbiIsInN0YXJ0IjozMTIzOSwiZW5kIjozMTI0
M30sIjVieE5teUhmdjVORVJiZXciOnsidGV4dCI6IkxlZnQgYW
5kIHJpZ2h0IHBhbmVscyBkZXBpY3QgY29kaW5nIGFuZCB0ZW1w
bGF0ZSBzdHJhbmRzIiwic3RhcnQiOjMxNDgxLCJlbmQiOjMxNT
M3fSwiNVd4T0dxVnc5NkpoRGtFTSI6eyJ0ZXh0IjoiYnkgTk1N
IHN0YWJpbGlzZWQgRzRzIGFwcGVhcnMiLCJzdGFydCI6MzE4MT
AsImVuZCI6MzE4Mzl9LCJvV1o5MWY4c0JvMkxvVjNDIjp7InRl
eHQiOiJUaGlzIHdhcyBzdXJwcmlzaW5nIGFzIGl0IGlzIGluIG
Rpc2FncmVlbWVudCB3aXRoIFBvbCBJSSBvY2N1cGFuY3kgcHJv
ZmlsZXMgaW7igKYiLCJzdGFydCI6MzI4NjcsImVuZCI6MzMwNT
B9LCI0Z21IQ0h0alF2SXV0Y1h5Ijp7InRleHQiOiJ3YXMgZ3Jl
YXRlciBQb2wgSUkgb2NjdXBhbmN5IGF0IHRoZSBUU1MgYW5kIG
luIHRoZSBUU1MgcHJveGltYWwgcGFydCBvZiB0aGUgZ2Vu4oCm
Iiwic3RhcnQiOjMzMzQ5LCJlbmQiOjMzNDM0fSwiT1A3b3RFdl
BZZnk1M0tDRyI6eyJ0ZXh0IjoibGFyZ2VzdCIsInN0YXJ0Ijoz
NTMxMiwiZW5kIjozNTMxOX0sIm9uVWhpZExhbjdEUW4wZEsiOn
sidGV4dCI6ImRhdGEgbm90IHNob3duIiwic3RhcnQiOjM1NDQ1
LCJlbmQiOjM1NDU5fSwiVlRSQjJVVVlpbEpwZTN2ViI6eyJ0ZX
h0IjoiMDA5Iiwic3RhcnQiOjM2MTY3LCJlbmQiOjM2MTcwfSwi
U3k3dUxjdWtuNmxRYVZCdiI6eyJ0ZXh0IjoicHJlc2VudCBjb2
RpbmcgYW5kIHRlbXBsYXRlIHN0cmFuZCIsInN0YXJ0IjozNzUw
NiwiZW5kIjozNzU0MH0sIkU2dmsxdlNPWnZWYmJhODEiOnsidG
V4dCI6IlNjYXR0ZXIgcGxvdCBzaG93aW5nIG1lYXN1cmVkIGV4
cHJlc3Npb24gaW4gbG9nMiBjb3VudHMgcGVyIG1pbGxpb24gKG
xvZ0NQTSkgZm/igKYiLCJzdGFydCI6MzcwNDksImVuZCI6Mzcy
NjV9LCI5QWczTVU3U0lvMjRjSHV1Ijp7InRleHQiOiIwLjkpIi
wic3RhcnQiOjM3NDAxLCJlbmQiOjM3NDA2fSwiRWdRZGJUSkFt
T0V1b1pmVSI6eyJ0ZXh0IjoiKmQpKiogMycgVVRSLCIsInN0YX
J0Ijo0Mjk5MSwiZW5kIjo0MzAwNH0sIlBhZWFweVBMdEl6ZFU3
MWoiOnsidGV4dCI6InN0cm9uZyIsInN0YXJ0Ijo0Mzk0MCwiZW
5kIjo0Mzk0Nn0sImkwVnY3NnkyWEFUY3JrelEiOnsidGV4dCI6
IipOTU0gZXhwcmVzc2lvbiBjaGFuZ2Ugb2YgQXJhYmlkb3BzaX
MgYW5kICpaIiwic3RhcnQiOjQ2ODY2LCJlbmQiOjQ2OTEwfSwi
UVZRQlNnZ3ZLS1Nyd2tTZiI6eyJ0ZXh0IjoiY2F1c2VkIiwic3
RhcnQiOjMzMTMsImVuZCI6MzMxOX0sInRITnhGNllPSkVlcWFi
ZWwiOnsidGV4dCI6IiMjIyBQbGFudCBHcm93dGggQ29uZGl0aW
9ucyIsInN0YXJ0Ijo0NTgxLCJlbmQiOjQ2MDh9LCJTaFVsTVlq
aFdtNVkxclN6Ijp7InRleHQiOiJzdXJmYWNlIHN0ZXJpbGlzZW
QsIiwic3RhcnQiOjQ4NTQsImVuZCI6NDg3M30sIndmMW5Ccm42
MWRRWDRUd1MiOnsidGV4dCI6InN0cmF0aWZpZWQgZm9yIDItMy
BkYXlzIiwic3RhcnQiOjQ4NzQsImVuZCI6NDg5N30sImV2Sjh0
Z0tiM2l6NDB3MFkiOnsidGV4dCI6Ik11cmFzaGlnZSAmIFNrb2
9nIChNUykiLCJzdGFydCI6NDk0NCwiZW5kIjo0OTY2fSwiTVFH
bzVRa1R4RDdrRnIyeCI6eyJ0ZXh0IjoiTVMgbGlxdWlkIG1lZG
lhIiwic3RhcnQiOjUwODEsImVuZCI6NTA4MX0sImR0SXlNUUk2
YWJzcGpWTHEiOnsidGV4dCI6Ik5NTSBhIiwiZW5kIjo1MDgxLC
JzdGFydCI6NTA4MX0sIklhUFRYOWJYNUdwdVFFYjciOnsidGV4
dCI6IldoaXRlIGFuZCBLYXBlciIsInN0YXJ0Ijo1MjgxLCJlbm
QiOjUyOTZ9LCJieVFvUGRhdUJvV1JWU3A1Ijp7InRleHQiOiJz
dGVyaWxlIHdhdGVyIiwic3RhcnQiOjUzNTMsImVuZCI6NTM2Nn
0sImp2VWxNck50bGFla0RsT3oiOnsidGV4dCI6IkFyYWJpZG9w
c2lzIEdlbmUgMS4wIFNUIiwic3RhcnQiOjU2NjMsImVuZCI6NT
Y4Nn0sIk96SUhsQ0NNVHFHTzl0anIiOnsidGV4dCI6IlJOQSBp
bnRlZ3JpdHkgYW5kIGFidW5kYW5jZSB3YXMgbWVhc3VyZWQgdX
NpbmcgYW4gQWxpZ2VudCBCaW9hbmFseXNlciAyMTAwIiwic3Rh
cnQiOjU3MDMsImVuZCI6NTc3N30sIjZTMWFjV1h4WkFGdGFvVn
UiOnsidGV4dCI6Ik5hbm9Ecm9wIDEwMDAgU3BlY3Ryb3Bob3Rv
bWV0ZXIiLCJzdGFydCI6NTQ0MywiZW5kIjo1NDc0fSwiaDBkdV
BOMVhuVXdVN2lBNyI6eyJ0ZXh0IjoiTWljcm9hcnJheSBBbmFs
eXNpcyIsInN0YXJ0Ijo1OTIwLCJlbmQiOjU5Mzl9LCJHdEw0aG
1lMXVQeGNDRFpyIjp7InRleHQiOiJBbmFseXNpcyBvZiBwcmV2
aW91c2x5IHB1Ymxpc2hlZCBtaWNyb2FycmF5IGRhdGEiLCJzdG
FydCI6NzY1MCwiZW5kIjo3Njk4fSwiOFd6QkxEd2RrYVIxcEto
ZSI6eyJ0ZXh0IjoiYG9saWdvYCBhbmQgYHB1bWFgIiwic3Rhcn
QiOjYwMTcsImVuZCI6NjA1MH0sIjlnOW1JREdKbjNWMzRjT0Ii
OnsidGV4dCI6IkVuc2VtYmwgYW5ub3RhdGlvbnMiLCJlbmQiOj
c1MDEsInN0YXJ0Ijo3NTAxfSwiRjQ5V21hVTdPZzV1RGd4MiI6
eyJ0ZXh0IjoiQXJhcG9ydDExIGdlbm9tZSBhbm5vdGF0aW9uIi
wic3RhcnQiOjEwMDk5LCJlbmQiOjEwMTI2fSwiTW9jbkNuTEdi
bjBEM2pTdSI6eyJ0ZXh0IjoiTk1NIHRyZWF0bWVudCBhcHBlYX
JlZCB0byBjYXVzZSBsYXJnZSBjb25zaXN0ZW50IGNoYW5nZXMg
aW4gZ2VuZSBleHByZXNzaW9uIHdoaeKApiIsInN0YXJ0Ijo2MT
g2LCJlbmQiOjYzMjh9LCJnaGx1TDZDZDlnelNmMVFRIjp7InRl
eHQiOiJSZXN1bHRzIHdlcmUgZmlsdGVyZWQgdXNpbmcgYSBkeW
5hbWljIHByb2dyYW1taW5nIGFwcHJvYWNoLCBjb21tb25seSB1
c2VkIGluIGlu4oCmIiwic3RhcnQiOjExMDIyLCJlbmQiOjExMT
E4fSwiREFraWJWQk9pU0pKbXdGUyI6eyJ0ZXh0IjoiU2VsZiBP
cmdhbmlzaW5nIE1hcCBBbmFseXNpcyIsInN0YXJ0IjoxMzYyOC
wiZW5kIjoxMzY1Nn0sImNYbWJGd1RXY0txUjk2aWoiOnsidGV4
dCI6IlBvbCBJSSBDaElQLXRpbGluZyBhcnJheSBhbmFseXNpcy
IsInN0YXJ0IjoxNDk0NSwiZW5kIjoxNDk3OH0sIk1kcWlTUWNi
bDlJamNpMzYiOnsidGV4dCI6Ikdsb2JhbCBSdW4gT24gKEdSTy
kgYW5kIFJOQSBzZXF1ZW5jaW5nIGRhdGEgZnJvbSBHRU8gQWNj
ZXNzaW9uIEdTRTgzMSIsInN0YXJ0IjoxNjg1OCwiZW5kIjoxNz
AxMX0sInZjMGdkMlJFTGVNc0lVM20iOnsidGV4dCI6Im1pY3Jv
YXJyYXkgYW5hbHlzaXMgdXNpbmcgUk5BIGZyb20gNyBkYXkgb2
xkIHNlZWRsaW5ncyB0cmVhdGVkIHdpdGggTk1NIGZvciA2IGji
gKYiLCJzdGFydCI6MjAyMjIsImVuZCI6MjAzMDZ9LCJuOEU3TF
pUZmZET0dFRHVqIjp7InRleHQiOiIzODgyIGRvd25yZWd1bGF0
ZWQsIDE5MzAgdXByZWd1bGF0ZWQpIiwic3RhcnQiOjIwNzQ4LC
JlbmQiOjIwNzg1fSwidVRCaGhkUnV2UDlNRWczRyI6eyJ0ZXh0
IjoiTk1NIGNhdXNlcyBnbG9iYWwgY2hhbmdlIGluIGdlbmUgZX
hwcmVzc2lvbiIsInN0YXJ0IjoyMDA3OSwiZW5kIjoyMDEyMn0s
IkJieWVkU0l6QW1zaHJPc3QiOnsidGV4dCI6InByb21vdGVyIH
JlZ2lvbiIsInN0YXJ0Ijo0ODYwNSwiZW5kIjo0ODYyMH19LCJj
b21tZW50cyI6eyJ6d3pHdlR3VEo1T3FZbnNxIjp7ImRpc2N1c3
Npb25JZCI6InF6Rms5Tzl5aEVMZjMwNHkiLCJzdWIiOiJnbzox
MDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoid2hvIGlzIF
wid2VcIi4gU2F5IGlmIHRoaXMgaXMgeW91IGFyZSBzb21lb25l
IGVsc2UuIiwiY3JlYXRlZCI6MTUzNTEyNzY0Mjg2Mn0sImZCVk
lNa3dzRnRLMGN0eDkiOnsiZGlzY3Vzc2lvbklkIjoiNXNYSFg3
dXczNk9OZTlpNyIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MT
AxMDY3NyIsInRleHQiOiJZb3UgbmVlZCB0byBleHBsYWluIHdo
YXQgUFBMUiBpcyBhbmQgd2h5IGl0IGlzIG1vcmUgYXBwcm9wcm
lhdGUgdGhhbiBub3JtYWwgbWV0cmljcyBoZXJlLiIsImNyZWF0
ZWQiOjE1MzUxMjc2OTk0Mzl9LCJCQWEzR2dkell4b0VvRDM4Ij
p7ImRpc2N1c3Npb25JZCI6ImlvaEszQVpTWGd4RVpoWWEiLCJz
dWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ij
oiSG93IGRvIHlvdSBrbm93IHRoaXMgaXNuJ3QgYSBub3JtYWxp
emF0aW9uIGFydGlmYWN0LCBwYXJ0aWN1bGFybHkgYXMgeW91IG
FyZSB1c2luZyBhIHdlaXJkIG5vcm1hbGlzYXRpb24gdGVjaG5p
cXVlPyIsImNyZWF0ZWQiOjE1MzUxMjc3MzI3MjZ9LCJWUGNWY3
NZM1IwaGhVeFIwIjp7ImRpc2N1c3Npb25JZCI6IklreXNlWlps
NTZRZ1h0T3YiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMT
A2NzciLCJ0ZXh0IjoiRXhwYW5kIGFuZCBleHBsYWluLiIsImNy
ZWF0ZWQiOjE1MzUxMjc3NjAyODV9LCJLZU45SmRIeGV5Mk40U2
dtIjp7ImRpc2N1c3Npb25JZCI6InNyMlJ4ZHphTDFRa1Y2Sm0i
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiVGhlIHN0cmFpZ2h0IGxpbmUgYXQgdGhlIGxlZnQgaGFu
ZCBlbmQgb2YgdGhpcyBwbG90IGlzIGEgYml0IHdlaXJkLiBUaG
lzIGlzIHByZXN1bWFibHkgZ2VuZXMgd2hlcmUgdGhlIGV4cHJl
c3Npb24gaW4gb24gY29uZGl0aW9uIGlzIDAgYW5kIGlzIFwibm
90IHplcm9cIiBpbiB0aGUgb3RoZXIgY29uZGl0aW9uLiIsImNy
ZWF0ZWQiOjE1MzUxMjc4NTY3ODF9LCJNNDZ3UHJKbmhJMXVuUn
RFIjp7ImRpc2N1c3Npb25JZCI6IjFMYk5wQ0R2dm54ZnZNM1Ai
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiYSkgdGhpcyBkb2Vzbid0IGxvb2sgdHJ1ZSB0byBtZSBm
cm9tIGxvb2tpbmcgYXQgdGhlIGZpZ3VyZS5cblxuYikgVGhpcy
BpcyBpbnRlcnByZXRhdGlvbiBhbmQgc2hvdWxkIGJlIGluIHRo
ZSB0ZXh0LCBub3QgdGhlIGZpZ3VyZSBsZWdlbmQuIiwiY3JlYX
RlZCI6MTUzNTEyNzk2MDA5MX0sImxDTVFoUHpXZWd6QzVLeXMi
OnsiZGlzY3Vzc2lvbklkIjoiekNXSVdXM3RGaGNrb0QxdyIsIn
N1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQi
OiJsYWJlbCBvbiBwbG90IiwiY3JlYXRlZCI6MTUzNTEyNzk3Mz
cwMn0sInN6SFNzUnc4VXZ5UG9HS2QiOnsiZGlzY3Vzc2lvbklk
IjoieGpNbFFISkdsMldHUUtRSCIsInN1YiI6ImdvOjEwMjIwNT
c5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaHkgaXMgdGhpcyBh
bnkgbW9yZSBcInN0cmlraW5nXCIgdGhhbiBvbiB0aGUgY29kaW
5nIHN0cmFuZCwgd2hlcmUgdGhlIGRpZmZlcmVuY2UgaXMgMTEl
PyIsImNyZWF0ZWQiOjE1MzUxMjgwNjQ2MTF9LCJMbTlsdUx2NF
N4OGlUbDQyIjp7ImRpc2N1c3Npb25JZCI6ImdhUkZVTlNmYm9N
WWVNODUiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2Nz
ciLCJ0ZXh0IjoiTGFiZWwgZnVsbC1nZW5lL2Nkcy81JyB1dHIv
IDMnIHV0ci9pbnRyb25zLiBMYWJlbCBjb2RpbmcvdGVtcGxhdG
UuXG5cblAtdmFsdWVzPyBIYXZlIHlvdSB0cmllZCBNYW5uLVdo
aXRuZXktVSB0ZXN0cyBpZiB5b3UgdGhpbmsgdC10ZXN0cyBhcm
VuJ3QgYXBwcm9wcmlhdGU/IiwiY3JlYXRlZCI6MTUzNTEyODIw
ODI4N30sIkx4Mkp6cGFObUNUeXpBZTEiOnsiZGlzY3Vzc2lvbk
lkIjoiSXV2Qno0ZkF4Vk13MmU4dyIsInN1YiI6ImdvOjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJOb3QgY2xlYXIgd2
hhdCBpcyBoYXBwZW5pbmcgaGVyZS4iLCJjcmVhdGVkIjoxNTM1
MTI4MzMzNjQ1fSwiQUNuYk8wNWJKeWNVcEJjNSI6eyJkaXNjdX
NzaW9uSWQiOiI2ckFQeFJhRTFjQmNuVExvIiwic3ViIjoiZ286
MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IlJlZmVyZW
5jZSB0byBmaWd1cmUgaW4gaW50cm9kdWN0aW9uIHRoYXQgZXhw
bGFpbnMgd2hhdCB0aGlzIGlzLiIsImNyZWF0ZWQiOjE1MzUxMj
gzNDk4MDN9LCJaM0M4NzNqUGl6dncySlFUIjp7ImRpc2N1c3Np
b25JZCI6IlE5a0dreXdnaHZ6bXdNSEQiLCJzdWIiOiJnbzoxMD
IyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiSSdtIG5vdCBz
dXJlIHlvdSBzaG93IHRoaXMuIiwiY3JlYXRlZCI6MTUzNTEyOD
M4MzM4OX0sIk0wMkRtUjBZM2Ewa0lHYnAiOnsiZGlzY3Vzc2lv
bklkIjoiOFF2bURQeDJMWkxFakdCTiIsInN1YiI6ImdvOjEwMj
IwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaGF0IGlzIHRo
ZSByZWxhdGlvbnNoaXAgYmV0d2VlbiBnZW5lcyB3aXRoIGxhcm
dlIGhpZ2hlc3QgZGVuc2l0eSBhbmQgZ2VuZXMgd2l0aCBqdXN0
IGEgbGFyZ2UgYW1vdW50IG9mIEc0cyBpbiBnZW5lcmFsPyIsIm
NyZWF0ZWQiOjE1MzUxMjg0NTc2MzB9LCJXTHlGUkduazBtaFRX
MEVFIjp7ImRpc2N1c3Npb25JZCI6InBsREgxbklhbEx3b3oyOH
EiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0
ZXh0IjoieS1heGlzIGxhYmVsIHNob3VsZCBzYXkgbWVhbi4iLC
JjcmVhdGVkIjoxNTM1MTI4NTE5MDg5fSwiYW5IUEladEtQTUo5
TFI1MCI6eyJkaXNjdXNzaW9uSWQiOiI1YnhObXlIZnY1TkVSYm
V3Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3Iiwi
dGV4dCI6ImxhYmVsIG9uIHBsb3QuIiwiY3JlYXRlZCI6MTUzNT
EyODUzMTc4OX0sIkVQcWZBZTJad29aS1N5Z0oiOnsiZGlzY3Vz
c2lvbklkIjoiNVd4T0dxVnc5NkpoRGtFTSIsInN1YiI6ImdvOj
EwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJncmFtbWVy
IiwiY3JlYXRlZCI6MTUzNTEyODYzODU0Nn0sIkpvN3J2aEc0OH
RyZ3hsSFMiOnsiZGlzY3Vzc2lvbklkIjoib1daOTFmOHNCbzJM
b1YzQyIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3Ny
IsInRleHQiOiJIdW1hbiBwcm9maWxlcyBhY3RhdWxseSB0ZW5k
IHRvIHNob3cgYSBwYXVzZSBhdCBib3RoIHRoZSBzdGFydCBhbm
QgdGhlIGZpbmlzaC4iLCJjcmVhdGVkIjoxNTM1MTI4NjY0NTY5
fSwiWXFsUndrZkFmTHQ0anE0YiI6eyJkaXNjdXNzaW9uSWQiOi
I0Z21IQ0h0alF2SXV0Y1h5Iiwic3ViIjoiZ286MTAyMjA1Nzk3
Mjc2OTQxMDEwNjc3IiwidGV4dCI6IkhvdyBhcmUgdGhlc2UgcG
xvdHMgbm9ybWFsaXNlZD8gRG9lcyB0aGUgYXJlYSB1bmRlciBl
YWNoIGN1cnZlIGhhdmUgdG8gc3VtIHRvIDA/IFRodXMgY291bG
QgdGhlIGFwcGFyZW50IGluY3Jhc2UgYXQgdGhlIDUnIGVuZCBi
ZSBkdWUgdG8gYSBkZWNyYXNlIGF0IHRoZSAzJyBlbmQ/IiwiY3
JlYXRlZCI6MTUzNTEyODc5MjMyMH0sIk52SXlmNVhBWWM3VTR3
RkwiOnsiZGlzY3Vzc2lvbklkIjoiT1A3b3RFdlBZZnk1M0tDRy
IsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRl
eHQiOiJsYXJnZXN0PyBXaGF0IG90aGVyIHRoaW5ncyBkaWQgeW
91IHRyeT8iLCJjcmVhdGVkIjoxNTM1MTMwNzc0Mzc3fSwiTVUx
VTl6aWZmZnBRNTJDTyI6eyJkaXNjdXNzaW9uSWQiOiJvblVoaW
RMYW43RFFuMGRLIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQx
MDEwNjc3IiwidGV4dCI6IldoeSBub3Q/IFRoaXMgc2VlbXMgbG
lrZSBhbiBpbXBvcnRhbnQgcG9zaXRpdmUgY29udHJvbC4gRXZl
cnlvbmUga25vd3MgdGhhdCBuY1JOQXMgYXJlIHVuc3RhYmxlLi
IsImNyZWF0ZWQiOjE1MzUxMzA4MDE4MDB9LCJxTGRNczFhQ2ph
eXVwaU9qIjp7ImRpc2N1c3Npb25JZCI6IlZUUkIyVVVZaWxKcG
UzdlYiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2Nzci
LCJ0ZXh0IjoidGVzdD8iLCJjcmVhdGVkIjoxNTM1MTMwODIwOD
AyfSwiYnBTbDRQcEZvUnpTaGMycCI6eyJkaXNjdXNzaW9uSWQi
OiJTeTd1TGN1a242bFFhVkJ2Iiwic3ViIjoiZ286MTAyMjA1Nz
k3Mjc2OTQxMDEwNjc3IiwidGV4dCI6ImxhYmVsIiwiY3JlYXRl
ZCI6MTUzNTEzMDg1NjM5NX0sInRYREZ5RGM0QkRuakxZZjAiOn
siZGlzY3Vzc2lvbklkIjoiRTZ2azF2U09adlZiYmE4MSIsInN1
YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOi
JJcyBpdCB3b3J0aCBzaG93aW5nIGNvbnRybCBhbmQgTk1NIFJO
QXNlcSBleHByZXNzaW9uIHNlcGVyYXRlbHkuIERvZXMgdGhlIE
dSTy9STkFzZXEgcmF0aW9uIGNoYW5nZSBmb3IgRzQgZ2VuZXM/
IiwiY3JlYXRlZCI6MTUzNTEzMDkyNTE3NX0sImpGWU1rVEVQMj
BmeVRLRlYiOnsiZGlzY3Vzc2lvbklkIjoiOUFnM01VN1NJbzI0
Y0h1dSIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3Ny
IsInRleHQiOiJpbiB0ZXh0IHRoaXMgc2F5cyAwLjk1IiwiY3Jl
YXRlZCI6MTUzNTEzMDk1NzY4N30sIkQ2Qm5vaTJ5ekpmc1Nydm
0iOnsiZGlzY3Vzc2lvbklkIjoiRWdRZGJUSkFtT0V1b1pmVSIs
InN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleH
QiOiJJZiB0aGUgZXJyb3IgYmFycyBhcmUgNjglIENJLCBJJ20g
Z29pZ24gdG8gc2F5IHRoYXQgdGhlc2UgZGlmZmVyZW5jZSBkb2
4ndCBsb29rIHBhcnRpY3VhbHJseSBzaWduaWZpY2FudC4iLCJj
cmVhdGVkIjoxNTM1MTMxMDMzMDg3fSwiS3VPYlhXdHVqVWJudW
1vdCI6eyJkaXNjdXNzaW9uSWQiOiJQYWVhcHlQTHRJemRVNzFq
Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidG
V4dCI6IklzIHRoaXMgYW55IHN0cm9uZ2VyIHRoYW4gdGhlIG92
ZXJsYXAgYmV0d2VlbiBETkEgZGFtYWdlIGFuZCBOTU0/IFF1YW
50aWZ5IHlvdXIgb3ZlcmxhcHMgd2l0aCBmb2xkIGVucmljaG1l
bnRzIG9yIE9kZHMgUmF0aW9zLiIsImNyZWF0ZWQiOjE1MzUxMz
EwODk2NDF9LCJQM29SbmljM3hpUHNmVzNlIjp7ImRpc2N1c3Np
b25JZCI6ImkwVnY3NnkyWEFUY3JrelEiLCJzdWIiOiJnbzoxMD
IyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiUm93IGFuZCBj
b2x1bW4gbGFiZWxzLiBcblxuTVdVIHRlc3RzLiIsImNyZWF0ZW
QiOjE1MzUxMzExMzIyNzN9LCJReUkyOVFsMlFXSDhoNnNIIjp7
ImRpc2N1c3Npb25JZCI6IlFWUUJTZ2d2S0tTcndrU2YiLCJzdW
IiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoi
dG9vIG1hbnkgY2F1c2VkcyIsImNyZWF0ZWQiOjE1MzUzNTU5MT
k3OTN9LCI0QlBGa3pGMTdBS0VqWFFkIjp7ImRpc2N1c3Npb25J
ZCI6InRITnhGNllPSkVlcWFiZWwiLCJzdWIiOiJnbzoxMDIyMD
U3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiRGlkIHlvdSBkbyBh
bGwgdGhpcyBsYWIgd29yaywgb3IgZGlkIHNvbWVvbmUgZWxzZS
4gSWYgeW91IGRpZCBpdCwgdGhlbiB0aGVzZSBtZXRob2RzIG5l
ZWQgdG8gYmUgbXVjaCBtb3JlIGRldGFpbGVkLiBJZiBzb21lb2
5lIGVsc2UgZGlkIGl0IHRoZW4geW91IG5lZWQgdG8gc2F5IHNv
LiIsImNyZWF0ZWQiOjE1MzUzNTU5OTkwMjV9LCIxQWI5U2QzZD
U1WUdtaW9LIjp7ImRpc2N1c3Npb25JZCI6IlNoVWxNWWpoV201
WTFyU3oiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2Nz
ciLCJ0ZXh0IjoiSG93PyIsImNyZWF0ZWQiOjE1MzUzNTYwMTE1
MjB9LCI5RHA4V0liWHg0ZE1wT1JDIjp7ImRpc2N1c3Npb25JZC
I6IndmMW5Ccm42MWRRWDRUd1MiLCJzdWIiOiJnbzoxMDIyMDU3
OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiSG93PyIsImNyZWF0ZW
QiOjE1MzUzNTYwMjQ0Mjd9LCJSTmhua3pER3NiUVdNT1c2Ijp7
ImRpc2N1c3Npb25JZCI6ImV2Sjh0Z0tiM2l6NDB3MFkiLCJzdW
IiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoi
UmVjaXBlPyIsImNyZWF0ZWQiOjE1MzUzNTYwNDA1OTZ9LCJVQV
dBTkVwbW9Xdk5IMlExIjp7ImRpc2N1c3Npb25JZCI6Ik1RR281
UWtUeEQ3a0ZyMngiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5ND
EwMTA2NzciLCJ0ZXh0IjoiUmVjaXBlPyIsImNyZWF0ZWQiOjE1
MzUzNTYwNTc1NzZ9LCJJZ2ptcGJ2Wk9Ubkhwcm1aIjp7ImRpc2
N1c3Npb25JZCI6ImR0SXlNUUk2YWJzcGpWTHEiLCJzdWIiOiJn
bzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiQ29uY2
VudHJhdGlvbj8gQWRkZWQgd2hlbiB0byBtZWRpYSBhdCB3aGF0
IHBvaW50cz8iLCJjcmVhdGVkIjoxNTM1MzU2MTUwMDMxfSwic0
V4Tzhwd1lVTnFPa1h2SSI6eyJkaXNjdXNzaW9uSWQiOiJJYVBU
WDliWDVHcHVRRWI3Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OT
QxMDEwNjc3IiwidGV4dCI6IklmIHlvdSBkaWQgdGhpcywgdGhl
biBcImFzIGRlc2NyaWJlZCBpblwiIGlzbid0IGdvb2QgZW5vdW
doLiIsImNyZWF0ZWQiOjE1MzUzNTYyMjEwMzV9LCJpU3BROHNY
eHg2bHhEc2FtIjp7ImRpc2N1c3Npb25JZCI6ImJ5UW9QZGF1Qm
9XUlZTcDUiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2
NzciLCJ0ZXh0IjoiSG93IG11Y2g/IiwiY3JlYXRlZCI6MTUzNT
M1NjIzMjQ3MH0sIlNKRE1lc1hQakh0RmRLVmsiOnsiZGlzY3Vz
c2lvbklkIjoianZVbE1yTnRsYWVrRGxPeiIsInN1YiI6ImdvOj
EwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJTdXBwbGll
cnMsIHBhcnQgbnVtYmVyIiwiY3JlYXRlZCI6MTUzNTM1NjI1NT
U4Mn0sIjRJWjRweTlhTHJWNkhla2giOnsiZGlzY3Vzc2lvbklk
IjoiT3pJSGxDQ01UcUdPOXRqciIsInN1YiI6ImdvOjEwMjIwNT
c5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJQcm90b2NvbCIsImNy
ZWF0ZWQiOjE1MzUzNTYyNzIyMDl9LCI5UWNVZWxSSWg3Y2FYSm
5JIjp7ImRpc2N1c3Npb25JZCI6IjZTMWFjV1h4WkFGdGFvVnUi
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiSG93IG11Y2ggd2FzIHVzZWQgZm9yIHF1YW50aWZpY2F0
aW9uPyIsImNyZWF0ZWQiOjE1MzUzNTYyODk3Mjh9LCI3bXZia1
g0cTNKNENJWGtQIjp7ImRpc2N1c3Npb25JZCI6ImgwZHVQTjFY
blV3VTdpQTciLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMT
A2NzciLCJ0ZXh0IjoiV2hlcmUgYXJlIHRoZSBhbmFseXNpcyBz
Y3JpcHRzPyIsImNyZWF0ZWQiOjE1MzUzNTYzMzU2Nzl9LCJlSz
E2cllQZWlaTG9oN1dnIjp7ImRpc2N1c3Npb25JZCI6Ikd0TDRo
bWUxdVB4Y0NEWnIiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5ND
EwMTA2NzciLCJ0ZXh0IjoiRGl0dG8iLCJjcmVhdGVkIjoxNTM1
MzU2MzQ0ODU5fSwiMjluSk9kdXBtdjRFeXRHSiI6eyJkaXNjdX
NzaW9uSWQiOiI4V3pCTER3ZGthUjFwS2hlIiwic3ViIjoiZ286
MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IlZlcnNpb2
5zPyIsImNyZWF0ZWQiOjE1MzUzNTYzNTg3Nzd9LCI3MVQ3WHJx
UHBCejFWczNBIjp7ImRpc2N1c3Npb25JZCI6IjlnOW1JREdKbj
NWMzRjT0IiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2
NzciLCJ0ZXh0IjoiRW5zZW1ibCB2ZXJzaW9uPyIsImNyZWF0ZW
QiOjE1MzUzNTYzODAyNzJ9LCJqWWpNN0hDVGJUYVBkR0lCIjp7
ImRpc2N1c3Npb25JZCI6IkY0OVdtYVU3T2c1dURneDIiLCJzdW
IiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoi
WW91IHNhaWQgYWJvdmUgdGhhdCBFbnNlbWJsIHdhcyB1c2VkPy
BJcyBFbnNlbWJsIGJhc2VkIG9uIEFyYXBvcnQ/IiwiY3JlYXRl
ZCI6MTUzNTM1NjQzMzc4Nn0sInFXRnRqUURreUhaZHhiU1AiOn
siZGlzY3Vzc2lvbklkIjoiTW9jbkNuTEdibjBEM2pTdSIsInN1
YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOi
JIb3cgZGlkIHlvdSBkZXRlcm1pbmUgdGhpcz8iLCJjcmVhdGVk
IjoxNTM1MzU2NDcxMjkyfSwieTFGdHdleHFuRmR3VmlsNyI6ey
JkaXNjdXNzaW9uSWQiOiJnaGx1TDZDZDlnelNmMVFRIiwic3Vi
IjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ik
Rlc2NyaWJlZCwgb3IgYXQgbGVhc3QgbmFtZSBhbmQgcmVmZXJl
bmNlIHRoZSBhbGdvLiIsImNyZWF0ZWQiOjE1MzUzNTY1MjMwMj
d9LCJDVlpZb3JHTnlSQWtmZ1dmIjp7ImRpc2N1c3Npb25JZCI6
IkRBa2liVkJPaVNKSm13RlMiLCJzdWIiOiJnbzoxMDIyMDU3OT
cyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiU2NyaXB0cz8iLCJjcmVh
dGVkIjoxNTM1MzU2NTY1ODQ1fSwidXBTVmJTRTRIa3JGemphai
I6eyJkaXNjdXNzaW9uSWQiOiJjWG1iRndUV2NLcVI5NmlqIiwi
c3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dC
I6IlBhY2thZ2UgdmVyc2lvbnMgYW5kIHNjcmlwdHMuIiwiY3Jl
YXRlZCI6MTUzNTM1NjcxOTAzMX0sIjhFakhuTFNGSTV5QXA0SV
EiOnsiZGlzY3Vzc2lvbklkIjoiTWRxaVNRY2JsOUlqY2kzNiIs
InN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleH
QiOiJJbmNvbnNpc3RlbnQuIFlvdSBsaXN0ZWQgdGhlIGdvd3Ro
IGNvbmRpdGlvbnMgYWJvdmUsIGJ1dCBub3QgaGVyZS4iLCJjcm
VhdGVkIjoxNTM1MzU2NzY1NzU4fSwiYmdtUGtGNmRsOVVPZkxp
cCI6eyJkaXNjdXNzaW9uSWQiOiJNZHFpU1FjYmw5SWpjaTM2Ii
wic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4
dCI6Ik5lZWRzIHBhY2thZ2UgdmVyc2lvbnMgYW5kIHNjcmlwdH
MgdXNlZC4iLCJjcmVhdGVkIjoxNTM1MzU2NzkyMTQ4fSwidmtK
alR3Z2hCdFZ0YlJFRiI6eyJkaXNjdXNzaW9uSWQiOiJ2YzBnZD
JSRUxlTXNJVTNtIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQx
MDEwNjc3IiwidGV4dCI6IkhvdyBtYW55IHJlcHM/IFdoZXJlIG
lzIHRoZSBRQz8iLCJjcmVhdGVkIjoxNTM1MzU2ODUxMzk5fSwi
R1dmWWtBSUxlejM3YnVIMCI6eyJkaXNjdXNzaW9uSWQiOiJuOE
U3TFpUZmZET0dFRHVqIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2
OTQxMDEwNjc3IiwidGV4dCI6IlRoZXNlIGFyZSB0aGUgb3RoZX
Igd2F5IGFyb3VuZCB0byBhYm92ZS4iLCJjcmVhdGVkIjoxNTM1
MzU2ODcyMTU5fSwiNmVndEdqS01FTEhiYlR4cSI6eyJkaXNjdX
NzaW9uSWQiOiJ1VEJoaGRSdXZQOU1FZzNHIiwic3ViIjoiZ286
MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IllvdSBuZW
VkIHRvIHRhbGsgc29tZXdoZXJlIGFib3V0IHRoZSBkZWNpc2lv
biB0byBub3JtYWxpc2UgdXNpbmcgUFVNQSByYXRoZXIgdGhhbi
BsaW1tYS4gSG93IGRvIHlvdSBrbm93IHRoYXQgUFVNQSBpc24n
dCBqdXN0IGdldHRpbmcgdGhlIG5vcm1hbGlzYWl0b24gd3Jvbm
c/IiwiY3JlYXRlZCI6MTUzNTM1NjkxNzYwMH0sIjdJTHdoTFNJ
T1g3SmZma1IiOnsiZGlzY3Vzc2lvbklkIjoiQmJ5ZWRTSXpBbX
Nock9zdCIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3
NyIsInRleHQiOiJEaWQgeW91IGV2ZXIgY2hlY2sgcHJvbW90b3
IgcmVnaW9ucz8iLCJjcmVhdGVkIjoxNTM1MzU3MDA0ODgxfX0s
Imhpc3RvcnkiOlsxMzY0OTQ0ODksLTIwNTQ2MTYxNTMsLTM0Mz
g2MzkzNywyMDIzNDcyNDUxLDgxMjY0NTA5OCwtMjA3MjU2ODY0
MSwtMTAwNzQxNTA4NiwxODg5ODU0NDY1LDY2ODM2ODk3MywxMz
QyNTM4MTg0LC0xOTE0NTc4ODkzLDIyNjc2NDc5LC01Mjc3NDYx
NjksLTM3OTEyMjEzNywtODY4MDM1MDU5LC0xMTExMTEyMDkwLD
EwMzM2MTUzNjQsLTE2OTk0ODIzOCwxMjg0NDA3MTA5LDMzNDY2
NTQ4MF19
-->
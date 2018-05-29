# Chapter 5: Effect of G-Quadruplexes on expression and splicing of the Extensin gene family

\newpage

## Introduction

\newpage

## Methods

### Plant growth conditions and drug treatments

* Growth conditions & drug treatments for qPCR
* Growth conditions & drug treatments for RNAseq
* Growth conditions for sanger sequencing

### RNA extraction for qPCR and RNAseq

Total nucleic acid isolation protocol was carried out by phenol-chloroform extraction as described by White and Kaper (1989). The resulting pellets were resuspended in sterile water and stored at -80C. The RNA concentration and quality was checked using the NanoDrop 1000 Spectrophotometer (ThermoScientific).

### RNAseq analysis

RNA was sent to Sheffield Children's Hospital Genomics Facility for library preparation and sequencing. Polyadenylated RNA was enriched using NEBNext Poly(A) mRNA Magnetic Isolation, and libraries were produced using the NEBNext Ultra II Directional RNA kit for Illumina. The chemical fragmentation step was adjusted to increase estimated insert size to the 400-450bp range. Paired end sequencing was conducted across two lanes of an Illumina HiSeq 2500 in rapid run mode, with 220bp read length. The run returned 157 million pass filter reads in lane 1, and 154 million in lane 2. Initial read preprocessing and adaptor trimming was conducted by the Genomics Facility.

Read quality was assessed locally using `FastQC`. Differential expression analysis was conducted by using `salmon` pseudoalignment to estimate transcript abundance against Araport11 cDNA and ncRNA. Mean insert size for each sample was assessed to be in the 400-500bp range. Gene level abundance was then aggregated from this using `tximport` and differential expression testing was conducted using `edgeR` linear modelling. Normalised log2 counts per million (CPM) were calculated and used for plotting. P values were adjusted using the Benjamini Hochberg multiple testing correction. Barplots of NMM and DMSO expression were produced using `seaborn`.

Spliced reads were mapped to the TAIR10 genome using `STAR`. The parameters for spliced mapping were adjusted to increase precision. A minimum of 8bp overhang was used for unannotated splice junctions, and 5bp for annotated splice junctions (following ENCODE guidelines). Minimum intron size was set to 60bp and maximum intron size to 10000bp. Output BAM files were sorted and indexed using `samtools`.

### Gene Ontology analysis

For PG4 enrichment Ontology analysis, two tetrad Quadparser PG4s were predicted in the TAIR10 genome using `g4predict`. The number of PG4s overlapping each strand of the flattened exon models for each gene in Araport11 was calculated using `bedtools intersect`. To calculate enrichment, a permutation test experiment was conducted, where each gene was assigned a weighting proportional to the total length of its exonic sequence. PG4s were then shuffled randomly amongst all genes. For each Ontology group, the number of PG4s observed in genes of that group was compared to the expected numbers if PG4s were distributed randomly across the transcriptome. 10000 permuations were used for testing, and two tailed P values were calculated as $max(min(\frac{\sum_{i=0}^nexp_i < obs}{n}, \frac{\sum_{i=0}^nexp_i > obs}{n}), \frac{1}{n})$, where $n$ is the total number of permutations and $exp_i$ is the expected value from the $i$th permuation. P values were adjusted using the Benjamini Hochberg multiple testing correction.

For Gene Ontology analysis of differentially expressed genes, `GOseq` was used. Up and downregulated genesets were produced using a log2 fold change threshold of 1 and an FDR of 0.05. Weighting factors used were the median transcript length for each gene. P values for enrichment were produced by `GOseq` using the Wallenius approximation method and were corrected for multiple testing using the Benjamini Hochberg method.

Tables of enriched GO terms were generated using `pandas` and formatted with `inkscape`.

### Extensin gene family total PG4 estimation

For estimated PG4 numbers in the table in Fig \ref{ext_table}, PG4s were predicted using three different methods. All instances of the dinucleotide `GG` were identified in each Extensin gene, and then a graph was built (using `networkx`) where each `GG` was a node and nodes were connected by an edge if dinucleotides were less than 7bp apart from each other. The number of overlapping PG4 conformations was then calculated as the number of subgraphs in the graph with exactly four members, whilst the number of merged PG4s was calculated as the number of unconnected subgraphs with four or more members. To identify the number of non-overlapping PG4s, a dynamic programming method was used. Overlapping PG4s were grouped and scored by inverse length, then filtered for the maximum number of high scoring, non-overlapping PG4s. Code for these calculations is available in Appendix XXXX.

### Quantitative PCR experiments

### Analysis of public RNAseq data

Root RNAseq from Li et al 2016 was downloaded in FASTQ format from ENA. Quality assessment was performed with `FastQC` and `fastq-screen`, and adapter contamination was removed using `Cutadapt`. The data was mapped using `STAR` with default settings, except than max intron length was set to 10000bp. Output BAM files were sorted and indexed using `samtools`.

Normalised gene expression estimates were generated using `featureCounts` to get raw counts of alignments overlapping each gene in Araport11, and `edgeR` to perform estimation of log2 counts per million.

### Splice junction analyses

Splice junction sites were extracted from aligned reads using `pysam`. For all analyses, reads were filtered to produce a set of unique donor/acceptor site pairs. Scatter plots of spliced read percentages and frequency plots of in frame exitrons were produced using `matplotlib` and `seaborn`.

### Sequence logo generation

Consensus sequence logos were generated using an in-house Python module `matplotlib_logo`. Unique splice junction pairs from spliced reads were identified, and the corresponding sequence information (8bp up and downstream of donor and acceptor) was extracted from the TAIR10 genome using `pysam`. Position frequency matrices were generated from these sequences, and entropy score in bits was calculated and plotted.

### RNAseq read simulation and bootstrapping experiment

For RNAseq read simulation, read counts generated using `featureCounts` for each file in the Li et al. 2016 dataset were used to inform a `polyester` simulation of Illumina 125bp paired end reads, using the Araport11 annotation. Fragment length was simulated using a normal distribution with mean 250bp and standard deviation of 25bp. Error rate was simulated at a uniform 0.5% across samples, reads, and nucleotides.

For the bootstrapping experiment, one or more pairs of real and simulated samples were randomly sampled from the dataset and unique splice donor/acceptor pairs in EXT9 were identified. Only splice sites from primary alignments, with overhangs greater than 20bp, and with more than 2 suppporting reads were kept. For each splice site, the exonic overhang sequence 20bp upstream of the donor site and 20bp downstream of the acceptor site were extracted from the read and concatenated. Any splice site whose concatenated sequence was present as a contiguous kmer in the reference EXT9 sequence was assumed to be a mapping error and filtered out. Finally, splice sites were deduplicated by directional clustering of sequences with edit distance of one or less using `umi_tools` directional clusterer. 500 iterations were used for bootstrapping. 67% condifence intervals were produced and plotted using `seaborn`.

### Mappability analyses

Mappability scores were generated for the TAIR10 genome using `gem-mappability` with a kmer size of 75bp, and converted to `BigWig` format using `gem-2-wig` and `wigToBigWig`. Minimum mappability scores for each Extensin gene were extracted using `pyBigWig` and plotted against spliced fraction using `matplotlib`.

### Sanger sequencing analysis

### Differential Splicing Analysis

For differential splicing analysis, we identified novel exitrons in our RNAseq dataset to augment the existing Araport11 annotation. Gene models from Araport11 were flattened into "exon chunks", defined as ranges of bases which are all present in the same set of transcripts. Exitrons were identified by extracting reads with spliced segments which fell within a single exon chunk (i.e. were flanked by exon on both sides). Only those exitrons with 20 or more reads supporting them were kept for downstream analyses. Exitrons were then added into the "exon chunk" gene models, which were then used for differential exon usage analyses.

Read counts per exon chunk were calculated using `featureCounts` at exon level. Strand specific counting was used, and only single mapping concordantly mapped reads were counted. Reads that overlapped more than one exon chunk were counted once for each chunk. Differential exon usage between NMM and DMSO treatments was conducted in R using `DEXSeq` with default settings, and P values were adjusted using Benjamini Hochberg correction.

\newpage

## Results

### Gene Ontology shows plant cell wall specific genes are enriched in PG4s and downregulated by NMM

To identify gene ontology groups which are specifically enriched with exonic PG4s, we compared the observed levels of PG4s per gene to expected levels if PG4s were randomly distributed across all genes (weighted by gene length). These observed and expected levels were summarised for each gene ontology group. Sorting the results for groups with the greatest positive observed/expected ratio of PG4s on the template strand, we discovered that gene ontology groups involved with functions at the cell periphery, particularly in the plasma membrane and cell wall, had strong enrichments (Fig \ref{go_table}). The log2 fold enrichment in `GO:0005199`, which contains structural cell wall genes, was +4.4 (FDR < 4.8e-4). This corresponded to an observed number of 992 PG4s in only 32 genes (the average expection under the null hypothesis was 46 PG4s). These PG4 dense gene ontology groups were also strongly enriched in the set of genes which are significantly downregulated by NMM in our RNAseq dataset (Fig\ref{go_table}, 50% of expressed genes in `GO:0005199` were downregulated by NMM, FDR = 9.6e-7).

\newpage

![**Gene Ontology groups enriched in template stranded putative G-Quadruplexes** Table showing the top ten Gene Ontology groups most enriched for exonic PG4s compared to null distribution. The top two groups, both containing genes involved in cell wall structure and organisation, are also enriched for genes downregulated by NMM. \label{go_table}](figures/g4_nmm_gene_ontology_table.svg)

\newpage

### The proline rich Extensin gene family contain large numbers of hardcoded PG4s

We discovered that the `GO:0005199` geneset was primarily made up of genes from the Extensin cell family (29/32 genes, 90.6%), including classical SP4/SP5 Extensin genes and chimeric Leucine Rich Repeat/Extensin (LRX) genes. These genes were found to be extremely PG4 rich on the template strand, with many genes containing greater than 10 PG4s per kilobase of exon (Fig \ref{ext_genes}a). Upon visualisation of these genes, we noted that in the majority of cases the PG4s were regularly spaced along the gene, and were contained solely within the coding region (CDS) of the gene (Fig \ref{ext_genes}b).

\newpage

![**Expression and PG4 density of genes in the Cell Wall Structural Ontology group `GO:0005199`** **a)** Panels showing gene expression (top panel) and PG4 density (bottom panel) for genes in the `GO:0005199` group. Expression in DMSO (blue) and NMM (orange) conditions is shown at log2 counts per million (from root RNAseq dataset). Errorbars show standard deviation. Genes which are differentially expressed with FDR < 0.05 are labelled with asterisks. In PG4 panel, the exonic PG4 density per kilobase is shown separately for coding (blue) and template (orange) strands of the gene. **b)** Gene tracks showing the location of predicted two tetrad PG4s in orange for (from top to bottom) LRX1, EXT13, and EXT9. Gene models from Araport11 are shown in blue. In gene models, thin boxes represent untranslated regions (UTRs), fat boxes represent coding regions (CDS), and connecting lines represent intronic regions. \label{ext_genes}](figures/extensin_gene_ontology_group_expression_g4s.svg)

\newpage

From a search of the literature, we discovered that Extensin genes are highly repetitive, proline rich proteins. These proteins polymerise to function as a structural matrix in the protein component of the plant cell wall. In particular, we noted that these proteins are characterised by large numbers of the SP3-5 repeat, which is made up of the sequence $Ser(Hyp)_{3-5}$, where Hyp is Hydroxyproline, a proline derivative. Since the codon for proline is `CCN`, the DNA which encodes SP4 and SP5 motifs will conform to the two tetrad Quadparser motif on the template strand of the gene (Fig \ref{cd_spec}a). This is the source of the PG4 density of Extensin genes, and PG4 counts in these genes are well correlated with SP3-5 repeats (Fig \ref{ext_table}). Since the SP4 motif is required for the function of the protein, and is restricted by the codon for proline, these PG4s are "hardcoded" into the body of the gene.

\newpage

![**The Extensin gene family contains large numbers of hardcoded PG4s** Table showing extended Extensin gene family, their expression patterns, SP4 motif counts, PG4 counts and expression during NMM treatment. Adapted from Showalter et al. 2010 \label{ext_table}](figures/ext_family_showalter_table.svg){height=750px}

\newpage

To demonstrate that the PG4 from Extensin genes could form a G4 structure in vitro we used circular dichroism spectroscopy (CD). We performed these experiments at physiologically relevant temperatures for Arabidopsis. An oligo representative of the SP4 repeat was designed (`AGAGGTGGTGGTGGTATG`) using 3bp flanks upstream and downstream of the PG4. CD showed the G4 oligo had peak absorbance at 260nm and trough at 240nm, indicative of a parallel G4 structure (Fig \ref{cd_spec}b). Removing the PG4 by mutating the sequence (mutated sequence: `AGAGGTGATGGTGGTATG`) or removing potassium ions from the buffer abolished this absorbance profile.

\newpage

![**The Extensin SP4 motif forms a G-Quadruplex *in vitro*.** **a)** Schematic showing how the Extensin SP4 protein motif hardcodes a two tetrad PG4 into the template strand of the gene body. **b)** CD spectroscopy of an Extensin repeat sequence (left) and a mutated control which does not conform the the Quadparser motif (right) show that the Extensin repeat forms a G4 *in vitro*. This is indicated by the peak in ellipticity at 260nm and the trough at 240nm, which are characteristic of a parallel G4. \label{cd_spec}](figures/cd_spectroscopy1.svg){height=750px}

\newpage

### Extensins are strongly downregulated by NMM and Berberine

To confirm that the Extensin genes are downregulated by NMM, we performed RNA extraction and quantitative RT-PCR (qPCR) on root tissue from 7 day old Arabidopsis seedlings treated for 6 hours with NMM at varying concentrations. EXT13 and LRX1 were chosen as representative classical and chimeric Extensins, respectively. The change in expression of both genes upon treatment was negatively correlated with the concentration of NMM applied (Fig \ref{nmm_berb_qpcr}a). Treatment with the G4 intercalating drug, Berberine, also caused strong downregulation of EXT13 and LRX1 (Fig \ref{nmm_berb_qpcr}b). Since NMM and berberine are very different drugs which stabilise G4s through different methods, taken together our results suggest downregulation of Extensins is caused by G4 stabilisation.

### Downregulation of Extensins by NMM is translation independent

To confirm whether downregulation of EXT13 and LRX1 by NMM was direct, or the result of a perturbation the levels of a transcription factor, we conducted qPCR experiments with combinatorial treatment of NMM and Cyclohexamide (CHX). CHX is an inhibitor of translation which is commonly used to determine whether interactions by transcription factors on a gene are direct. If effects are indirect (i.e. if the transcription factor of interest regulates transcription of some intermediate transcriptional factor, which regulates the gene of interest), then treatment with CHX will prevent regulation, since any intermediate factors will not be able to be translated. In the case of NMM treatment, this was used to see whether NMM acts directly on EXT13 and LRX1 through G4 stabilisation, or through other changes in the transcriptome. Seedlings were pretreated with CHX for two hours, before NMM was added and treatment was continued for another 6 hours. This experiment showed that that Extensin downregulation by NMM still occurs even when translation is blocked, suggesting that NMM acts directly upon the Extensin genes.

\newpage

![**Expression of EXT13 and LRX1 during treatment with G4-binding ligands** Scatter/strip plots showing qPCR results for EXT13 (left panels) and LRX1 (right panels). Log2 fold change in expression (ΔΔCT) of Extensin genes decreases with increasing concentrations of **a)** NMM and **b)** Berberine. **c)** NMM downreguation of EXT13 and LRX1 is not affected by concurrent Cyclohexamide treatment, suggesting a mechanism independent of translation. For all panels, each point is a single technical replicate, and colours represent different biological replicates. A small amount of jitter has been added to the X axis for better visualisation of results. \label{nmm_berb_qpcr}](figures/ext13_lrx1_qpcr.svg){height=750px}

\newpage

### RNAseq suggests Extensin genes contain exitronic splice sites

We noted from studying *de novo* assembled splice isoforms from a root specific RNAseq dataset (Li et al. 2016) that many of the Extensin genes had large numbers of novel spliced isoforms. EXT9 was found to have the most novel spliced forms of any gene in the dataset. These were not present in the annotation. In fact, the majority of Extensin domains are annotated in both TAIR10 and Araport11 as intronless. These splice isoforms are presumably therefore a product of "exitronic" splicing (Marquez et al. 2015), where sections of constitutive exons, flanked on both sides by exonic sequence, are spliced out of a gene. We hypothesised that these unusual exitrons could be a result of slow Pol II elongation allowing splicing to occur at weak splice sites.

A hallmark of most true splice junctions is the GT/AG intron motif, which is the conserved canonical sequence in all higher eukaryotes, including Arabidopsis (Fig \ref{splice_junct}a). To determine whether Extensin exitrons had canonical splice motifs, we produced splice junction sequence logos for predicted introns from the dataset produced by Li et al., for EXT9 and LRX3, both of which were highly spliced. We found that splice junctions in these genes had near universal GT/AG motifs (Fig \ref{splice_junct}b). Upon inspection of the methods for Li et al., however we discovered that CuffLinks was used for *de novo* transcript assembly. Since the RNAseq dataset is unstranded, Cufflinks requires the upstream mapping tool (here, STAR) to annotate the orientation of spliced reads using the intron motif (i.e. positive strand for GT/AG and negative strand for CT/AC). This setting means reads which do not conform to the intron motif are discarded, leading to serious bias.

To remove this bias, we remapped reads from the Li et al. dataset using STAR without filtering by intron motif. Since assemblers like CuffLinks and StringTie require strandedness information derived from the intron motif, transcript assembly was not possible. Instead, we simply extracted spliced reads aligning to EXT3 and LRX3 and identified all unique splice site starts and ends. The corresponding sequences were then used to produce sequence logos (Fig \ref{splice_junct}c). These logos showed only a weak enrichment for the GT/AG motif in EXT9, and CT/AC in LRX1.

\newpage

![**Splice junction motifs for EXT9 and LRX3** Sequence logo plots showing consensus splice site sequences around donor (left panels) and acceptor (right panels) splice sites. Putative intronic sequences are shown on shaded blue background. **a)** Splice junction consensus sequence logo for Arabidopsis, calculated from junctions in the Araport11 annotation. **b)** Splice junction consensus sequence logo produced from **de novo** assembled transcripts (Li et al. 2016) for EXT9 and LRX3. **c)** Splice junction consensus sequence logo produced from unique donor/acceptor pairs identified from spliced reads on EXT9 and LRX3. \label{splice_junct}](figures/splice_site_sequence_logo.svg)

\newpage

Since the Extensin exitrons appear in coding regions of DNA, if the spliced out region is not a multiple of three, then the resulting mRNA would be frameshifted, producing truncated and potentially deleterious proteins. We therefore tested whether the spliced reads in EXT9 and LRX3 contained gaps which were multiples of three or not. For both genes, almost all of the unique splice junction pairs were multiples of three (Fig \ref{splice_frame}a-b). This could be evidence that these exitrons are genuine and produce function gene products. On the other hand, we noted that splicing tended to occur between regions with high protein and DNA sequence level homology. EXT9, for example, is an incredibly repetitive gene with high self-homology (Fig \ref{ext_mapp}a). These in frame splice junctions could therefore simply be the result of mapping errors from the spliced aligner STAR, which utilises heuristics which may result in some reads from contiguous parts of the genome being mapped as spliced. If the homologous regions which could cause mapping errors within a gene also have homology at the protein level, as is the case in EXT9, then it is probable that erroneously spliced reads would be a multiple of three in intron length.

If spliced reads mapping to EXT9 were the result of some systematic error in mapping, one might not expect to see much variation in the percent of reads mapping to a gene being spliced. We therefore correlated the expression of EXT9 in each sample from the root RNAseq dataset (measured in log2 counts per million or needs to be!) with the percent of reads which mapped with a splice site. We found a slight positive correlation between expression and splicing (Fig \ref{splice_frame}c).

As a further precaution against these erroneous spliced mappings, we performed read simulation for each sample in the root RNAseq dataset. Put simply, the expression of each gene was quantified for each sample by counting the number of mapped reads, then Polyester (an Illumina sequencing read simulator) was run to generate reads from the reference transcriptome with the same read counts. These simulated reads were then remapped with STAR using the same parameters as the original mapping. We then performed a bootstrap analysis for EXT9 where we sampled one or more real/simulated sample pairs, and counted the number of unique splice donor/acceptor pairs that occurred in each. Junctions with the same exonic flanking sequence (using 20bp overhangs) or with edit distance of only one base were collapsed. Any junctions with flanking sequence that appeared as a contiguous kmer in the reference sequence of EXT9 were also removed. Despite this, we saw a consistently larger number of unique donor/acceptor splice pairs in the real data than in the simulated data (Fig \ref{splice_frame}d).

\newpage

![**Splice junction motifs for EXT9 and LRX3** **a** & **b)** Frequency barplots showing number of In/Out of frame splice junctions for EXT9 and LRX3 respectively. **c)** Scatterplot showing percentage of reads with splicing versus log2 counts per million for EXT9. **d)** Bootstrapped splicing simulation showing number of unique EXT9 splice junctions discovered with increasing numbers of samples for real root RNAseq data versus paired simulated RNAseq data. Errorbars are 67% confidence intervals. \label{splice_frame}](figures/splice_site_frame_and_simulation.svg)

\newpage

Since the Extensin genes are highly repetitive, this reduces the ability of read aligners to map to them. This "mappability" can be quantified using tools such as GEM which measure, for each genomic position, how often the sequence kmer that is found there occurs in the rest of the genome. We utilised GEM to score the mappability of the Arabidopsis genome, and compared the median mappability score for extensin genes to the percent of spliced reads for each gene. Only genes which were annotated as intronless (i.e. no constitutive introns, genes with exitrons were allowed) were included. We found a clear negative correlation between mappability and the number of mapped spliced reads (Fig \ref{ext_mapp}b). For EXT9, the regions with lowest mappability are clearly those with the most annotated splice sites, including in the Araport11 reference (Fig \ref{ext_mapp}c). 

\newpage

![**Extensin genes with greater spliced mapped reads have low mappability** **a)** Dotplot showing self homology of the EXT9 cDNA (unspliced isoform). Positional identity was calculated using 15bp windows across the gene. Positions with identity less than 75% were filtered to remove noise from the plot. **b)** Scatter plot showing the minimum mappability score of unspliced Extensin genes  against the percentage of spliced reads for that gene in mature root RNAseq samples (Li et al. 2016). Errorbars are standard deviation. **c)** Gene track showing the mappability score across EXT9 (orange). Most splice forms cross these low mappability regions, including in the reference annotation Araport11 (shown in blue). \label{ext_mapp}](figures/ext9_dotplot_mappability.svg)

\newpage

### Sanger sequences identifies LRX1 and EXT9 splice variants

To experimentally confirm whether Extensin gene exitron splicing exists, we performed RT-PCR of LRX1 and EXT9 mRNAs. PCR products of both genes showed multiple products, characteristic of several spliced forms (data not shown). PCR products which did not correspond to the full length of the unspliced mRNA were gel extracted, cloned and sanger sequenced to identify their origin. We identified a number of mRNA fragments originating from the LRX1 and EXT9 genes. Alignment of these products using BLAT identified a number of spliced isoforms in both genes. To identify whether these isoforms contained canonical splice sites, we produced sequence logos. Neither gene showed a clear pattern conforming to the canonical intron motif GT/AG, though the products from  LRX1 showed the reverse complement of this pattern, CT/AC.

\newpage

![**Sanger sequencing of LRX1 and EXT9 cDNA identifies spliced forms** **a)** Gene track showing aligned sanger sequencing products for **a)** LRX1 and **b)** EXT9. Products aligned to the forward strand are shown in green, and products aligned to the negative strand are shown in orange. Gene models are from the Araport11 annotation. **c)** Sequence logos for sanger product splice junctions for LRX1 (top panel) and EXT9 (lower panel). \label{sanger}](figures/sanger_splice_variants.svg)

\newpage

### NMM treated plants do not have increased splicing of Extensin genes.

\newpage

## Discussion

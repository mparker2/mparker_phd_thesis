# Effect of G-Quadruplexes on expression and splicing of the Extensin gene family
\label{chap:extensins}

## Introduction

\label{sec:extensin_intro}

It is well characterised that G4s have the ability to stall polymerases, including both DNA and RNA polymerases [@Han1999; @Siddiqui-Jain2002; @Dexheimer2006; @Cogoi2006; @Chambers2015; @Kwok2016]. This has been demonstrated through *in vitro* Polymerase stop assays as well as by identifying transcriptionally dependent DNA breaks during G4 ligand treatment [@Rodriguez2012]. Furthermore a number of helicases that are associated with the transcription initiation factor (TFIIH) and elongation, including BLM, WRN and XPD, have been shown to preferentially bind and unwind G4 DNA *in vitro* [@Sun1998; @Mohaghegh2001; @Gray2014]. XPD has also been associated with human TSSs containing PG4s by ChIP-seq [@Gray2014]. This evidence points to an effect of G4s on elongation by Pol II. We showed in \autoref{chap:global_nmm} that PG4 dense genes appeared to have slower elongation of Pol II and hence increased Pol II occupancy at the start of the gene. This reduction in Pol II speed may have knock-on effects on co-transcriptional modifications, such as splicing of mRNA. Furthermore, if G4 stabilisation can be regulated, this could constitute a new mechanism for the regulation of splicing.

It has been estimated that around 80% of all splicing occurs co-transcriptionally in higher eukaryotes [@CarrilloOesterreich2010; @Ameur2011; @Khodor2011; @Girard2012; @Tilgner2012; @Windhager2012]. This is likely due to the coupling of splicing to export and quality control mechanisms [@Straßer2002; @Reed2002; @Fan2017]. Which splice sites are used is highly dependent on Pol II speed, since changes in speed can alter how strong and weak splice junctions compete for assembly of spliceosomes [@DelaMata2010; @Jonkers2015]. The classic model for differential acceptor site usage is shown in (Fig \ref{speed_splice}a). When Pol II is elongating at high speed, acceptor sites which are more strongly canonical but further from the donor in sequence space are favoured. This is because on average, the greater strength of the canonical junctions outweigh the extra time that weaker junctions, which are closer in sequence space, have to be spliced. When Pol II elongates more slowly, however, weak acceptor sites which are more proximal may have much more time to be utilised before stronger distal acceptors are transcribed into the nascent RNA. This tips the balance towards utilisation of weak splice junctions, and can result in alternative splicing (AS) [@DelaMata2010; @Jonkers2015].

\newpage

![**Pol II elongation speed affects co-transcriptional splicing** **a)** Mechanism for effect of Pol II speed on co-transcriptional splicing: when elongation occurs rapidly (first row), the more canonical but distal acceptor site is used, resulting in the exclusion of the alternate exon (shown in green). When Pol II elongates more slowly (second row), there is more time for the weaker proximal site to be utilised, resulting in the inclusion of the alternate exon. **b)** Example mechanism for how G4s could affect splicing. When G4s are not present (top row), the constitutive splice acceptor is used and the green exonic chunk is excluded. When G4s are formed in the template strand of the DNA (second row), these slow down Pol II elongation, allowing the weaker proximal splice junction to be utilised, and including the alternate exon chunk. \label{speed_splice}](figures/polii_speed_splicing.svg)

\newpage

In mammalian systems, the most common form of AS appears to be exon skipping (>40% of AS exons) [@Ner-Gaon2004; @Kim2007]. Intron retention occurs at low rates in human transcriptomes (10% of AS exons), whilst in Arabidopsis, it is much more common (30-35% of AS exons) [@Ner-Gaon2004; @Kim2007]. Marquez et al. recently identified a novel form of intron retention, called the exitron, in human and Arabidopsis transcriptomes [@Marquez2012; @Marquez2015]. Exitrons are intronic regions which more closely resemble exonic sequence in their GC content, are multiples of three in length and do not contain in-frame stop codons [@Staiger2015; @Marquez2012; @Marquez2015]. These occur in coding regions of the mRNA, contributing to the protein sequence when spliced in, and causing truncations of the protein sequence when spliced out. Marquez et al. found exitrons to be relatively common in long exons of Arabidopsis genes [@Marquez2012; @Marquez2015]. The study also found evidence of protein products from these exitrons, indicating that unlike intron retention events, which are thought to be regulatory or caused by splicing errors, exitronic transcripts are functional at the protein level.

Here we identify a novel family of PG4 rich genes, which are downregulated by NMM, and appear to be heavily exitronically spliced. We characterise this splicing and identify whether it can be altered through G4 stabilisation by NMM.

\newpage

## Methods

\label{sec:extensin_methods}

### Plant growth conditions and drug treatments

All experiments used Arabidopsis thaliana Columbia (Col-0) ecotype and all mutant plants used were of Col-0 background. For all experiments seeds were surface sterilised, stratified for 2-3 days at 4°C then sown on vertical plates containing Murashige & Skoog (MS) agar with 1% sucrose and 0.8% agar and transferred to growth cabinets at constant light at 23°C (growth time refers to the time after transfer to growth cabinet). Seedlings used for qPCR and RNAseq were grown for 7 days on MS plates, treated for 6 hours by flooding the plate with MS liquid media containing the drug after which roots and shoots are harvested separately. Drugs used were N-Methyl Mesoporphyrin (Frontier Scientific, NMM580), Berberine (Sigma B3251) and Cyclohexamide (Sigma C7698). For combined Cyclohexamide and NMM treatment, plants were pretreated for 2 hours with Cyclohexamide before adding NMM and treating for a further 6 hours.

### RNA extraction for qPCR and RNAseq

Total nucleic acid isolation protocol was carried out by phenol-chloroform extraction as described by White and Kaper [@White1989]. The resulting pellets were resuspended in sterile water and stored at -80C. The RNA concentration and quality was checked using the NanoDrop 1000 Spectrophotometer (ThermoScientific). DNA contamination was removed using DNase I (Sigma AMPD1).

### RNAseq analysis

RNA was sent to Sheffield Children's Hospital Genomics Facility for library preparation and sequencing. Polyadenylated RNA was enriched using NEBNext Poly(A) mRNA Magnetic Isolation, and libraries were produced using the NEBNext Ultra II Directional RNA kit for Illumina. The chemical fragmentation step was adjusted to increase estimated insert size to the 400-450bp range. Paired end sequencing was conducted across two lanes of an Illumina HiSeq 2500 in rapid run mode, with 220bp read length. The run returned 157 million pass filter reads in lane 1, and 154 million in lane 2. Initial read preprocessing and adaptor trimming was conducted by the Genomics Facility.

Read quality was assessed locally using `FastQC` version 0.11.3 [@Andrews2010]. Differential expression analysis was conducted by using `salmon` version 0.8.2 pseudoalignment [@Patro2017] to estimate transcript abundance against Araport11 cDNA and ncRNA [@Cheng2017]. Library type was automatically determined as ISR (Inward facing pairs, stranded, read one on reverse strand). Mean insert size for each sample was assessed to be in the 400-500bp range. Gene level abundance was then aggregated from this using `tximport` version 1.6.0 [@Soneson2015] and differential expression testing was conducted using `edgeR` linear modelling (version 3.20.7) [@Robinson2010; @McCarthy2012]. Normalised log2 counts per million (CPM) were calculated and used for plotting. P values were adjusted using the Benjamini Hochberg multiple testing correction. Barplots of NMM and DMSO expression were produced using `seaborn` version 0.8.1 [@Waskom2014].

Spliced reads were mapped to the TAIR10 genome [@Initiative2000] using `STAR` version 2.4.2a [@Dobin2013]. The parameters for spliced mapping were adjusted to increase precision: a minimum of 8bp overhang was used for unannotated splice junctions, and 5bp for annotated splice junctions, following ENCODE guidelines [@Encode2016]. Minimum intron size was set to 60bp and maximum intron size to 10000bp. Output BAM files were sorted and indexed using `samtools` version 1.4.1 [@Li2009].

### Gene Ontology analysis

For PG4 enrichment Ontology analysis, two tetrad Quadparser PG4s were predicted in the TAIR10 genome using `g4predict` (available on GitHub https://www.github.com/mparker2/g4predict). The number of PG4s overlapping each strand of the flattened exon supertranscript models [@Davidson2017] for each gene in Araport11 [@Cheng2017] was calculated using `bedtools intersect` version 2.27.1 [@Quinlan2010]. TAIR10 Gene Ontology annotations were downloaded from arabidopsis.org. To calculate enrichment, a permutation test experiment was conducted, where each gene was assigned a weighting proportional to the total length of its exonic sequence. PG4s were then shuffled randomly amongst all genes. For each Ontology group, the number of PG4s observed in genes of that group was compared to the expected numbers if PG4s were distributed randomly across the transcriptome. 10000 permutations were used for testing, and two tailed P values were calculated as $max(min(\frac{\sum_{i=0}^nexp_i < obs}{n}, \frac{\sum_{i=0}^nexp_i > obs}{n}), \frac{1}{n})$, where $n$ is the total number of permutations and $exp_i$ is the expected value from the $i$th permutation. P values were adjusted using the Benjamini Hochberg multiple testing correction.

For Gene Ontology analysis of differentially expressed genes, `GOseq` version 1.30.0 was used [@Young2010]. Up and downregulated genesets were produced using a log2 fold change threshold of 1 and an FDR of 0.05. Weighting factors used were the median transcript length for each gene. P values for enrichment were produced by `GOseq` using the Wallenius approximation method and were corrected for multiple testing using the Benjamini Hochberg method.

Tables of enriched GO terms were generated using `pandas` [@Mckinney2011] and formatted with `inkscape`.

### Extensin gene family total PG4 estimation

For estimated PG4 numbers in the table in Fig \ref{ext_table}, PG4s were predicted using three different methods. All instances of the dinucleotide `GG` were identified in each Extensin gene, and then a graph was built where each `GG` was a node and nodes were connected by an edge if dinucleotides were less than 7bp apart from each other. The number of overlapping PG4 conformations was then calculated as the number of subgraphs in the graph with exactly four members, whilst the number of merged PG4s was calculated as the number of unconnected subgraphs with four or more members. To identify the number of non-overlapping PG4s, a interval scheduling method was used [@Ray]. Overlapping PG4s were grouped and scored by inverse length, then filtered for the maximum number of high scoring, non-overlapping PG4s.

### Quantitative PCR experiments

Total RNA was reverse transcribed into cDNA using the High Capacity cDNA Reverse Transcription Kit with a PolyT primer (Invitrogen Cat. No. 4368813). qPCR was carried out using SYBR® Green JumpStart™ Taq ReadyMix™ (Sigma-Aldrich Cat. No. S4438) using the Mx3000P qPCR System (Agilent Technologies). Thermal cycling conditions consist of a denaturation step at 94 °C for 2 minutes, 40 cycles of 15-second denaturation at 94 °C and 1-minute extension at 60 °C, then a final dissociation step of 2 minutes at 94°C, 1 minute at 60 °C and 2 minutes at 94 °C.

### Analysis of public RNAseq data

Root RNAseq from PRJNA323955 [@Li2016] was downloaded in FASTQ format from ENA. Quality assessment was performed with `FastQC` version 0.11.3 [@Andrews2010] and `fastq-screen` version 0.5.1 [@Wingett2011], and adapter contamination was removed using `Cutadapt` version 1.9.1 [@Martin2011]. The data was mapped using `STAR` version 2.4.2a [@Dobin2013] with default settings, except than max intron length was set to 10000bp. Output BAM files were sorted and indexed using `samtools` version 1.4.1 [@Li2009].

Normalised gene expression estimates were generated using `featureCounts` [@Liao2013] to get raw counts of alignments overlapping each gene in Araport11 [@Cheng2017], and `edgeR` [@Robinson2010; @McCarthy2012] to perform estimation of log2 counts per million.

### Splice junction analyses

Splice junction sites were extracted from aligned reads using `pysam` [@Heger2014]. For all analyses, reads were filtered to produce a set of unique donor/acceptor site pairs. Scatter plots of spliced read percentages and frequency plots of in frame exitrons were produced using `matplotlib` and `seaborn` [@Hunter2007; @Waskom2014].

### Sequence logo generation

Consensus sequence logos were generated using an in-house Python module `matplotlib_logo` (available on GitHub https://www.github.com/mparker2/matplotlib_logo). Unique splice junction donor/acceptor pairs from spliced reads were identified, and the corresponding sequence information (8bp up and downstream of donor and acceptor) was extracted from the TAIR10 genome [@Initiative2000] using `pysam` [@Heger2014]. Position frequency matrices were generated from these sequences, and entropy score in bits was calculated and plotted.

### RNAseq read simulation and bootstrapping experiment

For RNAseq read simulation, read counts generated using `featureCounts` [@Liao2013] for each file in the Li et al. 2016 dataset were used to inform a `polyester` [@Frazee2015] simulation of Illumina 125bp paired end reads, using the Araport11 annotation [@Cheng2017]. Fragment length was simulated using a normal distribution with mean 250bp and standard deviation of 25bp. Error rate was simulated at a uniform 0.5% across samples, reads, and nucleotides.

For the bootstrapping experiment, one or more pairs of real and simulated samples were randomly sampled from the dataset and unique splice donor/acceptor pairs in EXT9 were identified. Only splice sites from primary alignments, with overhangs greater than 20bp, and with more than 2 supporting reads were kept. For each splice site, the exonic overhang sequence 20bp upstream of the donor site and 20bp downstream of the acceptor site were extracted from the read and concatenated. Any splice site whose concatenated sequence was present as a contiguous kmer in the reference EXT9 sequence was assumed to be a mapping error and filtered out. Finally, splice sites were deduplicated by directional clustering of sequences with edit distance of one or less using `umi_tools` directional clusterer [@Smith2017]. 500 iterations were used for bootstrapping. 67% confidence intervals were produced and plotted using `seaborn` [@Waskom2014].

### Mappability analyses

Mappability scores were generated for the TAIR10 genome using `gem-mappability` [@Derrien2012] with a kmer size of 75bp, and converted to `BigWig` format using `gem-2-wig` and `wigToBigWig` [@Kent2010]. Minimum mappability scores for each Extensin gene were extracted using `pyBigWig` [@Ryan2018] and plotted against spliced fraction using `matplotlib` [@Hunter2007].

### Sanger sequencing analysis

For Sanger sequencing, cDNA was produced using the Invitrogen High Capacity cDNA Reverse Transcription Kit with a PolyT primer. Gene specific primers were used to amplify transcripts of interest by PCR. Products were then cloned using the Thermo Scientific CloneJET PCR cloning kit, and transformed into DH5alpha competent cells. Colonies were grown overnight on ampicillin and tested by colony PCR. Picked colonies were then grown up for a further 24 hours in liquid media before miniprepping with Qiagen QIAprep Spin Miniprep Kit. Plasmids were sent, with pJet1.2 forward and reverse sequencing primers, for Sanger sequencing at the Sheffield Hallamshire Hospital Core Genomics Facility. Sequences were aligned to the TAIR10 genome using `BLAT` [@Kent2002].

### Differential Splicing Analysis

For differential splicing analysis, we identified novel splice junctions in our RNAseq dataset to augment the existing Araport11 annotation. Splice junctions for each gene were extracted from reads mapped to correct strand of the annotated gene body using Python and `pysam` [@Heger2014]. Junctions were kept if they were supported by at least 20 reads across all samples. Differential splice junction usage between NMM and DMSO treatments was then conducted in R using `limma-voom` and `DiffSplice` [@Ritchie2015; @Law2014], and P values were adjusted using Benjamini Hochberg correction. Differentially utilised junctions were identified using an absolute log2 fold change in expression of 0.5 and an FDR of 0.2. Unstrigent thresholds were used due to the low read counts associated with splice junctions, to capture as many true positives as possible.

To produce junction categories based on relation to reference annotation, Araport11 GTF annotation [@Cheng2017] was flattened using the python module `CGAT.GTF` [@Sims2014] to produce a single model for each gene, which was converted to bed12 format. This was then used to identify junctions which shared one or both of their donor and acceptor sites with reference introns. These were labelled "alternate" and "constitutive" junctions, respectively. Constitutive junctions which spanned an internal exon where labelled "skipping" junctions. Junctions which were contained wholly within a single contiguous exonic region were labelled as "retained intronic / exitronic" (the distinction between these is whether retention or splicing of the region is more common). Finally, junctions which do not contain a donor or acceptor present in the reference, and which span a mixture of exonic and intronic regions, were labelled "other" junctions. Violin plots of distribution of junction type expression, and stacked barplots of class distribution amongst differentially utilised junctions, were produced using `seaborn` and `matplotlib` [@Hunter2007; @Waskom2014].

For barplots of spliced read percentages in Extensin genes, `pysam` [@Heger2014] was used to extract all uniquely mapped and properly paired reads covering each gene. Each read was counted separately (i.e. pairs were not counted as one fragment). Reads that were mapped across an exitronic splice junction were counted and divided by the total number of reads to get a percentage of exitronic reads. These percentages were compared between NMM and DMSO treatments.

\newpage

### Primer Sequences used

\DTLloaddb{primers}{chapter_6/primers.csv}
\begin{center}
\begin{longtable}{ll}\toprule
    \textbf{Name} & \textbf{Sequence}%
    \DTLforeach*{primers}{\Name=name,\Sequence=sequence}{%
        \DTLiffirstrow{\\\cmidrule{1-2}}{\\}%
        \Name & \texttt{\Sequence}
    }%
\end{longtable}
\end{center}

\newpage

## Results

\label{ssec:extensin_res}

### Gene Ontology shows plant cell wall specific genes are enriched in PG4s and downregulated by NMM

\label{ssec:extensin_go}

To identify gene ontology groups which are specifically enriched with exonic PG4s, we compared the observed levels of PG4s per gene to expected levels if PG4s were randomly distributed across all genes (weighted by gene length). These observed and expected levels were summarised for each gene ontology group. Sorting the results for groups with the greatest positive observed/expected ratio of PG4s on the template strand, we discovered that gene ontology groups involved with functions at the cell periphery, particularly in the plasma membrane and cell wall, had strong enrichments (Fig \ref{go_table}). The log2 fold enrichment in `GO:0005199`, which contains structural cell wall genes, was +4.4 (FDR < 4.8e-4). This corresponded to an observed number of 992 PG4s in only 32 genes (the average expectation under the null hypothesis was 46 PG4s). These PG4 dense gene ontology groups were also strongly enriched in the set of genes which are significantly downregulated by NMM in our RNAseq dataset (Fig\ref{go_table}, 50% of expressed genes in `GO:0005199` were downregulated by NMM, FDR = 9.6e-7).

\newpage

![**Gene Ontology groups enriched in template stranded PG4s** Table showing the top ten Gene Ontology groups most enriched for exonic PG4s compared to null distribution. The top two groups, both containing genes involved in cell wall structure and organisation, are also enriched for genes downregulated by NMM. \label{go_table}](figures/g4_nmm_gene_ontology_table.svg)

\newpage

### The proline rich Extensin gene family contain large numbers of hardcoded PG4s

\label{ssec:extensin_hardcoded}

We discovered that the `GO:0005199` geneset was primarily made up of genes from the Extensin cell family (29/32 genes, 90.6%), including classical SP4/SP5 Extensin genes and chimeric Leucine Rich Repeat/Extensin (LRX) genes [@Showalter2010; @Liu2016]. These genes were found to be extremely PG4 rich on the template strand, with many genes containing greater than 10 PG4s per kilobase of exon (Fig \ref{ext_genes}a). Upon visualisation of these genes, we noted that in the majority of cases the PG4s were regularly spaced along the gene, and were contained solely within the coding region (CDS) of the gene (Fig \ref{ext_genes}b).

\newpage

![**Expression and PG4 density of genes in the Cell Wall Structural Ontology group `GO:0005199`** **a)** Panels showing gene expression (top panel) and PG4 density (bottom panel) for genes in the `GO:0005199` group. Expression in DMSO (blue) and NMM (orange) conditions is shown at log2 counts per million (from root RNAseq dataset). Errorbars are standard deviation of three biological replicates. Genes which are differentially expressed with FDR < 0.05 are labelled with asterisks. In PG4 panel, the exonic PG4 density per kilobase is shown separately for coding (blue) and template (orange) strands of the gene. **b)** Gene tracks showing the location of predicted two tetrad PG4s in orange for (from top to bottom) LRX1, EXT13, and EXT9. Gene models from Araport11 are shown in blue. In gene models, thin boxes represent untranslated regions (UTRs), fat boxes represent coding regions (CDS), and connecting lines represent intronic regions. \label{ext_genes}](figures/extensin_gene_ontology_group_expression_g4s.svg)

\newpage

From a search of the literature, we discovered that Extensin genes are highly repetitive, proline rich proteins [@Kieliszewski1994; @Showalter2010; @Liu2016]. These proteins polymerise to function as a structural matrix in the protein component of the plant cell wall [@Lamport1960; @Lamport1965; @Showalter1993]. In particular, we noted that these proteins are characterised by large numbers of the SP3-5 repeat, which is made up of the sequence $Ser(Hyp)_{3-5}$, where Hyp is Hydroxy-proline, a proline derivative [@Showalter2010]. Since the codon for proline is `CCN`, the DNA which encodes SP4 and SP5 motifs will conform to the two tetrad Quadparser motif on the template strand of the gene (Fig \ref{cd_spec}a). This is the source of the PG4 density of Extensin genes, and PG4 counts in these genes are well correlated with SP3-5 repeats (Fig \ref{ext_table}). Since the SP4 motif is required for the function of the protein, and is restricted by the codon for proline, these PG4s are "hardcoded" into the body of the gene.

\newpage

![**The Extensin gene family contains large numbers of hardcoded PG4s** Table showing extended Extensin gene family, their expression patterns, SP4 motif counts, PG4 counts and expression during NMM treatment. Adapted from Showalter et al. 2010 \label{ext_table}](figures/ext_family_showalter_table.svg){height=750px}

\newpage

To demonstrate that the PG4 from Extensin genes could form a G4 structure in vitro we used circular dichroism spectroscopy (CD). We performed these experiments at physiologically relevant temperatures for Arabidopsis. An oligo representative of the SP4 repeat was designed (`AGAGGTGGTGGTGGTATG`) using 3bp flanks upstream and downstream of the PG4. CD showed the G4 oligo had peak absorbance at 260nm and trough at 240nm, indicative of a parallel G4 structure (Fig \ref{cd_spec}b). Removing the PG4 by mutating the sequence (mutated sequence: `AGAGGTGATGGTGGTATG`) or removing potassium ions from the buffer abolished this absorbance profile.

\newpage

![**The Extensin SP4 motif forms a G-Quadruplex *in vitro*.** **a)** Schematic showing how the Extensin SP4 protein motif hardcodes a two tetrad PG4 into the template strand of the gene body. **b)** CD spectroscopy of an Extensin repeat sequence (left) and a mutated control which does not conform the the Quadparser motif (right) show that the Extensin repeat forms a G4 *in vitro*. This is indicated by the peak in ellipticity at 260nm and the trough at 240nm, which are characteristic of a parallel G4. \label{cd_spec}](figures/cd_spectroscopy1.svg){height=750px}

\newpage

### Extensins are strongly downregulated by NMM and Berberine

\label{ssec:extensin_downreg}

To confirm that the Extensin genes are downregulated by NMM, we performed RNA extraction and quantitative RT-PCR (qPCR) on root tissue from 7 day old Arabidopsis seedlings treated for 6 hours with NMM at varying concentrations. EXT13 and LRX1 were chosen as representative classical and chimeric Extensins, respectively [@Showalter2010]. The change in expression of both genes upon treatment was negatively correlated with the concentration of NMM applied (Fig \ref{nmm_berb_qpcr}a). Treatment with the G4 intercalating drug, Berberine, also caused strong downregulation of EXT13 and LRX1 (Fig \ref{nmm_berb_qpcr}b). Since NMM and berberine are very different drugs which stabilise G4s through different methods, taken together our results suggest downregulation of Extensins is caused by G4 stabilisation.

### Downregulation of Extensins by NMM is translation independent

\label{sec:extensin_cyclohex}

To confirm whether downregulation of EXT13 and LRX1 by NMM was direct, or the result of a perturbation the levels of a transcription factor, we conducted qPCR experiments with combinatorial treatment of NMM and Cyclohexamide (CHX). CHX is an inhibitor of translation which is commonly used to determine whether interactions by transcription factors on a gene are direct [@William2003; @Taverner2004]. If effects are indirect (i.e. if the transcription factor of interest regulates transcription of some intermediate transcriptional factor, which regulates the gene of interest), then treatment with CHX will prevent regulation, since any intermediate factors will not be able to be translated. In the case of NMM treatment, this was used to see whether NMM acts directly on EXT13 and LRX1 through G4 stabilisation, or through other changes in the transcriptome. Seedlings were pretreated with CHX for two hours, before NMM was added and treatment was continued for another 6 hours. This experiment showed that that Extensin downregulation by NMM still occurs even when translation is blocked, suggesting that NMM acts directly upon the Extensin genes.

\newpage

![**Expression of EXT13 and LRX1 during treatment with G4-binding ligands** Scatter/strip plots showing qPCR results for EXT13 (left panels) and LRX1 (right panels). Log2 fold change in expression (ΔΔCT) of Extensin genes decreases with increasing concentrations of **a)** NMM and **b)** Berberine. **c)** NMM downreguation of EXT13 and LRX1 is not affected by concurrent Cyclohexamide treatment, suggesting a mechanism independent of translation. For all panels, each point is a single technical replicate, and colours represent different biological replicates. A small amount of jitter has been added to the X axis for better visualisation of results. \label{nmm_berb_qpcr}](figures/ext13_lrx1_qpcr.svg){height=750px}

\newpage


### RNAseq suggests Extensin genes contain exitronic splice sites

\label{ssec:extensin_splice_sites}

We noted from studying *de novo* assembled splice isoforms from a root specific RNAseq dataset [@Li2016] that many of the Extensin genes had large numbers of novel spliced isoforms. EXT9 was found to have the most novel spliced forms of any gene in the dataset. These were not present in the annotation. In fact, the majority of Extensin domains are annotated in both TAIR10 and Araport11 as intronless. These splice isoforms are presumably therefore a product of "exitronic" splicing [@Marquez2012; @Marquez2015], where sections of constitutive exons, flanked on both sides by exonic sequence, are spliced out of a gene. We hypothesised that these unusual exitrons could be a result of slow Pol II elongation through PG4 dense regions, allowing splicing to occur at weak splice sites.

A hallmark of most true splice junctions is the GT/AG intron motif, which is the conserved canonical sequence in all higher eukaryotes, including Arabidopsis (Fig \ref{splice_junct}a) [@Mount1982; @Shapiro1987]. To determine whether Extensin exitrons had canonical splice motifs, we produced splice junction sequence logos for predicted introns from the dataset produced by Li et al., for EXT9 and LRX3, both of which were highly spliced. We found that splice junctions in these genes had near universal GT/AG motifs (Fig \ref{splice_junct}b). Upon inspection of the methods for Li et al., however we discovered that `CuffLinks` was used for *de novo* transcript assembly [@Trapnell2012]. Since the RNAseq dataset is unstranded, `Cufflinks` requires the upstream mapping tool (here, `STAR` [@Dobin2013]) to annotate the orientation of spliced reads using the intron motif (i.e. positive strand for GT/AG and negative strand for CT/AC). This setting means reads which do not conform to the intron motif are discarded, leading to serious bias.

To remove this bias, we remapped reads from the Li et al. dataset using STAR without filtering by intron motif. Since assemblers like `CuffLinks` and `StringTie` [@Pertea2015] require strandedness information derived from the intron motif, transcript assembly was not possible. Instead, we simply extracted spliced reads aligning to EXT3 and LRX3 and identified all unique splice site starts and ends. The corresponding sequences were then used to produce sequence logos (Fig \ref{splice_junct}c). These logos showed only a weak enrichment for the GT/AG motif in EXT9, and CT/AC in LRX1.

\newpage

![**Splice junction motifs for EXT9 and LRX3** Sequence logo plots showing consensus splice site sequences around donor (left panels) and acceptor (right panels) splice sites. Putative intronic sequences are shown on shaded blue background. **a)** Splice junction consensus sequence logo for Arabidopsis, calculated from junctions in the Araport11 annotation. **b)** Splice junction consensus sequence logo produced from **de novo** assembled transcripts (Li et al. 2016) for EXT9 and LRX3. **c)** Splice junction consensus sequence logo produced from unique donor/acceptor pairs identified from spliced reads on EXT9 and LRX3. \label{splice_junct}](figures/splice_site_sequence_logo.svg)

\newpage

Since the Extensin exitrons appear in coding regions of DNA, if the spliced out region is not a multiple of three, then the resulting mRNA would be frameshifted, producing truncated and potentially deleterious proteins. We therefore tested whether the spliced reads in EXT9 and LRX3 contained gaps which were multiples of three or not. For both genes, almost all of the unique splice junction pairs were multiples of three (Fig \ref{splice_frame}a-b). This could be evidence that these exitrons are genuine and produce function gene products. On the other hand, we noted that splicing tended to occur between regions with high protein and DNA sequence level homology. EXT9, for example, is an incredibly repetitive gene with high self-homology (Fig \ref{ext_mapp}a). These in frame splice junctions could therefore simply be the result of mapping errors from the spliced aligner `STAR` [@Dobin2013], which utilises heuristics which may result in some reads from contiguous parts of the genome being mapped as spliced. If the homologous regions which could cause mapping errors within a gene also have homology at the protein level, as is the case in EXT9, then it is probable that erroneously spliced reads would be a multiple of three in intron length.

If spliced reads mapping to EXT9 were the result of some systematic error in mapping, one might not expect to see much variation in the percent of reads mapping to a gene being spliced. We therefore correlated the expression of EXT9 in each sample from the root RNAseq dataset (measured in log2 counts per million or needs to be!) with the percent of reads which mapped with a splice site. We found a slight positive correlation between expression and splicing (Fig \ref{splice_frame}c).

As a further precaution against these erroneous spliced mappings, we performed read simulation for each sample in the root RNAseq dataset. Put simply, the expression of each gene was quantified for each sample by counting the number of mapped reads, then `Polyester` (an Illumina sequencing read simulator) was run to generate reads from the reference transcriptome with the same read counts [@Frazee2015]. These simulated reads were then remapped with `STAR` using the same parameters as the original mapping [@Dobin2013]. We then performed a bootstrap analysis for EXT9 where we sampled one or more real/simulated sample pairs, and counted the number of unique splice donor/acceptor pairs that occurred in each. Junctions with the same exonic flanking sequence (using 20bp overhangs) or with edit distance of only one base were collapsed. Any junctions with flanking sequence that appeared as a contiguous kmer in the reference sequence of EXT9 were also removed. Despite this, we saw a consistently larger number of unique donor/acceptor splice pairs in the real data than in the simulated data (Fig \ref{splice_frame}d).

\newpage

![**Splice junction motifs for EXT9 and LRX3** **a** & **b)** Frequency barplots showing number of In/Out of frame splice junctions for EXT9 and LRX3 respectively. **c)** Scatterplot showing percentage of reads with splicing versus log2 counts per million for EXT9. **d)** Bootstrapped splicing simulation showing number of unique EXT9 splice junctions discovered with increasing numbers of samples for real root RNAseq data versus paired simulated RNAseq data. Errorbars are 67% confidence intervals. \label{splice_frame}](figures/splice_site_frame_and_simulation.svg)

\newpage

Since the Extensin genes are highly repetitive, this reduces the ability of read aligners to map to them. This "mappability" can be quantified using tools such as `GEM` which measure, for each genomic position, how often the sequence kmer that is found there occurs in the rest of the genome [@Derrien2012]. We utilised `GEM` to score the mappability of the Arabidopsis genome, and compared the median mappability score for extensin genes to the percent of spliced reads for each gene. Only genes which were annotated as intronless (i.e. no constitutive introns, genes with exitrons were allowed) were included. We found a clear negative correlation between mappability and the number of mapped spliced reads (Fig \ref{ext_mapp}b). For EXT9, the regions with lowest mappability are clearly those with the most annotated splice sites, including in the Araport11 reference (Fig \ref{ext_mapp}c). 

\newpage

![**Extensin genes with greater spliced mapped reads have low mappability** **a)** Dotplot showing self homology of the EXT9 cDNA (unspliced isoform). Positional identity was calculated using 15bp windows across the gene. Positions with identity less than 75% were filtered to remove noise from the plot. **b)** Scatter plot showing the minimum mappability score of unspliced Extensin genes  against the percentage of spliced reads for that gene in mature root RNAseq samples (Li et al. 2016). Errorbars are standard deviation of three biological replicates. **c)** Gene track showing the mappability score across EXT9 (orange). Most splice forms cross these low mappability regions, including in the reference annotation Araport11 (shown in blue). \label{ext_mapp}](figures/ext9_dotplot_mappability.svg)

\newpage

### Sanger sequences identifies LRX1 and EXT9 splice variants

\label{ssec:extensin_sanger}

To experimentally confirm whether Extensin gene exitron splicing exists, we performed RT-PCR of LRX1 and EXT9 mRNAs. PCR products of both genes showed multiple products, characteristic of several spliced forms (data not shown). PCR products which did not correspond to the full length of the unspliced mRNA were gel extracted, cloned and sanger sequenced to identify their origin. We identified a number of mRNA fragments originating from the LRX1 and EXT9 genes. Alignment of these products using `BLAT` identified a number of spliced isoforms in both genes [@Kent2002]. To identify whether these isoforms contained canonical splice sites, we produced sequence logos. Neither gene showed a clear pattern conforming to the canonical intron motif GT/AG, though the products from  LRX1 showed the reverse complement of this pattern, CT/AC.

\newpage

![**Sanger sequencing of LRX1 and EXT9 cDNA identifies spliced forms** **a)** Gene track showing aligned sanger sequencing products for **a)** LRX1 and **b)** EXT9. Products aligned to the forward strand are shown in green, and products aligned to the negative strand are shown in orange. Gene models are from the Araport11 annotation. **c)** Sequence logos for sanger product splice junctions for LRX1 (top panel) and EXT9 (lower panel). \label{sanger}](figures/sanger_splice_variants.svg)

\newpage

### NMM treated plants do not have increased splicing of Extensin genes

\label{ssec:extensin_nmm_splice}

We hypothesise that G4s cause the exitronic splicing of Extensin genes, by slowing down Pol II elongation. To test whether increased stabilisation of G4s changes this splicing pattern, we performed RNAseq of root tissue from plants treated with NMM using 220bp paired reads, to identify novel splicing isoforms. Mapping parameters for `STAR` [@Dobin2013] were made more stringent than defaults in an attempt to increase the precision of mapping over Extensin genes without attenuating recall of splice junctions too strongly. A common method for conducting differential splicing analysis is to use differential exon usage methods such as `DEXseq` [@Anders2012]. This is sometimes conducted on exon "chunks" which are the contiguous genomic ranges which each appear in a distinct set of transcripts of a gene. Since there are many overlapping exitrons in the Extensin genes, these chunks would be extremely short and be complex to interpret. Furthermore, if spliced transcripts are in lower abundance than full length transcripts, a large change in the use of a particular exitron may only lead to a small change in the expression of the exon chunk which is being spliced out. We therefore opted for counting the number of reads which support each unique junction in a gene (junction counts), and performing differential junction usage on this. The downside to this analysis is that the number of reads supporting each junction may be lower than the number per exon, leading to reduced power to detect differential usage.

In order to perform junction level differential splicing analysis, we identified all spliced reads with the correct first-in-pair strand orientation overlapping each gene. Only splice junctions with at least 20 supporting reads total across the 6 samples were kept for analysis. Splice junctions were categorised into five types based on how they related to the flattened reference annotation in Araport11: constitutive, alternative, retained intron/exitronic, exon skipping, or other. See Fig \ref{diff_junc}a legend for an explanation of these different categories.

We performed differential junction usage identification using `limma-voom` and `limma-diffSplice` [@Ritchie2015; @Law2014]. Unstrigent thresholds were used due to the low read counts associated with splice junctions, to capture as many true positives as possible. Using an FDR threshold of 0.2 and an absolute fold change threshold of 0.5, we identified 338 junctions in 302 genes with increased use during NMM treatment, and 189 junctions in 162 genes with decreased usage. 27 genes contained junctions which showed both increased and decreased usages, i.e. some type of switching. None of the junction isoforms identified in the Extensin genes were significantly differentially utilised, however.

\newpage

![**Differential Junction Usage during NMM treatment** **a)** Schematic showing the five classes used to categorise splice junctions based on flattened reference annotation. Junctions which are present as introns in the reference are labelled constitutive junctions. Constitutive junctions which cause skipping of an exon are labelled skipping junctions. Junctions which share a donor or acceptor site with one in the reference are labelled alternate junctions. Junctions which are wholly contained within a single exon of the reference are labelled retained/exitronic junctions. Finally, junctions which share no donor or acceptor with the reference and span a mixture of exonic and intronic sequence are labelled "other". **b)** Violin plot showing distribution of average expression in log2 counts per million for each junction class. **c)** Proportion of junction in each class for all splice junctions vs. those with significantly increased or decreased usage during NMM treatment (Absolute logFC > 0.5, FDR < 0.2). **d)** Categorisation of detected junctions in Extensin genes. \label{diff_junc}](figures/diff_splice_results.svg)

\newpage

The linear model fit by `diffSplice` is looking for differences in the usage of a junction relative to the usage of the rest of the junctions in that gene. It is possible therefore, if utilisation of all splice junctions in the Extensin genes were changed by a relatively similar amount, that `diffSplice` would not detect any differentially utilised junctions. We therefore examined the total percentage of reads mapping to each Extensin gene which were exitronic, and looked to see if this percentage changed during NMM treatment. Despite large gene-level effects on many Extensin genes with Extronic splicing (Fig \ref{perc_spliced}a), there was no strong effect on the overall level of spliced reads mapping to these genes (Fig \ref{perc_spliced}b), suggesting that either NMM does not affect the splicing of Extensins, only the expression, or that spliced mapping to these genes is a systematic mapping error that occurs at approximately the same rate per gene regardless of the read count.

\newpage

![**NMM does not cause significant changes in the amount of Exitronic Splicing of Extensins** Barplots showing **a)** gene expression in log2 counts per million, and **b)** percentage of reads exhibiting Exitronic splicing, for expressed Extensin genes with an average of at least 5% spliced reads. Levels during mock (DMSO) treatment are shown in blue, and levels during NMM treatment are shown in orange. Ordering of genes (by average percent spliced) in upper and lower panels is matched. Errorbars are standard deviation of three biological replicates. \label{perc_spliced}](figures/extensin_percent_spliced.svg)

\newpage

## Discussion

\label{sec:extensin_discuss}

We have identified a new family of G4 regulated genes: the Extensins. This is an extremely unusual and difficult to study family of genes, due to high levels of self homology and paralogy. The repetition which makes Extensins less tractable is also what makes them interesting to us: it is caused by large numbers of poly-proline rich sequences, called SP4 motifs, which give rise to template stranded two tetrad PG4s. Our *in vitro* CD spectroscopy work suggests that these PG4 sequences indeed form G4s which may be stable in the physiological temperature range that Arabidopsis lives at.

SP4 motifs are found throughout the coding regions of Extensin domains, usually in regularly spaced patterns. Since they contain large numbers of the motif, Extensins are some of the most PG4 dense genes in Arabidopsis, potentially in any organism. The density of PG4 sequences may also lead to large combinatorial increases in the G4 topologies that can form: an analysis of all the possible overlapping PG4 structures in Extensins suggests that some have more than 500 potential conformations per kilobase. Greater numbers of potential G4 conformations of these sequences may increase the entropy and therefore the thermodynamic stability of the folded population. This seems likely to have biological implications in transcription or regulation of the genes.

The majority of Extensin family genes exhibit a sharp reduction in gene expression when plants are treated with the G4-stabilising ligand NMM. We have demonstrated this robustly through independent analyses of microarray, RNAseq and qPCR data. Furthermore, treatment with another G4 binding ligand, Berberine, also downregulated selected Extensins, shown by qPCR. Berberine has a very different structure and binds G4s through a totally different mechanism to NMM, suggesting that this is likely to be a direct G4 effect. As we showed in \autoref{chap:global_nmm}, NMM appears to have a strong global effect on the expression of genes with template stranded PG4s. Whilst this effect may be in part attributable to the Extensin genes, it is not solely due to them, since 50% more NMM downregulated genes than might be expected by chance contain template stranded PG4s. Our hypothesis is that template stranded PG4s cause blockages which slow the progression of Pol II along the gene during transcription. When NMM is added, this may cause increased slowing or even premature termination of transcription, resulting in decreased expression. We demonstrated that the effect of NMM is translation independent (i.e. not caused by transcription factor activity) by treating plants with both Cyclohexamide and NMM at the same time, and showing that NMM still caused downregulation.

Examination of *de novo* assembled splice isoforms collated from many RNAseq samples in Li et al 2016 led us to speculate about splicing of Extensin genes [@Li2016]. EXT9 was found to have the most different splice isoforms of any gene in this dataset. These splice sites tended to occur within exonic sequence and remove regions which on the template strand were PG4 rich. Analysis of how these splicing patterns might affect the protein encoded by the mRNA showed that they were multiples of three and would not cause frameshifts. Instead they would simply cause greater variation in the length of protein products. This is an intriguing idea since the Extensins are structural components of the cell wall protein matrix, and it is possible that different lengths of building block might result in differences in flexibility of the wall. Furthermore, previous work has shown that truncated versions of LRX1 with fewer SP4 repeats still function and are able to rescue an *lrx1* null mutation [@Baumberger2001]. It became clear, however, that the number and type of splice junctions which could be discovered in Extensin genes was highly dependent on parameters used in mapping, suggesting that at least some of these splice junctions are spurious and caused by technical error.

To try to identify true splice isoforms of Extensin genes, we performed RNAseq, using longer 220bp paired end reads to attempt to capture more splice junctions and reduce the number of multimapping reads. Since splicing occurs co-transcriptionally, and Pol II elongation speed is linked to utilisation of weak splice sites, we hypothesised that splicing of Extensins might be controlled by formation of G4s which slow down transcription. We therefore also performed RNAseq with NMM treated plants to see how G4 stabilisation affects splicing. Our RNAseq dataset identified large numbers of potential splice variants, however again many of these did not have the hallmarks of canonical splice junctions, and their abundances were highly sensitive to changes in mapping parameters, suggesting the possibility of some technical defect. Furthermore we did not see any strong or consistent change in the percentage of spliced reads identified when plants were treated with NMM, despite some very strong changes in the expression of Extensin genes as a whole. Despite these negative results, we were able to identify PCR products from cDNA which appeared to be truncated forms of both EXT9 and LRX1. Again, these products did not have any consistent traits of splice variants. One explanation for these results, and the results of the RNAseq, could be artefacts introduced during PCR amplification. These can occur in repetitive regions due to incorrect annealing of single stranded DNA and can result in deletions or expansions. To overcome these issues in the future, the ideal techniques for detecting true splice forms would be Northern blots on a gene by gene basis, or using direct RNA Nanopore sequencing for global identification.

\newpage
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJydkluUXlVNFlRblg5V3g5Ijp7In
RleHQiOiJYUEQgaGFzIGFsc28gYmVlbiBhc3NvY2lhdGVkIHdp
dGggaHVtYW4gVFNTcyBjb250YWluaW5nIFBHNHMgYnkgQ2hJUC
1zZXEiLCJzdGFydCI6NzgyLCJlbmQiOjg1NH0sInM1Mmt0VWtC
dFZOT1IxM3YiOnsidGV4dCI6ImhpZ2hlciIsInN0YXJ0IjoxND
I5LCJlbmQiOjE0MzV9LCJrWTVvVE1ranpNQmQ4UUJ6Ijp7InRl
eHQiOiJtb3JlIHN0cm9uZ2x5IGNhbm9uaWNhbCIsInN0YXJ0Ij
oyMDMxLCJlbmQiOjIwNTR9LCJEcTlyUzBMaXBKcEhMVVgxIjp7
InRleHQiOiJnICg+NDAlIG9mIEFTIGV4b25zKS4iLCJzdGFydC
I6MzYzNiwiZW5kIjozNjgzfSwieUlaNVVCTWdhMlNIbjRHWSI6
eyJ0ZXh0IjoiIyMgTWV0aG9kcyIsInN0YXJ0Ijo0OTk5LCJlbm
QiOjUwMDl9LCIwNUJtR2xpWGRkMm5oNW82Ijp7InRleHQiOiJE
cnVncyB1c2VkIHdlcmUgTi1NZXRoeWwgTWVzb3BvcnBoeXJpbi
AoRnJvbnRpZXIgU2NpZW50aWZpYywgTk1NNTgwKSwgQmVyYmVy
aW5l4oCmIiwic3RhcnQiOjU3MjcsImVuZCI6NTg1Mn0sIk9rbW
1oRDRoaktuZzVkZFgiOnsidGV4dCI6IiMjIyBSTkFzZXEgYW5h
bHlzaXMiLCJzdGFydCI6NjQyNSwiZW5kIjo2NDQ0fSwiQ09xcG
RkYlNyZXlxUlBQUiI6eyJ0ZXh0IjoiUmVhZCBxdWFsaXR5IHdh
cyBhc3Nlc3NlZCBsb2NhbGx5IHVzaW5nIGBGYXN0UUNgIiwic3
RhcnQiOjcxMTcsImVuZCI6NzE2NX0sIlBBVEVaZ0FtS09aUzF1
c2wiOnsidGV4dCI6InNpbmcgYHNhbG1vbmAgcHNldWRvYWxpZ2
5tZW50IFtAUGF0cm8yMDE3XSIsInN0YXJ0Ijo3MjQ4LCJlbmQi
OjczMDR9LCIzTXV1RUJvSHVORnZsRDdjIjp7InRleHQiOiJNZW
FuIGluc2VydCBzaXplIGZvciBlYWNoIHNhbXBsZSB3YXMgYXNz
ZXNzZWQgdG8gYmUgaW4gdGhlIDQwMC01MDBicCByYW5nZS4iLC
JzdGFydCI6NzQ5NSwiZW5kIjo3NTcwfSwiY3hXbGdQMWltZUEz
c2lRMCI6eyJ0ZXh0IjoiZW5lIGxldmVsIGFidW5kYW5jZSB3YX
MgdGhlbiBhZ2dyZWdhdGVkIGZyb20gdGhpcyB1c2luZyBgdHhp
bXBvcnRgIFtAU29uZXNvbjIwMeKApiIsInN0YXJ0Ijo3NTcyLC
JlbmQiOjc4NzZ9LCJnV3RBaXRzMFBPcWVZRGs0Ijp7InRleHQi
OiJUaGUgcGFyYW1ldGVycyBmb3Igc3BsaWNlZCBtYXBwaW5nIH
dlcmUgYWRqdXN0ZWQgdG8gaW5jcmVhc2UgcHJlY2lzaW9uIiwi
c3RhcnQiOjgxNjIsImVuZCI6ODIzMn0sIm5FZDRZbVZHSk1RTz
E3OEIiOnsidGV4dCI6ImBTVEFSYCIsInN0YXJ0Ijo4MTI2LCJl
bmQiOjgxMzJ9LCJJYWc4ZkRyM0NBUEo4Tm9VIjp7InRleHQiOi
Jgc2FtdG9vbHNgIiwic3RhcnQiOjg1MDUsImVuZCI6ODUxNX0s
ImxEV2l6YVZGTEE1TUtCVngiOnsidGV4dCI6ImZsYXR0ZW5lZC
BleG9uIG1vZGVscyIsInN0YXJ0Ijo4ODA0LCJlbmQiOjg4NDF9
LCJia0JQdFNxU21GbEVCcGFUIjp7InRleHQiOiJGb3IgZWFjaC
BPbnRvbG9neSBncm91cCwiLCJzdGFydCI6OTI2MiwiZW5kIjo5
Mjg2fSwiRFFialN3UTc3cElheW5nbCI6eyJ0ZXh0IjoiVG8gY2
FsY3VsYXRlIGVucmljaG1lbnQsIGEgcGVybXV0YXRpb24gdGVz
dCBleHBlcmltZW50IHdhcyBjb25kdWN0ZWQsIHdoZXJlIGVhY+
KApiIsInN0YXJ0Ijo5MDQyLCJlbmQiOjk4MDN9LCIxakRyT1dz
bWZNS0NnQnhTIjp7InRleHQiOiJHT3NlcWAiLCJzdGFydCI6OT
g2OCwiZW5kIjo5ODc0fSwiS2FJMnpqUGZZclBJNlhhMiI6eyJ0
ZXh0IjoibWVkaWFuIHRyYW5zY3JpcHQgbGVuZ3RoIGZvciBlYW
NoIGdlbmUuIiwic3RhcnQiOjEwMDQ5LCJlbmQiOjEwMDg4fSwi
UU01NUFpa1NWQ3ZBYXoyTCI6eyJ0ZXh0IjoiYSBkeW5hbWljIH
Byb2dyYW1taW5nIG1ldGhvZCB3YXMgdXNlZCIsInN0YXJ0Ijox
MTA0NCwiZW5kIjoxMTA4MX0sIkhKem0xZG1ITXVqVHVqeUEiOn
sidGV4dCI6ImBDdXRhZGFwdGAiLCJzdGFydCI6MTIxMDIsImVu
ZCI6MTIxMTJ9LCI0TVl5QzdnV0ljMFBkNnhYIjp7InRleHQiOi
JTcGxpY2UganVuY3Rpb24gc2l0ZXMgd2VyZSBleHRyYWN0ZWQg
ZnJvbSBhbGlnbmVkIHJlYWRzIHVzaW5nIGBweXNhbWAiLCJzdG
FydCI6MTI2NTEsImVuZCI6MTI3MjB9LCJzUjl2c3c1RXZNNTJF
WFVFIjp7InRleHQiOiJvciBhbGwgYW5hbHlzZXMsIHJlYWRzIH
dlcmUgZmlsdGVyZWQgdG8gcHJvZHVjZSBhIHNldCBvZiB1bmlx
dWUgZG9ub3IvYWNjZXB0b3LigKYiLCJzdGFydCI6MTI3MzYsIm
VuZCI6MTI4MjR9LCJTUUxlZTlNTE5ST3NNYjFiIjp7InRleHQi
OiJgbWF0cGxvdGxpYl9sb2dvYCIsInN0YXJ0IjoxMzA4OCwiZW
5kIjoxMzEwNX0sIlNXZ0c1bkttN0xKSEFIbEoiOnsidGV4dCI6
IlVuaXF1ZSBzcGxpY2UganVuY3Rpb24gcGFpcnMiLCJzdGFydC
I6MTMxNzcsImVuZCI6MTMyMjB9LCJlcGpJSGMwaXpWYkpVZkZE
Ijp7InRleHQiOiJudHJvcHkgc2NvcmUgaW4gYml0cyIsInN0YX
J0IjoxMzQ5OSwiZW5kIjoxMzUxOX0sIkJjakJGWWtGSjZnZHZD
eDMiOnsidGV4dCI6IlJOQXNlcSByZWFkIHNpbXVsYXRpb24gYW
5kIGJvb3RzdHJhcHBpbmcgZXhwZXJpbWVudCIsInN0YXJ0Ijox
MzU1MywiZW5kIjoxMzYwNH0sImNieFZGb1dVZXQ1VVJGVEgiOn
sidGV4dCI6IkZvciB0aGUgYm9vdHN0cmFwcGluZyBleHBlcmlt
ZW50LCBvbmUgb3IgbW9yZSBwYWlycyBvZiByZWFsIGFuZCBzaW
11bGF0ZWQgc2FtcGzigKYiLCJzdGFydCI6MTQwNjQsImVuZCI6
MTQ1MzF9LCJSOFdRVmpZT296Z3QwSWViIjp7InRleHQiOiJgZ2
VtLW1hcHBhYmlsaXR5YCIsInN0YXJ0IjoxNTA4OCwiZW5kIjox
NTEwNX0sIk9TWm5mWUdFQ1hGYW5uTlciOnsidGV4dCI6ImBweU
JpZ1dpZ2AiLCJzdGFydCI6MTUzMDAsImVuZCI6MTUzMTB9LCI2
eXdaMUgxUzY1bExNcWJmIjp7InRleHQiOiJgQkxBVGAiLCJzdG
FydCI6MTYxNjcsImVuZCI6MTYxNzN9LCJLeXdVOHNTUDhyYzlR
TGpTIjp7InRleHQiOiJhdCBsZWFzdCAyMCByZWFkcyIsInN0YX
J0IjoxNjU1OSwiZW5kIjoxNjU3Nn0sIndiOEhUa2V1YjV3c0xj
UTYiOnsidGV4dCI6IkRpZmZlcmVudGlhbCBzcGxpY2UganVuY3
Rpb24gdXNhZ2UgYmV0d2VlbiBOTU0gYW5kIERNU08gdHJlYXRt
ZW50cyB3YXMgdGhlbiBjb27igKYiLCJzdGFydCI6MTY1OTcsIm
VuZCI6MTY3NDh9LCJ3eGJReU1tbjZPUTgxTTNWIjp7InRleHQi
OiJGRFIgb2YgMC4yLiIsInN0YXJ0IjoxNjkyNywiZW5kIjoxNj
kzOH0sInh5eGwxbDBERnlaV0VLRDQiOnsidGV4dCI6InNpbmds
ZSBjb250aWd1b3VzIGV4b25pYyByZWdpb24iLCJzdGFydCI6MT
c2NjksImVuZCI6MTc3MDB9LCJ6MEo5Y2Jkdnk2OEo2cG1DIjp7
InRleHQiOiJ1bmlxdWVseSBtYXBwZWQiLCJzdGFydCI6MTgzND
YsImVuZCI6MTgzNjF9LCJxOTNiVFFhd2IyVlNEWHRaIjp7InRl
eHQiOiJQcmltZXIgU2VxdWVuY2VzIHVzZWQiLCJzdGFydCI6MT
g3MjEsImVuZCI6MTg3NDJ9LCJwRGhoU1ZnNTVYdGlHSElSIjp7
InRleHQiOiJvdXIgUk5Bc2VxIGRhdGFzZXQiLCJzdGFydCI6Mj
AyMDksImVuZCI6MjAyMjd9LCJkRm5JbFBpTUt0b1djSW9kIjp7
InRleHQiOiJUbyBkZW1vbnN0cmF0ZSB0aGF0IHRoZSBQRzQgZn
JvbSBFeHRlbnNpbiBnZW5lcyBjb3VsZCBmb3JtIGEgRzQgc3Ry
dWN0dXJlIGluIHZp4oCmIiwic3RhcnQiOjIzODgyLCJlbmQiOj
I0MDA5fSwiZGEwUjc5TW5qemkxaEFOTiI6eyJ0ZXh0IjoiVG8g
Y29uZmlybSB0aGF0IHRoZSBFeHRlbnNpbiBnZW5lcyBhcmUgZG
93bnJlZ3VsYXRlZCBieSBOTU0sIHdlIHBlcmZvcm1lZCBSTkEg
ZeKApiIsInN0YXJ0IjoyNTIwOCwiZW5kIjoyNTMyMH0sImJYU0
huMlV1TXhXSVJKYWkiOnsidGV4dCI6ImMpKiogTk1NIGRvd25y
ZWd1YXRpb24gb2YgRVhUMTMgYW5kIExSWDEgaXMgbm90IGFmZm
VjdGVkIGJ5IGNvbmN1cnJlbnQgQ3ljbG9oZXjigKYiLCJzdGFy
dCI6Mjc1NzgsImVuZCI6Mjc2Njd9LCJlSGh2V1NXRmNGY2JpT2
wxIjp7InRleHQiOiJoYXQgbWFueSBvZiB0aGUgRXh0ZW5zaW4g
Z2VuZXMgaGFkIGxhcmdlIG51bWJlcnMgb2Ygbm92ZWwgc3BsaW
NlZCBpc29mb3Jtcy4iLCJzdGFydCI6MjgyMTgsImVuZCI6Mjgy
OTN9LCJxZ1FvS0pYWGRRYWFXZmpQIjp7InRleHQiOiJkZXJpdm
VkIGZyb20gdGhlIGludHJvbiBtb3RpZiIsInN0YXJ0IjozMDE1
MSwiZW5kIjozMDE4MH0sIlJBV3Z1cTRiakZtWmtBU0giOnsidG
V4dCI6Im9yIG5lZWRzIHRvIGJlISIsInN0YXJ0IjozMjgxMCwi
ZW5kIjozMjgyNX0sIk1DZW9lQ24yTEVwcjVKSUwiOnsidGV4dC
I6IkNUL0FDLiIsInN0YXJ0IjozNzM0OCwiZW5kIjozNzM1NH0s
IlQyQUlyRW5sdjF5R1BYY0oiOnsidGV4dCI6IiFbKipTYW5nZX
Igc2VxdWVuY2luZyBvZiBMUlgxIGFuZCBFWFQ5IGNETkEgaWRl
bnRpZmllcyBzcGxpY2VkIGZvcm1zKiogKiphKSoqIEfigKYiLC
JzdGFydCI6MzczNjYsImVuZCI6Mzc4NTZ9LCJHQ1dlWW9zUFJN
eEozZldZIjp7InRleHQiOiJPbmx5IHNwbGljZSBqdW5jdGlvbn
Mgd2l0aCBhdCBsZWFzdCAyMCBzdXBwb3J0aW5nIHJlYWRzIHRv
dGFsIGFjcm9zcyB0aGUgNiBzYW1w4oCmIiwic3RhcnQiOjM5Nj
Y5LCJlbmQiOjM5Nzc1fSwiakw0dW1BdWdiWUdzUUhwcSI6eyJ0
ZXh0IjoiZyBgbGltbWEtdm9vbWAgYW5kIGBsaW1tYS1kaWZmU3
BsaWNlYCIsInN0YXJ0Ijo0MDEyNiwiZW5kIjo0MDE2M30sIkZw
QjhnZUs0RFpNMlVPengiOnsidGV4dCI6ImFuIEZEUiB0aHJlc2
hvbGQgb2YgMC4yIiwic3RhcnQiOjQwMzM2LCJlbmQiOjQwMzU5
fSwiem5YVGdMVlJPTWpyOWp6diI6eyJ0ZXh0IjoicGxpY2VkIG
1hcHBpbmcgdG8gdGhlc2UgZ2VuZXMgaXMgYSBzeXN0ZW1hdGlj
IG1hcHBpbmcgZXJyb3IgdGhhdCBvY2N1cnMgYXQgYXBwcuKApi
IsInN0YXJ0Ijo0Mjc3NywiZW5kIjo0MjkxOH0sImY5aEZaaERo
WmF3djNWeUkiOnsidGV4dCI6IkRlc3BpdGUgdGhlc2UgbmVnYX
RpdmUgcmVzdWx0cywgd2Ugd2VyZSBhYmxlIHRvIGlkZW50aWZ5
IFBDUiBwcm9kdWN0cyBmcm9tIGNETkHigKYiLCJzdGFydCI6ND
g2NjksImVuZCI6NDg4MDd9fSwiY29tbWVudHMiOnsiRHF3Y1Z6
ZTZPYXNXOE1RVyI6eyJkaXNjdXNzaW9uSWQiOiJydkluUXlVNF
lRblg5V3g5Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEw
Njc3IiwidGV4dCI6ImVucmljaGVkIGF0ID8iLCJjcmVhdGVkIj
oxNTM2NjYzOTM3NDM2fSwiUnFZRnozaWVuNDk0UEVlNSI6eyJk
aXNjdXNzaW9uSWQiOiJzNTJrdFVrQnRWTk9SMTN2Iiwic3ViIj
oiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ildo
YXQgYXJlIFwiaGlnaGVyXCIgZXVrYXJ5b3Rlcz8gRG8geW91IG
1lYW4gcGxhbnRzIGFuZCBtZXRhem9hPyIsImNyZWF0ZWQiOjE1
MzY2NjQ4NTQxNzV9LCJ4ckF3Y1Q4UjdWRnlZa21pIjp7ImRpc2
N1c3Npb25JZCI6ImtZNW9UTWtqek1CZDhRQnoiLCJzdWIiOiJn
bzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiRGVmaW
5hdGUgXCJtb3JlIHN0cm9uZ2x5IGNhbm9uaWNhbFwiPyIsImNy
ZWF0ZWQiOjE1MzY2NjQ5NzM2NjF9LCI1czQyVlZPR29nTmowaH
lnIjp7ImRpc2N1c3Npb25JZCI6IkRxOXJTMExpcEpwSExVWDEi
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiSXMgdGhlIHJlZmVyZW5jZSBmb3IgdGhpcyB0aGUgb25l
IGluIHRoZSBuZXh0IHNlbnRlbmNlPyIsImNyZWF0ZWQiOjE1Mz
Y2NjUwNTkzNDZ9LCJJeGhkR25BYXJGSVkzTnA2Ijp7ImRpc2N1
c3Npb25JZCI6InlJWjVVQk1nYTJTSG40R1kiLCJzdWIiOiJnbz
oxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiV2UndmUg
YWxyZWFkeSB0YWxrZWQgYWJvdXQgYmVpbmcgbW9yZSBleHBsaW
NpdCBhYm91dCB3aG8gZGlkIHdoYXQuIFlvdSBuZWVkIG1vcmUg
ZGV0YWlsIHdoZXJlIHlvdSBkaWQgaXQsIGFuIHRvIHNheSB3aG
VuIHlvdSBkaWRuJ3QuIiwiY3JlYXRlZCI6MTUzNjY2NTExOTc3
OH0sIlRZMmxPWFgwUFZrdGNENlMiOnsiZGlzY3Vzc2lvbklkIj
oiMDVCbUdsaVhkZDJuaDVvNiIsInN1YiI6ImdvOjEwMjIwNTc5
NzI3Njk0MTAxMDY3NyIsInRleHQiOiJDb25jZW50cmF0aW9ucy
BuZWVkZWQgaGVyZSwgZXZlbiBpZiB5b3UgZGlkbid0IGRvIGl0
LiIsImNyZWF0ZWQiOjE1MzY2NjUxNjkxNzl9LCJwYmhPVjY3Qn
BkSWNBMTZtIjp7ImRpc2N1c3Npb25JZCI6Ik9rbW1oRDRoaktu
ZzVkZFgiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2Nz
ciLCJ0ZXh0IjoiRXhwZXJpbWV0YWwgZGVzaWduLCBudW1iZXIg
b2YgcmVwbGljYXRlcywgd2hhdCB3ZXJlIHRoZSBjb25kaXRpb2
5zIGV0Yy4iLCJjcmVhdGVkIjoxNTM2NjY1MjEwNzUzfSwiUlAx
cVJ3Z2FxTll2dklUYyI6eyJkaXNjdXNzaW9uSWQiOiJDT3FwZG
RiU3JleXFSUFBSIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQx
MDEwNjc3IiwidGV4dCI6IldoYXQgd2VyZSB5b3VyIGNyaXRlcm
lhPyAgRG8gd2UganVzdCBhY2NlcHQgdGhhdCB0aGUgc2FtcGxl
cyBwYXNzZWQ/IiwiY3JlYXRlZCI6MTUzNjY2NTI0OTY1Mn0sIl
o5QkRPd0lwVGNDTXMyM1QiOnsiZGlzY3Vzc2lvbklkIjoiUEFU
RVpnQW1LT1pTMXVzbCIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Nj
k0MTAxMDY3NyIsInRleHQiOiJXaGF0IHBhcmFtZXRlcnM/Iiwi
Y3JlYXRlZCI6MTUzNjY2NTI2MzU4MH0sInFyR3FtWFVmOHRRbj
FJS0UiOnsiZGlzY3Vzc2lvbklkIjoiM011dUVCb0h1TkZ2bEQ3
YyIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsIn
RleHQiOiJUYWJsZSB3aXRoIGluc2VydCBzaXplPyIsImNyZWF0
ZWQiOjE1MzY2NjUyODUwNTh9LCJ0aGpaUmxCa2ZGT1ZrSHY2Ij
p7ImRpc2N1c3Npb25JZCI6ImN4V2xnUDFpbWVBM3NpUTAiLCJz
dWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ij
oiUGFja2FnZSB2ZXJzaW9ucyBuZWVkZWQuIElkZWFsbHkgYWxz
byBzY3JpcHRzLiIsImNyZWF0ZWQiOjE1MzY2NjUzMTUyNzJ9LC
JPM05YZVdaWHVzSTJsTmNwIjp7ImRpc2N1c3Npb25JZCI6ImdX
dEFpdHMwUE9xZVlEazQiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNz
Y5NDEwMTA2NzciLCJ0ZXh0IjoiV2hhdCB3ZXJlIHRoZSBuZXcg
cGFyYW1ldGVycz8gUGVyaGFwcyBwcm92aWRlIHRoZSBjb21tYW
5kbGluZS4iLCJjcmVhdGVkIjoxNTM2NjY1MzU4NDM3fSwib0o5
elZkMmJodEl0bTVBVyI6eyJkaXNjdXNzaW9uSWQiOiJuRWQ0WW
1WR0pNUU8xNzhCIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQx
MDEwNjc3IiwidGV4dCI6InZlcnNpb24iLCJjcmVhdGVkIjoxNT
M2NjY1MzcyNTk4fSwiNVYyQ3diTUczUUJoVlVuZyI6eyJkaXNj
dXNzaW9uSWQiOiJJYWc4ZkRyM0NBUEo4Tm9VIiwic3ViIjoiZ2
86MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6InZlcnNp
b24iLCJjcmVhdGVkIjoxNTM2NjY1MzgwMjk0fSwiNU1RZ29PbE
hIQ1V4WDdyaSI6eyJkaXNjdXNzaW9uSWQiOiJsRFdpemFWRkxB
NU1LQlZ4Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNj
c3IiwidGV4dCI6ImV4cGxhaW4gXCJmbGF0dGVuZWQgZXhvbiBt
b2RlbHNcIiIsImNyZWF0ZWQiOjE1MzY2NjU0MDg4NzV9LCJhOU
xHSFQyRDdXQXc1QmtQIjp7ImRpc2N1c3Npb25JZCI6ImJrQlB0
U3FTbUZsRUJwYVQiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5ND
EwMTA2NzciLCJ0ZXh0IjoiV2hlcmUgYXJlIHlvdSBnZXR0ZXIg
eW91ciBnZW5lOm9udG9sb2d5IG1hcHBpbmcgZnJvbT8iLCJjcm
VhdGVkIjoxNTM2NjY1NDU5Njc1fSwiOGk0cVJsVWtlWmRpOVU5
NCI6eyJkaXNjdXNzaW9uSWQiOiJEUWJqU3dRNzdwSWF5bmdsIi
wic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4
dCI6InNjcmlwdHMgZm9yIGFsbCB0aGlzPyIsImNyZWF0ZWQiOj
E1MzY2NjU0OTU2NDl9LCJoWU81V3p6ak1FYWFBelhlIjp7ImRp
c2N1c3Npb25JZCI6IjFqRHJPV3NtZk1LQ2dCeFMiLCJzdWIiOi
JnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiVmVy
c2lvbiwgYW5ub3RhdGlvbiBzb3VyY2U/IiwiY3JlYXRlZCI6MT
UzNjY2NTUxMzc2M30sImFFejV2aGtDc0R1VEVQRHkiOnsiZGlz
Y3Vzc2lvbklkIjoiS2FJMnpqUGZZclBJNlhhMiIsInN1YiI6Im
dvOjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJzY3Jp
cHQiLCJjcmVhdGVkIjoxNTM2NjY1NTMwMDcwfSwienhDWnE2WF
A3TnozOVZ6YSI6eyJkaXNjdXNzaW9uSWQiOiJRTTU1QWlrU1ZD
dkFhejJMIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNj
c3IiwidGV4dCI6ImRlc2NyaWJlIG9yIHJlZmVyZW5jZS4iLCJj
cmVhdGVkIjoxNTM2NjY1NTQ3NzEwfSwiQVJpVmdXUzI3UjlsV3
hRUSI6eyJkaXNjdXNzaW9uSWQiOiJISnptMWRtSE11alR1anlB
Iiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidG
V4dCI6InBhcmFtZXRlcnMvb3B0aW9ucy4iLCJjcmVhdGVkIjox
NTM2NjY1NTg3MzE4fSwibmhrYzcwRUxScmZwQnNZYSI6eyJkaX
NjdXNzaW9uSWQiOiI0TVl5QzdnV0ljMFBkNnhYIiwic3ViIjoi
Z286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ikhvdz
8iLCJjcmVhdGVkIjoxNTM2NjY1NjA0MjI2fSwiQW1hczU1d0t3
UExBWWV6cSI6eyJkaXNjdXNzaW9uSWQiOiJzUjl2c3c1RXZNNT
JFWFVFIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3
IiwidGV4dCI6Ikhvdz8iLCJjcmVhdGVkIjoxNTM2NjY1NjEyND
kzfSwibTRFc0U2ZGZCU01tbVdwdSI6eyJkaXNjdXNzaW9uSWQi
OiJTUUxlZTlNTE5ST3NNYjFiIiwic3ViIjoiZ286MTAyMjA1Nz
k3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IldoZXJlIGRvIEkgZmlu
ZCB0aGlzPyIsImNyZWF0ZWQiOjE1MzY2NjU2NjI2NjZ9LCJTVk
JDd0hxVTBybk54RHRqIjp7ImRpc2N1c3Npb25JZCI6IlNXZ0c1
bkttN0xKSEFIbEoiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5ND
EwMTA2NzciLCJ0ZXh0IjoiV2hhdCBkbyB5b3UgbWVhbiBieSBc
InVuaXF1ZSBzcGxpY2UganVuY3Rpb24gcGFpcnNcIj8iLCJjcm
VhdGVkIjoxNTM2NjY1Njg1OTA1fSwiNnFRVkNOVGJMRzJKUGZh
RSI6eyJkaXNjdXNzaW9uSWQiOiJlcGpJSGMwaXpWYkpVZkZEIi
wic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4
dCI6IkNhbGN1bGF0ZWQgaG93PyIsImNyZWF0ZWQiOjE1MzY2Nj
U5MzQ1NzF9LCJJZHVnYlBTTGVnZVpuOUg0Ijp7ImRpc2N1c3Np
b25JZCI6IkJjakJGWWtGSjZnZHZDeDMiLCJzdWIiOiJnbzoxMD
IyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiU2NyaXB0cyBu
ZWVkZWQiLCJjcmVhdGVkIjoxNTM2NjY1OTUyMTg2fSwiNndWOW
FGdFM2V2Q4NXhObiI6eyJkaXNjdXNzaW9uSWQiOiJjYnhWRm9X
VWV0NVVSRlRIIiwic3ViIjoiZ286MTAyMjA1Nzk3Mjc2OTQxMD
EwNjc3IiwidGV4dCI6IkhvdyB3YXMgdGhpcyBkb25lPyBweXNh
bT8gd2hlcmUgaXMgdGhlIGNvZGU/IiwiY3JlYXRlZCI6MTUzNj
Y2NTk3NzE4Nn0sIlF3bTVrdFBSZWYyYjF4TkQiOnsiZGlzY3Vz
c2lvbklkIjoiUjhXUVZqWU9vemd0MEllYiIsInN1YiI6ImdvOj
EwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJ2ZXJzaW9u
IiwiY3JlYXRlZCI6MTUzNjY2NjA0NTQ3Nn0sIkNrVk83TDVNUE
NHdXlVcE8iOnsiZGlzY3Vzc2lvbklkIjoiT1NabmZZR0VDWEZh
bm5OVyIsInN1YiI6ImdvOjEwMjIwNTc5NzI3Njk0MTAxMDY3Ny
IsInRleHQiOiJzY3JpcHQvY29kZSIsImNyZWF0ZWQiOjE1MzY2
NjYwNTYzMDh9LCJDUGwyQkZ3Y29obHR1dFZaIjp7ImRpc2N1c3
Npb25JZCI6IjZ5d1oxSDFTNjVsTE1xYmYiLCJzdWIiOiJnbzox
MDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoic2V0dGluZ3
MvcGFyYW1ldGVycz8gT3Igd2Vic2l0ZSBhZGRyZXNzPyIsImNy
ZWF0ZWQiOjE1MzY2NjYwODk3MDZ9LCJiTXRWdTd3UWtKa1pHMl
NTIjp7ImRpc2N1c3Npb25JZCI6Ikt5d1U4c1NQOHJjOVFMalMi
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0Ijoid2h5IDIwPyIsImNyZWF0ZWQiOjE1MzY2NjYxMTM3MzB9
LCJGcTNJQUlqcEZlWjh5UVg2Ijp7ImRpc2N1c3Npb25JZCI6In
diOEhUa2V1YjV3c0xjUTYiLCJzdWIiOiJnbzoxMDIyMDU3OTcy
NzY5NDEwMTA2NzciLCJ0ZXh0IjoiSG93IGNvdW50ZWQ/IFZlcn
Npb25zL3NjcmlwdHM/IiwiY3JlYXRlZCI6MTUzNjY2NjE2Mjk4
M30sIjdDdWJtUWlXREVBRWRQaXciOnsiZGlzY3Vzc2lvbklkIj
oid3hiUXlNbW42T1E4MU0zViIsInN1YiI6ImdvOjEwMjIwNTc5
NzI3Njk0MTAxMDY3NyIsInRleHQiOiJZb3UgbmVlZCB0byBjb2
1tZW50IHNvbWV3aGVyZSBhYm91dCB0aGUgbGF4bmVzcyBvZiB0
aGlzIHRocmVzaG9sZC4gTWF5YmUgaW4gdGhlIHJlc3VsdHMgb3
IgdGhlIGRpc2N1c3Npb24uIiwiY3JlYXRlZCI6MTUzNjY2Njg2
OTY1M30sIkpMNEZJQjNoTFhWVlFDa1EiOnsiZGlzY3Vzc2lvbk
lkIjoieHl4bDFsMERGeVpXRUtENCIsInN1YiI6ImdvOjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJDb250aWd1b3VzIG
luIHRoZSBmbGF0dGVuZWQgbW9kZWw/IiwiY3JlYXRlZCI6MTUz
NjY2NjkwNjU4OH0sIkZYeGtlTTFtSVM1WkVLTUQiOnsiZGlzY3
Vzc2lvbklkIjoiejBKOWNiZHZ5NjhKNnBtQyIsInN1YiI6Imdv
OjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJpZGVudG
lmaWVkIGhvdz8gRnJvbSB0aGUgTkggdGFnPyBGcm9tIHRoZSBt
YXBwaW5nIHF1YWxpdHk/IiwiY3JlYXRlZCI6MTUzNjY2Njk1MT
IxMH0sInVDcXhETjFaWndYekNuRnciOnsiZGlzY3Vzc2lvbklk
IjoicTkzYlRRYXdiMlZTRFh0WiIsInN1YiI6ImdvOjEwMjIwNT
c5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJNYXliZSB5b3UgY291
bGQgaGF2ZSBhbiBhcHBlbmRpeCB3aXRoIFByb2dyYW0gYW5kIH
BhY2thZ2UgdmVyc2lvbi4gU2NyaXB0IG5hbWVzLCBsb2NhdGlv
bnMgYW5kIHdoZXJlIHRoZXkgd2VyZSB1c2VkLiIsImNyZWF0ZW
QiOjE1MzY2NjcwMTAzNDd9LCJwaGRCVmZ4VUJRRGJPbzh6Ijp7
ImRpc2N1c3Npb25JZCI6InBEaGhTVmc1NVh0aUdISVIiLCJzdW
IiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ijoi
V2hhdCBSTkFzZXEgZGF0YXNldD8gWW91IGhhdm4ndCBpbnRyb2
R1Y2VkIHRoaXM/IFdoYXQgd2FzIHRoZSBkZXNpZ24/IERvZXMg
aXQgcmVjYXB0aXVsYXRlIHRoZSBwYXR0ZXJucyB5b3UgZGlzY3
Vzc2VkIGluIHRoZSBwcmV2aW91cyBjaGFwdGVyPyIsImNyZWF0
ZWQiOjE1MzY2NjcxNzk0MjZ9LCJneFZwWkhCZlJXeXpXNnBYIj
p7ImRpc2N1c3Npb25JZCI6ImRGbklsUGlNS3RvV2NJb2QiLCJz
dWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0Ij
oiTm90IGluIG1ldGhvZHM/IiwiY3JlYXRlZCI6MTUzNjY2OTc1
MzM2M30sIjZHbDl5eUFHbXc4dDlzOTciOnsiZGlzY3Vzc2lvbk
lkIjoiZGEwUjc5TW5qemkxaEFOTiIsInN1YiI6ImdvOjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaGF0IGhhcHBlbn
MgdG8gdGhlaXIgZXhwcmVzc2lvbiBpbiB0aGUgZGF0YXNldHMg
ZnJvbSB0aGUgcHJldmlvdXMgY2hhcHRlcj8iLCJjcmVhdGVkIj
oxNTM2NjY5Nzg3MDc1fSwiZG8zVlZuZVhHclZFTlhYMyI6eyJk
aXNjdXNzaW9uSWQiOiJiWFNIbjJVdU14V0lSSmFpIiwic3ViIj
oiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Iklu
dGVyc3RpbmcgdGhhdCB0byBzb21lIGV4dGVudCB0aGUgZWZmZW
N0IGluIENIWCBpcyBldmVuIGhpZ2hlci4iLCJjcmVhdGVkIjox
NTM2NjY5ODIyNDE0fSwibUhSTzlMc3RpZTB2WnlSSiI6eyJkaX
NjdXNzaW9uSWQiOiJlSGh2V1NXRmNGY2JpT2wxIiwic3ViIjoi
Z286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IklzIH
RoaXMgaW4gdGhhdCBwYXBlciwgb3IgaXMgaXQgc29tZXRoaW5n
IHlvdSBleHRyYWN0ZWQgZnJvbSB0aGVpciByZXN1bHRzLiIsIm
NyZWF0ZWQiOjE1MzY2Njk4NTQ4MzZ9LCJqQVg0eW9EZHJUTlYw
NkcyIjp7ImRpc2N1c3Npb25JZCI6InFnUW9LSlhYZFFhYVdmal
AiLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0
ZXh0IjoiT3IgZnJvbSBzdHJhbmRlZCBzZXF1ZW5jaW5nLiIsIm
NyZWF0ZWQiOjE1MzY2Njk4NzgwNTl9LCJhYVVpeURiMU1TamdI
TDhEIjp7ImRpc2N1c3Npb25JZCI6IlJBV3Z1cTRiakZtWmtBU0
giLCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0
ZXh0IjoiPz8hISE/IiwiY3JlYXRlZCI6MTUzNjY2OTkwMDQ2MH
0sIkZOVEU2a0NNUWdKb0N1SHIiOnsiZGlzY3Vzc2lvbklkIjoi
TUNlb2VDbjJMRXByNUpJTCIsInN1YiI6ImdvOjEwMjIwNTc5Nz
I3Njk0MTAxMDY3NyIsInRleHQiOiJ3ZWlyZCIsImNyZWF0ZWQi
OjE1MzY2NzIyNjY1MzF9LCJwOUhXNFFyYWpOR3ZzTVQwIjp7Im
Rpc2N1c3Npb25JZCI6IlQyQUlyRW5sdjF5R1BYY0oiLCJzdWIi
OiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiTW
FyayBvbiBzZXF1ZW5jaW5nL2Nsb25pbmcgcHJpbWVycyIsImNy
ZWF0ZWQiOjE1MzY2NzI5NTk2Nzd9LCJNbXk1R3VyVUhWOFpzNk
14Ijp7ImRpc2N1c3Npb25JZCI6IkdDV2VZb3NQUk14SjNmV1ki
LCJzdWIiOiJnbzoxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiV2h5IDIwPyBXaGF0IGZyYWN0aW9uIG9mIGp1bmN0aW9u
cyB3ZXJlIGtlcHQ/IFdoYXQgZnJhY3Rpb24gb2YganVuY3Rpb2
5zIGluIHRoZSBFWFQgZ2VuZXMuIiwiY3JlYXRlZCI6MTUzNjY3
MzAxMTk4N30sIlJzb1U3eXpUM0JTVnB3UHciOnsiZGlzY3Vzc2
lvbklkIjoiakw0dW1BdWdiWUdzUUhwcSIsInN1YiI6ImdvOjEw
MjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaHkgdGhpcy
BhbmQgbm90IERFWFNlcT8iLCJjcmVhdGVkIjoxNTM2NjczMDMx
NDY3fSwiSXY4bE9lU3A3c3FLUWpORyI6eyJkaXNjdXNzaW9uSW
QiOiJGcEI4Z2VLNERaTTJVT3p4Iiwic3ViIjoiZ286MTAyMjA1
Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ik5lZWRzIGNvbW1lbm
QsIG90aGVyd2lzZSBpdCBsb29rcyBsaWtlIHlvdXIganVzdCB0
cnlpbmcgdG8gZ2V0IGF3YXkgd2l0aCBpdC4iLCJjcmVhdGVkIj
oxNTM2NjczMDUxOTE2fSwiSTNjZnRUN0o0bG5sbXV5RiI6eyJk
aXNjdXNzaW9uSWQiOiJ6blhUZ0xWUk9NanI5anp2Iiwic3ViIj
oiZ286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ildo
YXQgYWJvdXQgbGltaXRpbmcgdGhlIGFuYWxheXNpcyB0byBjYW
5vbmljYWwgc3BsaWNlIGp1bmN0aW9ucz8iLCJjcmVhdGVkIjox
NTM2NjczMDk3NDgxfSwid2NEd25NMTVSbVBxRXZZeCI6eyJkaX
NjdXNzaW9uSWQiOiJmOWhGWmhEaFphd3YzVnlJIiwic3ViIjoi
Z286MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IklzIH
RoZXJlIGEgY2hhbmdlIGluIHRoZSBkaXN0cmlidXRpb24gb2Yg
UENSIGJhbmRzICtOTU0/IiwiY3JlYXRlZCI6MTUzNjY3NDI2ND
Q4Mn19LCJoaXN0b3J5IjpbLTU5NzE0NTE5NiwtMTY1OTU0MDM3
NywtMTYxMDc3MzA1Myw3OTg4ODY4NCwtMjAxNjkxMTEwNywtMT
gyODM3ODA1OCw0NDcxMjc0MTksLTE1MDU3NjAwMDgsLTEyMDM4
NDQ5N119
-->
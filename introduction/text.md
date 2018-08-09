# Introduction

## The G-Quadruplex

\label{sec:intro_g4s}

### What is a Quadruplex?

\label{ssec:what_is_a_g4}

The genome is often referred to as the "blueprint" for life. This metaphor suggests a static set of instructions, which simply encodes data and does not change its form. In reality, however, chromosomes are highly dynamic structures which are able to undergo various covalent modifications to the DNA and proteins, as well as change topologies on a global and local scale [@Misteli2007]. At the global level, changes in how the chromatin is packed result in the closing off or opening up of specific regions, changing the level of transcription of the genes contained within them. Genes and regulatory sequences which are far apart in sequence space can be brought together through looping to allow co-regulation [@Feuerborn2015]. At the smaller scale, the DNA itself is able to fold into a number of different shapes, including various types of duplex, triplex, and quadruplex. Examples include B-DNA (the classic double helix), A-DNA (duplex), R-loops (three stre) [@Roberts1992; @Korzheva2000], i-motifs (quadruplex) [@Gehring1993; @Tua2000; @Zeraati2018] and G-Quadruplexes [@Gellert1962; @Sen1988; @Sundquist1989; @Hansel2016]. These structures have different relative stabilities depending on the conditions of their local environment, e.g. the local concentration of solutes, complementary RNAs, or stabilising proteins, the level of molecular crowding, or the intracellular pH [@Gehring1993; @Schultze1999; @Gaynutdinov2008; Rajendran2010; @Heddi2011; @Zhou2016]. Furthermore, some structures form only in DNA containing specific sequence attributes, such as high GC content [@Huppert2005]. Whilst non-B forms of DNA have been known to form *in vitro* for some time, it is only recently that evidence of their formation *in vivo*, and their effects on biological systems, has begun to accumulate.

One of the more well studied non-B DNA/RNA forms is the G-Quadruplex (G4), a guanine rich, four stranded DNA helix. The properties of guanine which allow G4s to form were first hinted at by German chemist Ivar Bang in 1910, more than four decades before Watson, Crick and Franklin deduced the structure of the double helix. Bang noted that guanosine nucleotides in concentrations of around 25 mg/ml will form a viscous gel [@Bang1910]. It was not until 1962, however, that Gellert et al. were able to use the technique of X-ray diffraction to identify the interactions which caused this property. They noted that guanine monomers were able to interact to form square, planar quartets, which then stacked to form a helical structure [@Gellert1962].

Many more recent studies have shown conclusively that far from being just an unusual property of monomeric guanine, square planar "G-quartets" can form from polymeric DNA, so long as it contains a high enough local concentration of guanines, in the correct arrangement. The guanines in a quartet interact through non-canonical Hoogsteen base pairing, with two hydrogen bonds occurring between each adjacent base (i.e. each base participates in four hydrogen bonds, for a total of eight bonds per quartet) (Fig \ref{g4_struct}a) [@Gellert1962; @Sen1988; @Sundquist1989]. These quartets can then stack on top of one another through hydrophobic interactions (Fig \ref{g4_struct}b). The number of stacked quartets in a G4 is usually referred to as the number of tetrads. The stability of G4s is generally proportional to the number of tetrads they contains. Mono or divalent cations, usually potassium, fit into the central channel of the G4. These sit equidistant from the guanines of each quartet, and between each pair of adjacent tetrads [@Rovnyak2000; @Wu2003]. The number of potassium cations is therefore one less than the number of tetrads (Fig \ref{g4_struct}b). Potassium strongly stabilises G4s, and so the stability of G4s is dependent upon potassium concentration [@Schultze1999; @Gaynutdinov2008].

The structure of G4s is highly polymorphic. Quadruplexes can form between multiple molecules (intermolecular) or from the same molecule (intramolecular). Intermolecular G4s are considered less common *in vivo* than *in vitro* however, as the effective concentration of G4 forming motifs *in vivo* is generally much lower than under experimental *in vitro* conditions [@Huppert2005]. On top of this, G4s can fold with DNA strands in the same orientation (parallel), in different orientations (anti-parallel) or in a mixed hybrid conformation (Fig \ref{g4_struct} c) [@Parkinson2002; @Phan2003; @Ambrus2006]. Which conformation is chosen can be dependent upon the sequence from which the G4 is formed, as well as environmental conditions. The loops which connect the G-rich "pillars" of the G4 can be connected laterally (resulting in anti-parallel conformation), diagonally (also anti-parallel) or through a "propeller" like fold (resulting in parallel strands). An anti-parallel G4 with only lateral loops is generally referred to as a "chair" like G4, whilst an anti-parallel G4 with a diagonal loop is referred to as "basket" like (Fig \ref{g4_struct} c). Many G4 forming sequences will fold into multiple conformations with different rates, resulting in two or more subpopulations of folded G4 from molecules of the same sequence. The human telomeric sequence has been found to adopt both parallel, anti-parallel and hybrid forms *in vitro*, suggesting a dynamic population [@Parkinson2002; @Phan2003; @Ambrus2006].

\newpage

![**Structure of a G-Quadruplex** **a)** The molecular structure of a G-quartet. Four Guanosines (only guanine base is shown, sugar-phosphate is represented as R) interact through Hoogsteen base pairing around a central monovalent cation. **b)** Cartoon showing the basic structure of a three tetrad G4. Three G-quartets are stacked through interactions between delocalised pi electrons. The structure is made up of four pillars of homopolymeric G-runs (shown in orange) joined by loop sequences. Potassium cations (shown in blue) sit between each tetrad. **c)** Cartoons showing how loop arrangement can contribute to the structural polymorphism of G4s. Loops can be lateral, diagonal or propeller like, resulting in anti-parallel, anti-parallel, and parallel G4s respectively. Anti-parallel G4s with all lateral loops are referred to as "chair"-like, whilst anti-parallel G4s with a diagonal loop are referred to as basket like. G4s can also contain a mix of parallel and anti-parallel strands. \label{g4_struct}](figures/g4_structure.svg)

\newpage

![**The Human Telomeric Sequence forms a G-Quadruplex** **a)** top and **b)** side view of a parallel G4 structure solved by X-ray crystallography (PDB: 1KF1, Parkinson et al. 2002). The sequence used is (TTAGGG)4, corresponding to four of the human telomeric repeat. Gs are coloured in yellow, As in green and Ts in orange. Potassium cations are coloured in purple.](figures/xray_structure.svg){height=750px}

\newpage

### G-Quadruplex Prediction from Sequence

\label{ssec:predict_g4s}

Genomic G4s form in sequences with high GC content; i.e. high ratio of G:C to A:T base pairs, and high GC skew; i.e. high ratio of G to C on one strand of the DNA. Because of the dependence of G4 structure on sequence, it is theoretically possible to computationally predict the genome-wide prevalence of G4s from sequence information alone, assuming that other conditions such as potassium concentration are held constant. The first attempt to characterise putative G4s (PG4s) at whole genome scale was conducted by Huppert and Balasubramanian in 2005. They formulated rules describing the general patterns that PG4 forming sequences tend to follow. Their first observation was that intermolecular G4s are unlikely to be common *in vivo* due to low strand concentration of the DNA. They also noted that the pillars of the G4 tended to be formed from contiguous guanine homopolymers, or G-runs: there have to be four such G-runs in close proximity to create a PG4. The length of the shortest G-run will determine the maximum number of stacked tetrads which can be formed. A minimum of three tetrads was suggested for prediction: whilst 2 tetrad G4s are possible, they are less stable. Finally, they suggested that to make folding of the G4 favourable, the length of the loop sequences connecting the G-runs should be relatively short. They suggested, using evidence from molecular modelling and CD spectroscopy, an upper limit of 7bp. Again, whilst loops of much longer length are possible [@Guedin2010], they were thought likely to be unstable. Their observations were combined to create the folding pattern $G_XN_{1-7}G_XN_{1-7}G_XN_{1-7}G_X$, where $X \geq 3$. This was named the Quadparser method [@Huppert2005], and can be applied to search genomes using simple regular expression machinery.

The Quadparser method has been successfully used to identify G4s in many organisms, however the it is not perfect and results in many false negatives as well as false positives. Various adjustments can and have be made to the pattern, including increasing loop length to a maximum of 12bp, allowing two tetrad PG4s, and allowing short bulges in G-runs [@Chambers2015]. Several tools have been released which attempt to include these sequences amongst matches [@Varizhuk2014; @Dhapola2016; @Hon2017]. These tend to increase the recall of the method but can also greatly increase the number of false positives. Other methods have been proposed, such as G4Hunter [@Bedrat2016], which allow PG4s to be given a numeric score based on the GC content and skew of the sequence. G4Hunter is generally performed using a sliding window between 20 and 40bp in length, and is evaluated for each window by the following method:

```python
score = 0
for base, run_length in run_length_encode(sequence):
    if base == 'G':
        score += min(run_length, 4) ** 2
    elif base == 'C':
        score -= min(run_length, 4) ** 2
    else:
        pass
score /= len(sequence)
```

Sequences which have high PG4 forming ability on the positive strand will therefore be given strong positive scores, whilst sequences with PG4 forming ability on the negative strand will be given strong negative scores. A threshold value is chosen below which to filter out non-PG4 forming sequences. Bedrat et al argued that this method was an improvement over the Quadparser technique because it was more flexible, however it also results in a much greater number of false positives when applied to a whole genome, since there are no constraints on how the G-runs are arranged in the windowed sequence. These methods will be discussed further in \autoref{chap:g4seeqer}.

\newpage

### Methods used in the characterisation of G-Quadruplexes

\label{ssec:biophys_char_g4s}

Gellert et al. first characterised G4 structure using X-ray diffraction of fibres from dehydrated guanine gels [@Gellert1962]. Since then, biophysical techniques have become key in the study of G4 structures *in vitro*. Since the advent of chemical DNA oligonucleotide synthesis in the 1980s, it has become relatively cheap to order high purity single stranded oligonucleotides for PG4 sequences, and produce micromolar concentration solutions which can be probed by EMSA, DMS footprinting, FRET, CD spectroscopy or NMR.

Circular dichroism (CD) spectroscopy utilises the difference in absorbance of circularly polarised light by molecules with chiral structures (i.e. with non-superimposable mirror images). Parallel and anti-parallel G4s both exhibit unique CD absorbance spectra which are distinct from the spectra of disordered single stranded DNA [@DelVillar-Guerra2017]. Solutions containing multiple subpopulations of different parallel and antiparallel G4s will produce spectra which are more complex to interpret, but are still clearly distinct from unordered DNA. Melting temperatures of G4s can be determined using CD or UV spectroscopy temperature gradients measured at 295nm [@DelVillar-Guerra2017].

Nuclear magnetic resonance spectroscopy (NMR), specifically Proton-exchange spectroscopy (1H-NMR), can also be used to identify G4 DNA [@daSilva2007]. NMR is conducted in deuterated water, since deuterium has a spin of 1 and therefore does not contributed to the NMR spectrum. In single or double stranded DNA, imino protons in the guanine nucleotides will be exchanged with the solvent on short timescales, resulting in loss of the imino proton signal as they are replaced with deuterated protons. In G4 DNA, on the other hand, imino protons are located centrally within the G4 structure, and are therefore protected from exchange. This means that the 1H spectra can be used to distinguish G4s from unordered ssDNA. It does not, however, identify whether the G4 has parallel or anti-parallel topology. Further characterisation can be conducted using Nuclear Overhauser effect spectroscopy (NOESY) to identify spatial relationships between protons in the G4 [@daSilva2007].

Finally, since folded G4s with short loops and flanking sequence are relatively globular, a number of G4 structures have been crystallised from oligomers. The structures of these crystals have then been solved by X-ray crystallography [@Parkinson2002; @Nicoludis2012].

Studies of the structure of G4s has led to the development of a variety of G4-binding ligands, which might have potential as anti-cancer agents [@Mergny1998; @Neidle2010; @Balasubramanian2011]. These have a wide range of structures and bind to the various G4 topologies with different strengths [@Monchaud2008]. The three major classes of G4 binding ligands are: 
* **external stacking ligands**, which have delocalised pi electron systems, and stack on the the hydrophobic surfaces of the outer G4 tetrads;
* **intercalating ligands** , which bind in the groove between stacked tetrads;
* **external groove binding ligands**, which insert into the groove of the G4 helix, or interact with loop sequences [@Monchaud2008]. 

Small molecules which bind G4s include:
- **N-methyl-mesoporphyrin (NMM)**, an external stacking porphyrin ligand [@Nicoludis2012]; TMPyP4, another porphyrin derived ligand with multiple binding modes including external stacking and intercalation [@Haq1999; @Phan2005] (Fig \ref{drugs}a).  NMM is known to prefer parallel G4s [@Nicoludis2012];
- **TMPyP4**, a porphyrin derived ligand with multiple binding modes including external stacking and intercalation [@Haq1999; @Phan2005]. Externally stacked TMPyP4 has been reported to prefer parallel G4s over anti-parallel [@Arora2008; @Nagesh2010];
-  **Pyridostatin**, an external stacking molecule specifically designed for G4 binding [@Rodriguez2008] (Fig \ref{drugs}b).  Pyridostatin is thought to bind to all G4s equally well [@Muller2010];
-  **Berberine**, a naturally occurring alkaloid and G4 external stacking agent (Fig \ref{drugs}c) [@Wu1999; @Franceschin2006]. Berberine has been reported to bind to both parallel and anti-parallel G4s [@Bazzicalupi2013; @Li2017].

These molecules have been used to study G4s *in vitro* and also *in vivo*. G4 binding ligands have been shown to have different binding affinities for G4s derived from different sequences, hinting at the potential of drugs targeting G4s in specific gene promoters [@Weldon2018].

\newpage

![**G-Quadruplex Stabilising Ligands** Structures and mode of action of Pyridostatin (a), NMM (b), and Berberine(c). Crystal structures of human telomeric DNA (parallel G4 form) in complex with NMM and Berberine are from PDB entries 4FXM and 3R6R respectively (Nicoludis et al. 2012, Bazzicalupi et al. 2013). Potassium ions and solvent molecules have been hidden for visualisation purposes. \label{drugs}](figures/drugs.svg)

\newpage

Whilst biophysical methods have led to a wealth of data on G4 formation *in vitro*, biological evidence of G4 formation has been a much later development. One common approach to identify cellular G4s, and study their biological implications, is to treat cells or organisms with G4 binding ligands such as NMM or Pyridostatin [@Rodriguez2012; @Hershman2008; @Nakagawa2012]. This has been shown to have various effects on replication, genome stability, transcription and development [REFs?]. Naturally fluorescent or fluorophore labelled small molecules are commonly used to localise G4s by microscopy, or study their *in vitro* folding by single molecule Förster resonance energy transfer (smFRET) [@Maleki2017, @Hou2017]. Biotinylated pyridostatin has also been used to pull down G4 DNA structures from human DNA [@Muller2010].

Moving beyond small molecules, Biffi et al developed an antibody using synthetic phage display technology that specifically recognises G4 DNA [@Biffi2013]. They used the antibody, named BG4, in immunofluorescence experiments to visualise G4s in human chromatin [@Biffi2013]. More recently, it was used in Chromatin Immunoprecipitation sequencing (ChIP-seq) experiments to identify the specific regions of human chromatin where folded G4s occur [@Hansel2016]. The biological findings of these experiments will be discussed in \autoref{sec:biological_roles}.

A number of techniques have been employed for whole genome or transcriptome mapping of G4s. In a method they named G4-seq, Chambers et al. introduced potassium or pyridostatin into the buffer of Illumina sequencing-by-synthesis reactions [@Chambers2015]. This resulted in G4 formation in the single stranded DNA fragments, which caused stalling of DNA polymerase, resulting in sequencing errors. They conducted experiments using this method on single stranded DNA derived from the human genome. Clusters of identical reads generated by bridge amplification in an Illumina flow cell were first sequenced normally [@Metzker2010]. The products of the first sequencing-by-synthesis reaction were then washed off, and the clusters were resequenced in the presence of the G4 binding agents. This yielded an initial high quality, mappable read, and an error prone read. When they mapped the initial reads to the genome, the number of sequencing errors in each position was considered an indicator of the G4 forming potential. The authors showed using this technique that only 30% of human genomic sequences found to form G4s *in vitro* conform to the Quadparser motif [@Chambers2015]. This does not, however, address the question of whether these motifs form G4s *in vivo*.

Yoshida et al. also used a similar method to identify G4 clusters in human genomic DNA, by PCR amplifying sequences in the presence and absence of the G4 binding ligand telomestatin [@Yoshida2018]. They showed that regions which contained G4s were amplified with lower efficiency in the presence of telomestatin, due to polymerase stalling events. This resulted in a quantifiable reduction in the number of reads mapping to G4 containing regions of the genome in the drug positive samples, relative to the drug negative samples. The authors identified clusters of motifs predicted to form multiple G4s were more effective at stalling polymerase than single G4s [@Yoshida2018].

G4s which form in RNA *in vitro* can also be mapped globally, using a technique called rG4-seq, developed by Kwok et al. This method again utilised stalling at G4s, in the presence of potassium or pyridostatin, this time by the RNA templated DNA polymerase Reverse transcriptase (RT) [@Kwok2016]. They identified positions in mRNAs where a reproducible drop in reads occurred in samples where RT mediated DNA synthesis was conducted in the G4 stabilising conditions, relative to unstabilised controls. Guo and Bartel elaborated on this technique to map RNA G4s *in vivo*, by treating cells with dimethyl sulphate (DMS) prior to *in vitro* RT stalling. DMS can penetrate into living cells and methylate the N7 position of guanines that are not involved Hoogsteen base paired structures, i.e. G4s [@Wells2000]. This methylation prevents the refolding of mRNA into G4s *in vitro*, thereby abolishing RT stalling compared to untreated mRNAs [@Guo2016]. Guo and Bartel identified that the majority of G4s in yeast and human mRNAs are maintained in an unfolded state *in vivo* [@Guo2016]. The biological implications of this are discussed in \autoref{ssec:mrna_stability}.

### G-Quadruplex stability prediction using Machine Learning

\label{ssec:g4_machine_learning}

A growing number of G4 forming oligonucleotide sequences have now been characterised *in vitro* by various methods, particularly by CD spectroscopy and UV melting, and the melting temperatures have been published online. The number of these is now great enough that several groups have utilised them to train machine learning models to predict G4 forming potential. Machine learning is the application of statistical techniques to identify patterns in large datasets, and produce rules which will generalise to new data [@Larranaga2006]. The use of such statistical models in bioinformatic tools has grown greatly since the advent of techniques such as high throughput sequencing, as they produce the large amounts of data required for training such methods [@Larranaga2006]. The use of machine learning in nucleic acid motif prediction will be discussed in more detail in \autoref{chap:g4seeqer}.

Stegle et al. implemented a Gaussian Process model incorporating extracted features from Quadparser conforming G4s [@Stegle2009]. These features were the number of tetrads, the length of each of the three loops, the total loop length, and the frequencies of adenine, cytosine and guanine in the sequence, as well as the raw sequence itself. A second kernel which incorporated features about the conditions the melting temperature was acquired under, i.e. the concentrations of potassium, sodium, ammonium, and magnesium ions, was also used in the model. They trained this model on a set of 260 DNA G4 melting temperatures which were acquired from a literature search. In a cross validation experiment using 100 random 50% hold out splits, the authors were able to achieve a good level of test set accuracy with an average of 80% of predictions within 5 degrees of the true melting temperature [@Stegle2009]. Furthermore, their model was interpretable, and they were able to identify tetrad number and the length of the central loop as the most important sequence features in PG4 stability. The authors employed active learning to identify candidates from human promoter sequences with high uncertainty in the model, and used CD spectroscopy to characterise them.

Whilst Stegle et al.'s model was successful on Quadparser conforming motifs, it is estimated that ~70% of G4s in the human genome do not conform to this motif [@Chambers2015]. More recently, Garant et al. published a method for predicting RNA G4s which was trained on a set of 368 experimentally determined sequences, 149 of which were G4 positive and 179 of which were G4 negative [@Garant2017]. From these sequences the trinucleotide contents were extracted, and used to train a densely connected multi-layer perceptron model, with a single hidden layer containing 35 nodes. This trinucleotide trained model had the advantage being more flexible for G4s that do not conform to the Quadparser motif. Their model achieved an average AUC score of 0.92 on hold out sets in a 5 fold cross validation experiment. When tested on the rG4-seq dataset of RT stalled RNA G4s [@Kwok2016], the method did not perform as well as G4Hunter [@Bedrat2016]. This indicates that some positional information is lost when sequences are converted to trinucleotide content features.

Finally, the more recent efforts to produce high-throughput methods for identifying genomic G4s, such as G4-seq developed by Chambers et al., have created much better in depth datasets for training machine learning models. Sahakyan et al. used the G4-seq dataset in their model. This was a extreme gradient boosted machine model developed using `xgboost` [@Chen2016], which regressed the percentage mismatch score of sequences from the G4-seq dataset which conform to the Quadparser method [@Sahakyan2017]. The authors extracted features from Quadparser conforming PG4s similar to those employed by Stegle et al., including tetrad number, loop length, and mono-, di- and triunucleotide contents of the PG4, and flanking regions. This model was very successful at identifying Quadparser conforming motifs which did or did not actually form G4s, achieving a root mean squared error score of 8.14 (units used were mismatch score in G4-seq experiment, in percentage format) [@Sahakyan2017]. The method could not identify non-Quadparser conforming G4 motifs, however, which make up a large proportion of experimentally characterised G4s [@Chambers2015].

\newpage

### Other G-rich Nucleic Acid Structures

G4s are not the only DNA or RNA structures which occur specifically in sequences with high GC content and skew. Another structure is the R loop, which can form when transcription of a C-rich template strand to a G-rich RNA molecule occurs [@Reaban1994; @Li2005; @Ginno2012; @Ginno2013]. This RNA molecule is complementary to the template strand, and can therefore form a DNA:RNA hybrid duplex, leaving the coding strand of the DNA in a single stranded form [@Skourti-Stathaki2014]. Once formed, these hybrids are more thermodynamically stable than normal DNA:DNA duplexes [@Roberts1992]. This could be partially explained to the formation of G4 structures in the G-rich single stranded DNA of the coding strand [@Duquette2004]. The biological implications of R loop formation and their interactions with G4 structures are discussed briefly in \autoref{ssec:rloop_csr}.

\newpage

## Biological Roles of G-Quadruplexes

\label{sec:biological_roles}

### Genome Stability & DNA Replication

\label{ssec:genome_stability}

The distribution of G4s in genomes has been predicted from sequence for a wide variety of organisms [@Huppert2005; @Hershman2008; @Du2008; @Mullen2010; @Garg2016], and has been experimentally determined by the techniques mentioned above for the human genome. There is conclusive evidence that G4s are not uniformly distributed throughout genomes, but tend to be clustered at functional locations. By far the strongest enrichment of G4s is seen at telomeres. G4s are also found more than would be expected in origins of replication, gene promoters, and inside gene bodies, particularly the 5' and 3' untranslated regions (UTRs). In these locations, it has been demonstrated that G4 formation has effects on the processes of DNA replication, genome stability, transcription, and translation.

Perhaps the best characterised biological role for is at the telomeres. Telomeres are the protein-DNA structures found at the ends of linear eukaryotic chromosomes. They consist of thousands of tandem repeats of a G-rich sequence [@Moyzis1988], with a single stranded overhang of around 100-200bp [@Makarov1997]. Due to functional limitations in templated DNA replication, the very ends of these cannot be duplicated during cell division. This means that without intervention, the chromosome will gradually shorten with each division. Telomeres therefore serve as protective caps that prevent the loss of important coding DNA from the genome. In humans, the telomeric repeat is (TTAGGG)n [@Moyzis1988]. This has been identified through various methods as a G4 forming sequence.

Telomeres are coated in architectural proteins, called telomere end binding proteins (TEBPs), which protect the DNA from recognition by DNA damage response pathways. Giraldo & Rhodes showed that a yeast TEBP, RAP1, induces formation of G4 structures in telomeres *in vitro* [@Giraldo1994]. More recently, immunofluorescence experiments using the BG4 antibody showed that G4 foci overlap with fluorescent foci produced by Fluorescent *In Situ* Hybridisation (FISH) of the telomeric repeat in human HEK 293T cells [@Moye2015]. These results suggest that telomeric sequences do form G4s *in vivo*, which may be bound and stabilised by TEBPs. Furthermore, Moye et al. identified that telomerase, the template-independent DNA polymerase which synthesises new telomeric repeats, is able to bind to and partially unwind parallel G4s, but not anti-parallel G4s [@Moye2015]. There is therefore the suggestion of a G4 regulated mechanism for telomere maintenance.

G4s seem to also play roles in other aspects of DNA replication. It is well documented that G4s are capable of stalling polymerases *in vitro*. There is also growing evidence that without the assistance of G4 unwinding helicases, G4s might cause polymerase stalling *in vivo*. Work by Rodriguez et al. demonstrated that treatment of human cancer cells with G4 stabiliser pyridostatin causes an increase in the DNA damage marker γH2AX, suggesting an increase in double strand breaks (DSB) [@Rodriguez2012]. This damage was ameliorated by treatment with a DNA replication inhibitor, suggesting the damage was caused during replication. This is likely to be the result of replication fork collapse at G4 DNA blockages. The DNA helicase FANCJ, which is a tumour-suppressor gene often mutated in breast and ovarian cancers, is involved in DSB repair and has been shown to preferentially bind and unwind G4s [@Wu2016]. Mutation of the FANCJ ortholog DOG1 in *C. elegans* results in genome instability, and accumulation of deletions upstream of G4s [@Kruisselbrink2008].

The Bloom gene (BLM) encodes a DNA helicase which, when mutated in humans, causes Bloom's Disease [@German1993]. This is characterised by a reduction in genome stability resulting from a large increase in sister chromatid exchange events (SCE). SCEs are caused by double strand breaks, which are repaired via homologous recombination [@Wu2007]. BLM has been shown to prevent SCEs by unwinding Holliday Junctions (a common form of branched duplex), and restarting collapsed replication forks [@Davies2007]. Work by Sun et al. first showed that BLM is also capable of binding and resolving G4s *in vitro* [@Sun1998]. Recently, van Wietmarschen et al. used a single cell sequencing technique to map the locations of SCEs in cells lacking function BLM, and found that many SCEs were at sites predicted to form G4s [@VanWietmarschen2016]. This suggests that direct action of helicases is required to prevent genome instability at G4 loci.

In humans, DNA replication occurs from tens to hundreds of thousands of origins of replication which are found at regular distances of around 10-100kb apart [@Huberman1968; @Besnard2012]. Using genome-wide mapping of replication origins by short nascent strand sequencing, Besnard et al. identified that the majority of human origins were in GC rich regions of DNA, and that 67% overlapped with motifs conforming to the Quadparser pattern [@Besnard2012]. 91% of origins were associated with G4s with extended loop lengths of up to 15bp. Furthermore, they found an association between the number of G4 motifs, and the strength of usage of origins, suggesting that G4s might be a recruiting factor for replication machinery.

\newpage

![**Role of G4s in DNA replication** Illustration showing the possible mechanisms of G4 involvement in DNA replication. Besnard et al. identified that a majority of human replication origins contain a PG4, suggesting that G4s may recruit replication machinery (Besnard et al., 2012). In the absence of G4 unwinding helicase FANCJ, Kruisselbrink et al. identified deletions downstream of PG4 loci, suggesting G4s cause replication fork collapse. \label{replic}](figures/replication.svg)

\newpage

### Transcription

\label{ssec:transcription}

Transcription is the process by which DNA is copied into messenger RNA (mRNA) or non-coding RNA (ncRNA), by a DNA-templated RNA polymerase. In eukaryotic systems, all mRNA is transcribed by RNA Polymerase II (Pol II). Initiation of transcription is often catalysed by general or specific transcription factors which bind to the promoter region, upstream of the transcriptional start site (TSS). Human promoter sequences are enriched for PG4 motifs conforming to the Quadparser motif [@Eddy2006; @Huppert2007]. Tumour suppressor genes have fewer promoter PG4s than might be expected by chance, whilst proto-oncogenes contain more than might be expected [@Eddy2006]. PG4s also overlap with regions of open chromatin, detected by methods such as DNase Hypersensitivity sequencing (DNase-seq) or Assay for Transposable-Accessible Chromatin by Sequencing (ATAC-seq), suggesting they are often actively transcribed [@Huppert2007]. Hänsel-Hertsch et al. used the BG4 antibody to perform ChIP-seq of G4 structures in conjunction with ATAC-seq and RNA-seq, in normal human keratinocytes and an immortalised cell line [@Hansel2016]. They found that BG4 peaks were indeed associated with open promoters, and were found upstream of expressed genes [@Hansel2016]. Interestingly, many more BG4 peaks were identified in the immortalised cells than in normal cells, despite having similar levels of open chromatin. Furthermore, genes which only had a promoter BG4 peak in the immortalised cells tended to be more highly expressed in those cells compared to the normal cells, suggesting that promoter G4s may increase gene expression. A potential mechanism for this increase in expression might be the result of recruitment of positive transcription factors.

Perhaps the most well-studied promoter G4 is the Nuclease Hypersensitive Element III (NHEIII) found in the promoter of the c-MYC oncogene. The NHEIII contains a number of G-rich tracts which have been shown to form G4s *in vitro* by a variety of methods [@Simonsson1998; @Siddiqui-Jain2002; @Seenisamy2004; @Ambrus2005]. Formation of a G4 by this region has a strongly repressive effect on gene expression. Siddiqui-Jain et al. found that treatment of cells with the G4 stabilising agent TMPyP4 led to repression of c-MYC, whilst G4 abolishing mutations in a c-MYC promoter luciferase assay caused a three-fold increase in expression [@Siddiqui-Jain2002]. Pull-down of NHEIII binding factors by González et al. identified Nucleolin as a potential G4 interacting partner [@Gonzalez2009]. Nucleolin is a multi-functional protein implicated in ribosome synthesis, transcription and chromatin remodelling. González et al. went on to show that Nucleolin binds to the c-MYC promoter *in vivo*, and that nucleolin overexpression results in downregulation of c-MYC gene expression [@Gonzalez2009], by promoting G4 formation [@Gonzalez2010].

G4s which form within the gene body may have differing effects on transcription, depending on the strand they occur in. Analysis of human gene expression data has suggested that genes containing coding strand G4s downstream of transcriptional start sites tend to have higher expression at the mRNA level than those that do not [@Du2008], even when other factors such as gene function are controlled for. It has been speculated that G4 formation competes with double stranded DNA, creating single stranded "bubbles" which promote Pol II binding and transcription [@Rhodes2015]. Since the coding strand is not directly used by Pol II, coding strand G4s will not cause polymerase stalling. G4s which form in the template strand, however, may form blockages which could slow or pause the progression of Pol II.

Transcription progresses by using the template strand as an antisense copy to replicate the coding strand sequence in mRNA. G4s which form in the template strand of gene body DNA will therefore need to be resolved before Pol II can move through them. Due to the relative stabilities of dsDNA and G4s, G4 formation may only occur in the gene body after a pioneering round of transcription, during which the DNA is in single stranded form [@Eddy2011]. Rodriguez et al. showed that some of the DNA damage caused by treating cells with pyridostatin was transcription-dependent, and could be ameliorated by treating cells additionally with an inhibitor of transcription [@Rodriguez2012].

\newpage

![**Possible Mechanisms for DNA and RNA G4 function in transcription and gene expression** **a)** Possible mechanism for function of G4s located in the coding strand. Since G4s form from single stranded DNA, G4s in the coding strand may promote melting of double stranded DNA, increasing transcription levels. G4s which form in the coding strand of the exonic DNA of a gene will also be present in the mRNA produced from that locus, and G4s which form in the coding strand of introns will be present in pre-mRNA. These RNA G4s might also influence gene expression through alteration of pre-mRNA splicing, mRNA stability, or translation. **b)** Possible mechanism for the function of G4s located in the template strand. The template strand of genes is scanned by RNA Polymerase II during transcription, and G4s which form in ahead of the transcription complex may cause slowing or stalling if they cannot be correctly resolved. RNA Polymerase II translocation speed is linked to a number of co-transcriptional processes, including splicing. Adapted from Figure 4. Rhodes and Lipps, 2015. \label{mech}](figures/transcription.svg)

\newpage

Methods for estimating Pol II elongation rates across genes, such as GRO-seq or BruDRB-seq, have associated changes in speed with various features such as specific histone modifications, exon density, and sequence features. Veloso et al. correlated elongation rates from BruDRB-seq data with GC content of genes, and found that genes with higher GC content tended to have slower elongation [@Veloso2014]. They hypothesised that this could be due to the greater stability of GC-rich duplexes, which have extra hydrogen bonds. It is also possible, however, that this effect could be partially due to greater numbers of G4s in GC-rich genes.

In human cells, profiling of Pol II occupancy by ChIP-seq has demonstrated that there is a large peak of paused polymerase in the first 30-60bp downstream of the TSS [@Jonkers2015]. This pausing is an important checkpoint ensuring Pol II is correctly modified before elongation begins. Genes which require large and rapid increases in expression in response to environmental stresses, such as heat shock proteins, also have large amounts of paused Pol II which can be activated quickly. During initiation of transcription, Pol II is recruited to the TSS by specific or general transcription factors, and transcribes for a short distance before becoming paused. Formation of paused Pol II, referred to as the Pre-Initiation Complex (PIC), is stabilised by the Negative elongation factor (NELF) and DRB-sensitivity-inducing factor (DSIF), as well as by phosphorylation of the carboxy-terminal domain (CTD) of Pol II at Serine 5. Productive elongation can then be restarted by the action of the positive transcription elongation factor-b (P-TEFb) complex, which phosphorylates NELF and DSIF, causing the former to be released from the PIC, and the latter to switch to becoming a positive elongation factor. P-TEFb also phosphorylates the CTD at Serine 2, which is considered a hallmark modification of active Pol II. How Pol II pausing is regulated is still not clear [@Jonkers2015; @Liu2015]. One hypothesis is that the sequence content of promoters and promoter proximal regions may be important for regulating pausing. A number of promoter motifs, such as the GAGA motif or the downstream promoter element (DPE) have been associated with promoters with high levels of stalling [@Hendrix2008].

It is well established that the promoter proximal regions of genes in many organisms have, on average, higher GC content than the rest of the gene body [@Veloso2014]. Eddy et al. identified that the first 200bp downstream of the TSS tends to be more GC-rich in genes with high levels of proximal pausing than in genes which do not exhibit pausing [@Eddy2011]. The G4-forming potential of these genes also tended to be greater on the coding strand, meaning that G4 structures might also form in the nascent mRNA. Eddy et al. hypothesised that these 5' mRNA G4s might signal back to the polymerase to produce pausing. A similar mechanism involving an RNA hairpin has been implicated in pausing of *E coli* RNA Polymerase [@Toulokhonov2007].

The enrichment of G4s in promoters and promoter proximal regions suggests that proteins involved in transcriptional complexes may bind specifically to these structures. The general transcription initiation factor complex, TFIIH, contains 11 subunits, and is required for both transcriptional processes and DNA repair through the Nucleotide Excision Repair (NER) pathway [@Compe2016]. The DNA helicases XPB and XPD are essential components of TFIIH, which catalyse the denaturation of DNA in promoters or around lesions [@Coin2007]. Through *in vitro* binding assays, Gray et al. identified XPB and XPD as G4 interacting proteins, which bind G4s preferentially over dsDNA [@Gray2014]. XPD was also found to unwind G4s *in vitro*. Gray et al. went on to perform ChIP-seq of XPD and XPB, and showed that it was enriched at TSS loci containing G4s. Approximately 40% of XPD/XPB peaks contained PG4s conforming to the Quadparser pattern (with loop lengths >= 12bp). Furthermore, Hänsel-Hertsch et al. also reported a strong overlap between these XPD/XPB peaks and G4 loci observed by BG4 ChIP-seq [@Hansel2016]. This suggests that TFIIH may be recruited to G4 containing promoters to initiate transcription.

\newpage

### mRNA Processing

\label{ssec:mrna_processing}

Nascent pre-mRNA which is newly transcribed by Pol II must undergo 5' capping, splicing, RNA modification, poly-adenylation and quality control before it can mature into mRNA which is exported to the cytoplasm. Many of these processes occur co-transcriptionally and are tightly co-ordinated to prevent mistakes. Multiple independent studies in different organisms estimate that between 75%-85% of splicing is conducted in a cotranscriptional manner [@CarrilloOesterreich2010; @Ameur2011; @Khodor2011; @Girard2012; @Tilgner2012; @Windhager2012]. Oesterreich et al. found that in yeast, 10% of intron splicing is complete when Pol II is only 26bp downstream of the intron acceptor site, and 50% complete when Pol II is 45bp downstream. [@CarrilloOesterreich2016]. By modifying Pol II to increase its speed 2.3x, they also showed that splicing could become rate limiting when elongation rate is greater. Furthermore, modifications to splice site sequences in a reporter reduced the rate of splicing, presumably by reducing the strength of recognition by snRNAs, and thereby the rate of spliceosome assembly [@CarrilloOesterreich2016]. This indicates that interplay of splice site strength and Pol II elongation speed determine the relative efficiency of splicing.

It has been estimated that greater than 90% of human genes undergo some form of alternative splicing [@Wang2008]. During this process, different donor and acceptor sites compete to be utilised. The most common form of alternative splicing in humans is exon skipping, where constitutive donor and acceptor sites are paired such that an intervening exon is removed from the mature mRNA [@Kim2007]. Other forms include alternate donor or acceptor usage, where the other site used is constitutive, or intron retention, where splice sites are simply not used at all [@Wang2015]. Regulation of alternative splicing can occur via protein splicing factors, as well as through changes in Pol II elongation speed [@DelaMata2010; @Jonkers2015]. Through this mechanism, alternative splicing may be linked to changes in Pol II speed as it transcribes through template stranded G4s or other DNA structures.

Aside from their effects on Pol II speed, G-rich sequences with the potential to form G4s have been implicated as important intron motifs for splicing. These PG4s are predicted in the coding strand, meaning that they could form G4s in either the DNA or the nascent pre-mRNA. Analysis of the first exon-intron boundary of human genes by Eddy & Maizels revealed that about 50% of boundaries contain PG4s within the first 100bp of intronic sequence [@Eddy2008]. They noted that a number of hnRNP family proteins, such as hnRNP A1 and hnRNP H, bind to G-rich motifs in RNA. hnRNP A1 has been called the "swiss army knife of gene expression", due to its ability to bind both chromatin and mRNA, and its putative roles in transcription, mRNA splicing, telomere maintenance, mRNA export, and translation. Interestingly, hnRNP A1 is capable of binding and unwinding DNA G4s, and has been demonstrated to increase expression of the KRAS and c-MYC oncogenes by resolving repressive G4s in their promoters [@Chu2016]. Furthermore, there is evidence that hnRNPs F and H are capable of binding to G-rich RNA sequences to regulate splicing [@Xiao2009]. Xiao et al. identified that G-rich sequences downstream of donor splice sites with intermediate levels of homology to the snRNA U1 were strongly conserved. They also noted that the expression of genes with these intermediate splice sites and G-runs was sensitive to the knock-down of hnRNP H [@Xiao2009].

A model gene for the study of this G-rich motif dependent splicing is Bcl-X, a regulator of cell death which has two major spliced forms. The dominant isoform is the longer, Bcl-XL, which is anti-apoptotic. A switch in splicing leads to formation of the shorter form, Bcl-XS, which is pro-apoptotic [@Boise1993]. This switch involves differential donor site usage in the splicing of the second intron. Garneau et al. showed that alternative splicing of Bcl-X is mediated by hnRNP F/H binding to two exonic G-rich regions, one upstream of the Bcl-XL donor site, and another downstream of the Bcl-XS donor site [@Garneau2005]. Mutation of these G-runs abolished hnRNP binding and removed the effect of recombinant hnRNP F treatment on Bcl-X splicing. Weldon et al. recently demonstrated that both of the G-rich regions are capable of forming G4s *in vitro*, and that treatment of *in vitro* splicing assays with ellipticine derived G4 binding agents was able to alter the ratios of the spliced forms [@Weldon2017; @Weldon2018]. The RNA recognition domain of hnRNP F binds to single stranded DNA, suggesting that G4s modulate splicing by preventing hnRNP binding [@Dominguez2010; @Samatanga2013]. Furthermore, Weldon et al. modelled RNA structure constraining the nucleotides which form the G4s, and suggested that G4 formation near the Bcl-XS donor site may abolish a long stem loop structure in the pre-mRNA [@Weldon2018]. Stem loop structures are thought to inihibit donor-site usage [@Eperon1988; @Nasim2002], meaning that G4 formation at this site might promote Bcl-XS production. Weldon et al. suggested that  G4 upstream of the Bcl-XL donor site might overlap with the donor splice site itself, explaining how G4 formation at this site has an inhibitory effect on its usage, through blocking of U1 snRNP binding [@Weldon2018].

* 3'UTR enconucleolytic
* APA

\newpage

![**G-Quadruplexes control the splicing of Bcl-X pre-mRNA** **a)** Diagram showing the splice isoforms Bcl-XS and Bcl-XL, which are derived from alternative splicing of exon 2. **b)** RNA structure of Bcl-X-681 (a synthetic transcript derived from Bcl-X) showing positions of G4s and splice donor/acceptor sites. Formation of a G4 near Bcl-XS donor site promotes usage of the XS splice donor. Formation of a G4 near Bcl-XL donor site inhibits usage of the XL splice donor. Adapted from Weldon et al. 2018](./figures/bclx_splicing.svg)

\newpage

### mRNA Stability and Translation

\label{ssec:mrna_stability}

G4s which form in mRNA have the potential to be more stable than those formed in DNA, as they may not have to compete with double stranded forms. Furthermore, structural studies have suggested that the extra hydroxyl groups in ribose compared to deoxy-ribose allow RNA G4s to form more hydrogen bonds within the Quadruplex structure [@Collie2010]. This increases the enthalpic favourability of the RNA G4 whilst also reducing the entropic cost of hydrogen bonds with ordered water molecules. Analysis of yeast and human mRNAs using DMS-seq and SHAPE-seq have shown that eukaryotic mRNAs are much less structured *in vivo* than predicted from *in silico* and *in vitro* analyses [@Rouskin2014; @Ding2014; @Loughrey2014; @Spitale2015]. Rouskin et al. showed that this reduction in structure was ATP-dependent, suggesting that active unwinding of mRNA secondary structure occurs regularly *in vivo* [@Rouskin2014]. Kwok et al. have used rG4-seq to identify mRNA G4s that form *in vitro* in the presence of potassium or pyridostatin, by their ability to stall RT. They found a total of 3383 G4s in mRNAs, of which 62% were contained in the 3' UTR, 16% in the 5' UTR and 21% were in the CDS [@Kwok2016]. Simulation of RNA folding with the RNA G4 region constrained yielded very different structures on average to folding without constraints, suggesting that G4 formation might act as a molecular switch, changing the overall structure of the mRNA [@Kwok2016]. RNA G4s have been implicated in a number of regulatory roles, including translational regulation in the 5' UTR, and mRNA stability in the 3' UTR. Guo and Bartel used a similar DMS-seq method to identify G4s which fold *in vivo*, but found that the majority were in an unfolded state [@Guo2016]. This unfolding was not dependent on the cellular levels of ATP, suggesting that RNA G4s may be kept unfolded *in vivo* through some passive process, e.g. through protein binding to ssDNA [@Guo2016]. Biffi et al. have used the BG4 antibody to visualise RNase sensitive G4 structures in the cytoplasm of human cells, however, suggesting that at least a fraction of RNA G4s are formed *in vivo* [@Biffi2014].

Translational initiation in eukaryotes begins with the binding of a 43S preinitiation complex (PIC) to the 5' cap of the linear mRNA [@Hinnebusch2014]. The PIC is pre-loaded with a Methionine aminoacyl-tRNA, and scans along the 5' UTR of the mRNA until the first methionine codon AUG is identified. This identification catalyses the recruitment of the large 60S subunit of the ribosome, to produce a translationally active complex [@Hinnebusch2014]. The 5' UTR of the mRNA often contains structures such as hairpins which must be resolved by RNA helicases to allow the passage of the PIC. G4s which form in the 5' UTR can have similar consequences, impeding the scanning of the PIC, or, if close enough to the m7G cap, preventing the loading of the complex onto the mRNA [@Bugaut2012].

A number of oncogenes, such as NRAS, BCL-2, FGF-2, and VEGF have been identified as containing 5' UTR G4s [@Kumari2007; @Shahid2010; @Morris2010]. When folded, these often to have repressive effects on translation. The KRAS 5' UTR, for example, causes reduced expression when attached to a luciferase reporter, compared to the same 5' UTR with mutations that abolish G4 folding [@Kumari2007]. Wolfe et al. showed that treatment of cells with silvestrol, an inhibitor of the RNA helicase EIF4A, caused ribosome stalling on the 5' UTRs of genes containing PG4s [@Wolfe2014]. EIF4A has previously been identified as a required helicase for the translation of transcripts with long and structurally complex 5' UTRs [@Parsyan2011]. This evidence suggests that 5' UTR G4s can act as negative regulations of translation initiation.

Interestingly, studies of the FGF-2 and VEGF 5' UTR have identified G4s which act as positive regulators of translation [@Bonnal2003; @Morris2010]. FGF-2 is a human growth factor encoding gene, whose mRNA produces five different protein isoforms through differential start codon usage. One of these isoforms is produced through cap-dependent translation, whilst the other four are translated in a cap-independent manner from an internal ribosome entry site (IRES) [@Bonnal2003]. IRESs are sequence-dependent structures in the 5' UTR which can control cap independent binding and initiation of ribosome scanning. They are often found in very long and GC-rich 5'UTRs which are more highly structured, and maintain their translation during stresses such as hypoxyia [@Jackson2013]. Bonnal et al. found that the IRES of FGF-2 contains a RNA G4 forming motif, which was able to stall RT *in vitro* [@Bonnal2003]. VEGF has an extremely long 5' UTR of greater than 1000 bases, which contains two distinct IRES [@Morris2010]. Morris et al. showed that a G4 structure located at IRES-A was essential for internal ribosome initiation [@Morris2010]. These studies indicates that 5' UTR G4s can have complex roles in translational regulation, acting as both positive and negative regulators, or promoting translation of alternative isoforms.

G4s in 3' UTRs have the potential to act as *cis* regulatory elements controlling translation efficiency or mRNA stability. One hypothesised method of action is through competition between G4 structures and unstructured UTRs which can bind microRNAs (miRNAs). miRNAs are short (20-25bp) regulatory RNAs which act through complementarity to their target mRNA. Hybridisation of miRNAs to their targets results in changes in mRNA degradation through P-body localisation, or translational repression [@Ambros2004]. The local structure of 3' UTRs has been identified as a key modulator of their binding efficiency [@Long2007]. Rouleau et al. performed a characterisation of the overlap between predicted miRNA binding sites and PG4 sequences in human 3' UTRs [@Rouleau2017]. The authors found that 54% of PG4s in 3' UTRs overlapped with at least one predicted miRNA target site, though the enrichment or statistical significance of this overlap were not calculated. A candidate mRNA, FADS2, which is regulated by mir331-3p, was tested for G4 formation. Rouleau et al. showed that a 3' UTR sequence which binds mir331-3p also forms a G4, and that G4 formation prevents miRNA binding *in vitro* [@Rouleau2017].

\newpage

![**G4 function in translation and mRNA stability** Schematic showing how 5' and 3' UTR G4s can effect translation of mRNAs. **a)** G4s in the 5' UTR of mRNAs can cause a reduction in gene expression by preventing scanning of the PIC along the UTR. **b)** and **d)** G4 formation in the 3'UTR could block the recognition of miRNA binding sites by Argonaute proteins. **c)** G4s which form in long and structured 5' UTRs can contribute to internal ribosome entry site (IRES) which promote cap independent initiation of translation](figures/mRNA_stability_translation.svg)

\newpage

### Epigenetics and Chromatin Remodelling

\label{ssec:chromatin}

Chromatin conformational dynamics are important regulatory mechanisms for gene expression. The major protein component of chromatin is made up of histones, which assemble into the nucleosomes around which DNA is wrapped. G4 formation is highly dependent upon chromatin status, and G4s isolated by BG4 ChIPseq have been shown to form in regions of open chromatin [@Hansel2016]. Hänsel-Hertsch et al. showed that treatment of cells with the histone deacetylase inhibitors, which cause a relaxation of the chromatin, also increased the number of G4 loci identified by ChIPseq [@Hansel2016].

G4s have been shown to cause transient and non-heritable changes to epigenetic status through the causing of DSB, which are known to induce phosphorylation of H2AX (gamma-H2AX) via an ATM-dependent pathway [@Rodriguez2013; @Shanbhag2010; @Ayrapetov2014]. A number of authors have also reported heritable epigenetic changes to the BU-1 and beta-globin gene loci caused by G4s in chicken cell lines [@Guilbaud2017; @Sarkies2010; @Sarkies2011; @Schiavone2014; @Schiavone2016]. These changes occur in cells lacking G4 resolving helicases, such as REV1, FANCJ, WRN and BLM [@Sarkies2010; @Sarkies2012; @Schiavone2014]; polymerases, such as PrimPol [@Schivone2016]; or in cells treated with G4 stabilising agents, such as NMM, PhenDC3 and PDS [@Guilbaud2017]. None of these conditions caused DSB or mutations in the BU-1 gene, but did cause loss of active chromatin marker H3K4me3, and gain of repressive H3K9me3, as well as an increase in repressive DNA methylation. There was a concommitant decrease in gene expression accompanying these changes [@Sarkies2010; @Sarkies2011; @Schiavone2014; @Schiavone2016; @Guilbaud2017]. During replication, parental histones from upstream of the replication fork are removed; and recycled into the newly synthesised DNA behind the replication machinery [@DeKoning2007]. This recycling maintains the epigenetic status of replicated DNA. Sarkies et al. proposed that in the absence of G4 unwinding helicases, or in the presence of a stabilising drug, leading strand synthesis is stalled at the G4, and uncoupled from histone recycling. This leads to the incorporation of newly synthesised H3/4 without the correct epigenetic markers, rather than parental histones [@Sarkies2010; @Sarkies2011; @Schiavone2014; @Schiavone2016; @Guilbaud2017]. Abolishment of the G4 motif through mutation removed this heritable change in epigenetic status. Furthermore, switching of the strand of the G4, which was suggested to switch the replication strand of the G4 from the leading to the lagging strand, also abolished the heritable change in expression [@Sarkies2010; @Sarkies2012; @Schiavone2014; @Schiavone2016].

\newpage

![**G4 stalling causes epigenetic reprogramming** Illustration demonstrating the mechanism by which G4 stalling is thought to cause epigenetic reprogramming. Histone recycling from upstream of the replication fork into newly synthesised DNA is tightly linked to DNA replication. This is important to maintain the correct histone modifications in the daughter cells. When G4s form in the leading strand of the DNA, they cause blockages that decouple DNA replication from histone removal, resulting in the insertion of newly synthesised histones into the DNA. This causes the loss of histone modifications from this region of the DNA.](figures/histone_recycling.svg)

\newpage

Chromatin remodelling is the process by which modifications are made to histones. These modifications include switching of histone isoforms, covalent modifications, and histone sliding. The latter involves the use of energy from ATP hydrolysis to translocate histones along the DNA polymer without ever removing them [@Langst2015]. One class of remodellers is SWI/SNF complexes, which are able to slide or eject nucleosomes. ATRX is an X-linked SWI/SNF family member, which is mutated in ATR-X syndrome, an intellectual disability disorder. The ATRX protein has been shown to locate at heterochromatic and pericentromeric regions, as well as at telomeres [@McDowell1999], regions which are enriched in G-rich PG4 structures. During ATR-X syndrome, global DNA methylation patterns at these regions are disregulated [@Gibbons2000]. Knockdown of ATRX during the S-phase of the cell cycle, when DNA replication occurs, results in an increase in the DSB response marker gamma-H2AX at telomeric regions [@Wong2010]. Heaphy et al. showed that human pancreatic cancers which had lost the expression of ATRX or its partner protein DAXX had elongated telomeres that were characteristic of the Alternative Lengthening of Telomeres pathway. This pathway is telomerase independent, suggesting that ATRX is important for normal telomere repair [@Heaphy2011].

ATRX mutation also effects gene expression. Gibbons et al. identified gene expression changes in the alpha globin gene cluster in response to ATRX mutation [@Gibbons1991]. More recently, Law et al. performed ChIP-seq using an ATRX specific antibody to identify other gene targets. They found that ATRX binds specifically to G-rich tandem repeats which form G4s *in vitro*. Furthermore, the strength of the effect of ATRX mutation on gene expression was proportional to the number of tandem G-rich repeats. This suggests that the effect is increased by the number or stability of G4s [@Law2010].

\newpage

### R loop Biology

\label{ssec:rloop_csr}

R loops have been characterised *in vitro* in a number of organisms, from initial studies in *E. coli* [@Drolet1995] to genome wide profiling in *S. cerevisiae* [@Chan2014] and *H. sapiens* cell lines [@Chen2017]. R loops have been implicated in maintaining chromatin status and gene expression [@Powell2013; @Chen2015], as well as termination of transcription [@Skourti-Stathaki2011; @Skourti-Stathaki2011]. Chen et al. used an inactive RNase H enzyme, which binds specifically to RNA:DNA heteroduplexes, to immunoprecipitate R-loops in human chromatin [@Chen2017]. They identified a large number of R loop peaks located at promoter proximal regions, suggesting R loop formation early in transcription. The authors found that R loops were more common on genes that showed higher levels of proximal pausing, suggesting that R loop formation may be associated with the pausing of Pol II [@Chen2017]. Furthermore, regions containing R loops were also strongly associated with coding strand G4s, with 54% of R loop peaks containing a G4 observed by G4seq [@Chen2017; @Chambers2015]. G4 formation in the coding strand, which is in a single stranded state in the R loop, may help to stabilise the RNA:DNA hybrid [@Duquette2004]. 

R loop formation has been identified as a key component of class switch recombination (CSR) of immunoglobulin genes [@Daniels1995; @Reaban1990; @Yu2003]. CSR is the process by which the constant region of the immunoglobulin heavy chain (IgH) gene is switched, changing the properties of the resultant antibody [@Stavenezer1996]. This process is mediated by the activation-induced cytidine deaminase (AID), an enzyme which deaminates cytosines in ssDNA, within the switch regions of the IgH locus. This deamination results in DSB, causing the recombination of the gene. Transcription of the IgH gene produces a long non-coding RNA (lncRNA) with a constitutive first exon containing the variable antigen binding domain. Splicing to an alternative second exon is determined by the pattern of cytokines that the cell is exposed to. The intron lariat removed by this splicing is debranched and re-linearised by the enzyme DBR1 [@Zheng2015]. This intron lariat then forms an R-loop with the switch region of the IgH, resulting in ssDNA on the coding strand which becomes the substrate for AID. Zheng et al. showed that AID is able to bind specifically to G4 DNA and RNA, and that targeting of AID to the switch region is via binding to the debranched RNA lariat [@Zheng2015]. Furthermore, Ribeiro de Almeida et al. identified that resolution of G4s in the intron derived RNA is required for R-loop formation and CSR. This G4 unwinding was found to occur via the RNA Helicase DDX1 [@RibeirodeAlmeida2018].

\newpage

![**An R-loop and G4 dependent mechanism for Class Switch Recombination** **a)** Schematic showing the formation of an R-loop from G-rich RNA and C-rich template DNA, behind transcribing RNA Polymerase II. The G-rich coding strand is left in a single stranded form, and may fold into G4s. **b)** Proposed mechanism for Class Switch Recombination of IgH gene through R-loop formation. Transcription and splicing of the IgH locus results in a G-rich intron lariat which is debranched by DBR1. This folds into RNA G4s which recruit AID. Resolution of G4s by DDX1 promotes R loop formation and targets AID to the switch region, causing deamination and double strand breaks.](figures/rloops.svg)

\newpage

## Role of G-Quadruplexes *in planta*

\label{sec:g4s_in_plants}

### G-Quadruplex Distribution

\label{ssec:plant_distribution}

Whilst the majority of G4 studies have been performed in mammalian systems, particularly *H. sapiens*, there is a growing body of evidence for functions of G4s in plant species including *Arabidopsis thaliana*, *Zea mays*, and *Oryza sativa* [@Mullen2010; @Andorf2014; @Wang2015; @Garg2016]. Monocotyledonous flowering plants such as *O. sativa* and *Z. mays* generally have higher GC content than dicotyledonous flowering plants and non-flowering plants [@Smarda2014]. Monocots also have higher PG4 content, presumably as a result of this [@Garg2016]. Three tetrad PG4 densities of plants tend to be lower than those of the well studied *H. sapiens* and *M. musculus*, however. *Arabidopsis thaliana*, the model plant, has a small genome and low PG4 density, with only 1200 three tetrad PG4s [@Mullen2010]. This represented a greater than two fold depletion of PG4s over what was expected in a windowed markov chain modelled genome [@Mullen2010]. Mullen et al. noted that Arabidopsis contains 43000 two tetrad PG4s, however, which they suggested might be stable at the temperature ranges that Arabidopsis lives at [@Mullen2010; @Mullen2012]. Garg et al. identified 1331 genes which have conserved two tetrad or greater PG4s in all dicot species tested, suggesting functional conservation of these motifs [@Garg2016].

Arabidopsis telomeres are made up of the heptameric sequence TTTAGGG [@Richard1988]. This is slightly different sequence repeat from that found in most vertebrates, including *H. sapiens* [@Kim2018]. A number of plants species have also been reported to contain the same TTAGGG repeat as vertebrates, however [@Kim2018]. Both Arabidopsis and *O. sativa* have been shown to be sensitive to the telomerase inhibitor telomestatin [@Zhang2006]. Telomestatin is a G4 binding agent which has been shown to interact with the human telomeric repeat, and inhibit telomerase causing telomere shortening *in vivo* [@Kim2002]. Zhang et al. showed that sensitivity of cultured Arabidopsis and *O. sativa* cells to telomestatin occurs via the same method of telomere shortening, suggesting that plant telomeres may also form G4s *in vivo* [@Zhang2006].

### Translational Regulation

\label{ssec:plant_translation}

An enrichment of PG4s in 5' UTRs has been identified in a number of plant species, including *A. thaliana* [@Mullen2010; @Garg2016], *Z. mays* [@Andorf2014], and *O. sativa* [@Wang2015]. Several of these have been shown to form in RNA [@Kwok2015; @Cho2018]. Kwok et al. identified a PG4 motif on the coding strand of the ATR gene in Arabidopsis, which forms in the 5' UTR of the mRNA [@Kwok2015]. They showed using a transient reporter gene assay that this G4 inhibits translation, but not transcription, presumably by preventing scanning of the UTR by the ribosomal PIC [@Kwok2015].

Cho et al. also identified an RNA G4 which forms in the 5'UTR of SMXL4/5, a gene involved in the differentiation of the phloem [@Cho2018]. They found that this G4 was bound and stabilised by the G4 binding protein JULGI. Using SELEX, they showed that JULGI was able to bind a number of PG4 forming sequences [@Cho2018]. G4 formation in the SMXL4/5 5'UTR resulted in a reduction of the translation of mRNA, and a restriction of the phloem, whilst knockout of JULGI strongly promoted phloem development [@Cho2018]. A PG4 motif was identified in the 5'UTR of SMXL4/5 in all vascular plants analysed, suggesting that this mechanism was evolved early in the development of vascular plants [@Cho2018]. 

### Plant Development

\label{ssec:plant_dev}


It is plausible that molecular mechanisms involving G4s are central to many plant development pathways [@Cho2018; @Nakagawa2012]. Nakagawa et al. analysed the effects of the G4 binding agents NMM and Berberine on Arabidopsis development, and found that both drugs caused defects in leaf polarity [@Nakagawa2012]. Furthermore, plants with double mutations in the genes *ASYMMETRIC LEAVES 1* and *2* were hypersensitive to G4 binding agents [@Nakagawa2012].

### Stress Response

\label{ssec:plant_stress}

Analysis by Mullen et al. identified that the greatest enrichment of two tetrads in the Arabidopsis genome was inside genic regions, and that these PG4s might form in mRNA [@Mullen2010]. Furthermore, they demonstrated that many two tetrad PG4s are only strongly stable in high potassium concentrations, and postulated that these might act as molecular switches sensing environmental salt concentrations [@Mullen2012]. Since cellular potassium concentrations are increased from around 100mM to as much as 600mM during osmotic stresses (e.g. drought stress) [@Walker1996; @Hummel2010], they suggested that drought stress responsive genes could contain these molecular switches. Analysis of the *Z. mays* genome by Andorf et al. has shown that PG4s are indeed enriched in a number of gene ontology classes involved in stress response, including hypoxia response and nutrient deprivation responsive genes [@Andorf2014].

\newpage

## Summary

\label{ssec:intro_summary}


There is now clear evidence for the *in vivo* function of DNA and RNA G4s in a wide range of molecular processes, including genome stability, replication, transcription, mRNA processing, and translation. In this thesis, we will introduce a new method for the detection of DNA and RNA G4s, using information learned from high throughput sequencing datasets [@Chambers2015; @Kwok2016]. We will then apply this and other methods to the detection of PG4s in the Arabidopsis genome, and identify their enrichment in various genic features including 5'UTRs and CDSs. Finally, we will attempt to identify the effect of G4 stabilisation on Arabidopsis gene expression by treating plants with the G4 binding agent NMM. We will attempt to shed some light on whether the G4 rich genes affected by NMM are naturally regulated through G4 dependent mechanisms.

\newpage
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJWY0NnREpSMjNXMHFJN2VlIjp7In
RleHQiOiJhcyB3ZWxsIGFzIGNoYW5nZSB0b3BvbG9naWVzIG9u
IGEiLCJzdGFydCI6Mzc5LCJlbmQiOjQxMn0sIjVTWVhPNnU1Yn
laUDJZa2giOnsidGV4dCI6IiFbKipTdHJ1Y3R1cmUgb2YgYSBH
LVF1YWRydXBsZXgqKiAqKmEpKiogVGhlIG1vbGVjdWxhciBzdH
J1Y3R1cmUgb2YgYSBHLXF1YXJ0ZXTigKYiLCJzdGFydCI6NTQ1
NywiZW5kIjo2NDk4fSwiVVN3ZWJmOUh6QW5CcWs2UiI6eyJ0ZX
h0IjoiIyMjIEctUXVhZHJ1cGxleCBQcmVkaWN0aW9uIGZyb20g
U2VxdWVuY2UiLCJzdGFydCI6Njk4MiwiZW5kIjo3MDIzfSwiRW
s0SWo2UUVFSVJOdFRvQyI6eyJ0ZXh0IjoiZiBkaXNvcmRlcmVk
IHNpbmdsZSBzdHJhbmRlZCBETkEiLCJzdGFydCI6MTE2MDYsIm
VuZCI6MTE2Mzh9LCJScVFYUnFDWUs5ODdGRWpHIjp7InRleHQi
OiIhWyoqRy1RdWFkcnVwbGV4IFN0YWJpbGlzaW5nIExpZ2FuZH
MqKiBTdHJ1Y3R1cmVzIGFuZCBtb2RlIG9mIGFjdGlvbiBvZiBQ
eXJpZG9z4oCmIiwic3RhcnQiOjE1Mjg1LCJlbmQiOjE1NzEwfS
wiSWpMVVJSYU92YUFDZUh3biI6eyJ0ZXh0IjoiT25lIGNvbW1v
biBhcHByb2FjaCB0byBzdHVkeWluZyB0aGUgZWZmZWN0IG9mIE
c0IHN0YWJpbGlzYXRpb24gb24gYmlvbG9naWNhbCBwcuKApiIs
InN0YXJ0IjoxNTg3NywiZW5kIjoxNjA0Nn0sIktiSUJOUUtJME
9zaUJET1AiOnsidGV4dCI6Ik5hdHVyYWxseSBmbHVvcmVzY2Vu
dCIsInN0YXJ0IjoxNjIxMiwiZW5kIjoxNjIzM30sIm9PMkRzT2
tqQUxrR3lUT00iOnsidGV4dCI6Ii4gVGhlIGF1dGhvcnMgc2hv
d2VkIHVzaW5nIHRoaXMgdGVjaG5pcXVlIHRoYXQgb25seSAzMC
Ugb2YgZXhwZXJpbWVudGFsbHkgb2JzZXLigKYiLCJzdGFydCI6
MTgxMjYsImVuZCI6MTgyODZ9LCJ0R3Jycnd5ekZqbzA1TXl6Ij
p7InRleHQiOiJ0aGF0IHJlZ2lvbnMiLCJzdGFydCI6MTg5Mzcs
ImVuZCI6MTg5Mzd9LCJId3ZHb2lPOW9MZkh0dUNGIjp7InRleH
QiOiJUaGlzIG1vZGVsIHdhcyB2ZXJ5IHN1Y2Nlc3NmdWwgYXQg
aWRlbnRpZnlpbmcgUXVhZHBhcnNlciBjb25mb3JtaW5nIG1vdG
lmcyB3aGlj4oCmIiwic3RhcnQiOjI0MjczLCJlbmQiOjI0NTMw
fSwic0hNSzA5NHIyc0FqRlphZiI6eyJ0ZXh0IjoiRy1RdWFkcn
VwbGV4IHN0YWJpbGl0eSBwcmVkaWN0aW9uIHVzaW5nIE1hY2hp
bmUgTGVhcm5pbmciLCJzdGFydCI6MjAyMjcsImVuZCI6MjAyOD
N9LCJDdFdWQmI3RkI5MWxMaU54Ijp7InRleHQiOiJHNHMgd2l0
aCBsb29wIGxlbmd0aCBvZiB1cCB0byAxNWJwIiwic3RhcnQiOj
MwNzgyLCJlbmQiOjMwODI2fSwicjhpOWhzQ3RRVkNBcm5rMCI6
eyJ0ZXh0IjoiKEVkZHkgJiBNYWl6ZWxzIDIwMDYpIiwic3Rhcn
QiOjMyMjIzLCJlbmQiOjMyMjMxfSwiMlVmcTE0cTFqRThRYmpP
SiI6eyJ0ZXh0IjoiVGhlIE5IRUlJSSBjb250YWlucyBhIG51bW
JlciBvZiBHLXJpY2ggdHJhY3RzIHdoaWNoIGhhdmUgYmVlbiBz
aG93biB0byBmb3JtIEc0c+KApiIsInN0YXJ0IjozMzQ1MSwiZW
5kIjozMzU2NX0sInBPR2c0YWY1aThkM3FRTWsiOnsidGV4dCI6
ImV2ZW4gd2hlbiBvdGhlciBmYWN0b3JzIHN1Y2ggYXMgZ2VuZS
BmdW5jdGlvbiBhcmUgY29udHJvbGxlZCBmb3IuIEl0IGhhcyBi
ZWVuIHPigKYiLCJzdGFydCI6MzQ3ODQsImVuZCI6MzUwMDd9LC
JMUTN4SEJ0N3R0TzFHSEplIjp7InRleHQiOiJUcmFuc2NyaXB0
aW9uIHByb2dyZXNzZXMgYnkgdXNpbmcgdGhlIHRlbXBsYXRlIH
N0cmFuZCBhcyBhbiBhbnRpc2Vuc2UgY29weSB0byBy4oCmIiwi
c3RhcnQiOjM1MjU0LCJlbmQiOjM1OTM2fSwieUlzSGdxWWxIb3
dIdmdROCI6eyJ0ZXh0IjoiSXQgaXMgYWxzbyBwb3NzaWJsZSwg
aG93ZXZlciwgdGhhdCB0aGlzIGVmZmVjdCBjb3VsZCBiZSBwYX
J0aWFsbHkgZHVlIHRvIGdyZWF0ZeKApiIsInN0YXJ0IjozNzYx
MiwiZW5kIjozNzcyNX0sIjVFOUZ0V2VlV3c0VHdLenQiOnsidG
V4dCI6IlRoZSBHNC1mb3JtaW5nIHBvdGVudGlhbCBvZiB0aGVz
ZSBnZW5lcyBhbHNvIHRlbmRlZCB0byBiZSBncmVhdGVyIG9uIH
RoZSBjb2RpbmfigKYiLCJzdGFydCI6Mzk3ODEsImVuZCI6Mzk4
Njd9LCJ2OURyQXkyeDlaVnBSNGlIIjp7InRleHQiOiIhWyoqRy
1RdWFkcnVwbGV4ZXMgY29udHJvbCB0aGUgc3BsaWNpbmcgb2Yg
QmNsLVggcHJlLW1STkEqKiAqKmEpKiogRGlhZ3JhbSBzaG934o
CmIiwic3RhcnQiOjQ2OTAwLCJlbmQiOjQ3NDM0fSwiTXA1OTF5
clFhenB1bVQwQSI6eyJ0ZXh0IjoiYmVjYXVzZSBSTkFzIHRlbm
QgdG8gZm9ybSBtb3JlIGNvbXBsZXggc3RydWN0dXJlcyIsInN0
YXJ0Ijo0NzU5NywiZW5kIjo0NzU5OH0sIngxaUpsd29wOG9ubk
R3YTIiOnsidGV4dCI6IlNIQVBFLXNlcSIsInN0YXJ0Ijo0ODA1
NSwiZW5kIjo0ODA2NH0sImloTUxXQ0RVN25EMWt2eVUiOnsidG
V4dCI6Ikc0IGZvcm1hdGlvbiBtaWdodCBhY3QgYXMgYSBtb2xl
Y3VsYXIgc3dpdGNoIiwic3RhcnQiOjQ4ODY0LCJlbmQiOjQ4OT
A4fSwiVGZJS0ZBWTRqWUlrWUVQOCI6eyJ0ZXh0IjoiTW9ub2Nv
dHMgYWxzbyBoYXZlIGhpZ2hlciBQRzQgY29udGVudCwiLCJzdG
FydCI6NjM5MTcsImVuZCI6NjM5NTV9LCJOaHRmQXRIUEpTalZQ
MzR3Ijp7InRleHQiOiJtaWdodCBiZSBzdGFibGUgYXQgdGhlIH
RlbXBlcmF0dXJlIHJhbmdlcyIsInN0YXJ0Ijo2NDQ5NCwiZW5k
Ijo2NDUzNX0sIkJFQWtpQmJ0OENGbjMyb3UiOnsidGV4dCI6Ik
EgUEc0IG1vdGlmIiwic3RhcnQiOjY2NzY2LCJlbmQiOjY2Nzc3
fSwicURJbmZUMnRpcHBRZTREMyI6eyJ0ZXh0IjoiQW5hbHlzaX
MgYnkgTXVsbGVuIGV0IGFsLiBpZGVudGlmaWVkIHRoYXQgdGhl
IGdyZWF0ZXN0IGVucmljaG1lbnQgb2YgdHdvIHRldHJhZOKApi
IsInN0YXJ0Ijo2NzUwNCwiZW5kIjo2ODQxOX0sImdaUjlIUG03
Z3Y1aEJjV28iOnsidGV4dCI6IiMjIFJvbGUgb2YgRy1RdWFkcn
VwbGV4ZXMgKmluIHBsYW50YSoiLCJzdGFydCI6NjMzMTcsImVu
ZCI6NjMzNTR9LCJvajQxNldaUDU5YzBpSFMwIjp7InRleHQiOi
JJdCBpcyBwbGF1c2libGUgdGhhdCIsInN0YXJ0Ijo2Njk5OSwi
ZW5kIjo2NzAxOX19LCJjb21tZW50cyI6eyJyb2Q2QUdUM2ZJdj
l3VUtGIjp7ImRpc2N1c3Npb25JZCI6IlZjQ2dESlIyM1cwcUk3
ZWUiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiVG8gd2hhdCBleHRlbnQgaXMgYSBxdWFkcnVwbGV4IGEg
XCJjaGFuZ2VcIiB0byBzb21lIGJhc2UgbGluZSBzdGF0ZSwgb3
IgaXMgaXQgYW4gaW5uYXRlIHByb3BlcnR5IG9mIHRoZSBzZXF1
ZW5jZT8iLCJjcmVhdGVkIjoxNTMxODM3NTMyOTY5fSwiekhTY1
dpTEpkMVc2SWNORyI6eyJkaXNjdXNzaW9uSWQiOiI1U1lYTzZ1
NWJ5WlAyWWtoIiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEwNj
c3IiwidGV4dCI6IklzIHRoaXMgZmlndXJlIGVudGlyZWx5IHlv
dXIgb3duIGNyZWF0aW9uPyBOb3QgZXZlbiBpbnNwaXJlZCBieS
BhIHJlZmVyZW5jZT8iLCJjcmVhdGVkIjoxNTMxODQ1MDUyMzA0
fSwiVjJmNWJNdzBOcFVHeFgwQyI6eyJkaXNjdXNzaW9uSWQiOi
JVU3dlYmY5SHpBbkJxazZSIiwic3ViIjoiMTAyMjA1Nzk3Mjc2
OTQxMDEwNjc3IiwidGV4dCI6IkkgZmVlbCBsaWtlIHRoaXMgc2
VjdGlvbiBjb3VsZCBiZW5lZml0IGZyb20gc29tZSBhY3RhdWwg
bnVtYmVycyBvbiB0aGUgcGVyZm9ybWFuY2Ugb2YgdGhlc2UgYW
xnb3JpdGhtcz8gV2lsbCB5b3UgaW5jbHVkZSB0aGVzZSBpbiB0
aGUgY2hhcHRlciBvbiB5b3VyIG5ldyBtZXRob2Q/IGF0IHRoZS
Btb21lbnQgc3RhdGVtZW50cyBsaWtlIFwibWFueSBmYWxzZSBw
b3NpdGl2ZXNcIiBzZWVtIHVuc3VwcG9ydGVkLiIsImNyZWF0ZW
QiOjE1MzIzNDM3NTU4MjZ9LCJza2QwTnpCNkhLRklhcnBKIjp7
ImRpc2N1c3Npb25JZCI6IkVrNElqNlFFRUlSTnRUb0MiLCJzdW
IiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiSXMg
aXQgZ2VuZXJhbGx5IHRoZSBjYXNlIHRoYXQgcGVvcGxlIGFzay
BpZiB0aGUgc2luZ2xlIHN0cmFuZGVkIG9saWdvIHdpbGwgZm9y
bSBhIEc0PyBXaGF0IGFib3V0IHRoZSBkb3VibGUgaGVsaXg/IE
lzbid0IHRoZSBrZXkgcXVlc3Rpb24gd2hldGhlciBhIGRvdWJs
ZSBzdHJhbmRlZCBtb2xlY3VsZSB3aWxsIGZhdm91ciBhIEc0IG
luIG9uZSBvZiB0aGUgc3RyYW5kcyBvdmVyIHRoZSBoZWxpeD8g
T3IgaXMgaXQgYXNzdW1lZCB0aGF0IGluIGFsbCBjYXNlcyBzb2
1ldGhpbmcgZWxzZSBpcyBob2xkaW5nIHRoZSBETkEgaW4gc2lu
Z2xlIHN0cmFuZGVkIHN0YXRlPyIsImNyZWF0ZWQiOjE1MzIzND
M5NTkxMzJ9LCJvalpNZEZReWRQbTRZc1lqIjp7ImRpc2N1c3Np
b25JZCI6IlJxUVhScUNZSzk4N0ZFakciLCJzdWIiOiIxMDIyMD
U3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiSnVzdCB1c2UgdGhl
IGZpZ3VyZSBsZWdlbmQgdG8gZGVzY3JpYmUgdGhlIGZpZ3VyZS
wgcmF0aGVyIHRoYW4gdG8gY29udmV5IG5ldyBpbmZvcm1hdGlv
bi4gU29tZSBvZiB0aGlzIHN0dWZmIHNob3VsZCBkZWZpbmF0ZW
x5IGJlIGluIHRoZSBtYWluIHRleHQuIEFsc28gYnJlYWsgdXAg
aW50byBwYW5uZWxzLiIsImNyZWF0ZWQiOjE1MzIzNDU3NTI0Mj
R9LCI1TzlqU3pzcFljSGdheXJBIjp7ImRpc2N1c3Npb25JZCI6
IklqTFVSUmFPdmFBQ2VId24iLCJzdWIiOiIxMDIyMDU3OTcyNz
Y5NDEwMTA2NzciLCJ0ZXh0IjoiVGhpcyByZWFkcyBhcyB0aG91
Z2ggaXQgaXMgYSBkaXJlY3QgZm9sbG93IG9uIGZyb20gdGhlIH
ByZXZpb3VzIHNlbnRlbmNlLCBidXQgaW4gdGhhdCBzZW50ZW5j
ZSB5b3UgYXJlIHRhbGtpbmcgYWJvdXQgdGhlICpmb3JtYXRpb2
4qIG9mIEc0cyBhbmQgaW4gdGhpcyBzZW50ZW5jZSB5b3Ugc2Vl
bSB0byBiZSB0YWxraW5nIGFib3V0IHRoZSBiaW9sb2dpY2FsIG
NvbnNlcXVlbmNlcyBvZiB0aGVpciBmb3JtYXRpb24uIiwiY3Jl
YXRlZCI6MTUzMjM0NTg1NTU1NH0sIkh6RU5JVTF0Nk94NnZXQ0
YiOnsiZGlzY3Vzc2lvbklkIjoiS2JJQk5RS0kwT3NpQkRPUCIs
InN1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOi
JBcmUgc29tZSBvZiB0aGUgRzQgYmluZGluZyBtb2xlY3VsZXMg
eW91IHRhbGtlZCBhYm91dCBiZWZvcmUgbmF0dXJhbGx5IGZsdW
9yZXNjZW50PyIsImNyZWF0ZWQiOjE1MzIzNDYxMDMzNTR9LCJD
c25iZm5XeWpWbG1SblQyIjp7ImRpc2N1c3Npb25JZCI6Im9PMk
RzT2tqQUxrR3lUT00iLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEw
MTA2NzciLCJ0ZXh0IjoiU2hvdWxkIG1ha2UgdGhlIHBvaW50IH
RoYXQgdGhpcyBpcyBvbmx5IHNlcXVlbmNlIHdpdGggdGhlICpj
YXBhY2l0eSogdG8gZm9ybSBHNHMsIG5vdCBzZXF1ZW5jZSB0aG
F0IGlzIGluIGEgRzQgY29uZmlybWF0aW9uICppbiB2aXZvKi4g
WW91IHNob3VsZCBkaXNjdXNzIGhvdyB3ZWxsIGFubm90YXRpb2
5zIGJ5IENEIHNwZWMsIEJHNCBjaGlwc2VxLCBHNC1zZXEgZXRj
IGFncmVlIHdpdGggZWFjaCBvdGhlci4iLCJjcmVhdGVkIjoxNT
MyMzQ2MTk0NjM2fSwia3Q1Vk5keHFaR3pNdVZEaCI6eyJkaXNj
dXNzaW9uSWQiOiJ0R3Jycnd5ekZqbzA1TXl6Iiwic3ViIjoiMT
AyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6ImRvbid0IGxp
a2UgdGhpcyB3b3JkaW5nLiIsImNyZWF0ZWQiOjE1MzIzNDYyNT
U5MTN9LCIxbm1uWmFaTHFOd2lxUmxmIjp7ImRpc2N1c3Npb25J
ZCI6Ikh3dkdvaU85b0xmSHR1Q0YiLCJzdWIiOiIxMDIyMDU3OT
cyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiU2hvdWxkIHlvdSBtYWtl
IGl0IGV4cGxpY3QgdGhhdCBpdCBjYW5ub3QgZmluZCBub24tY2
9uZm9ybWluZyBHNHMuIiwiY3JlYXRlZCI6MTUzMjM1NjUyNDA1
NH0sIkZTSGlNMUJRQlk4dUdEaEMiOnsiZGlzY3Vzc2lvbklkIj
oic0hNSzA5NHIyc0FqRlphZiIsInN1YiI6IjEwMjIwNTc5NzI3
Njk0MTAxMDY3NyIsInRleHQiOiJUaGlzIHNlY3Rpb24gaXMgYW
xsIHByZXR0eSBmYWN0dWFsLiBQZXJoYXBzIHlvdSBzaG91bGQg
YWRkIHNvbWUgY3JpdGljYWwgZXZhbHVhdGlvbiBvbiB0b3AuIF
doaWNoIG1vZGVscyBhcmUgdGhlIGJlc3QuIFdoYXRzIHdyb25n
IHdpdGggdGhlIG90aGVycy4uLiIsImNyZWF0ZWQiOjE1MzIzNT
Y1NzE2NDV9LCJNMzhLYWFYa2wxellVR2t2Ijp7ImRpc2N1c3Np
b25JZCI6IkN0V1ZCYjdGQjkxbExpTngiLCJzdWIiOiIxMDIyMD
U3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiUmVtaW5kIG1lIGlm
IHRoaXMgaXMgbG9uZyBvciBzaG9ydD8iLCJjcmVhdGVkIjoxNT
MyMzU3MDMzMDA5fSwiRFRhZmdIWmhqYXpNVFU4TyI6eyJkaXNj
dXNzaW9uSWQiOiJyOGk5aHNDdFFWQ0FybmswIiwic3ViIjoiMT
AyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IlVuIGZvcm1h
dHRlZCBjaXRhdGlvbi4iLCJjcmVhdGVkIjoxNTMyMzYwODU2Nj
kzfSwiWGpiVnZyQXBVVVFScjE4diI6eyJkaXNjdXNzaW9uSWQi
OiIyVWZxMTRxMWpFOFFiak9KIiwic3ViIjoiMTAyMjA1Nzk3Mj
c2OTQxMDEwNjc3IiwidGV4dCI6IldoYXQgZGlkIEJHNCBjaGlw
IHNheSBhYm91dCB0aGlzIHJlZ2lvbj8iLCJjcmVhdGVkIjoxNT
MyMzYxMTEyMjY3fSwidWlWZ1VyMlZZYjdVbmlqQiI6eyJkaXNj
dXNzaW9uSWQiOiJwT0dnNGFmNWk4ZDNxUU1rIiwic3ViIjoiMT
AyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ik9yIHRoYXQg
aGlnaGx5IHRyYW5zY3JpYmVkIGdlbmVzIHNwZW5kIG1vcmUgb2
YgdGhlcmUgdGltZSBpbiBhIHNpbmdsZSBzdHJhbmRlZCBzdGF0
ZSwgZmF2b3VyaW5nIEc0IGZvcm1hdGlvbj8iLCJjcmVhdGVkIj
oxNTMyMzYyNTU2OTk5fSwiNVRZRmhpaEdpdW05V3JtUyI6eyJk
aXNjdXNzaW9uSWQiOiJMUTN4SEJ0N3R0TzFHSEplIiwic3ViIj
oiMTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ik1pZGRs
ZSBzZW50ZW5jZSBpbiB0aGlzIHBhcmEgZG9lc24ndCBzZWVtIH
RvIGJlbG9uZy4gTGFzdCBzZW50ZW5jZSBvZiBsYXN0IHBhcmEg
ZmVlbHMgbGlrZSBpdCBkb2Vzbi4gXG5cbk5pZXRoZXIgZXZpZG
VuY2Ugbm9yIGNpdGF0aW9uIGluIG1pZGRsZSBzZW50ZW5jZS4i
LCJjcmVhdGVkIjoxNTMyMzYyNzIwNjI2fSwiaGwzVExqdmdzT2
93UXJHMCI6eyJkaXNjdXNzaW9uSWQiOiJ5SXNIZ3FZbEhvd0h2
Z1E4Iiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidG
V4dCI6IkNhbiB0aGlzIGJlIHRlc3RlZCBieSBjb21wYXJpbmcg
c2ltaWxhcmx5IEctcmljaCBzZXF1ZW5jZXMgdGhhdCBlaXRoZX
IgZG8gb3IgZG8gbm90IG1hdGNoIHRoZSBRdWFkcGFyc2VyIG1v
dGlmPyIsImNyZWF0ZWQiOjE1MzIzNjQwNjQzMDR9LCJRUWI4Y0
pTdTZ3VjVINlhvIjp7ImRpc2N1c3Npb25JZCI6IjVFOUZ0V2Vl
V3c0VHdLenQiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2Nz
ciLCJ0ZXh0IjoiQXMgaW4gdGhlIEc0IGZvcm1pbmcgcG90ZW50
aWFsIGlzIGhpZ2ggQU5EIHRoYXQgcG90ZW50aWFsIHRlbmRzIH
RvIGJlIGdyZWF0ZXIgb24gdGhlIHRoZSBjb2Rpbmcgc3RyYW5k
IHRoYW4gdGhlIHRlbXBsYXRlIHN0YXJ0P1xuXG5PUlxuXG5UaG
UgRzQgZm9ybWluZyBwb3RlbnRpYWwgd2FzIGhpZ2ggb24gdGhl
IHRlbXBsYXRlIHN0cmFuZCBhbmQgYWxzbyB0aGUgY29kaW5nIH
N0cmFuZD8iLCJjcmVhdGVkIjoxNTMyMzY0Mjc1NTI3fSwiWlVn
U3M3eERlNHN4VkVyaCI6eyJkaXNjdXNzaW9uSWQiOiJ2OURyQX
kyeDlaVnBSNGlIIiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEw
Njc3IiwidGV4dCI6IklzIHRoaXMgeW91cidzIG9yIGlzIGl0IG
FkYXB0ZWQgZnJvbSBzb21ld2hlcmU/IiwiY3JlYXRlZCI6MTUz
MjM2NTE1ODI3Mn0sIkRzVlZvMDFkcDFNbE53MDciOnsiZGlzY3
Vzc2lvbklkIjoiTXA1OTF5clFhenB1bVQwQSIsInN1YiI6IjEw
MjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJSTkEgZm9ybX
MgY29tcGxleCBzdHJ1Y3R1cmUgYmVjYXVzZSBSTkEgZm9ybXMg
Y29tcGxleCBzdHJ1Y3R1cmU/IiwiY3JlYXRlZCI6MTUzMjM2Nz
k3NTMyOH0sImwzSjQxS0lGMjZLSExZVUwiOnsiZGlzY3Vzc2lv
bklkIjoieDFpSmx3b3A4b25uRHdhMiIsInN1YiI6IjEwMjIwNT
c5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaGF0IHdvdWxkIEc0
cyBsb29rIGxpa2UgaW4gU0hBUEUtc2VxPyIsImNyZWF0ZWQiOj
E1MzIzNjgwNDc2NTF9LCJvbzNFZmpWNm84UzR0UGJ1Ijp7ImRp
c2N1c3Npb25JZCI6ImloTUxXQ0RVN25EMWt2eVUiLCJzdWIiOi
IxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiT25seSBp
ZiBmb3JtYXRpb24gb2YgRzQgd2FzIGl0c2VsZiByZWd1bGF0aW
9uPyIsImNyZWF0ZWQiOjE1MzIzNjgxMTcxMDd9LCJpQ3dQallU
WjllWHZhYWdpIjp7ImRpc2N1c3Npb25JZCI6IlRmSUtGQVk0al
lJa1lFUDgiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2Nzci
LCJ0ZXh0IjoiaXMgdGhpcyBkZXRlcm1pbmVkIGp1c3QgdGhyb3
VnaCBRdWFkcGFyc2VyIG1hdGNoZXM/IiwiY3JlYXRlZCI6MTUz
MjM3NjM5ODgxOX0sImtKOW5acWhjSkVobXdYQUsiOnsiZGlzY3
Vzc2lvbklkIjoiTmh0ZkF0SFBKU2pWUDM0dyIsInN1YiI6IjEw
MjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJUaGlzIGlzIH
RoZSBmaXJzdCB0aW1lIHlvdSd2ZSBtZW50aW9uZWQgdGhlIGRp
ZmZlcmVuY2VzIGluIHRlbXBlcmF0dXJlIGFuZCB0aGVpciBlZm
ZlY3Qgb24gc3RhYmlsaXRpZXMgb2YgRzRzLiBQZXJoYXBzIHlv
dSBzaG91bGQgbWFrZSBtb3JlIG9mIHRoaXMgc29tZXdoZXJlPy
IsImNyZWF0ZWQiOjE1MzIzNzY0OTY2Njh9LCJwZ1ZYMGhaU1Ns
YnNvSXBlIjp7ImRpc2N1c3Npb25JZCI6IkJFQWtpQmJ0OENGbj
Myb3UiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0
ZXh0IjoiSXMgaXQgdGhlIHNhbWUgbW90aWYgb3IgZG9lcyBpdC
Btb3ZlIGFyb3VuZD8iLCJjcmVhdGVkIjoxNTMyMzc2NjIyMTQz
fSwiQlBHN1NGQmdHUTAySm1ldiI6eyJkaXNjdXNzaW9uSWQiOi
JxREluZlQydGlwcFFlNEQzIiwic3ViIjoiMTAyMjA1Nzk3Mjc2
OTQxMDEwNjc3IiwidGV4dCI6IkhvdyBkb2VzIHRoaXMgaW50ZX
JhY3Qgd2l0aCBjb2RvbiBiaWFzPyBDb3VsZCBpdCBleHBhbGlu
LCBvciBiZSBleHBsYWluZWQgYnkgY29kb24gYmlhcz8gSSBrbm
93IHRoYXQgeW91IHdpbGwgY29tZSBiYWNrIHRvIHRoaXMgcXVl
c2l0b24sIGJ1dCBpdCBtaWdodCBiZSB3b3J0aCBwb2ludGluZy
B0aGlzIG91dC4iLCJjcmVhdGVkIjoxNTMyMzc2Nzg0MTIxfSwi
QTJGSXBjNXEzT2duYUVyciI6eyJkaXNjdXNzaW9uSWQiOiJnWl
I5SFBtN2d2NWhCY1dvIiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQx
MDEwNjc3IiwidGV4dCI6IklzIHRoaXMgcmVhbGx5IGV2ZXJ5IH
BhcGVyIG9uIEc0cyBpbiBwbGFudHMhISEgVmlyZ2luIHRlcnJp
dG9yeSBpbmRlZWQhIiwiY3JlYXRlZCI6MTUzMjM3Njg1MjMxNX
0sIkkxZHJvZ1NBb3k2Y0cweUYiOnsiZGlzY3Vzc2lvbklkIjoi
b2o0MTZXWlA1OWMwaUhTMCIsInN1YiI6IjEwMjIwNTc5NzI3Nj
k0MTAxMDY3NyIsInRleHQiOiJJIGZlZWwgbGlrZSB0aGlzIGlz
IHNvbWV0aGluZyB5b3Ugd291bGQgc2F5IGFmdGVyIHlvdSBoYW
QgZmluaXNoZWQgcmV2aWV3aW5nIHRoZSBldmlkZW5jZS5cblxu
U2luY2UgeW91IGFscmVhZHkgdGFsa2VkIGFib3V0IHBobG9lbS
BkZXZlbG9wbWVudCBpbiB0aGUgbGFzdCBwYXJhLCB3aHkgbm90
IHN0YXJ0IHdpdGggc29tZSBsaWtlIGxpa2UgXCJmdXJ0aGVyIG
V2aWRlbmNlIGZvciBhIHJvbGUgb2YgRzRzIGluIGRldmVsb3Bt
ZW50IGNvbWVzIGZyb20uLi5cIlxuXG5Zb3UgY291bGQgZmluaX
NoIHdpdGggYSBqdWRnZW1lbnQgb24gdGhlIGV2aWRlbmNlIHRv
IHRoZSBlZmZlY3QgdGhhdCBpdCBpcyBwbGF1c2libGUuIiwiY3
JlYXRlZCI6MTUzMjM3NzA1MTg0NH0sInlKdmR5MXpCSElveUpO
d0IiOnsiZGlzY3Vzc2lvbklkIjoiVmNDZ0RKUjIzVzBxSTdlZS
IsInN1YiI6IjEwODUyMDAyOTMwMjI5NDY1MDQxNyIsInRleHQi
OiJJIHRoaW5rIHdoYXQgaSdtIHRyeWluZyB0byBnZXQgYXQgaX
MgdGhhdCB0aGUgbW9sZWN1bGUgaXMgcGxhc3RpYyBhbmQgY2Fu
IGhhdmUgbWFueSBkaWZmZXJlbnQgc3RhdGVzLCB3aGljaCBtaW
dodCBoYXZlIGRpZmZlcmVudCBiaW9sb2dpY2FsIGltcGxpY2F0
aW9ucyIsImNyZWF0ZWQiOjE1MzIzODI1MzU4OTR9LCJEOWRFYV
A0OWxQNTJ2b2Z3Ijp7ImRpc2N1c3Npb25JZCI6IjVTWVhPNnU1
YnlaUDJZa2giLCJzdWIiOiIxMDg1MjAwMjkzMDIyOTQ2NTA0MT
ciLCJ0ZXh0IjoidGhlcmUgaXMgYSBmaWd1cmUgbGlrZSB0aGlz
IGluIGV2ZXJ5IEc0IHBhcGVyIHNpbmNlIGFib3V0IDIwMDUsIE
kgY291bGQgY2l0ZSBhIHBvcHVsYXIgcmV2aWV3PyIsImNyZWF0
ZWQiOjE1MzIzODI2NTQxODV9LCJueUJnaXZzdVhQTnRSOEp0Ij
p7ImRpc2N1c3Npb25JZCI6Im9PMkRzT2tqQUxrR3lUT00iLCJz
dWIiOiIxMDg1MjAwMjkzMDIyOTQ2NTA0MTciLCJ0ZXh0IjoiaX
MgdGhpcyBhbiBpbXByb3ZlbWVudD8iLCJjcmVhdGVkIjoxNTMy
MzgyOTkzMDEwfSwiWjBVZFhvSmY5VVV5NjV0eSI6eyJkaXNjdX
NzaW9uSWQiOiI1U1lYTzZ1NWJ5WlAyWWtoIiwic3ViIjoiMTAy
MjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ik5haCwgYXMgbG
9uZyBhcyB0aGUgYXN0aGV0aWMgb2YgaXQgaXNuJ3QgYWN0dWFs
bHkgY29waWVkIGZyb20gc29tZXdoZXJlLiIsImNyZWF0ZWQiOj
E1MzI0NDkzODE0ODJ9LCJNOU0zUEtyWlZ5Wk5tUDVqIjp7ImRp
c2N1c3Npb25JZCI6Im9PMkRzT2tqQUxrR3lUT00iLCJzdWIiOi
IxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0IjoiUGVyaGFw
cyB5b3UgY291ZGwgc2F5IHNvbWV0aGluZyBsaWtlIFwiLCBvZi
Bjb3Vyc2UgdGhpcyBkb2Vzbid0IGFkZHJlc3MgdGhlIHF1ZXN0
aW9uIG9mIHdodGhlciB0aGVzZSBzZXF1ZW5jZXMgZm9ybSBHNH
MgKmluIHZpdm8qXCIiLCJjcmVhdGVkIjoxNTMyNDQ5NjQ5MTk0
fX0sImhpc3RvcnkiOls4MDAyNzYyMTMsMTgxOTMzNjg0NSwtMT
I3NzU5NjUyMiwxNjM2OTc2MzE1LC0xMDQ5NjM1ODMxLDEyNzg4
MTMyNTEsLTk4MTg0MjQ0NSwxNDI2NTI5MzAxLDE1OTc5MjYwMD
YsMTQyNjUyOTMwMSwxOTI2MzY2MDMzLDYxNDc3NDQ1MCw1NzY2
NzAzNzksMjEzMTgxNjQ2Nyw5MjMyMDQyNzQsMTU0ODc5NDk4MS
w0OTc2Njg2NDMsMTU0Mzk4NzU1MSwtNzQ2NTQ4MjUwLC05NzI2
NTk4ODRdfQ==
-->
# Discussion

We have presented an in-depth analysis of the distribution of genic G4s in *Arabidopsis thaliana*. As demonstrated by other groups, we found the levels of potential three tetrad G4s in Arabidopsis are low compared to the number of two tetrad PG4s [@Mullen2010]. Our analysis also shows that Arabidopsis genes have high levels of template stranded PG4 sequences, both in the 5' UTR and in the start codon proximal end of the CDS. These levels constitute an enrichment over what would be expected in random sequence with identical GC content and dinucleotide frequencies, and random sequences which could encode the same protein. Furthermore, we show that many G4s are fully or partially hardcoded by protein sequence, suggesting that in some cases protein sequence might be influenced by the G4 forming potential of the DNA. Poly-glycine and poly-proline motifs, which induce hardcoded PG4s on the coding and template strands respectively, were found to account for a large number of Arabidopsis CDS PG4s. We also showed, however, that the percentage of non-hardcoded template strand PG4s is greater at the start codon proximal end of the CDS, indicating codon selection to increase G4 forming potential in this part of the gene body. We hypothesised that G4s may arise in the start of CDSs when 5'UTRs are too short to contain them, however we did not find any correlation between the presence of G4s in the first 100bp of the CDS, and the length of the 5' UTR, suggesting that this is not the case.

Our microarray and RNAseq datasets have shown reproducibly that treatment of Arabidopsis with the G4 stabilising ligand NMM causes widespread changes in gene expression. Genes which contain template stranded two tetrad G4s in their 5' UTRs and CDS regions tended to be downregulated by NMM treatment. Genes with coding strand G4s in the gene body were not as strongly affected, indicating a potential mechanism involving the template strand. Since the template, or non-transcribed, strand is the one scanned by Pol II during transcription, we hypothesised that stabilised G4s on the template strand may form blockages that prevent the elongation of Pol II, whilst coding strand G4s will not. These blockages may result in premature termination of transcription. We found that large numbers of clustered PG4s in a 200bp window were also indicative of reduced expression, suggesting that clusters of stabilised G4s are more problematic to Pol II than single G4s.

We analysed publicly available Pol II ChIP-chip data [@Chodvarupu2012] to show that PG4 dense genes have altered Pol II profiles, with greater Pol II density at the TSS proximal end of genes and lower Pol II density at the TSS distal end, compared to the profiles of most over genes. This could be considered indicative of a slowing or pausing of Pol II over G4 dense regions under normal conditions, i.e. not in the presence of G4 stabilising agents. We did not see higher ratios of reads originating from nascent to mature mRNAs on PG4 dense genes, however, suggesting that this slowing does not naturally result in premature termination or degradation. Despite this, it is well documented that Pol II elongation speed affects co-transcriptional processes such as splicing, and so Pol II slowing by G4s may affect these processes.

Work by Mullen et al. and others previously identified drought stress and hypoxia stress genes as containing greater levels of genic PG4s [@Mullen2010; @Mullen2012; @Yadav2017; @Andorf2014]. They hypothesised that increased cellular potassium concentrations during drought might increase the stability of G4s. G4s might therefore act as a molecular switch to sense drought. Building on this, we found that drought stress downregulated genes contained higher levels of template stranded PG4s in their 5'UTRs. This geneset also overlapped with genes downregulated by NMM. This corroborates the evidence found by Mullen et al. and Yadav et al., and also suggests that the mechanism for drought responsiveness in these genes might be Pol II stalling at G4s.

A family of poly-proline rich protein-coding genes, called the Extensins, were identified as strongly downregulated by NMM. The Extensins were found to be extremely PG4 rich on the template strand, with a number of genes containing as many as 200 different overlapping PG4 registers. The Extensin repeat was analysed by CD Spectroscopy and found to form a G4 structure *in vivo*.  We also found that Extensin genes are also downregulated by Berberine, another G4 binding agent, suggesting that downregulation is indeed caused by G4 stabilisation and not by any off-target effects. Furthermore, downregulation by NMM treatment was not affected by pre-treatment with the translation inhibitor cyclohexamide. Whilst this does not rule out regulation through post-translational effects such as phosphorylation of existing transcription factors, this is good evidence that NMM has a direct effect on Extensin gene expression, and is not simply affecting the expression of a transcription factor.

The Extensins have relatively high levels of spliced reads mapping to them, despite their annotated transcript models having only single exons. These exitronic splice junctions tend to be over the PG4 rich regions of the gene. We hypothesised that these spliced reads may be the result of Pol II slowing at G4 rich regions of the gene, which allows co-transcriptional splicing to occur at weak splice junctions. Since the majority of exitrons were found to be a multiple of three in length, the resultant truncated transcripts would also result in shortened proteins. The diversity in truncated proteins, if regulated correctly, could potentially be beneficial to the plant in changing the properties of the plant cell wall (e.g. the number of glycosylation sites, the flexibility of the protein matrix). This hypothesis is supported by the work of Baumberger et al., who showed that truncated forms of the LRX1 gene, containing fewer Extensin repeats, were able to partially complement an *lrx1* mutant [@Baumberger2001].

Because of the repetitiveness of the Extensin genes, it is possible that some of these spliced reads might result from mapping errors, however we did not find as many unique splice junctions in simulated datasets as in the real RNAseq data. Furthermore, were able to isolate, clone, and sanger sequence some of these spliced transcripts. None of these truncated transcripts appeared to have canonical splice junction sequences, however, possibly indicating that they may be PCR artefacts caused by mispriming. In future we could use Nanopore direct RNA sequencing, which has no PCR step and is able to sequence whole mRNA molecules regardless of repetitiveness, to identify whether the Extensin exitrons are real or artefactual. Regardless of whether the Extensin splice variants are real or not, the levels of splicing do not appear to be responsive to NMM. This could also be considered indicative of the splicing being artefactual, or alternatively that splicing is not linked to G4 formation.

Overall, we have found good evidence that NMM is able to affect the expression of genes with template stranded G4s in their gene bodies. Whilst the mechanism for this action is not fully clear, there is some evidence that impairment of Pol II elongation is involved. This fits with the current literature which suggests that stabilised G4s are able to stall polymerases both *in vitro* and *in vivo* [@Han1999; @Siddiqui-Jain2002; @Dexheimer2006; @Cogoi2006; @Chambers2015; @Kwok2016; @Rodriguez2012]. Furthermore, the altered Pol II profile over PG4 dense genes in the absence of any G4 stabilising ligands indicates that slow Pol II progression at G4s may occur naturally in the plant. One condition during which this may occur is drought stress, when the cellular potassium concentration is greater [@Mullen2010; @Mullen2012; @Yadav2017]. The effects of G4 location (i.e. position and strand) should therefore be considered in future when modelling transcription efficiency and gene expression in Arabidopsis.
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJrOTk2RUJTbWQxdEFYaU11Ijp7In
RleHQiOiJYIiwic3RhcnQiOjQzMTIsImVuZCI6NDMxMn0sIkhT
QmtYM2V0QWtBY2dTM1ciOnsidGV4dCI6IkZ1cnRoZXJtb3JlLC
Bkb3ducmVndWxhdGlvbiBieSBOTU0gdHJlYXRtZW50IHdhcyBu
b3QgYWZmZWN0ZWQgYnkgcHJlLXRyZWF0bWVudOKApiIsInN0YX
J0Ijo0NjQ3LCJlbmQiOjQ3NzF9LCJ4SDhoeGFJVXlyYWFORTRm
Ijp7InN0YXJ0Ijo2NDIwLCJlbmQiOjY1MDgsInRleHQiOiJOb2
5lIG9mIHRoZXNlIHRydW5jYXRlZCB0cmFuc2NyaXB0cyBhcHBl
YXJlZCB0byBoYXZlIGNhbm9uaWNhbCBzcGxpY2UganVuY3Rpb2
7igKYifSwiRTdHblFQRVFaZEpqWmJSUiI6eyJzdGFydCI6Njc5
NSwiZW5kIjo2ODEwLCJ0ZXh0Ijoib3IgYXJ0ZWZhY3R1YWwuIn
19LCJjb21tZW50cyI6eyJZVHA0VWszMVZlU0ZBVkFKIjp7ImRp
c2N1c3Npb25JZCI6Ims5OTZFQlNtZDF0QVhpTXUiLCJzdWIiOi
JnaDo1MzkxNzU4IiwidGV4dCI6IkRvbid0IGZvcmdldCB0aGlz
IiwiY3JlYXRlZCI6MTUzNzA4NTY5NTQxMH0sIldRY0hCdVIxdD
dMY3RmUHMiOnsiZGlzY3Vzc2lvbklkIjoiSFNCa1gzZXRBa0Fj
Z1MzVyIsInN1YiI6ImdoOjUzOTE3NTgiLCJ0ZXh0IjoiRG9lc2
4ndCBydWxlIG91dCBpbmRpcmVjdCAgZWZmZWN0cyB2aWEgcG9z
dC10cmFuc2xhdGlvbmFsIHJlZ3VsYXRpb24uIiwiY3JlYXRlZC
I6MTUzNzA4NTc4MzgxNH0sIlhRTFZBMGpUTjliWTN3bDYiOnsi
ZGlzY3Vzc2lvbklkIjoieEg4aHhhSVV5cmFhTkU0ZiIsInN1Yi
I6ImdoOjUzOTE3NTgiLCJ0ZXh0IjoiSXMgdGhpcyB0cnVlIG9m
IHRoZSBjbG9uZWQgdmFyaWFudHMsIEkgb25seSByZW1lbWJlci
B0aGlzIGFubGF5c2lzIGZvciB0aGUgYXNzZW1ibGVkIHRyYW5z
Y3JpcHRzLCB0aG91Z2ggSSBtaWdodCBiZSB3cm9uZy4iLCJjcm
VhdGVkIjoxNTM3MDg1OTc5MzU0fSwiOE0zblV3aHFxQ2pyUms1
MSI6eyJkaXNjdXNzaW9uSWQiOiJFN0duUVBFUVpkSmpaYlJSIi
wic3ViIjoiZ2g6NTM5MTc1OCIsInRleHQiOiJZb3UgZG9uJ3Qg
bWVudGlvbiBhbnl0aGluZyBhYm91dCB0aGUgbGFjayBvZiByZX
Nwb25zZSB0byBOTU0uIiwiY3JlYXRlZCI6MTUzNzA4NjA4OTgz
M319LCJoaXN0b3J5IjpbMTYwOTcxMjgwLC00NTA0NTcyMjMsMT
A3OTk4MDc5NiwtMTExNTgxMTA5MywtMTk2NzkwMjIzMSw5MDcy
NTQ0OTMsMTE1NzI1ODk4XX0=
-->
# A Recurrent Neural Network to Predict G Quadruplex Structure
\label{chap:g4seeqer}

## Introduction

Because of the dependence of G4 structure largely on sequence information, it is possible to make predictions about the propensity of specific sequences to form G4s using pattern matching analyses. The initial rule which was employed for putative G4 (PG4) detection in the human genome was the Quadparser method [@Huppert2005]. This is a simple regular expression following the pattern $G_{X}N_{1-7}G_{X}N_{1-7}G_{X}N_{1-7}G_{X}$ where $X\geqslant 3$. The Quadparser method was chosen based upon early Circular Dichroism and UV melting data, which suggested that G4s tended to be more stable with three or more tetrads, relatively short loop lengths, and no bulges in tetrads. These fairly stringent rules for formation mean that in general, the Quadparser method is considered to be fairly conservative, and misses a lot of sequences with high G4 forming potential.

More recently there has been a large increase in the number of available methods for predicting G4s [@Bedrat2016; @Hon2017; @Garant2017; @Sahakyan2017]. The contribution from Bedrat et al., named G4Hunter, is a scoring method based on a run length encoding the input sequence. Runs of Gs score positively whilst runs of Cs score negatively. It can be used with a sliding window approach to score an entire genome, with thresholding to identify high scoring PG4s. Whilst this approach is much more tolerant of imperfections which violate the Quadparser method, it is arguably too tolerant, producing many false positives. Furthermore, it does not take into account flanking A and T sequences which may contribute to the stability of the G4, e.g. through reducing the favourability of double stranded DNA.

Some middle ground is required to improve existing methods. The key to finding this middle ground may come from the fields of data science and machine learning. There is a current global interest in machine learning approaches, due to increases in computational power which allow larger datasets to be analysed. These techniques have been adopted by bioinformaticians for making global predictions of functional sequences such as splice sites or transcription factor binding motifs. For example, Mapleson et al. have developed a tool called `portcullis` which is able to filter out spurious splice sites caused by mapping errors from RNAseq data. This is done by training a random forest model on the fly using features derived from high and low confidence splice sites in the dataset. The trained model is then used to classify all other ambiguous splice sites in the dataset as either true or false positives [@Mapleson2017]. Neural networks have also been used to identify sequence motifs. Alipanahi et al. developed a convolutional neural network called `DeepBind` which they applied to identify transcription factor binding motifs and other regulatory sequences, using data from Protein binding microarrays (PBMs), SELEX or ChIP-seq for training [@Alipanahi2015]. Quang et al. were able to improve on this model by using a hybrid convolutional and recurrent neural network [@Quang2016]. This neural network architecture will be explained in detail in \autoref{ssec:results_model_choice}.

Machine Learning approaches have already been brought to bear on the problem of G4 prediction. For example, Garant et al. developed G4RNA screener, a densely connected neural network which is trained on the trinucleotide contents of input sequences [@Garant2017]. This model was trained by using a set of melting temperatures of RNA G4s obtained from a literature search. Currently this database contains only 368 sequences, however. Given the almost inconceivable number of potential G4 forming sequences (there are more than $10^{12}$ sequences which could match the original Quadparser pattern, not including flanking sequences), it is probable that this dataset does not capture all of the variety of possible RNA G4 forming sequences. To do this, a more high throughput method for measuring G4 forming propensity is required.

A new method for sequencing of genomic G4 structures, termed G4Seq, may provide this level of throughput [@Chambers2015]. To create this dataset, Chamber et al. sequenced human genomic DNA in an Illumina sequencing-by-synthesis machine, in the presence or absence of G4 stabilising potassium cations, or G4 binding ligand Pyridostatin. The G4 structures caused stalling of DNA polymerase, resulting in a large number of errors in the resulting read. When the authors then mapped the resultant reads, the mismatch rate which occurred at each position in the genome could be counted to create a map of G4-forming loci. Another high throughput method for sequencing G4s, this time in human mRNAs, was developed by Kwok et al. This protocol utilises the ability of potassium or pyridostatin stabilised RNA G4s to stall reverse transcriptase, resulting in a drop off of reads in the 3'->5' direction in RNAseq data [@Kwok2016].

Data from the G4Seq experiment was leveraged by Sahakyan et al. to build an extreme gradient boosted machine model, named Quadron, to predict G4 formation [@Sahakyan2017]. Quadron initially predicts PG4s using the Quadparser method, but with extended loop lengths of 1-12. Various derived features, such as tetrad number, loop length and mono-, di- and trinucleotide content, are then extracted from these patterns and used to train a model predicting G4Seq mismatch rate. The primary drawback of this method is that the Quadparser method is still used to make initial predictions. This means that Quadron is only able to improve the precision of the existing Quadparser method by rejecting false positives that do not form G4s. Since the Quadparser method is already quite conservative, a lot of potential G4 forming sequences are still missed.

Here we present a new method for G4 prediction, which builds on the work of Bedrat et al. and Sahakyan et al. We use the G4Seq dataset as training data, however a convolutional and recurrent neural network is used to process input sequences directly, meaning fewer prior assumptions are required about what constitutes a G4 forming sequence. Our new method, which we name G4Seeqer, performs better on both the G4Seq datasets and other datasets, including G4s immunoprecipitated from chromatin [@Hansel2016]. We also use transfer learning to apply the model learned from G4Seq data to RNA G4 prediction. Finally, we use G4Seeqer to characterise unnoticed features of the G4Seq dataset, and compare our scores to data from UV melting experiments.

\newpage

## Materials and Methods

### G4Hunter Algorithm

The G4hunter algorithm from Bedrat et al. was reimplemented in Cython (Python superset which can be compiled to C) with some alterations [@Bedrat2016]. Input sequences were run length encoded, and each run of Gs was scored as the square of length of the run, with a maximum score per run of 16. These scores were summed to give a positive strand total score. Runs of Cs were scored equivalently but in a separate negative strand score. Scores were divided by the length of the input sequence to get a normalised score.

### Training Data Preprocessing

The modified G4Hunter method was run on the hg19 genome using a window size of 50bp, a step size of 5 and a threshold of 0.75 to generate a total of 7484506 candidate G4 intervals, which were output in bed format. Intervals were increased in size by 39bp in each direction using `bedtools slop` to introduce flanking sequence information for classification [@Quinlan2010]. Overlapping intervals were filtered to yield the greatest number of non-overlapping intervals. Intervals were weighted by their G4Hunter score, such that the maximum number of high scoring non-overlapping intervals was yielded.

To ground truth score these sequences, `bedtools map` was used to intersect them with the G4Seq dataset [@Chambers2015], which was downloaded from GEO (GSE63874). Bedgraph files of this data contained percentage mismatch scores for each position of the human genome at 15bp resolution. The G4Seq dataset generated in the presence of potassium was chosen as it was deemed more likely to be of biological relevance than the dataset generated in the presence of Pyridostatin, a G4-binding drug.

Intervals files and corresponding mismatch scores were read into Python using `pandas` and histograms of log transformed mismatch scores were plotted using `matplotlib` [@VanRossum1995; @Mckinney2011; @Hunter2007]. The threshold of approximately 3 for separating G4-forming and non-G4 forming sequences was chosen using `scipy` to determine the local minimum in the histogram [@Jones2001]. Joint plots of percentage mismatch score against G4Hunter score were plotted using `seaborn` [@Waskom2014].

Since positive training examples were outweighed by negative ones in this dataset by a factor of 10:1, random under-sampling of negative examples was conducted using `imblanced-learn` to attain a ratio of 2:1 [@Lemaitre2001]. This filtered dataset was shuffled and written to disk in bed format, and `bedtools getfasta` was used to extract sequences for each interval from hg19 [@Quinlan2010]. Sequences were then one hot encoded, i.e. represented as binary matrices of size 128x4, and loaded into HDF5 format for training using `h5py` [@Collette2013]. For training of models on trinucleotide content, trinucleotide content statistics were extracted and loaded into HDF5 format.

### Model Training and Validation

All models were trained in Python using `Keras` with `TensorFlow` backend [@Chollet2018; @Abadi2016]. The trinucleotide Multi-Layer Perceptron (MLP) model contained three hidden layers with 16 units per layer. These were trained using the ADAM optimiser and binary crossentropy loss function, with a dropout rate of 0.2 on all layers. The convolutional portion of G4Seeqer was made up of two convolutional layers with 8 filters and kernel size of 3 and ReLu activation, followed by a maximum pooling layer with step size of 2. This was connected to a bidirectional Long Short Term Memory layer with 8 units and Tanh activation. The final hidden layer was a fully connected layer with 16 units, ReLu activation and a dropout rate of 0.5. G4Seeqer was trained using the RMSprop optimiser and binary crossentropy loss function.

All models were trained on 80% of the training data with 10% used for validation. Training was conducted for a maximum of 30 epochs, but with early stopping when the change in validation loss was less than 0.0005 for more than 3 Epochs. The densely connected network (also known as a Multi Layer Perceptron or MLP) trained on trinucleotide contents converged to this minimum change after 8 Epochs, whilst G4Seeqer converged after 15 Epochs.

Models were validated on 10% of the total data held out for testing purposes. Receiver Operator Characteristic (ROC) and Precision Recall (PR) curves were generated using `scikit-learn` and plotted with `matplotlib` [@Pedregosa2011; @Hunter2007]. ROC/PR curves for G4Hunter are produced using the modified method. For comparison against Quadron, the Quadron source code was downloaded from GitHub and installed [@Sahakyan2017]. Since the flanking sequences required for Quadron are longer than those used for G4Seeqer, test set sequences were increased in size by 50bp in each direction to 228bp. Sequences were extracted using `bedtools getfasta` and run through Quadron [@Quinlan2010]. For intervals which contained multiple Quadron scoring motifs, the highest score was used. For intervals which had no motifs scored by Quadron, a score of zero was assigned.

### BG4 Analysis

NarrowPeak BED files of BG4 ChIP-seq peaks were downloaded from GEO accession GSE76688 [@Hansel2016]. To accommodate Quadron's flanking sequence requirements, the size of the BG4 intervals was increased by 50bp in each direction using `awk` [@Aho1988]. A BG4-negative peak set was generated using `bedtools shuffle` [@Quinlan2010]. Shuffling was performed such that an equal number of simliarly sized intervals were selected that excluded gaps in the genome or BG4-positive peaks. Positive and negative peaks were concatenated and sequences were extracted using `bedtools getfasta` [@Quinlan2010]. Predictions were made on these sequences using G4Seeqer/G4Hunter/Quadron, and the maximum scoring interval per peak was assigned as the overall score of the peak. Where a model did not make any predictions in a peak, it was assigned a score of zero. Receiver Operator Characteristic (ROC, false positive rate plotted against true positive rate) and Precision Recall (PR, precision plotted against recall) curves were generated using `scikit-learn` and plotted with `matplotlib` [@Pedregosa2011; @Hunter2007].

### rG4seq Training Data Preprocessing

To produce training data for rG4seeqer, G4hunter windows were predicted in hg19 using a window size of 50bp and a threshold of 0.75. These were intersected with human exons using `bedtools` to get a set of 186279 putative RNA G4 forming sequences. RNA G4s identified by rG4seq in the presence of potassium [@Kwok2016] were downloaded from GSE77282 (`GSE77282_K_hits.bed.gz`), and intersected with the G4hunter windows to identify RNA G4 positive and negative examples. Of the 3383 identified RNA G4s in the rG4seq dataset, 2811 (83%) overlapped with G4hunter windows. Since the ratio of negative to positive examples was extremely high, negative examples were undersampled to a ratio of 2:1 using `imbalanced-learn` [@Lemaitre2001]. These were then shuffled and written to disk in bed format, and `bedtools getfasta` was used to extract sequences for each interval from hg19 [@Quinlan2010]. Sequences were then one hot encoded and loaded into HDF5 format for training using `h5py` [@Collette2013].

### rG4seq Transfer Learning

Model weights trained on the G4Seq dataset were reloaded in the same architecture for training on the rG4seq dataset, using `Keras` and `tensorflow` [@Chollet2018; @Abadi2016]. Weights from the initial convolutional layers were fixed (i.e. made untrainable), and only LSTM weights and final dense weights were trained. The model was trained on 6014 samples (80%), validated on 752 samples (10%), and tested on 752 samples (10%). Training was conducted as for G4Seeqer, but for a maximum of 200 epochs, with early stopping after 25 epochs if validation loss did not improve. Initial learning rate was set at 0.001 but reduced by a factor of 1/3 after 15 epochs when no reduction in validation loss was seen.

### Comparison to G4RNA Screener

Supplementary data containing 368 RNA sequences, their experimentally determined G4 forming status, and their predicted score from G4RNA Screener were downloaded [@Garant2017]. Sequences greater than 128bp were filtered to give a total of 347 sequences. Average rG4Seeqer and G4Seeqer predictions for these sequences were made by generating 1000 randomly paddings for each sequence, one hot encoding, and performing forward pass through the network. G4RNA Screener scores used in the ROC curve were taken directly from the supplemental material [@Garant2017].

### Mutation Mapping analysis

Mutation mapping was applied to cell type independent human promoter regions from the ENSEMBL regulatory build [@Zerbino2015], which was originally generated using `ChromHMM` [@Ernst2017]. Promoter sequences were extracted from hg38. Mutation mapping was implemented as in Alipanahi et al. 2015: candidate sequences were edited at each position to each nucleotide, and the resulting sequences were scored for G4 formation using G4Seeqer [@Alipanahi2015]. Heatmaps were generated in Python using `seaborn` [@Waskom2014].

### G-Triplex and Hairpin analyses

G-triplex motifs were predicted in the hg19 genome using an in house script, using the pattern $G_XN_{1-4}G_XN_{1-4}G_X$ where $3 \leqslant X \leqslant 6$. Candidate triplexes which overlapped with or were contained within Quadparser motifs (pattern $G_{X}N_{1-7}G_{X}N_{1-7}G_{X}N_{1-7}G_{X}$ where $3 \leqslant X \leqslant 6$) were subtracted to produce only triplex motifs which could not form "classical" G4s. To assign mismatch scores to these sequences, they were increased in size by 50bp in each direction using `bedtools slop` [@Quinlan2010], and mismatch scores from the G4Seq dataset were mapped using `bedtools map`. Distances to next G-run were measured in python using `pyfaidx` [@Shirley2015].

G-hairpin motifs were predicted in the hg19 genome using the same script, and the pattern $G_XN_{1-4}G_X$ where $4 \leqslant X \leqslant 6$. Candidate hairpins which overlapped with or were contained within Quadparser motifs or G-triplex motifs were filtered. Intervals were increased in size by 50bp in each direction using `bedtools slop` [@Quinlan2010], and mismatch scores from the G4Seq dataset were mapped using `bedtools map`. Distances to next G-hairpin were measured using `bedtools closest`. Triplex and hairpin histograms and boxplots were generated in Python using `matplotlib` and `seaborn` [@Hunter2007; @Waskom2014].

### Median model score experiments

For G4Seeqer scoring of experimentally validated G4s from Guédin et al 2010., sequences were recreated from the information in Table 1, and padded using uniformly sampled sequences to 128bp in length [@Guedin2010]. Left and right padding lengths were also varied at random. 1000 randomly padded sequences were generated and scored per input sequence. Scatter plots of UV melting temperature vs. median G4Seeqer score were produced using `matplotlib` [@Hunter2007]. Errorbars are 68% confidence intervals calculated from the variation in score for the same sequence with different random padding.

Synthetic sequences used for loop length and G-register experiments were generated by combining uniformly sampled padding and loops with GGG trinucleotides to create 3 tetrad Quadparser conforming sequences. Each sample contained 5000 randomly generated sequences. Left and right padding sizes and nucleotide contents for loops were varied within random samples. Loop nucleotide contents were also varied. For G-register experiments, the extra G per run was randomly assigned to either the left of right side of the G-run. Loop lengths of 3 were used. Line plots for loop length experiments were generated using `matplotlib`. Errorbars are 68% confidence intervals showing the variation in G4Seeqer score caused by left and right padding size and padding/loop nucleotide content, when holding loop lengths constant [@Hunter2007]. Boxplots were generated using `seaborn` [@Waskom2014].

### G-register Experiments

For G-register mismatch experiments, 3 tetrad Quadparser conforming G4s from the hg19 genome were identified and the corresponding mismatch score was extracted using `bedtools map` [@Quinlan2010]. The number of tetrads with G-register was counted. Boxplots were generated using `seaborn` [@Waskom2014].

### Human and Mouse G4 Subpopulation Analyses

Example Human and Mouse G4 populations were predicted in hg19 and mm10 using the Quadparser pattern with loop lengths of 3. All possible G4s conforming to this pattern were generated using python `itertools`. Venn diagrams were generated using `matplotlib_venn`. P-values were produced using hypergeometric tests. Dinucleotide complexity was defined as the total number of unique dinucleotides contained in the motif. Histograms and kernel density estimate plots were produced using `seaborn` [@Waskom2014]. For visualisation of PG4 distributions, a sample of 50000 motifs were randomly selected from the total population of all possible Quadparser motifs with loop lengths of three. These were transformed into two components using UMAP dimensionality reduction, with Hamming distance as the distance metric [@McInnes2018]. Sequences which appear in the human and mouse genome were extracted from the sample. 2D Kernel Density Estimate plots for the full sample, and hg19 and mm10 subsets, were generated using `seaborn`.

\newpage

## Results and Discussion

### Candidate G4 proposal

One major drawback to any machine learning method for G4 prediction over existing regular expression or pattern based methods is the relative expense of computation. It is therefore not sensible to train and classify a neural network model on all possible input sequences from a genome. Instead we decided to use an existing method, the G4hunter algorithm proposed by Bedrat et al. [@Bedrat2016], to produce candidate regions which could then be labelled as true positive G4s or non-G4s, and used as input for model training.

We reimplemented the G4hunter algorithm with some minor modifications [@Bedrat2016]. Bedrat et al's method run length encodes the sequence of interest and scored G-runs as the square of the run length, and C-runs as the negative of the square of the run length. These scores are then summed to give an overall score for the sequence. This method was chosen as it was assumed, based on earlier work, that a high G-content on both strands would make the G-Quadruplex unfavourable compared to double stranded DNA. Since we wished to make as few assumptions as possible, and given that the G4seq dataset stems from G4 formation in *in vitro* single stranded DNA, we altered the method to produce two scores. G runs score positively on the positive strand and C runs score positively on the negative strand. This means the sequence d(GGGCCC) would yield a high score on both strands rather than a single score of zero.

To produce candidate regions for model training, we ran the modified G4hunter method on the human genome (hg19) using a window size of 50bp, a step size of 5bp and a threshold of 0.75. Unstringent values were chosen to produce a high recall, i.e. capture as many true positive G4 structures as possible. These settings produced intervals which overlapped with all PG4 sequences predicted by the Quadparser method using maximum loop lengths of 12. It also produced significantly more sequences that did not conform to the Quadparser method, some of which are likely to form G4 structures.

### Training Data preprocessing

To created training sequences, the 50bp candidate regions from the G4hunter method were increased in size by 39bp in each direction to produce intervals of 128bp in length, since previous work has suggested that flanking regions are an important determinant of G4 stability. Clusters of overlapping intervals were filtered to produce only the interval with the highest G4Hunter score (in cases of ties, a random highest scoring interval was selected). This produced a total of 6,237,943 candidate intervals. Each candidate interval was then scored by mapping the value of the maximum scoring overlapping window from the G4seq dataset [@Chambers2015], which contains percentage mismatches (%mm) from sequencing in the presence of potassium vs. absence of potassium, in 12bp windows. Regions of high %mm on the positive strand indicate a G4 structure on the negative strand, and vice versa.

Plotting the distribution of the log of the %mm scores produced a bimodal distribution with a peak around 1 (corresponding to 2-3% mismatch) and another around 3.5 (corresponding to approximately 30% mismatch) (Fig. \ref{training_data}a). We determined the local minimum in the histogram between the two peaks to be approximately 3 (around 20% mismatch), therefore we used this value to split the data into G4 positive and G4 negative subsets. This yielded more than 10 times as many G4 negative sequences than positive, however (5,809,719 negative to 428,224 positive). Since maintaining such an imbalance in the training data would produce a poor classifier, we undersampled the G4 negative class to a ratio of 2:1, yielding a total of 1,284,672 sequences.

\newpage

![**Mismatch scores of candidate sequences identified by G4Hunter:** **a)** Histogram of log percentage mismatch score from the G4Seq dataset, for the 50bp sequences identified by G4Hunter (threshold of 0.75). Dashed line shows the threshold chosen to delimit G4 positive and G4 negative sequences. This corresponded to around 20% mismatch score. **b)** Joint plot of log mismatch score against G4Hunter score for 10000 randomly sampled sequences. Orange line shows lowess curve fit. \label{training_data}](figures/training_data_dist.svg)

\newpage

### Model selection and training

\label{ssec:results_model_choice}

Previously published G quadruplex prediction methods which utilise machine learning techniques have used derived features such as trinucleotide content to feed to models [@Garant2017; @Sahakyan2017]. These features result in the loss of some spatial information about the sequence, however. For example, the sequence `GGTGGTGGTGGGGGG` has the same trinucleotide composition as `GGGTGGGTGGGTGGG`, but is unlikely to have equivalent G4 forming propensity. Furthermore, Quadron derived features require input sequences to conform to the QuadParser regular expression [@Huppert2005], meaning that Quadron is only able to improve the precision of the QuadParser method, and not the recall. We opted for a neural network involving convolutional layers (those often used for image classification) that could make predictions directly from the sequence itself, without any derived features whatsoever. This allows us to make no assumptions about potential G4 patterns in the dataset. The overall architecture selected was a convolutional-recurrent neural network (Fig. \ref{architecture}), which has previously been used to identify regulatory motifs in DNA [@Quang2016]. The first layers of the architecture consists of two one dimensional convolutional layers with a kernel size of 3. Convolutional layers contain a number of separate filters which are convolved with the data to generate an output function. The values (aka weights) in these filters are determined by training. Each filter captures different local features in the sequence.  A maximum pooling layer then reduces the size of the output feature space by half. These features are then fed to a bidirectional Long Short Term Memory (LSTM) layer. LSTMs are recurrent layers, meaning that the nodes of the layer form a linear directed graph over the input sequence. The internal state of the node is dependent on the output from the previous node, meaning that the model is able to "remember" features it has previously seen in a  which is able to learn and recognise long distance relationships between features in the sequence. The model outputs a single value between zero and one of the probability of G4 formation. The dataset was split into three for training and testing: 1332565 sequences (~80%) were used for training, 166571 (~10%) for in-training model validation, and 166576 (~10%) for post-training model testing.

\newpage

![**G4Seeqer architecture** Adapted from Quang & Xie 2016. **a)** Sequences are one hot encoded to produce a matrix which can be processed by the neural network. **b)** Input matrices are passed through a convolutional layer. Each layer contains 8 filters which are trained to recognise local patterns on the scale of 3-6bp in size. **c)** Convolutional features are passed through a bidirectional Long Short Term Memory (LSTM) layer. This layer recognises long distance interactions between features which might combine to produce G4s. **d)** Finally, features are passed through a fully connected layer. Output from the model is a single probability of whether the sequence forms a G4. \label{architecture}](figures/architecture.svg)

\newpage

### Comparison to existing methods

We benchmarked our technique (hereafter referred to as G4Seeqer) using the 10% of the data reserved for testing. The model was compared to our modified G4Hunter method, as well as a multi-layer perceptron model trained on trinucleotide frequencies derived from the same dataset as was used to train G4Seeqer. This model allows us to compare the methodology of G4RNA Screener [@Garant2017] to our own method, since G4RNA Screener was originally trained on a database of RNA G Quadruplexes, and may not perform as well on a dataset derived from DNA. The performance of the methods were calculated using the Receiver Operating Characteristic area under curve (ROC AUC). We found that neither G4Hunter [@Bedrat2016] nor the G4RNA-like method performed as well as G4Seeqer on the test dataset (AUC G4Hunter 0.82, G4RNA Screener-like 0.90, G4Seeqer 0.94) (Fig. \ref{roc}a-b). This is likely due to the loss of sequence spatial information in the former methods.

To benchmark our method against the other G4Seq trained machine learning G4 prediction package, Quadron, we downloaded and installed the Quadron source code [@Sahakyan2017]. Quadron was used to score sequences from the held out test set and compared to the performance of other methods. Quadron requires larger flanking regions of 50bp, so test set intervals and sequences were increased in length by 50bp in both directions to produce test sequences of length 228. This was deemed the best way to compare the two methods, however was still not ideal since Quadron may have been trained on some or all of the test sequences. Regions which had predictions associated with them by G4Hunter and G4Seeqer but had no associated predictions from Quadron (due to not conforming to the Quadparser regular expression) were given a Quadron score of zero. The ROC curve of Quadron (AUC 0.7) ((Fig. \ref{roc}a-b)) had interesting properties: the model is capable of producing a true positive rate of around 30% with a false positive rate of less than 1%, however shortly beyond this point in the ROC curve the curve becomes linear. This is because Quadron is only capable of scoring sequences which conform to the Quadparser method, with all other sequences being scored zero. Since these sequences only account for a 25% of all the potential PG4 forming sequences in the test set, the model does not perform well on the full dataset. Because it makes no assumptions about the input sequence, G4Seeqer is able to make accurate predictions for all forms of PG4 in the test data, as reflected in the ROC curve.

We were interested in how G4Seeqer compared against Quadron on only sequences which conformed to the Quadparser motif. We therefore filtered our dataset for intervals on which Quadron had made a prediction, and replotted the ROC curves and Precision Recall curves for the filtered data (Fig. \ref{roc}c-d). As expected the AUC for Quadron was much better on this filtered set (AUC 0.93), however it was still outperformed by G4Seeqer (AUC 0.94), suggesting that G4Seeqer captures the same or more explanatory information directly from the input sequence.

\newpage

![**Validation curves for G4Seeqer method** **a)** Receiver Operator Characteristic (ROC) curve showing the performance of G4Seeqer, Multi-Layer Perceptron (MLP) trained on trinucleotide contents, Quadron, a Gradient Boosted Machine model (Sahakyan et al. 2017) and the G4Hunter method (Bedrat et al. 2016), on a held out test set of the G4Seq dataset (10% of total dataset). **b)** Precision-recall curves showing the performance of G4Seeqer, trinucleotide MLP, Quadron, and G4Hunter on the same dataset. **c)** ROC curve and **d)** Precision Recall curve showing the performance on sequences from the test set conforming to the Quadparser motif. \label{roc}](figures/test_set_roc_pr.svg)

\newpage

### BG4 ChIP-seq data evaluation

Further model validation was performed on G4s experimentally validated by an entirely different technique, namely G4-chromatin immunoprecipitation (BG4) [@Hansel2016]. The BG4 dataset is arguably more biologically relevant than G4seq since G4 structures are not induced by addition of potassium, and are captured from native chromatin. BG4 peak intervals were shuffled to produce a set of G4 negative sequences, and then the performance of the models was evaluated on the real and shuffled peaks. For each BG4 interval, the highest scoring overlapping prediction for each model was assigned. Any intervals with no overlapping predictions were scored zero for that model. We tested the G4Seeqer, G4Hunter and Quadron methods. As with the G4Seq test dataset, we found that Quadron performed reasonably for Quadparser conforming BG4 peaks, but was unable to identify most of the true positive BG4 peaks due to its restriction to the pattern. G4Seeqer performed better on all BG4 peaks, with an AUC of 0.71, however was only marginally better than the G4Hunter technique (AUC 0.7) (Fig. \ref{bg4}). These results suggest that the information within the G4Seq dataset, when captured by a suitable model, is predictive of G4s in an *in vivo* setting.

\newpage

![**Detection by BG4 peak sequences using G4Seeqer** **a)** Receiver Operator Characteristic (ROC) curve showing the performance of G4Seeqer, Quadron, and the G4Hunter method, on BG4 ChIP-seq peaks and randomly shuffled negative sequences. \label{bg4}](figures/bg4_roc.svg)

\newpage

### Transfer learning on RNA G4Seq (rG4seq) Dataset

PG4s with the same sequence are likely to have slightly different G4 forming potentials in DNA and RNA, due to the chemical differences in these molecules. The sugars which make up the backbone of RNA are riboses, which have an extra 2' hydroxyl group compared to the deoxyribose found in DNA. This extra hydroxyl group is thought to have a number of implications for G4 formation: it increases the number of backbone hydrogen bonds in the G4, increasing its enthalpic favourability and its entropic favourability (by reducing the number of coordinated water molecules) [@Collie2010]. Furthermore, the 2' hydroxyl introduces steric constraints which make parallel RNA G4s much more favourable than anti-parallel ones [@Collie2010]. Given these differences, it is likely therefore that a model specifically trained on DNA G4 sequences will not perform optimally on RNA G4s.

To address this issue, we decided to retrain G4Seeqer using the rG4seq dataset produced by Kwok et al. 2016, to create an rG4Seeqer model [@Kwok2016]. Data was prepared similarly to the data for G4Seeqer: candidate regions were selected from human exonic sequences using G4hunter with window size of 50bp and a threshold of 0.75. This yielded a set of 186279 putative RNA G4 forming sequences, which were increased by 39bp in each direction to get flanking sequences. These were then intersected with rG4seq hits collected under potassium stabilising conditions (Kwok et al. 2016). Of the 3383 identified RNA G4s in the rG4seq dataset, 2811 (83%) overlapped with G4hunter windows. rG4seq negative examples were undersampled with a ratio of 2:1 to yield 7518 training samples. 80% of these were used for training, with 10% for validation and 10% held out for testing.

Because of the significantly smaller size of the rG4seq derived training set, we found that the method for training which yielded most optimal results was transfer learning from the G4Seeqer model. Weights of the initial convolutional feature extraction layers were therefore held constant, and only the weights of the LSTM layers (which find long range interactions) and dense output layers were retrained.

Testing was first conducted on the held out set of 752 sequences using G4Hunter, G4Seeqer and the newly trained rG4Seeqer. G4Seeqer significantly outperformed G4Hunter (AUC 0.9 vs 0.83 respectively), suggesting that the information extracted from the G4Seq dataset is applicable to the rG4seq dataset (Figure \ref{rG4seeqer_test}). rG4Seeqer outperformed both methods, however (AUC 0.95), demonstrating that domain specific information is better for predicting RNA G4s in the rG4seq dataset (Figure \ref{rG4seeqer_test}).

\newpage

![**Validation curves for rG4Seeqer method** **a)** Receiver Operator Characteristic (ROC) curves showing the performance of rG4Seeqer, G4Seeqer and the G4Hunter method (Bedrat et al. 2016), on a held out test set of the rG4Seq dataset (10% of total dataset). **b)** Precision-recall curves showing the performance of rG4Seeqer, G4Seeqer, and G4Hunter on the same dataset. \label{rG4seeqer_test}](figures/rg4seq_test_set_roc_pr.svg)

\newpage

We sought to test rG4Seeqer on G4s identified by a variety of physical methods, using the set of RNA G4s curated by Garant et al. for their model, G4RNA Screener [@Garant2017]. We used 347 sequences from this dataset, of which 169 sequences were G4 positive and 178 were G4 negative. G4Seeqer and rG4Seeqer predictions were calculated by padding with random sequences to a length of 128bp before one hot encoding. This was conducted 1000 times for each sequence and the mean score was taken. G4RNA screener scores were taken directly from the supplemental information of Garant et al. 2017. G4RNA screener was found to perform best on the dataset (AUC 0.91), perhaps unsurprisingly since it was trained directly on the sequences. Perhaps more importantly, rG4Seeqer significantly outperformed G4Seeqer on the dataset (AUC 0.89 vs 0.82), showing that the rG4seq trained model generalises better to RNA G4 sequences than the G4seq trained model.

\newpage

![**Validation of rG4Seeqer on *in vitro* experimentally categorised RNA sequences** **a)** Receiver Operator Characteristic (ROC) curves showing the performance of rG4Seeqer, G4Seeqer and the G4RNA Screener method (Garant et al. 2017), on the G4RNA dataset curated by Garant et al. **b)** Precision-recall curves showing the performance of rG4Seeqer, G4Seeqer, and G4RNA Screener on the same dataset. \label{rG4seeqer_test}](figures/g4rna_roc_pr.svg)

\newpage

### Interpreting G4Seeqer output using Mutation Mapping

One common complaint about neural network techniques is that the complexity of the models they produce make them "black boxes" which are impossible for humans to understand or extract useful knowledge or rules from. It is possible, however, to visualise some of the output of a neural network through various means. One commonly used method for interpretation is the "saliency" of the network, which can produce heatmaps showing the attention of the network to specific regions of the input image or sequence. This can be used to determine the important aspects of the input in classification. For biological sequences, previous studies have used similar methods, called "Mutation Maps", to analyse the importance of individual nucleotide positions on convolutional neural network model predictions [@Alipanahi2015; @Quang2016]. Simply, the importance of each particular nucleotide is evaluated by replacing it with each of the other three bases, and calculating the change in model score. This is then used to build a heatmap which can be used to visualise the importance of each position.

Previous studies have highlighted a possible role for G4 forming sequences in promoters, with G4 formation tending to have a positive effect on expression, particularly in proto-oncogenes [@Eddy2006; @Hansel2016]. We therefore decided to use the mutation map approach to characterise PG4s in promoter regions extracted from the ENSEMBL regulatory build [@Zerbino2015]. Promoter sequences were screened using G4Seeqer and single base mutation maps were created for each candidate PG4, including those regions where the neural network score was low. Unsurprisingly, all of the most deleterious single base substitutions (causing a score reduction of more than 0.9) predicted by G4Seeqer mutation mapping were G->H changes.

We identified PG4 sequences scoring more than 0.9 for which a single G->H change resulted in a reduction of as much as 0.9 in score (i.e. switched the score from strongly PG4 positive to strongly PG4 negative) (Fig. \ref{triplex}a). Analysis of the mutation maps for these sequences showed that the majority of contained regions containing three G-runs with short connecting loops. These tended to have a long final loop to the next homopolymeric G-run, or no final G-run within the window size. Any G->H mutation in these G-dense regions strongly affected the predicted G4 forming ability. Recent work by Hou et al. has shown that formation of G4 structures in human telomeric sequences occurs via a stable G-triplex structure [@Hou2017]. We believe these results suggest that the G4seq dataset is either capturing mismatches caused by G-triplex structures, or by G-quadruplexes formed from short range G-triplex interaction with more long range single G-runs. To further illustrate this we predicted all short looped (1-4bp) G-triplexes in the human genome which did not overlap with a Quadparser predicted PG4 (loop lengths 1-12bp). We found that 35% of these were associated with a %mm score greater than 20%, suggesting the formation of some secondary structure (Fig. \ref{triplex}b). To see if medium range G-run interactions might stabilise these, we then measured, for each G-triplex, the distance to the next run of at least three Gs on the same strand, and compared this to the %mm score. We found a weak negative correlation between distance and %mm score (Spearmans rho -0.2) for G triplexes with a G-run less than 100bp away, suggesting that these longer range interactions do occur but become weaker with distance (Fig. \ref{triplex}c). This could be due to a reduction in G4 stability with loop length, but could equally be explained by a reduction in the likelihood of the next G-run being contained in the same sequenced fragment. No clear difference was observed in the correlation of %mm score with upstream or downstream G-runs (Spearmans rho -0.17 and -0.15 respectively) suggesting there is no preference for the long loop region to be at the 5' or 3' of the G triplex.

\newpage

![**Identification of G-triplex structures by G4Seeqer Mutation Maps** **a)** Mutation map showing the a high scoring (0.99) G4Seeqer motif which may form a G-triplex. Mutation of any base in the central G-run of the motif is sufficient to reduce the score by up to 80%. **b)** Histogram of log percentage mismatch score for motifs conforming to a G-triplex like pattern. The bimodal distribution suggests that many of these motifs form structures which disrupt polymerase in the presence of potassium. **c)** Boxplot showing the relationship between %mm score and distance to next G-run in G-triplex structures. The negative correlation suggests G-triplexes might recruit distant G-runs to form G4s. \label{triplex}](figures/g_triplex.svg)

\newpage

We also noted that G4Seeqer was positively labelling certain sequences which contained only two G-runs, usually of greater than 4 bases in length (Fig. \ref{hairpin}a). These scores were also very sensitive to G->H mutations. Based on the folding dynamics work by Hou et al, we hypothesised that these sequences might form G-hairpins, which could associate with other nearby hairpins to form G4 structures [@Hou2017]. To test this, G-hairpins with G-run lengths greater than four were predicted and filtered to remove any overlaps with predicted QuadParser G4s (loop length 1-12) or G-triplexes (loop length 1-4). 27% of these sequences had %mm scores greater than 20% (Fig. \ref{hairpin}b). For each predicted G-hairpin, the distance to the nearest G-hairpin was calculated and correlated with %mm. Again, distance was found to correlate negatively with %mm score (Spearmans rho -0.18), suggesting that these hairpins may associate with each other to form G4s (Fig. \ref{hairpin}c).

\newpage

![**Identification of G-hairpin structures by G4Seeqer Mutation Maps** **a)** Mutation map showing the a high scoring (0.99) G4Seeqer motif which may form a G-hairpin. Mutation of any base in the core motif is sufficient to reduce the score by up to 80%. **b)** Histogram of log percentage mismatch score for motifs conforming to a G-hairpin like pattern. The bimodal distribution suggests that many of these motifs form structures which disrupt polymerase in the presence of potassium. **c)** Boxplot showing the relationship between %mm score and distance to next G-hairpin for G-hairpin structures. The negative correlation suggests G-hairpins might interact with other relatively distant hairpins to form G4s. \label{hairpin}](figures/g_hairpin.svg)

\newpage

In order to study how interactions between pairs mutations might affect predicted G4 stability, we developed a pairwise mutation map, in which pairs of Gs in each sequence were combinatorially mutated to Ts. We then analysed the resultant mutation maps to identify pairs of G->T transversions which interact to reduce predicted G4 forming potential more strongly than each individual mutation. Perhaps unsurprisingly, we found that in sequences which had more than four G-runs, or had G4s containing features which might form G-triplexes or G-hairpins, mutations which disrupted peripheral G-runs did not have a strong effect on predicted stability. Combinations of mutations which disrupt multiple G-runs had a much stronger effect on stability, however.

### Loop Length and G4 stability

To determine whether G4Seeqer probability scores supported previous work on the relationship between G4 loop length and stability, we first downloaded UV melting temperatures for three tetrad G4 sequences from Table 1 of Guédin et al. 2010. The majority of these G4 sequences contain only runs of Gs and Ts. In each experiment, one or two of the three T loops were held at a constant length of either 1 or 3, and the other loops varied from 1 up to 15 bp in length [@Guedin2010]. For each sequence, we produced 1000 sequences padded to 128bp (the input length of G4Seeqer) using randomly generated bases, and used G4Seeqer to predict the stability. We found a very strong correlation (Spearmans rho 0.93, p = 1.1e-35) between empirically determined melting temperature in potassium, and G4Seeqer score (Fig. \ref{tm}), suggesting that G4Seeqer is successfully capturing information about G4 structure which is transferable between conditions. The midpoint of the curve appears to suggest that a melting temperature of around 65 degrees Celsius is required for significant mismatches to occur in the G4Seq dataset. We also noted that G4seeqer output was more variable for sequences with lower melting temperatures, suggesting that sequence context may be more important for these G4s.

\newpage

![**G4Seeqer scores correlate with experimentally determined melting temperatures** Scatter plot showing median G4Seeqer scores vs. UV melting temperature for sequences from Guédin et al. 2010. Error bars are 68% confidence intervals generated from 1000 iterations of prediction with randomly generated flanking sequences. \label{tm}](figures/tm_vs_score.svg)

\newpage

We next performed a similar *in silico* experiment to that of Guédin et al., whereby we generated random QuadParser conforming G4 sequences with two loop lengths held at 1bp and the third varied from 1-60bp. Loop regions were constructed by randomly selecting from A C or T. The G4Seeqer score for these sequences was then generated. Unsurprisingly, we found that G4Seeqer score reduced with increasing loop length. This effect was strongest for the central loop of the G4, presumably because varying this loop has a greater effect on the ability to form stable G-triplex intermediates (Fig. \ref{loop_len}a). We then set the length of the non-varying loops to 3bp and re-ran the analysis for third loop length 1-60bp. For these PG4s, we found that the probability of G4 prediction was much more sensitive to longer loop lengths (Fig \ref{loop_len}c). There was also no longer a strong difference in prediction when varying the central loop, compared to either loops 1 or 3, possibly suggesting that triplex formation in these sequences is less common.

\newpage

![**Effect of increasing loop length on G4Seeqer score** **a)** Effect of loop length on predicted G4 stability when other loops are held at a constant length of **a)** 1bp, **b)** 2bp or **c)** 3bp. Median value and 68% confidence intervals are produced using 5000 randomly generated sequences (including flanking regions and loop contents) for each pattern analysed. \label{loop_len}](figures/loop_length.svg)

\newpage

### Effect of G-register on G4 stability

Work by Harkness and Mittermaier has indicated that extra Guanines in some G-runs of a G4 forming sequence might increase the G4 forming potential of the sequence, by allowing exchange between different G4 conformations [@Harkness2016]. They termed this effect G-register. We analysed the G4seq dataset and the G4Seeqer model outputs to determine whether there was evidence of a relationship between stability and G-register. Firstly, for all PG4s in the human genome conforming to the three tetrad Quadparser motif, we counted the number of G-runs which contained four Gs rather than three. These motifs were then intersected with the G4seq dataset to identify the %mm score. We found that G4seq motifs with greater G-register indeed tend to have higher mismatch scores (Fig. \ref{register}b). To test whether G4seeqer had successfully identified this pattern, we then randomly generated three tetrad PG4 sequences with loop lengths of 3bp, and introduced addition guanines to increase the G-register. Higher G-register strongly increased the G4Seeqer score of the sequence, showing that the model has successfully learned this feature of G4s (Fig. \ref{register}a).

\newpage

![**Scoring of G-register effects by G4Seeqer** **a)** Boxplot showing G4Seeqer scores for randomly generated 3 tetrad Quadparser conforming sequences with zero to three additional Guanines per run, referred to by Harkness and Mittermaier as G-register. **b)** Boxplot showing relationship between G-register and %mm score in the G4Seq dataset for Quadparser conforming G4s. \label{register}](figures/g_register.svg)

\newpage

### Applicability of the model to other genomes

Whilst the Human genome contains a large number of G4 forming sequences, this is not even close to saturating the population of all potential G4s. Indeed, simply considering the Quadparser motif with loop lengths up to 12, there are $1.1\times10^{22}$ different conforming sequences, many orders of magnitude more than there are bases in the human genome. It is probable that the human genome contains a biased subpopulation of all G4s. This might mean that the G4Seeqer model does not generalise well to other genomes which contain different subpopulations. As an example, we analysed 3 tetrad Quadparser motifs from hg19, where all three loops were of length 3. There are $4^9$ (262,144) possible sequences fitting this description, of which only 1.4% (3,748) appear in hg19 (Fig. \ref{g4_space}a). The motifs that did appear tended to be those with lower complexity, as measured by the number of distinct dinucleotides in the sequence, than the total possible population (Fig. \ref{g4_space}b). The PG4 space of the *M. musculus* genome was also measured, and found to contain 1.7% (4330) of all possible 3 tetrad PG4s with loop length 3. There was a strong overlap of 27% (p < 2.2e-308) between sequences in the human and mouse genomes, however, suggesting that at least for this pattern, the PG4 populations of these genomes are comparable. The more complex the PG4 motifs become, however, the more likely it is that these subpopulations will be very different. There are clear differences in dinucleotide content between different genomes, which are often a result of differences in amino acid composition of proteins, or other environmental factors such as temperature. For genomes whose last common ancestor with *Homo sapiens* was longer ago, this divergence may be much greater. These systemic differences between may result in patterns to which the model has not previously been exposed, and reduce the performance of the model. G4Seeqer, or any other models which are trained on sequences from a single genome, should therefore be used with caution on others.

\newpage

![**Human and mouse genomes contain different G4 subpopulations** **a)** Venn diagram showing overlap of 3 tetrad Quadparser motifs populations with loop lengths of 3bp in the human (hg19) and mouse (mm10) genomes, compared to all possible sequences. **b)** Histogram and kernel density estimate of dinucleotide complexity for human, mouse and all possible 3 tetrad Quadparser motifs with loop length of 3bp. **c)** 2D Kernel density estimate plot showing distribution of all possible  tetrad Quadparser motifs (left), those found in the human genome (centre) and those found in the mouse genome. Dimensionality reduction was conducted using UMAP with Hamming distance as the metric. \label{g4_space}](figures/g4_space.svg)

\newpage

### Conclusion

We present G4Seeqer, the first Convolutional/Recurrent Neural Network model for prediction of G Quadruplex forming structures. G4Seeqer is implemented in Python, using Cython for speed-up of G4Hunter candidate region proposal, and `Keras` with `Tensorflow` backend for neural network prediction [@Chollet2018; @Abadi2016]. Weights have been trained on the G4Seq dataset [@Chambers2015] and transfer learned to the rG4Seq dataset [@Kwok2016], to produce models tailored for DNA and RNA G4s, respectively. It is able to process the whole human genome in approximately 1 hour on a 8 core i7 desktop computer with 16GB RAM. Because G4Seeqer is trained directly upon sequences from the human genome, rather than on derived sequence features, it is able to identify patterns in the G4Seq dataset that have not previously been reported, as well as removing false positive sequences which are flagged by pattern matching techniques. This greatly improves the accuracy of the model on various *in vitro* and *in vivo* datasets, from stabilities determined by UV melting to genomic regions identified by BG4 ChIP-seq [@Hansel2016].

\newpage
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJHdWtoM3JqY3h3dWhXSHRpIjp7In
RleHQiOiJwcm9kdWNpbmcgbWFueSBmYWxzZSBwb3NpdGl2ZXMi
LCJzdGFydCI6MTU2MCwiZW5kIjoxNTkwfSwiV0RDM1lpQXdsak
VxVERkVyI6eyJ0ZXh0IjoiaG93ZXZlciBhIGNvbnZvbHV0aW9u
YWwgYW5kIHJlY3VycmVudCBuZXVyYWwgbmV0d29yayIsInN0YX
J0Ijo2MDI2LCJlbmQiOjYwNzh9LCJtanZjN1ZHRUg3MTc3VG5t
Ijp7InRleHQiOiJUaGUgY29udm9sdXRpb25hbCBwb3J0aW9uIG
9mIEc0U2VlcWVyIHdhcyBtYWRlIHVwIG9mIHR3byBjb252b2x1
dGlvbmFsIGxheWVycyB34oCmIiwic3RhcnQiOjk4NzYsImVuZC
I6MTAxNjh9LCJ1MGxIYjc5V2FuSGdMUkxuIjp7InRleHQiOiJN
TFAiLCJzdGFydCI6MTA3MTYsImVuZCI6MTA3MTZ9LCJyOFZpZU
83Z1BDOXVITWIwIjp7InRleHQiOiJpbnRlcnZhbCBwZXIiLCJz
dGFydCI6MTIzODYsImVuZCI6MTIzOTh9LCJpWUZpeFpPZjBiS3
JyUmROIjp7InRleHQiOiJPZiB0aGUgMzM4MyBpZGVudGlmaWVk
IFJOQSBHNHMgaW4gdGhlIHJHNHNlcSBkYXRhc2V0LCAyODExIC
g4MyUpIG92ZXJsYXBwZWQgd2l04oCmIiwic3RhcnQiOjEzMzA3
LCJlbmQiOjEzNDA1fSwiQmhuTzVPMFBweUpSRkVzRSI6eyJ0ZX
h0Ijoib25lIGhvdCBlbmNvZGVkIiwic3RhcnQiOjkyNDAsImVu
ZCI6OTI1NX0sIlNCUGZOb0V0S3ljQXBMYkYiOnsidGV4dCI6Im
FwcGxpZWQgdG8gaHVtYW4gcHJvbW90ZXIgcmVnaW9ucyBmcm9t
IHRoZSBFTlNFTUJMIHJlZ3VsYXRvcnkgYnVpbGQiLCJzdGFydC
I6MTUyMjIsImVuZCI6MTUzMTF9LCJvb1Y3R3BhMWkxR3VGYnly
Ijp7InRleHQiOiJQeXRob24gc2NyaXB0cy4iLCJzdGFydCI6MT
Y0NjUsImVuZCI6MTY0NjZ9LCJSTnJaMk1USVlvTlZYYm96Ijp7
InRleHQiOiJwYWRkZWQgdXNpbmcgcmFuZG9tbHkgZ2VuZXJhdG
VkIHNlcXVlbmNlcyB0byAxMjhicCIsInN0YXJ0IjoxNzI3Nywi
ZW5kIjoxNzMyNn0sIlBxUGJ5TTg5WHNnRnoxalUiOnsidGV4dC
I6IkVycm9yYmFycyBhcmUgNjglIGNvbmZpZGVuY2UgaW50ZXJ2
YWxzLiIsInN0YXJ0IjoxNzYwMiwiZW5kIjoxNzczMn0sImJMUG
JNd3lzWXpXT1dlTGoiOnsidGV4dCI6IjUwMDAgcmFuZG9tbHkg
Z2VuZXJhdGVkIHNlcXVlbmNlcy4iLCJzdGFydCI6MTc5NjQsIm
VuZCI6MTc5OTh9LCJCNmJ3SEtDVEs0ZldIS3laIjp7InRleHQi
OiJycyBhcmUgNjglIGNvbmZpZGVuY2UgaW50ZXJ2YWxzIFsiLC
JzdGFydCI6MTgzNjQsImVuZCI6MTg1NTB9LCJuTjFvN3BMeFM5
aktEM0tyIjp7InRleHQiOiJHLXJlZ2lzdGVyIHdhcyBjb3VudG
VkIiwic3RhcnQiOjE4ODcyLCJlbmQiOjE4ODk0fSwidXUwclJW
N2sxcE5Tb1lzRiI6eyJ0ZXh0IjoiSW5zdGVhZCB3ZSBkZWNpZG
VkIHRvIHVzZSBhbiBleGlzdGluZyBtZXRob2QsIHRoZSBHNGh1
bnRlciBhbGdvcml0aG0gcHJvcG9zZWQgYuKApiIsInN0YXJ0Ij
oyMDM3NCwiZW5kIjoyMDYxMn0sInh6U1pLMDBBcUpUMXNCVHoi
OnsidGV4dCI6IlNpbmNlIHdlIHdpc2hlZCB0byBtYWtlIGFzIG
ZldyBhc3N1bXB0aW9ucyBhcyBwb3NzaWJsZSwgYW5kIGdpdmVu
IHRoYXQgdGhlIEc0c2XigKYiLCJzdGFydCI6MjExMjUsImVuZC
I6MjEzMTd9LCJCcG1sUkdhbHdLR3JQU2doIjp7InRleHQiOiJ0
IGFsc28gcHJvZHVjZWQgc2lnbmlmaWNhbnRseSBtb3JlIHNlcX
VlbmNlcyB0aGF0IGRpZCBub3QgY29uZm9ybSB0byB0aGUgUXVh
ZHBh4oCmIiwic3RhcnQiOjIxOTc3LCJlbmQiOjIyMTE2fSwiQU
lxbXN3RlZCenlEWEV4RiI6eyJ0ZXh0IjoiNiwyMzcsOTQzIiwi
c3RhcnQiOjIyNjI4LCJlbmQiOjIyNjM3fSwiRGhCcWJSU25iNH
NmbjNYMiI6eyJ0ZXh0IjoiU2luY2UgbWFpbnRhaW5pbmcgc3Vj
aCBhbiBpbWJhbGFuY2UgaW4gdGhlIHRyYWluaW5nIGRhdGEgd2
91bGQgcHJvZHVjZSBhIHBvb3IgY+KApiIsInN0YXJ0IjoyMzYx
MiwiZW5kIjoyMzc5OX0sIlNYbE1TajN6NjhxU2xMRWEiOnsidG
V4dCI6Ik9yYW5nZSBsaW5lIHNob3dzIGxvd2VzcyBjdXJ2ZSBm
aXQuIiwic3RhcnQiOjI0MjU5LCJlbmQiOjI0Mjk0fSwiZGdLVT
lld2dKWmZKaXNoWSI6eyJ0ZXh0IjoiTW9kZWwgc2VsZWN0aW9u
IGFuZCB0cmFpbmluZyIsInN0YXJ0IjoyNDM2NSwiZW5kIjoyND
M5M30sImwwY0d3Qk8yemZ3ZW9HV3kiOnsidGV4dCI6IldlIGZv
dW5kIHRoYXQgbmVpdGhlciBHNEh1bnRlciBbQEJlZHJhdDIwMT
ZdIG5vciB0aGUgRzRSTkEtbGlrZSBtZXRob2QgcGVyZm9ybWXi
gKYiLCJzdGFydCI6MjgyNzQsImVuZCI6Mjg0NzR9LCIzU2tTRE
JsbnY4VTk2TVI1Ijp7InRleHQiOiJTaW5jZSB0aGVzZSBzZXF1
ZW5jZXMgb25seSBhY2NvdW50IGZvciBhIDI1JSBvZiBhbGwgdG
hlIHBvdGVudGlhbCBQRzQgZm9ybWluZyBz4oCmIiwic3RhcnQi
OjI5ODI5LCJlbmQiOjI5OTIxfSwiY1hoUlZ1N1NZeUs2bzlXZS
I6eyJ0ZXh0IjoiRm9yIGVhY2ggQkc0IGludGVydmFsLCB0aGUg
aGlnaGVzdCBzY29yaW5nIG92ZXJsYXBwaW5nIHByZWRpY3Rpb2
4gZm9yIGVhY2ggbW9kZeKApiIsInN0YXJ0IjozMTk2MCwiZW5k
IjozMjA1Nn0sIk9kZGxvRzhkeEQzUEZ0NXYiOnsidGV4dCI6Il
RoZXNlIHJlc3VsdHMgc3VnZ2VzdCB0aGF0IHRoZSBpbmZvcm1h
dGlvbiB3aXRoaW4gdGhlIEc0U2VxIGRhdGFzZXQsIHdoZW4gY2
FwdHXigKYiLCJzdGFydCI6MzI1NTgsImVuZCI6MzI3MDd9LCJW
TnppS2tXbUg4V3M4ejNhIjp7InRleHQiOiJUaGVzZSByZXN1bH
RzIHN1Z2dlc3QgdGhhdCB0aGUgaW5mb3JtYXRpb24gd2l0aGlu
IHRoZSBHNFNlcSBkYXRhc2V0LCB3aGVuIGNhcHR14oCmIiwic3
RhcnQiOjMyNTU4LCJlbmQiOjMyNzA3fSwiYVk4bkJJVlBHa0xR
VHJlUyI6eyJ0ZXh0IjoiUmVjZWl2ZXIgT3BlcmF0b3IgQ2hhcm
FjdGVyaXN0aWMgKFJPQykgY3VydmUgc2hvd2luZyB0aGUgcGVy
Zm9ybWFuY2Ugb2YgRzRTZWVxZeKApiIsInN0YXJ0IjozMjc3OS
wiZW5kIjozMjk4OX0sIkF3VWV3TXJXRkhLSE9SRm0iOnsidGV4
dCI6InJHNFNlZXFlciwgRzRTZWVxZXIgYW5kIHRoZSBHNFJOQS
IsInN0YXJ0IjozNzMwNSwiZW5kIjozNzMzOH0sIm92eVY0dTNN
OFdPN0JHVm4iOnsidGV4dCI6Im5ldXJhbCBuZXR3b3JrIiwic3
RhcnQiOjM3Njg3LCJlbmQiOjM3NzAxfSwiOU90MUdXV1dlMFZ5
QnQ0diI6eyJ0ZXh0IjoiV2UgaWRlbnRpZmllZCBQRzQgc2VxdW
VuY2VzIHNjb3JpbmcgbW9yZSB0aGFuIDAuOSBmb3Igd2hpY2gg
YSBzaW5nbGUgRy0+SCBjaGFuZ+KApiIsInN0YXJ0IjozOTQ3NC
wiZW5kIjozOTY4M30sImRza3FTazl4WEVDUWh0TmgiOnsidGV4
dCI6IlRoaXMgY291bGQgYmUgZHVlIHRvIGEgcmVkdWN0aW9uIG
luIEc0IHN0YWJpbGl0eSB3aXRoIGxvb3AgbGVuZ3RoLCBidXQg
Y291bGQgZXHigKYiLCJzdGFydCI6NDEyMjUsImVuZCI6NDE0Mj
F9LCJncnNsQkZCMkZhNGlYakdmIjp7InRleHQiOiJBZ2Fpbiwg
ZGlzdGFuY2Ugd2FzIGZvdW5kIHRvIGNvcnJlbGF0ZSBuZWdhdG
l2ZWx5IHdpdGggJW1tIHNjb3JlIChTcGVhcm1hbnMgcmhv4oCm
Iiwic3RhcnQiOjQzMjI5LCJlbmQiOjQzNDEyfSwib2NYNmkxeW
RiVGR1elJEcSI6eyJ0ZXh0Ijoib3hwbG90IHNob3dpbmcgdGhl
IHJlbGF0aW9uc2hpcCBiZXR3ZWVuICVtbSBzY29yZSBhbmQgZG
lzdGFuY2UgdG8gbmV4dCBHLWhhaXJwaeKApiIsInN0YXJ0Ijo0
MzkxOSwiZW5kIjo0NDE3N30sInVtS2xGbmc4M1VYZWl2WnUiOn
sidGV4dCI6IlRoZSBQRzQgc3BhY2Ugb2YgdGhlICpNLiBtdXNj
dWx1cyogZ2Vub21lIHdhcyBhbHNvIG1lYXN1cmVkLCBhbmQgZm
91bmQgdG8gY29udGHigKYiLCJzdGFydCI6NTA4MjksImVuZCI6
NTEwMjN9LCJjUWRSbmZyVVpscUtWYlJuIjp7InRleHQiOiIyRC
BLZXJuZWwgZGVuc2l0eSBlc3RpbWF0ZSBwbG90IHNob3dpbmcg
ZGlzdHJpYnV0aW9uIG9mIGFsbCBwb3NzaWJsZSAgdGV0cmFkIF
F14oCmIiwic3RhcnQiOjUyMzMwLCJlbmQiOjUyNTk3fSwiaXlk
MWE3SlE1c1lXaWZpdSI6eyJ0ZXh0IjoiSXQgaXMgYWJsZSB0by
Bwcm9jZXNzIHRoZSB3aG9sZSBodW1hbiBnZW5vbWUgaW4gYXBw
cm94aW1hdGVseSAxIGhvdXIgb24gYSA4IGNvcuKApiIsInN0YX
J0Ijo1MzE2OSwiZW5kIjo1MzI4NH19LCJjb21tZW50cyI6eyI2
bHhaYmRLMEE1ejhxWnhrIjp7ImRpc2N1c3Npb25JZCI6Ikd1a2
gzcmpjeHd1aFdIdGkiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEw
MTA2NzciLCJ0ZXh0IjoiSG93IG1hbnkgZmFsc2UgcG9zaXRpdm
VzPyIsImNyZWF0ZWQiOjE1MzI1MTg1NzIxMTF9LCJxWng0bzY2
WHBLbkhWRVVXIjp7ImRpc2N1c3Npb25JZCI6IldEQzNZaUF3bG
pFcVREZFciLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2Nzci
LCJ0ZXh0IjoiV2hhdCBpcyBhIGNvbnZvbHV0aW9uYWwgYW5kIH
JlY3VycmVudCBuZXVyYWwgbmV0d29yaz8gSG93IGRvZXMgaXQg
d29yaz8gSG93IGlzIGl0IGRpZmZlcmVudCBmcm9tIG90aGVyIG
1hY2hpbmUgbGVhcm5pbmcgYWxnb3Mgb3Igb3RoZXIgTk5zPyBc
blxuV2h5IGRpZCB5b3UgY2hvb3NlIHRvIHVzZSBhIENSTk4/Ii
wiY3JlYXRlZCI6MTUzMjUxOTIxNTMxM30sImpITTl4bGJWM005
YjVEUGgiOnsiZGlzY3Vzc2lvbklkIjoiV0RDM1lpQXdsakVxVE
RkVyIsInN1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRl
eHQiOiJBbG1vc3QgY2VydGFpbmx5IHlvdXIgYXVkaWVuY2Ugd2
lsbCBub3QgYmUgYSBtYWNoaW5lIGxlYXJuaW5nIGV4cGVydC4g
WW91IG5lZWQgdG8gZXhwbGFpbiB0aGUgbWV0aG9kLCBhbmQgd2
h5IHlvdSBhcmUgZG9pbmcgaXQgdGhpcyB3YXkuIiwiY3JlYXRl
ZCI6MTUzMjUxOTI2NDkyNX0sIklMSnU3U1AwMGpqTmtlUTQiOn
siZGlzY3Vzc2lvbklkIjoibWp2YzdWR0VINzE3N1RubSIsInN1
YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJIb3
cgZGlkIHlvdSBkZWNpZGUgdG8gdXNlIHRoaXMgYXJjaGV0ZWN0
dXJlPyIsImNyZWF0ZWQiOjE1MzI4NzQ4MzY0NzZ9LCJ2SlcyYn
h2aDgyTjZrWHdJIjp7ImRpc2N1c3Npb25JZCI6InUwbEhiNzlX
YW5IZ0xSTG4iLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2Nz
ciLCJ0ZXh0IjoiRGVmaW5lIiwiY3JlYXRlZCI6MTUzMjg3NDg1
NzUxOX0sIlFmOE9LdEI0Tk5PWnNha2MiOnsiZGlzY3Vzc2lvbk
lkIjoicjhWaWVPN2dQQzl1SE1iMCIsInN1YiI6IjEwMjIwNTc5
NzI3Njk0MTAxMDY3NyIsInRleHQiOiJJcyB0aGlzIGJlY2F1c2
UgZWFjaCBhbGdvcml0aG0gcHJlZGljdHMgdGhlIHNjb3JlIGZv
ciBhbiBpbnRlcnZhbCwgYW5kIHRoZXNlIGludGVydmFscyBhcm
UgZ2VuZXJhbGx5IHNtYWxsZXIgdGhhbiB0aGUgcGVha3M/IEJl
IG1vcmUgZXhwbGljdC4iLCJjcmVhdGVkIjoxNTMyODc1MTAzNj
ExfSwieGlTMGZEVVh0TEpVVHNSUyI6eyJkaXNjdXNzaW9uSWQi
OiJpWUZpeFpPZjBiS3JyUmROIiwic3ViIjoiMTAyMjA1Nzk3Mj
c2OTQxMDEwNjc3IiwidGV4dCI6IldoYXQgd2FzIHRoaXMgbnVt
YmVyIGZvciB0aGUgRE5BIEc0c2VxIHNldD8iLCJjcmVhdGVkIj
oxNTMyODc1MzI0MzczfSwiWDVwNElVQjE5dVVNUGpzNyI6eyJk
aXNjdXNzaW9uSWQiOiJCaG5PNU8wUHB5SlJGRXNFIiwic3ViIj
oiMTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6ImRlZmlu
ZSIsImNyZWF0ZWQiOjE1MzI4NzU1MzQ5MjZ9LCJLOU1VSnVVZD
ZpbHV4VDJsIjp7ImRpc2N1c3Npb25JZCI6IlNCUGZOb0V0S3lj
QXBMYkYiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLC
J0ZXh0IjoiV2hpY2ggY2VsbCB0eXBlLCBvciB3YXMgaXQgcmVn
aW9ucyB0aGF0IHdlcmUgcHJvbW90b3JzIGluIGFueSBjZWxsIH
R5cGU/IiwiY3JlYXRlZCI6MTUzMjg3NTYzMTQwOH0sIkd5WWJt
Zm9uQ0VBNlBzY2siOnsiZGlzY3Vzc2lvbklkIjoib29WN0dwYT
FpMUd1RmJ5ciIsInN1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3
NyIsInRleHQiOiJXaGVyZSB3b3VsZCBJIGZpbmQgdGhlc2UgcH
l0aG9uIHNjcmlwdHM/IiwiY3JlYXRlZCI6MTUzMjg3NTY3ODgx
Mn0sInJ4eVBROEw2N0RZRG9BRlUiOnsiZGlzY3Vzc2lvbklkIj
oiUk5yWjJNVElZb05WWGJveiIsInN1YiI6IjEwMjIwNTc5NzI3
Njk0MTAxMDY3NyIsInRleHQiOiJVbmlmcm9tYWxseSBzYW1wbG
VkPyIsImNyZWF0ZWQiOjE1MzI4NzU3MTQ3ODJ9LCJha1JGRGVM
VFNZV2tzOWhtIjp7ImRpc2N1c3Npb25JZCI6IlBxUGJ5TTg5WH
NnRnoxalUiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2Nzci
LCJ0ZXh0IjoiQ2FsY3VsYXRlZCBob3c/IiwiY3JlYXRlZCI6MT
UzMjg3NTcyODE0Nn0sInNNYUpJYUNTMnFFc3VYN0ciOnsiZGlz
Y3Vzc2lvbklkIjoiYkxQYk13eXNZeldPV2VMaiIsInN1YiI6Ij
EwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJHZW5lcmF0
ZWQgaG93PyIsImNyZWF0ZWQiOjE1MzI4NzU3NDkxMTh9LCJPdU
RWc215d0xkUEg4MXFJIjp7ImRpc2N1c3Npb25JZCI6IkI2YndI
S0NUSzRmV0hLeVoiLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMT
A2NzciLCJ0ZXh0IjoiRGl0dG8uIiwiY3JlYXRlZCI6MTUzMjg3
NTc5NDc5NH0sIjU2UkFvTVB1S3dXOFU1Uk8iOnsiZGlzY3Vzc2
lvbklkIjoibk4xbzdwTHhTOWpLRDNLciIsInN1YiI6IjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJEZWZpbmUgb3IgWC
1yZWYgc2VjdGlvbiB3aXRoIGRlZmluaXRpb24uIiwiY3JlYXRl
ZCI6MTUzMjg3NTgzNzk3NX0sIkNCNzVDeWJodk9TblhlMjQiOn
siZGlzY3Vzc2lvbklkIjoidXUwclJWN2sxcE5Tb1lzRiIsInN1
YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJXaG
F0IHdhcyB0aGUgcnVuIHRpbWVzPyBXYXMgaXQgcmVhbGx5IG5v
dCBmZWFzaWxiZSB0byBydW4gdGhlIHdob2xlIGdlbm9tZSBvbi
BTSEFSQz8iLCJjcmVhdGVkIjoxNTMyODc1ODk3NTg5fSwiakpn
WXVQRmZXN1dRNnV3YyI6eyJkaXNjdXNzaW9uSWQiOiJ4elNaSz
AwQXFKVDFzQlR6Iiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEw
Njc3IiwidGV4dCI6IkJldHRlciAqaW4gdml0cm8qLCBidXQgd2
hhdCBhYm91dCAqaW4gdml2byo/IiwiY3JlYXRlZCI6MTUzMjg3
NTk1MzQ1NH0sInQ2Z3p0TVZTajdsR1JuWkQiOnsiZGlzY3Vzc2
lvbklkIjoiQnBtbFJHYWx3S0dyUFNnaCIsInN1YiI6IjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJIb3cgbWFueT8gV2
hhdCBwcm9wb3J0aW9uIG9mIHRoZXNlIG92ZXJsYXBwZWQgRzRT
ZXEgcHJlZGljdGlvbnM/IiwiY3JlYXRlZCI6MTUzMjg3NjA0MT
QyMH0sImJiZXJCRzNZZTRraG1pUGoiOnsiZGlzY3Vzc2lvbklk
IjoiQnBtbFJHYWx3S0dyUFNnaCIsInN1YiI6IjEwMjIwNTc5Nz
I3Njk0MTAxMDY3NyIsInRleHQiOiJBbmQgd2hhdCBmcmFjdGlv
biBvZiBHNFNlcSBwb3NpdGl2ZSByZWdpb25zIG92ZXJsYXBwZW
Qgb25lIG9mIHRoZXNlPyIsImNyZWF0ZWQiOjE1MzI4NzYwNjk0
Mzd9LCJLSWtWUkVQUHVtNkNUT29MIjp7ImRpc2N1c3Npb25JZC
I6IkFJcW1zd0ZWQnp5RFhFeEYiLCJzdWIiOiIxMDIyMDU3OTcy
NzY5NDEwMTA2NzciLCJ0ZXh0IjoiQ29tbWEgZm9ybWF0IGJpZy
BudW1iZXJzIGxpa2UgdGhlc2UuIiwiY3JlYXRlZCI6MTUzMjg3
NjA5MDk3NH0sIkJLMzUxM2J3Qno0eEtZRU0iOnsiZGlzY3Vzc2
lvbklkIjoiRGhCcWJSU25iNHNmbjNYMiIsInN1YiI6IjEwMjIw
NTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJIb3cgbWF5IEc0c2
VxIHFpbmRvd3Mgd2l0aCBhIHNjb3JlIGxpa2UgdGhpcyB3ZXJl
IG1pc3NlZCBieSB0aGUgRzRIdW50ZXIgc2Nhbj8iLCJjcmVhdG
VkIjoxNTMyODc2MTgyNjEyfSwiZ29yVW1aSHhVeWJ2cEtudSI6
eyJkaXNjdXNzaW9uSWQiOiJTWGxNU2ozejY4cVNsTEVhIiwic3
ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Ikxv
b2tzIGxpa2UgdHdvIHBvcHVsYXRpb25zIHJhdGhlciB0aGFuIG
EgY29udGludW91cyByZWxhdGlvbnNoaXAuIiwiY3JlYXRlZCI6
MTUzMjg3NjM0NzM1M30sImpEQzZoS0l4RUdKalFQdmMiOnsiZG
lzY3Vzc2lvbklkIjoiZGdLVTlld2dKWmZKaXNoWSIsInN1YiI6
IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJIb3cgbG
9uZyBkaWQgdHJhaW5pbmcgdGFrZSwgYW5kIG9uIHdoYXQgaGFy
ZHdhcmU/IiwiY3JlYXRlZCI6MTUzMjg3NjM3MTY0Mn0sImxxRU
FDekRENmtmUUcwREIiOnsiZGlzY3Vzc2lvbklkIjoibDBjR3dC
TzJ6Zndlb0dXeSIsInN1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMD
Y3NyIsInRleHQiOiJXaGF0IGFib3V0IHF1YWRwYXJzZXI/IEkg
dW5kZXJzdGFuZCB5b3UgY2FuJ3QgY2FsY3VsYXRlZCBBVUMsIG
J1dCBjb3VsZCBhZGQgYSBzaW5nbGUgcG9pbnQgb24gZWFjaCBw
bG90IHRvIHNob3cgcGVyZm9ybWFuY2UuIiwiY3JlYXRlZCI6MT
UzMjg3NjQ4MTQwMX0sImpTVE50MGhJa2hCelVQajAiOnsiZGlz
Y3Vzc2lvbklkIjoiM1NrU0RCbG52OFU5Nk1SNSIsInN1YiI6Ij
EwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJBcmUgeW91
IHNheWluZyB0aGF0IG9ubHkgMjUlIG9mIHlvdXIgdGVzdCBzZX
F1ZW5jZXMgY29uZm9ybSB0byBRdWFkcGFyc2VyPyBXaGF0IGlz
IHRoZSBzcGxpY3QgYmV0d2VlbiBHNHNlcSArLy0/IiwiY3JlYX
RlZCI6MTUzMjg3NjU0OTkyMH0sInFoOXgzR2s0cWJ6MUdYMjki
OnsiZGlzY3Vzc2lvbklkIjoiY1hoUlZ1N1NZeUs2bzlXZSIsIn
N1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJB
cmUgdGhlc2Ugc2VxdWVuY2VzIHBhc3NlZCB0aHJvdWdoIHRoZS
BHNEh1bnRlciB0aHJlc2hvbGQgZmlyc3Q/IiwiY3JlYXRlZCI6
MTUzMjg3NjYyMTIxOX0sIkZBNmxReXpZT290THpvbXIiOnsiZG
lzY3Vzc2lvbklkIjoiT2RkbG9HOGR4RDNQRnQ1diIsInN1YiI6
IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJEb2VzIG
l0PyBJZiBHNEh1bnRlciBwZXJmb3JtcyB0aGUgc2FtZT8gSG93
IGRvIHlvdSBrbm93IGl0IGhhc24ndCBqdXN0IGxlYXJudCB0by
BkbyB0aGUgc2FtZSBhcyBHNEh1bnRlcj8iLCJjcmVhdGVkIjox
NTMyODc2NjU3NjY2fSwid1RMRU8wNWNGTWxLVVFudSI6eyJkaX
NjdXNzaW9uSWQiOiJWTnppS2tXbUg4V3M4ejNhIiwic3ViIjoi
MTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6IkFsc28gaW
50ZXJlc3RpbmcgdGhhdCBmb3IgdGhlIGhpZ2hlc3QgcHJlY2lz
aW9ucywgRzRIdW50ZXIgaGFzIHRoZSBiZXN0IHJlY2FsbC4iLC
JjcmVhdGVkIjoxNTMyODc2NzM0NTQ2fSwiaUpWUTBkdWQ0NTJv
RDVIcCI6eyJkaXNjdXNzaW9uSWQiOiJhWThuQklWUEdrTFFUcm
VTIiwic3ViIjoiMTAyMjA1Nzk3Mjc2OTQxMDEwNjc3IiwidGV4
dCI6IkFyZSB5b3UgZ29pbmcgdG8gdGFsayBhYm91dCB0aGUgZG
lmZmVyZW5jZSBiZXR3ZWVuIGEgUk9DIGN1cnZlIGFuZCBhUFIg
Y3VydmUsIHdoZW4gb25lIGlzIHVzZWZ1bCBvdmVyIHRoZSBvdG
hlciBhbmQgd2hhdCB0aGUgdHdvIHRlbGwgeW91IGluIGNvbXBh
cmlzb24gdG8gZWFjaCBvdGhlciBpbiB0aGlzIHNwZWNpZmljIG
V4YW1wbGU/IiwiY3JlYXRlZCI6MTUzMjg3Njc3OTY4M30sIjVo
V3pjZWFMNjZWek5OV1IiOnsiZGlzY3Vzc2lvbklkIjoiQXdVZX
dNcldGSEtIT1JGbSIsInN1YiI6IjEwMjIwNTc5NzI3Njk0MTAx
MDY3NyIsInRleHQiOiJXaGF0IGFib3V0IEc0SHVudGVyPyIsIm
NyZWF0ZWQiOjE1MzI4NzY4MjIyNzR9LCJSMGJwekMzOTZpejdt
NHgyIjp7ImRpc2N1c3Npb25JZCI6Im92eVY0dTNNOFdPN0JHVm
4iLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZXh0
IjoiI25vdGFsbE1MIiwiY3JlYXRlZCI6MTUzMjg3Njg1NTg2OX
0sImtTTWZ1bUhCUUNNY252bVciOnsiZGlzY3Vzc2lvbklkIjoi
OU90MUdXV1dlMFZ5QnQ0diIsInN1YiI6IjEwMjIwNTc5NzI3Nj
k0MTAxMDY3NyIsInRleHQiOiJXaGF0IGZyYWN0aW9uIG9mIHRo
ZXNlIGNvbmZvcm1lZCB0byBRdWFkcGFyc2VyPyBIb3cgbWFueS
B3ZXJlIEc0SHVudGVyIGhpZ2g/IERpZCB0aGUgY2hhbmdlIGhh
dmUgYSBiaWcgZWZmZWN0IG9uIEc0SHVudGVyIHNjb3JlPyBcbl
xuV2hhdCBmcmFjdGlvbiBvZiBwb3NpdGl2ZSBwcmVkaWN0aW9u
cyBkaWQgdGhpcyBlcXVhdGUgdG8/IiwiY3JlYXRlZCI6MTUzMj
g3Njk1NDAzMn0sImRNbW1hSURmMUl1bkF6cWUiOnsiZGlzY3Vz
c2lvbklkIjoiZHNrcVNrOXhYRUNRaHROaCIsInN1YiI6IjEwMj
IwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJBcmUgdGhlcmUg
Qkc0IGNoaXAgcGVha3MgdGhhdCBjb250YWluIG5vIEc0U2VxIG
hpZ2ggaW50ZXJ2YWw/IENvdWxkIHRoZXNlIGJlIGV4cGxhaW5l
ZCBieSB0aGlzPyIsImNyZWF0ZWQiOjE1MzI4NzcwODA5MDR9LC
J2ZUpjRElOWVBkMHoyVXZmIjp7ImRpc2N1c3Npb25JZCI6Imdy
c2xCRkIyRmE0aVhqR2YiLCJzdWIiOiIxMDIyMDU3OTcyNzY5ND
EwMTA2NzciLCJ0ZXh0IjoiRGlkIHlvdSBjaGVjayB3aGV0aGVy
IHRoZXNlIGZlYXR1cmVzIGFwcGVhciBpbiBzZXF1ZW5jZXMgbm
90IHByZWRpY3RlZCB0byBiZSBHNHM/IiwiY3JlYXRlZCI6MTUz
Mjg3NzI1MTU2Mn0sImFWbWUxdUxNY29pRWFGa2MiOnsiZGlzY3
Vzc2lvbklkIjoib2NYNmkxeWRiVGR1elJEcSIsInN1YiI6IjEw
MjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOiJJbnB1dHMgdG
8gRzRTZWVrZXIgYXJlIDEyOGJwIHJpZ2h0PyBTaG91ZG4ndCBo
dGVyZSBiZSBhIHNpbmdsZSB3aW5kb3cgdGhhdCBjb250YWlucy
Bib3RoIHRoZSBjbHVzdGVyZWQgMnggb3IgM3ggRy1ydW5zIGFu
ZCB0aGUgdGhlIGRpc3RpbCBHLXJ1biwgYXQgbGVhc3QgZm9yIH
RoZSAxMy01MCBjYXRlZ29yeT8iLCJjcmVhdGVkIjoxNTMyODc3
MzMwNzQwfSwidm40dVVwQmpYcUNaRjJrSCI6eyJkaXNjdXNzaW
9uSWQiOiJ1bUtsRm5nODNVWGVpdlp1Iiwic3ViIjoiMTAyMjA1
Nzk3Mjc2OTQxMDEwNjc3IiwidGV4dCI6Iklzbid0IHRoaXMgY2
FsY3VsYXRpb24gZm9yIGFyYWJpZG9wc2lzIHJlbGV2YW50PyIs
ImNyZWF0ZWQiOjE1MzI4Nzc0NTg4OTJ9LCJzOEZSNVoyVFV6Yj
d5UUNpIjp7ImRpc2N1c3Npb25JZCI6ImNRZFJuZnJVWmxxS1Zi
Um4iLCJzdWIiOiIxMDIyMDU3OTcyNzY5NDEwMTA2NzciLCJ0ZX
h0IjoiSSBndWVzcyB5b3Ugc2hvdWxkIGVpdGhlciBsYWJlbCB0
aGUgYXhlcyAoZGltMSwgZGltMiBldGMpLCBvciBub3QgaGF2ZS
BudW1iZXJzIG9uIGF4aXMgKHBvc3NpYmxlIGJvdGgpIiwiY3Jl
YXRlZCI6MTUzMjg3NzUwNjEwM30sInA3UmZIMlliV1NSaHJUb1
EiOnsiZGlzY3Vzc2lvbklkIjoiaXlkMWE3SlE1c1lXaWZpdSIs
InN1YiI6IjEwMjIwNTc5NzI3Njk0MTAxMDY3NyIsInRleHQiOi
JUcmFpbm5nIG9yIHByZWRpY3Rpb24/IiwiY3JlYXRlZCI6MTUz
Mjg3NzUyNjQyNH0sImdlSEtHcHJhZnp1U1ZEZmsiOnsiZGlzY3
Vzc2lvbklkIjoib2NYNmkxeWRiVGR1elJEcSIsInN1YiI6IjEw
ODUyMDAyOTMwMjI5NDY1MDQxNyIsInRleHQiOiJ5ZXMsIGJ1dC
B0aGUgbW9kZWwgZG9lc24ndCBzZWVtIHRvIGNvbnNpZGVyIHRo
ZSBkaXN0aWwgRy1ydW4gaW1wb3J0YW50IGZvciBhIGhpZ2ggc2
NvcmUsIHByZXN1bWFibHkgY2F1c2UgdGhlcmUgYXJlIHF1aXRl
IGEgbG90IG9mIGNhc2VzIHdoZXJlIGl0IGNhbm5vdCBzZWUgaX
QuIEkgd2lsbCBtYWtlIHRoaXMgY2xlYXJlciIsImNyZWF0ZWQi
OjE1MzI5NDU1MTcxNTZ9LCJmQ1J4SElQUWFtckZZSjB4Ijp7Im
Rpc2N1c3Npb25JZCI6IlBxUGJ5TTg5WHNnRnoxalUiLCJzdWIi
OiIxMDg1MjAwMjkzMDIyOTQ2NTA0MTciLCJ0ZXh0IjoiaXMgdG
hpcyBhbnkgY2xlYXJlcj8iLCJjcmVhdGVkIjoxNTMzMTQzODI2
NDg1fSwiU2M5aG9Ld0ExcTNCNW0xUCI6eyJkaXNjdXNzaW9uSW
QiOiJiTFBiTXd5c1l6V09XZUxqIiwic3ViIjoiMTA4NTIwMDI5
MzAyMjk0NjUwNDE3IiwidGV4dCI6ImFueSBjbGVhcmVyPyIsIm
NyZWF0ZWQiOjE1MzMxNDM5Mjk5MTR9LCJYU3dVZENEdks3VEVV
bjJIIjp7ImRpc2N1c3Npb25JZCI6IkI2YndIS0NUSzRmV0hLeV
oiLCJzdWIiOiIxMDg1MjAwMjkzMDIyOTQ2NTA0MTciLCJ0ZXh0
IjoiYmV0dGVyPyIsImNyZWF0ZWQiOjE1MzMxNDQwNDM4ODd9fS
wiaGlzdG9yeSI6Wy0xODc0MDYyNjk4LC0xNjI4MjE0OTI0LDEz
NjIyOTg5MDMsMTQ2NDYyNDkyLC0xNjg1MzE5NzMxLC0xNDcwNz
YzMjAyLC0xNTU3NDIyMDQyLDE2ODE5ODQzODEsNTE1NTA3MDY5
LC00ODQxNzA0MDksMTQ1MTU5NzU5LDEwNjI4MjAyOTksMjU2MD
E3NTI5LDEwNjI4MjAyOTksLTYxNzk3OTc4NiwtMTQyNDc4MjU1
NiwtMjAzNjEwMzIsMTk1NTMyMTgzNSw0NTIyNTQwMzAsLTU0Nz
U2MDg5MF19
-->
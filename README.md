# PEKA
Positionally-enriched k-mer analysis (PEKA) is a software package for identifying enriched protein-RNA binding motifs from CLIP datasets. PEKA compares k-mer enrichment in proximity of high-confidence crosslink sites (tXn - thresholded crosslinks), located within crosslinking peaks and having a high cDNA count, relative to low-count crosslink sites located outside of peaks (oXn - outside crosslinks). This approach reduces the effects of technical biases, such as uridine-preference of UV crosslinking. Each k-mer is assigned a PEKA score, which is used to rank the k-mers from the most to the least enriched. Additionally, PEKA provides comprehensive visualisations of motif enrichment profiles around the high-confidence crosslink sites and clusters the motifs that display similar profiles. PEKA also enables motif discovery within specific transcriptomic regions, including or excluding repetitive elements.

To interactively explore PEKA applied to all ENCODE eCLIP data, visit [iMaps](https://imaps.goodwright.org/apps/peka/).

Author: aram.amalietti@gmail.com


**Dependencies** (due to a breaking change in pandas >=1 using the versions below is recommended):
```
python=3.7
matplotlib=3.1.2
numpy=1.17.4
pybedtools=0.8.0
scipy=1.3.1
seaborn=0.9.0
plumbum=1.6.8
scikit-learn=0.21.3
pandas=0.24.2
textdistance=4.1.3
```
**Usage**:

usage: peka.py [-h] -i INPUTPEAKS -x INPUTXLSITES -g GENOMEFASTA -gi GENOMEINDEX -r REGIONS [-k [{3,4,5,6,7}]] [-o [OUTPUTPATH]] [-w [WINDOW]] [-dw [DISTALWINDOW]] [-t [TOPN]] [-p [PERCENTILE]] [-c [CLUSTERS]] [-s [SMOOTHING]]
               [-re [{masked,unmasked,repeats_only}]] [-a [ALLOUTPUTS]] [-sr [{whole_gene,intron,UTR3,other_exon,UTR5,ncRNA,intergenic,genome}]] [-sub [SUBSAMPLE]]

Search for enriched motifs around thresholded crosslinks in CLIP data.

required arguments:
  -i INPUTPEAKS, --inputpeaks INPUTPEAKS
                        CLIP peaks (intervals of crosslinks) in BED file format
  -x INPUTXLSITES, --inputxlsites INPUTXLSITES
                        CLIP crosslinks in BED file format
  -g GENOMEFASTA, --genomefasta GENOMEFASTA
                        genome fasta file, ideally the same as was used for read alignment
  -gi GENOMEINDEX, --genomeindex GENOMEINDEX
                        genome fasta index file (.fai)
  -r REGIONS, --regions REGIONS
                        genome segmentation file produced as output of "iCount segment" function

optional arguments:
  -h, --help            show this help message and exit
  -k [{3,4,5,6,7}], --kmerlength [{3,4,5,6,7}]
                        kmer length [DEFAULT 4]
  -o [OUTPUTPATH], --outputpath [OUTPUTPATH]
                        output folder [DEFAULT current directory]
  -w [WINDOW], --window [WINDOW]
                        window around thresholded crosslinks for finding enriched kmers [DEFAULT 25]
  -dw [DISTALWINDOW], --distalwindow [DISTALWINDOW]
                        window around enriched kmers to calculate distribution [DEFAULT 150]
  -t [TOPN], --topn [TOPN]
                        number of kmers ranked by z-score in descending order for clustering and plotting [DEFAULT 20]
  -p [PERCENTILE], --percentile [PERCENTILE]
                        percentile for considering thresholded crosslinks eg. percentile 0.7 means that the top 70 percent of crosslinks will be considered thresholded [DEFAULT 0.7]
  -c [CLUSTERS], --clusters [CLUSTERS]
                        how many enriched kmers to cluster and plot [DEFAULT 5]
  -s [SMOOTHING], --smoothing [SMOOTHING]
                        window used for smoothing kmer positional distribution curves [DEFAULT 6]
  -re [{masked,unmasked,repeats_only}], --repeats [{masked,unmasked,repeats_only}]
                        how to treat repeating regions within genome (options: "masked", "unmasked", "repeats_only"). When applying any of the options with the exception of repeats == "unmasked", a genome with soft-masked
                        repeat sequences should be used for input, ie. repeats in lowercase letters.
  -a [ALLOUTPUTS], --alloutputs [ALLOUTPUTS]
                        window used for smoothing kmer positional distribution curves [DEFAULT True]
  -sr [{whole_gene,intron,UTR3,other_exon,UTR5,ncRNA,intergenic,genome}], --specificregion [{whole_gene,intron,UTR3,other_exon,UTR5,ncRNA,intergenic,genome}]
                        choose to run PEKA on a specific region only [DEFAULT None]
  -sub [SUBSAMPLE], --subsample [SUBSAMPLE]
                        window used for smoothing kmer positional distribution curves [DEFAULT True]
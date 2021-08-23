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
scipy=1.6.2
seaborn=0.9.0
plumbum=1.6.8
scikit-learn=0.21.3
pandas=0.24.2
textdistance=4.1.3
```

**Set up**:

We recommend running PEKA in a conda environment so all the dependencies are managed for you, to set this up run the following command from your PEKA directory:
```
conda env create -f environment.yml
```
Before you run PEKA, activate your environment:
```
conda activate peka
```

**Usage**:
```
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
                        kmer length [DEFAULT 5]
  -o [OUTPUTPATH], --outputpath [OUTPUTPATH]
                        output folder [DEFAULT current directory]. Make sure the specified folder exists before execution of the script.
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
  -re [{masked,unmasked,repeats_only}], --repeats [{masked,unmasked,repeats_only,remove_repeats}]
                        how to treat repeating regions within genome (options: "masked", "unmasked", "repeats_only", "remove_repeats"). When applying any of the options with the exception of repeats == "unmasked", a genome with soft-masked
                        repeat sequences should be used for input, ie. repeats in lowercase letters. [DEFAULT "unmasked"]
  -a [ALLOUTPUTS], --alloutputs [ALLOUTPUTS]
                        controls the number of outputs, can be True/False [DEFAULT False]
  -sr [{whole_gene,intron,UTR3,other_exon,UTR5,ncRNA,intergenic,genome}], --specificregion [{whole_gene,intron,UTR3,other_exon,UTR5,ncRNA,intergenic,genome}]
                        choose to run PEKA on a specific region only, to specify multiple regions enter them space separated [DEFAULT None]
  -sub [SUBSAMPLE], --subsample [SUBSAMPLE]
                        if the crosslinks file is large, they can be subsampled to reduce runtime, can be True/False [DEFAULT True]
```

**Common issues**

If you have one of the following errors, it is because the numpy/pandas versions you are running are incompatible with PEKA.
To ensure you are using the correct versions, we reccomend using our conda environment.
```
AttributeError: type object 'object' has no attribute 'dtype'
```
or
```
TypeError: sum() got an unexpected keyword argument 'axis'
```
The script needs writing permission in the staging directory to save results and make an environment variable `TMPDIR` for temporary files. 
If you get `KeyError: 'TMPDIR'` a solution would be to type `export TMPDIR=<path_to_folder>` in terminal where you want to run the script.

**Outputs**

The default outputs produced by PEKA for each specified genomic region are:
- A pdf file with graphs showing k-mer occurrence distributions around thresholded crosslink sites for top n most enriched k-mers. K-mers are clustered based on their sequence and distributions.
- A csv file with clusters of top n k-mers.
- A tsv file with summed occurrence distributions of k-mers within defined clusters. Distributions in this file correspond to the curves on the last plot in the .pdf file.
- A tsv file with calculated PEKA score and occurrence distribution for all possible k-mers in the window -48 to +50 around thresholded crosslinks. 
  - This file also contains several other values that are calculated for a particular k-mer during the analysis, such as:
    - mtxn :  position of occurrence distribution max peak
    - prtxn : positions used for calculating motif enrichment relative to sampled background sequences
    - DtXn : average k-mer occurence in a distal window around thresholded crosslink sites
    - artxn : average k-mer relative occurrence at prtxn positions around thresholded crosslink sites
    - aroxn : average k-mer relative occurrence at prtxn positions around background crosslink sites
    - etxn : log2 relative k-mer enrichment, calculated as log2(artxn/aroxn)

If the option ALLOUTPUTS is set to True, the following files are also outputted for each specified genomic region:
- a bed file of thresholded crosslink sites
- a bed file with background crosslink sites (oxn)
- a tsv file with relative occurrence distribution (rtxn) for all possible k-mers. Relative occurrences are obtained by dividing raw occurrences with the average k-mer occurrence in the distal region (DtXn).

Other outputs:
- A tsv file with saved run parameters

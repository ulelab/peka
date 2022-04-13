# PEKA
Positionally-enriched k-mer analysis (PEKA) is a software package for identifying enriched protein-RNA binding motifs from CLIP datasets. PEKA compares k-mer enrichment in proximity of high-confidence crosslink sites (tXn - thresholded crosslinks), located within crosslinking peaks and having a high cDNA count, relative to low-count crosslink sites located outside of peaks (oXn - outside crosslinks). This approach reduces the effects of technical biases, such as uridine-preference of UV crosslinking. Each k-mer is assigned a PEKA score, which is used to rank the k-mers from the most to the least enriched. Additionally, PEKA provides comprehensive visualisations of motif enrichment profiles around the high-confidence crosslink sites and clusters the motifs that display similar profiles. PEKA also enables motif discovery within specific transcriptomic regions, including or excluding repetitive elements.

To interactively explore PEKA applied to all ENCODE eCLIP data, visit [iMaps](https://imaps.goodwright.org/apps/peka/).
For a more detailed description of PEKA method please refer to our preprint https://www.biorxiv.org/content/10.1101/2021.12.07.471544v1.

Author: aram.amalietti@gmail.com


## Dependencies (due to a breaking change in pandas >=1 using the versions below is recommended):
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

## Set up

We recommend running PEKA in a conda environment so all the dependencies are managed for you, to set this up run the following command from your PEKA directory:
```
conda env create -f environment.yml
```
Before you run PEKA, activate your environment:
```
conda activate peka
```

You may then install PEKA to the environment with the command:
```
python -m pip install . -vv --ignore-installed --no-deps
```
This command is also useful for installing development versions of PEKA.

## Usage
```
usage: peka.py [-h] -i INPUTPEAKS -x INPUTXLSITES -g GENOMEFASTA -gi GENOMEINDEX -r REGIONS \
               [-k [{3,4,5,6,7}]] \
               [-o [OUTPUTPATH]] \
               [-w [WINDOW]] \
               [-dw [DISTALWINDOW]] \
               [-t [TOPN]] \
               [-p [PERCENTILE]] \
               [-c [CLUSTERS]] \
               [-s [SMOOTHING]] \
               [-re [{remove_repeats,masked,unmasked,repeats_only}]] \
               [-pos [RELEVANT_POS_THRESHOLD] | -relax {True,False}] \
               [-a {True,False}] \
               [-sr {genome,whole_gene,intron,UTR3,other_exon,ncRNA,intergenic}] \
               [-sub {True,False}] \
               [-seed {True,False}]
```

🔴 **IMPORTANT NOTE!** Make sure all the required inputs, i.e. bed files, genome in fasta format, genome index and regions file, follow the same naming convention for chromosome names. Either all files must use the UCSC (GENCODE) naming convention, which prefixes chromosome names with "chr" ("chr1", ..., "chrM") or all files should use Ensembl naming convention ("1", ..., "MT").

```
required arguments:
  -i INPUTPEAKS, --inputpeaks INPUTPEAKS
                        CLIP peaks (intervals of crosslinks) in BED file
                        format
  -x INPUTXLSITES, --inputxlsites INPUTXLSITES
                        CLIP crosslinks in BED file format
  -g GENOMEFASTA, --genomefasta GENOMEFASTA
                        genome fasta file, ideally the same as was used for
                        read alignment
  -gi GENOMEINDEX, --genomeindex GENOMEINDEX
                        genome fasta index file (.fai)
  -r REGIONS, --regions REGIONS
                        genome segmentation file produced as output of "iCount
                        segment" function

optional arguments:
  -h, --help            show this help message and exit
  -k [{3,4,5,6,7}], --kmerlength [{3,4,5,6,7}]
                        kmer length [DEFAULT 5]
  -o [OUTPUTPATH], --outputpath [OUTPUTPATH]
                        output folder [DEFAULT current directory]
  -w [WINDOW], --window [WINDOW]
                        window around thresholded crosslinks for finding
                        enriched kmers [DEFAULT 20]
  -dw [DISTALWINDOW], --distalwindow [DISTALWINDOW]
                        window around enriched kmers to calculate distribution
                        [DEFAULT 150]
  -t [TOPN], --topn [TOPN]
                        number of kmers ranked by z-score in descending order
                        for clustering and plotting [DEFAULT 20]
  -p [PERCENTILE], --percentile [PERCENTILE]
                        Percentile for considering thresholded crosslinks.
                        Accepts a value between 0 and 1 [0, 1). Percentile 0.7
                        means that a cDNA count threshold is determined at
                        which >=70 percent of the crosslink sites within the
                        region have a cDNA count equal or below the threshold.
                        Thresholded crosslinks have cDNA count above the
                        threshold. [DEFAULT 0.7]
  -c [CLUSTERS], --clusters [CLUSTERS]
                        how many enriched kmers to cluster and plot [DEFAULT
                        5]
  -s [SMOOTHING], --smoothing [SMOOTHING]
                        window used for smoothing kmer positional distribution
                        curves [DEFAULT 6]
  -re [{remove_repeats,masked,unmasked,repeats_only}], --repeats [{remove_repeats,masked,unmasked,repeats_only}]
                        how to treat repeating regions within genome (options:
                        "masked", "unmasked", "repeats_only",
                        "remove_repeats"). When applying any of the options
                        with the exception of repeats == "unmasked", a genome
                        with soft-masked repeat sequences should be used for
                        input, ie. repeats in lowercase letters. [DEFAULT
                        "unmasked"]
  -pos [RELEVANT_POS_THRESHOLD], --relevant_pos_threshold [RELEVANT_POS_THRESHOLD]
                        Percentile to set as threshold for relevant positions.
                        Accepted values are floats between 0 and 99 [0, 99].
                        If threshold is set to 0 then all positions within the
                        set window (-w, default 20 nt) will be considered for
                        enrichment calculation. If threshold is not zero, it
                        will be used to determine relevant positions for
                        enrichment calculation for each k-mer. If the -pos
                        option is not set, then the threshold will be
                        automatically assigned based on k-mer lengthand number
                        of crosslinks in region.
  -relax {True,False}, --relaxed_prtxn {True,False}
                        Whether to relax automatically calculated prtxn
                        threshold or not. Can be 'True' or 'False', default is
                        'True'. If 'True', more positions will be included for
                        PEKA-score calculation across k-mers. Using relaxed
                        threshold is recommended, unless you have data of very
                        high complexity and are using k-mer length <=5. This
                        argument can't be used together with -pos argument,
                        which sets a user-defined threshold for relevant
                        positions. [DEFAULT True]
  -a {True,False}, --alloutputs {True,False}
                        controls the number of outputs, can be True or False
                        [DEFAULT False]
  -sr {genome,whole_gene,intron,UTR3,other_exon,ncRNA,intergenic}, --specificregion {genome,whole_gene,intron,UTR3,other_exon,ncRNA,intergenic}
                        choose to run PEKA on a specific region only, to
                        specify multiple regions enter them space separated
                        [DEFAULT None]
  -sub {True,False}, --subsample {True,False}
                        if the crosslinks file is very large, they can be
                        subsampled to reduce runtime, can be True/False
                        [DEFAULT True]
  -seed {True,False}, --set_seeds {True,False}
                        If you want to ensure reproducibility of results the
                        option set_seeds must be set to True. Can be True or
                        False [DEFAULT True]. Note that setting seeds reduces
                        the randomness of background sampling.
```

## Common issues

If you have one of the following errors, it is because the numpy/pandas versions you are running are incompatible with PEKA.
To ensure you are using the correct versions, we recommend using our conda environment.
```
AttributeError: type object 'object' has no attribute 'dtype'
```
or
```
TypeError: sum() got an unexpected keyword argument 'axis'
```
The script needs writing permission in the staging directory to save results and make an environment variable `TMPDIR` for temporary files.
If you get `KeyError: 'TMPDIR'` a solution would be to type `export TMPDIR=<path_to_folder>` in terminal where you want to run the script.

## Outputs

The default outputs produced by PEKA for each specified genomic region are:
- A pdf file with graphs showing k-mer occurrence distributions around thresholded crosslink sites for top n most enriched k-mers in the region spannig 
-50...50 nt around thersholded crosslink sites. K-mers are clustered based on their sequence and distributions. 
    - **NOTE: Even though top k-mers are selected with respect to the user-defined window, the graphs of top k-mers are not scaled to the selected window and 
    are always plotted in the -50...50 nt range.**
    - y-axis on the plots denotes a % of thresholded crosslinks for which a particular k-mer occurs at a specified position (y-axis).
    - **NOTE:** If clustering could not be optimized to yield less or equal to user-defined number of clusters, the number of clusters will be determined
    automatically by the algorithm.
- A tsv file with summed occurrence distributions of k-mers within defined clusters. Distributions in this file correspond to the curves on the last plot in the .pdf file.
- A tsv file with calculated PEKA score and occurrence distribution for all possible k-mers in the window -48 to +50 around thresholded crosslinks. 
  - This file also contains several other values that are calculated for a particular k-mer during the analysis, such as:
    - mtxn :  position of occurrence distribution max peak
    - prtxn : positions used for calculating motif enrichment relative to sampled background sequences
    - DtXn : average k-mer occurence in a distal window around thresholded crosslink sites
    - artxn : average k-mer relative occurrence at prtxn positions around thresholded crosslink sites
    - aroxn : average k-mer relative occurrence at prtxn positions around background crosslink sites
    - etxn : log2 relative k-mer enrichment, calculated as log2(artxn/aroxn)
    - p-value : a p-value for each k-mer is obtained under assumption that a distribution of obtained PEKA-scores is normal. K-mers with no PEKA-score are **omitted** from distribution. A distribution of PEKA-scores is converted to distribution of z-scores and p-value is obtained as area under the resulting standard Gaussian probability density function. <br>
      🔴 **CAUTION! The distribution of PEKA-score is not necessarily normal. If using p-value to gauge significance, make sure that enough k-mers indeed get the PEKA-score and that they are normally distributed. This can be achieved by lowering k-mer length or by decreasing the threshold for relevant positions, using the -pos argument.**

If the option ALLOUTPUTS is set to True, the following files are also outputted for each specified genomic region:
- a bed file of thresholded crosslink sites
- a bed file with background crosslink sites (oxn)
- a tsv file with relative occurrence distribution (rtxn) for all possible k-mers. Relative occurrences are obtained by dividing raw occurrences with the average k-mer occurrence in the distal region (DtXn).
- A csv file with clusters of top n k-mers.

Other outputs:
- A tsv file with saved run parameters

## Plotting k-mer relative occurrences in heatmap format (preprint)
To produce k-mer relative occurrence heatmaps as shown in Figure 1g of our preprint https://www.biorxiv.org/content/10.1101/2021.12.07.471544v1, 
use the script plotRelativeOccurrenceHeatmap.py in peka environment as follows:
```
# Example of run command
python3 ./src/plotRelativeOccurrenceHeatmap.py \
PathToKmerOccurrenceTable \
PathToKmerRtxnTable \
PathToOutputDirectory \
WindowAroundCrosslinkSite\
NumberOfTopKmers
```
### Arguments
PathToKmerOccurrenceTable -  file with k-mer PEKA-scores *5mer_distribution_whole_gene.tsv* <br>
PathToKmerRtxnTable - file with k-mer relative occurrences *5mer_rtxn_whole_gene.tsv* <br>
PathToOutputDirectory - location to which files will be saved <br>
WindowAroundCrosslinkSite - n nucleotides up and downstream from crosslink to plot <br>
NumberOfTopKmers - n top k-mers to plot <br>

To run on test data head to _TestData/RelativeOccurrenceHeatmap_

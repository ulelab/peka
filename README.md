# peka
Find motifs enriched around prominent crosslinks

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
```
    python3 <path_to_script> <sites_file_path> <genome_path> <genome_fai_path> <regions_file_path> <window> <window_distal> <kmer_length_input> <output_path> <top_n> <percentile> <clusters> <smoothing> <repeats>
```

`peak_file`: *intervals of crosslinks in BED file format;*  
`sites_file`: *crosslinks in BED file format;*  
`genome`: *FASTA file format, preferably the same as was used for alignment;*  
`genome_fai`: *FASTA index file;*  
`regions_file`: *custom genome segmentation file;*  
`window`: *region around thresholded crosslinks where positional
distributions are obtained by counting kmers per position, usualy 25;*  
`window_distal`: *region considered for background distributio (recommended 150);*  
`kmer_length`: *length (in nucleotides) of kmers to be analysed between 3 and 7, higher is not tested);*  
`output_path`: *path to folder where the outputs will be saved;*
`top_n`: *number of kmers ranked by PEKA-score in descending order for clustering and plotting (recommended 20);*  
`percentile`: *used for thresholding crosslinks (recommended 0.7);*  
`clusters`: *number of clusters of kmers(recommended 5);*  
`smoothing`: *window used for smoothing kmer positional distribution curves (recommended 6);*  
`repeats`: *how to treat repeating regions within genome (options: 'masked', 'unmasked', 'repeats_only'). When applying any of the options with the exception of repeats == 'unmasked', a genome with masked repeat sequences should be used for input;*  
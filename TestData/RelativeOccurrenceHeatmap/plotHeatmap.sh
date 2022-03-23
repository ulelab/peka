#!/bin/bash

# The code plotRelativeOccurrenceHeatmap.py takes the following arguments:
# PathToKmerOccurrenceTable (*5mer_distribution_whole_gene.tsv*)
# PathToKmerRtxnTable (*5mer_rtxn_whole_gene.tsv*)
# PathToOutputDirectory
# WindowAroundCrosslinkSite (plot n nucleotides up and downstream from crosslink)
# NumberOfTopKmers (plot n top k-mers)


# Exmple of run command:
# python3 ../../src/plotRelativeOccurrenceHeatmap.py PathToKmerOccurrenceTable PathToKmerRtxnTable PathToOutputDirectory WindowAroundCrosslinkSite NumberOfTopKmers


# To run on test data:
python3 ../../src/plotRelativeOccurrenceHeatmap.py \
./inputs/HepG2-TIA1-merged_5mer_distribution_whole_gene.tsv \
./inputs/HepG2-TIA1-merged_5mer_rtxn_whole_gene.tsv \
./results \
20 \
40

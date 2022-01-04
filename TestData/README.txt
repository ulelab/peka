peka_run.sh contains an example command that runs PEKA on test data.

Attention!
1. Before running adjust paths in the script accordingly (TMPDIR and output directory).
2. Unzip genome before running the peka_run.sh script, using the gunzip command.

Inputs:
- EXAMPLE PEAK FILE: K562-TIA1-chr20.xl_peaks.bed.gz

- EXAMPLE CROSSLINK FILE: K562-TIA1-chr20.xl.bed.gz

- GENOME: GRCh38.p12.genome.masked.fa.gz (UNZIP BEFORE RUNNING)

- GENOME INDEX (.FAI, GENERATED WITH SAMTOOLS): GRCh38.p12.genome.masked.fa.fai

- iCount SEGMENTATION: sorted.regions.gtf.gz

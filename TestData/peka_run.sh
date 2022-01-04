export TMPDIR='/home/user/temp' #Adjust this path to your temporary directory
output_dir='./results' #Adjust this path to your output directory

#Running the command
python3 ../peka.py \
-i './inputs/K562-TIA1-chr20.xl_peaks.bed.gz' \
-x './inputs/K562-TIA1-chr20.xl.bed.gz' \
-g './inputs/GRCh38.p12.genome.masked.fa' \
-gi './inputs/GRCh38.p12.genome.masked.fa.fai' \
-r './inputs/sorted.regions.gtf.gz' \
-sr intron UTR3 whole_gene \
-o $output_dir \
-re unmasked \
-k 5 \
-s 6 \
-p 0.7 \
-t 20 \
-dw 150 \
-w 25
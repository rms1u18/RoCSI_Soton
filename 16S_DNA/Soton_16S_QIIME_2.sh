#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=01:30:00

# look at the trimmed-summary.qzv interactive quality plots to decide where to trim the forward and reverse read in dada2
#   aim to have atleast a medium quality score of >20
#   and overlap of >20 bases

# With amplicon length 464 - I chose F 267 bp R 228 bp, for a 31 bp overlap

# opening Qiime 2
module load qiime/2019.1
source activate /local/software/qiime/2019.1
cd $PBS_O_WORKDIR

# denoise data with dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-paired-end.qza \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trunc-len-f 253 \
--p-trim-left-r 0 \
--p-trunc-len-r 173 \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza \
--o-table table.qza

# create summaries for dada2
qiime feature-table summarize \
 --i-table table.qza \
 --o-visualization table.qzv \
 --m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv

qiime metadata tabulate \
 --m-input-file denoising-stats.qza \
 --o-visualization denoising-stats.qzv

# generate a tree for phylogentics diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
 --p-n-threads 0 \
 --i-sequences rep-seqs.qza \
 --o-alignment aligned-rep-seqs.qza \
 --o-masked-alignment masked-aligned-rep-seqs.qza \
 --o-tree unrooted-tree.qza \
 --o-rooted-tree rooted-tree.qza

 # use the table.qvz to decide the sampling depth for alpha and beta diversity analysis
 #   you need a compromise to get as greater sampling depth as possible without excluding too many samples
 # I chose 50923 - to only omit the blanks, or 61430 which is the lowest feature count from the time series samples

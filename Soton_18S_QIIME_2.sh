#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=01:30:00

# look at the trimmed-summary.qzv interactive quality plots to decide where to trim the forward and reverse read in dada2
#   aim to have atleast a medium quality score of >20
#   and overlap of >20 bases

# With total amplicon length ~260 (~220 -primers) - I chose F 238 bp R 212 bp, for a 190 bp overlap

# opening Qiime 2
conda activate qiime2-2021.4
cd Desktop/Soton/Soton_18S

# denoise data with dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-paired-end.qza \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trunc-len-f 178 \
--p-trim-left-r 0 \
--p-trunc-len-r 149 \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza \
--o-table table.qza

# create summaries for dada2
qiime feature-table summarize \
 --i-table table.qza \
 --o-visualization table.qzv \
 --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv

qiime metadata tabulate \
 --m-input-file denoising-stats.qza \
 --o-visualization denoising-stats.qzv

# generate a tree for phylogentics diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
 --i-sequences rep-seqs.qza \
 --o-alignment aligned-rep-seqs.qza \
 --o-masked-alignment masked-aligned-rep-seqs.qza \
 --o-tree unrooted-tree.qza \
 --o-rooted-tree rooted-tree.qza

 # use the table.qvz to decide the sampling depth for alpha and beta diversity analysis
 #   you need a compromise to get as greater sampling depth as possible without excluding too many samples
 # I chose 55 - to only omit the blanks, or 135 which is the lowest feature count from the time series samples

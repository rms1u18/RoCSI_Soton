#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=0:40:00

# To start you will need  a project directory containing:
#   this bin/bash file
#   a directory with the rawreads
#   metadata.txt/metadata.tsv
#   latest Silva 18S only .fasta (copy from .fna) See part 3 of this script for more details on obtaining the lastest Silva database
#   latest Silva 18S only taxonomy_7_levels.txt

# Pipeline check list:
#   check that primers are correct for the current project
#   check that ppn and core number match

# opening Qiime 2
# module load qiime/2019.1
# source activate /local/software/qiime/2019.1
# cd $PBS_O_WORKDIR
conda activate qiime2-2021.4
cd Desktop/Soton/Soton_18S

# import raw sequences
# I think this is the hardest part. If you're sequnces follow the casava format see here:
# then go ahead and use the code below
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path rawreads \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux-paired-end.qza

# if they don't follow that format then you can create a manifest.csv see here:
# use ="$PWD/rawreads/"&B2 to add the filepath in excell
# If you don't know the PHRED offset you can use the vsearch tool to guess it:
# vsearch --fastq_chars rawreads/samplename_r1.fq
# then use the following code with the input-format either Phred33 or Phred64:

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.txt \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

# import metadata file (check the metadata file is valid in google sheets with Keemei add-on and download as tsv)
qiime metadata tabulate \
--m-input-file metadata.tsv \
--o-visualization tabulated-sample-metadata.qzv

# show number of rawreads and quality
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux-summary.qzv

# trim primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences paired-end-demux.qza \
--p-cores 16 \
--o-trimmed-sequences trimmed-paired-end.qza \
--p-front-f GTACACACCGCCCGTC \
--p-front-r TGATCCTTCTGCAGGTTCACCTAC

# show number of trimmed reads and quality
qiime demux summarize \
  --i-data trimmed-paired-end.qza \
  --o-visualization trimmed-summary.qzv

# look at the trimmed-summary.qzv interactive quality plots to decide where to trim the forward and reverse read in dada2
#   aim to have atleast a medium quality score of >20
#   and overlap of >20 bases

# With total amplicon length ~260 (~220 -primers) - I chose F 178 bp R 149 bp, for a 67 bp overlap

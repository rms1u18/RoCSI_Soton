#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=03:00:00

# use the table.qvz to decide the sampling depth for alpha and beta diversity analysis
#   you need a compromise to get as greater sampling depth as possible without excluding too many samples
# I chose 2997 - to only omit the blanks and keep atleast 1 of 3 pcr replicates

# check which parameters you want to compare with PERMANOVA and adjust code/filenames accordingly

# choose a value for alpha-rarefaction plotting from the table.qzv - approx the median frequency usually works

# check that reference dataset file names are correct

# make sure the correct primers are used to extract reference reads and the min max length = amplicon read length +F +R primer length +- 10%

# opening Qiime 2
conda activate qiime2-2022.2
cd Desktop/Soton_16S_RNA

# alpha and beta diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 2997 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

# test for significant associations between categorical metadata columns and alpha diversity
  qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/faith-pd-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/evenness_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/evenness-group-significance.qzv

# PERMANOVA
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Preservative \
  --o-visualization core-metrics-results/unweighted-unifrac-Preservative-significance.qzv \
  --p-pairwise


  # Alpha rarefaction plotting - get max depth from table.qzv choose a value around median frequency
  qiime diversity alpha-rarefaction \
    --i-table table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 12926 \
    --m-metadata-file metadata.tsv \
    --o-visualization alpha-rarefaction.qzv

  # Obtaining and importing reference dataset (for most recent dataset visit https://www.arb-silva.de/no_cache/download/archive/qiime/) the
  # Download latest release
  #   get fna from rep_set/rep_set_18S_only/99/*.fna
  #   get taxonomy.txt from taxonomy/16S_only/99/
  # Copy and import reference dataset
# cp silva_132_99_16S.fna silva_132_99_16S.fasta
#  qiime tools import \
#  --type 'FeatureData[Sequence]' \
#  --input-path silva_132_99_16S.fasta \
#  --output-path 99_16S_otus.qza

#  qiime tools import \
#  --type 'FeatureData[Taxonomy]' \
#  --input-format HeaderlessTSVTaxonomyFormat \
#  --input-path taxonomy_7_levels.txt \
#  --output-path ref-taxonomy.qza

# the latest 138.1 database isn't qiime compatible yet so I will use RESCRIPt to make it compatible (see tutorial https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494 )
# install RESCRIPt (see here: https://github.com/bokulich-lab/RESCRIPt )

# easy way to download all the relevant silva files
qiime rescript get-silva-data \
--p-version '138.1' \
--p-target 'SSURef_NR99' \
--p-include-species-labels \
--o-silva-sequences silva-138.1-ssu-nr99-seqs.qza \
--o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

# Cull low quality sequences with cull-seq - removes sequences containing 5 or more ambiguous bases and any homopolymers that are 8 or more bases in length
qiime rescript cull-seqs \
--i-sequences silva-138.1-ssu-nr99-seqs.qza \
--o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

# filter sequence length by taxonomy
# removes rRNA gene sequences that do not meet the following criteria: Archaea (16S) >= 900 bp, Bacteria (16S) >= 1200 bp, and any Eukaryota (18S) >= 1400 bp
qiime rescript filter-seqs-length-by-taxon \
--i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
--i-taxonomy silva-138.1-ssu-nr99-tax.qza \
--p-labels Archaea Bacteria Eukaryota \
--p-min-lens 900 1200 1400 \
--o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
--o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

#Dereplicate sequences and taxonomy
# Using --p-mode 'uniq' will retain indentical sequence records that have differeing taxonomies - this allows you to use q2-clawback later on to retrain the classifier based on EMPO 3 habitat types (see https://forum.qiime2.org/t/using-q2-clawback-to-assemble-taxonomic-weights/5859 and Kaehler et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6789115/ in this case water (saline))
qiime rescript dereplicate \
--i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
--i-taxa silva-138.1-ssu-nr99-tax.qza \
--p-rank-handles 'silva' \
--p-mode 'uniq' \
--o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

# We can now make our clasifier for use on full-length SSU sequences
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads  silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep-uniq.qza \
--o-classifier silva-138.1-ssu-nr99-classifier.qza

# Now we make the amplicon specific classifier using our forward and reverse primers
# extract sequences with primers
qiime feature-classifier extract-reads \
--i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer CCGYCAATTYMTTTRAGTTT \
--p-n-jobs 2 \
--p-read-orientation 'forward' \
--o-reads silva-138.1-ssu-nr99-seqs-515f-926r.qza

#dereplicate extracted region
qiime rescript dereplicate \
--i-sequences silva-138.1-ssu-nr99-seqs-515f-926r.qza \
--i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
--p-rank-handles 'silva' \
--p-mode 'uniq' \
--o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-515f-926r-uniq.qza \
--o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-515f-926r-derep-uniq.qza

#Now we build the amplicon-region specific classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva-138.1-ssu-nr99-seqs-515f-926r-uniq.qza \
--i-reference-taxonomy silva-138.1-ssu-nr99-tax-515f-926r-derep-uniq.qza \
--o-classifier silva-138.1-ssu-nr99-515f-926r-classifier.qza


  # Taxonomic analysis (takes a long time but don't use p-n-jobs as the laptop restarts - I left this overnight)
  qiime feature-classifier classify-sklearn \
    --i-classifier silva-138.1-ssu-nr99-515f-926r-classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

  qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

# filter out missing features from table.qza (these next 2 steps are only necissary is qiime taxa barplot doesn't work
qiime feature-table filter-features \
--i-table table.qza \
--m-metadata-file taxonomy.qza \
--o-filtered-table id-filtered-table.qza

#create id-filtered-table summary (can compare table.qzv's in qiime2view)
qiime feature-table summarize \
--i-table id-filtered-table.qza \
--o-visualization id-filtered-table.qzv \
--m-sample-metadata-file metadata.txt

# Create taxa barplots
  qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv

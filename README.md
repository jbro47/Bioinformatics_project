# ------------------------------------------------------------------------------
# MICROBIOME PIPELINE README
# ------------------------------------------------------------------------------

# 1. INSTALL REQUIRED SOFTWARE

# Download Miniconda Installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install SRA Toolkit to download FASTQ files from EMP
conda install -c bioconda sra-tools

# Install QIIME2 (Amplicon Plugin) - adjust version as needed
wget https://data.qiime2.org/distro/core/qiime2-2024.2-py38-linux-conda.yml
conda env create -n qiime2-2024.2 --file qiime2-2024.2-py38-linux-conda.yml
conda activate qiime2-2024.2

# Install FastQC and MultiQC for quality assessment
conda install -c bioconda fastqc multiqc


# 2. DOWNLOAD FASTQ FILES FROM ACCESSIONS

# List of ERR accessions to download
ACCESSIONS=(
ERR1535953 ERR1535957 ERR1535960 ERR1535961 ERR1535962
ERR1535968 ERR1535976 ERR1535977 ERR1896636 ERR1896653
ERR1896667 ERR1896681 ERR1896688 ERR1896712 ERR1896717
ERR1896795 ERR1896796 ERR1896805 ERR1896819 ERR1896840
ERR1896854 ERR1896864 ERR1896999 ERR1897029
)

# Download and gzip FASTQ files using prefetch + fastq-dump
for acc in "${ACCESSIONS[@]}"; do
    prefetch $acc                       # Downloads SRA archive
    fastq-dump $acc                     # Extracts and compresses FASTQ (.fastq.gz) single-end files
done
gzip *.fastq                            # Use gzip zip all the fastq files

# 3. RUN QUALITY CONTROL

# Create a folder for FastQC output
mkdir fastqc_reports

# Run FastQC on all fastq.gz files
fastqc *.fastq.gz -o fastqc_reports

# Aggregate FastQC reports with MultiQC
multiqc fastqc_reports -o multiqc_report

# 4. IMPORT DATA INTO QIIME2

# Prepare a manifest TSV file: manifest.tsv (ensure absolute paths)
# Format:
# sample-id  absolute-filepath  direction
# sample1  /path/to/ERR1535953_1.fastq.gz  forward
# ...

# Import into QIIME2 as single-end data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33V2


# 5. DENOISE USING DADA2

# Remove low-quality reads, correct errors, and infer ASVs
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 140 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza


# 6. SUMMARIZE OUTPUTS

# Get overview of feature table and representative sequences
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv


# 7. GENERATE PHYLOGENETIC TREE

# Align sequences, mask, build tree (for phylogenetic diversity metrics)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-alignment.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


# 8. RUN CORE DIVERSITY ANALYSIS

# Compute alpha and beta diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 40000 \   # Adjust this based on rarefaction depth from feature table summary
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results

# 9. RUN FEATURE CLASSIFIER

# Download the pre-trained SILVA 138 Naive Bayes classifier for QIIME2 2024.2
wget https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza

# Use the classifier to assign taxonomy to your representative ASV sequences
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \  # input: pre-trained classifier
  --i-reads rep-seqs.qza \                         # input: ASV sequences from DADA2
  --o-classification taxonomy.qza                  # output: taxonomic assignments

# Create a visual summary of the taxonomy assignments
qiime metadata tabulate \
  --m-input-file taxonomy.qza \                    # input: taxonomy from classifier
  --o-visualization taxonomy.qzv                   # output: visual summary view in QIIME2 View

# Create interactive bar plots of taxonomy by sample
qiime taxa barplot \
  --i-table table.qza \                            # input: feature table (ASV counts)
  --i-taxonomy taxonomy.qza \                      # input: taxonomy assignments
  --m-metadata-file sample-metadata.tsv \          # input: sample metadata
  --o-visualization taxa-barplot.qzv               # output: taxonomy bar plot (interactive)

# 10. INSTALL AND ACTIVATE PICRUSt2 ENVIRONMENT
# ------------------------------------------------------------------------------
# Clone the PICRUSt2 repo and install the environment using conda

git clone https://github.com/picrust/picrust2.git
cd picrust2

# Install the PICRUSt2 conda environment
conda env create -f picrust2-env.yaml
conda activate picrust2

# (Optional) Run a test to verify successful installation
picrust2_pipeline.py -h

# Export QIIME2 artifacts to FASTA and BIOM
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path exported_rep_seqs

qiime tools export \
  --input-path table.qza \
  --output-path exported_table

# Subsample top 5000 ASVs using built-in PICRUSt2 script
# NOTE: Run this from within the PICRUSt2 directory
scripts/filter_top_n.py \
  --input_fasta exported_rep_seqs/dna-sequences.fasta \
  --input_table exported_table/feature-table.biom \
  --n 5000 \
  --output_fasta top5000.fasta \
  --output_table top5000.biom

# Run full PICRUSt2 pipeline
picrust2_pipeline.py \
  -s top5000.fasta \         # Subsampled representative sequences (FASTA)
  -i top5000.biom \          # Subsampled feature table (BIOM)
  -o picrust2_output \       # Output directory
  -p 4                       # Number of threads (adjust as needed)

# Functional Beta Diversity Analysis
# ---------------------------------
# This workflow calculates Bray-Curtis dissimilarities between samples based on predicted
# microbial functional profiles (e.g., from PICRUSt2), performs Principal Coordinates Analysis (PCoA),
# visualizes functional clustering by lake and continent, and assesses group differences using PERMANOVA.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova 

# ---------------------------------
# Step 1: Load the Bray-Curtis distance matrix (e.g., from QIIME2 export)
# The matrix contains pairwise functional distances between all samples.
dist_matrix = pd.read_csv('braycurtis_distance_matrix.csv', index_col=0)
dm = DistanceMatrix(np.ascontiguousarray(dist_matrix.values), ids=dist_matrix.index)

# Step 2: Load sample metadata
# This file contains sample IDs and associated metadata (e.g., lake, continent)
metadata = pd.read_csv('minimal_manifest.tsv', sep='\t')

# ---------------------------------
# Step 3: Run PCoA to reduce dimensionality and extract sample positions
# PCoA converts the distance matrix into principal coordinate space
pcoa_results = pcoa(dm)
pcoa_df = pcoa_results.samples

# Step 4: Merge PCoA coordinates with metadata for visualization
pcoa_df['sample-id'] = pcoa_df.index
pcoa_with_metadata = pcoa_df.merge(metadata, on='sample-id')

# ---------------------------------
# Step 5: PCoA Plot - Colored by Lake
# This figure shows how samples cluster based on lake origin.
plt.figure(figsize=(10,8))
sns.scatterplot(
    data=pcoa_with_metadata,
    x='PC1', y='PC2',
    hue='lake',
    style='lake',
    s=100
)
plt.title('PCoA of Functional Beta Diversity (Colored by Lake)')
plt.xlabel(f'PC1 ({pcoa_results.proportion_explained.iloc[0]*100:.1f}% variance)')
plt.ylabel(f'PC2 ({pcoa_results.proportion_explained.iloc[1]*100:.1f}% variance)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig('pcoa_by_lake.png', dpi=300)
plt.show()

# ---------------------------------
# Step 6: PCoA Plot - Colored by Continent
# This figure highlights broader geographic separation across continents.
plt.figure(figsize=(10,8))
sns.scatterplot(
    data=pcoa_with_metadata,
    x='PC1', y='PC2',
    hue='continent',
    style='continent',
    s=100
)
plt.title('PCoA of Functional Beta Diversity (Colored by Continent)')
plt.xlabel(f'PC1 ({pcoa_results.proportion_explained.iloc[0]*100:.1f}% variance)')
plt.ylabel(f'PC2 ({pcoa_results.proportion_explained.iloc[1]*100:.1f}% variance)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig('pcoa_by_continent.png', dpi=300)
plt.show()

# ---------------------------------
# Step 7: PERMANOVA Tests
# PERMANOVA statistically assesses whether grouping variables (lake or continent)
# significantly explain variation in functional composition.

# Create grouping vectors from metadata
grouping_by_lake = metadata.set_index('sample-id').loc[dist_matrix.index, 'lake']
grouping_by_continent = metadata.set_index('sample-id').loc[dist_matrix.index, 'continent']

# PERMANOVA test comparing differences between lakes
result_lake = permanova(dm, grouping_by_lake, permutations=999)
print("\nPERMANOVA results by Lake:")
print(result_lake)

# PERMANOVA test comparing differences between continents
result_continent = permanova(dm, grouping_by_continent, permutations=999)
print("\nPERMANOVA results by Continent:")
print(result_continent)


# Functional Beta Diversity Analysis: Extended Visualizations
# ------------------------------------------------------------
# This section generates visualizations and statistical tests to explore
# differences in predicted functional profiles across three freshwater lakes.
# It includes:
# - Heatmap of KO abundances
# - Barplots of top functions per lake
# - Shannon alpha diversity boxplots
# - Differential function testing (Mann-Whitney U and Kruskal-Wallis)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import entropy, mannwhitneyu, kruskal

# ------------------------------------------------------------
# Step 1: Load KO (KEGG Ortholog) abundance table and metadata
# - KO table: samples (columns) x functions (rows)
# - Metadata includes lake and continent info for each sample
ko = pd.read_csv("KO_normalized.tsv", sep='\t', index_col=0).T
meta = pd.read_csv("minimal_manifest.tsv", sep='\t')

# ------------------------------------------------------------
# Step 2: Heatmap of Top 50 Most Abundant KOs
# - Summarizes relative abundance of functional genes across samples
# - Clustering shows similarity patterns between samples
top_kos = ko.sum().sort_values(ascending=False).head(50).index
ko_top = ko[top_kos].copy()
ko_top['lake'] = meta.set_index('sample-id').loc[ko_top.index, 'lake']

sns.clustermap(
    ko_top.drop(columns='lake'),
    row_cluster=True,
    col_cluster=True,
    cmap='viridis',
    figsize=(15, 10)
)
plt.title('Top 50 KO Functional Profiles Across Samples')
plt.savefig('heatmap_top50_KOs.png', dpi=300)

# ------------------------------------------------------------
# Step 3: Barplot of Mean KO Abundance per Lake
# - Aggregates top KO functions by lake
# - Highlights lake-specific functional signatures
ko_top['lake'] = meta.set_index('sample-id').loc[ko_top.index, 'lake']
grouped = ko_top.groupby('lake').mean()

grouped[top_kos[:10]].T.plot(kind='bar', figsize=(12,6))
plt.ylabel("Mean Relative Abundance")
plt.title("Top 10 KO Functions by Lake")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('barplot_top10_KOs_per_lake.png', dpi=300)

# ------------------------------------------------------------
# Step 4: Functional Alpha Diversity - Shannon Index
# - Shannon diversity is computed from KO abundances
# - Boxplot compares diversity levels between lakes
ko['shannon'] = ko.drop(columns='lake', errors='ignore').apply(lambda x: entropy(x[x > 0]), axis=1)
ko['lake'] = meta.set_index('sample-id').loc[ko.index, 'lake']

plt.figure(figsize=(8,6))
sns.boxplot(data=ko, x='lake', y='shannon')
plt.ylabel('Shannon Diversity (Functional)')
plt.title('Functional Diversity by Lake')
plt.savefig('boxplot_shannon_diversity.png', dpi=300)

# ------------------------------------------------------------
# Step 5: Differential Function Test - Mann-Whitney U Test
# - Compares abundance of one specific KO (K00001) between two lakes
# - Useful for identifying KO functions present in one lake but absent in another
ko['group'] = meta.set_index('sample-id').loc[ko.index, 'lake']
group1 = ko[ko['group'] == 'Great Lakes']['K00001']
group2 = ko[ko['group'] == 'Tiefwaren Lake']['K00001']
stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
print(f"K00001: p-value = {p:.4f}")

# ------------------------------------------------------------
# Step 6: Kruskal-Wallis Test - Multiple Group Comparison
# - Tests if KO abundance differs across all three lakes
# - Applies to the top 20 most abundant KO functions
kruskal_results = []
for ko_id in top_kos[:20]:  # limit to top 20 to speed up
    groups = [ko[ko['group'] == lake][ko_id] for lake in ko['group'].unique()]
    H, pval = kruskal(*groups)
    kruskal_results.append((ko_id, pval))

# Show significantly differentially abundant KO functions
sig_kos = [(k, p) for k, p in kruskal_results if p < 0.05]
print("Significant KOs (p < 0.05):", sig_kos)
  



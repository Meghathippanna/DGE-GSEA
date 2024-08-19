# RNA-seq Differential Expression Analysis

This repository contains the R script for performing differential expression analysis using DESeq2.

## Installation

Before running the analysis, ensure that you have the necessary packages installed in R. You can install them using the following commands:

```r
# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install DESeq2
BiocManager::install("DESeq2")

# Install Org.Hs.eg.db for gene annotation
BiocManager::install("org.Hs.eg.db")

# Install additional required packages
install.packages("ggplot2")
install.packages("ComplexHeatmap")
BiocManager::install("EnhancedVolcano")
'''

##Files Required
1. samplesheet.csv
This file contains metadata about the samples used in the RNA-seq experiment. Each row represents a sample, and columns include sample identifiers and experimental conditions.

Example
Format:

SampleID	condition
Sample1	control
Sample2	control
Sample3	case
Sample4	case
SampleID: Unique identifier for each sample.
condition: The experimental condition (e.g., control or case).

2. count_file.csv
This file contains the raw count data from the RNA-seq experiment. Rows correspond to genes, and columns correspond to sample IDs. The first column contains gene identifiers, and subsequent columns contain raw counts for each sample.

Example Format:

GeneID	Sample1	Sample2	Sample3	Sample4
ENSG000001	100	150	200	120
ENSG000002	300	400	250	350
ENSG000003	50	75	60	70
GeneID: Unique identifier for each gene (e.g., ENSEMBL gene ID).
Sample1, Sample2, etc.: Raw counts for each gene in each sample.

Running the Analysis
Clone this repository to your local machine.
Place your samplesheet.csv and count_file.csv in the working directory.
Open the R script and update the setwd() function with your working directory path.
Run the R script to perform the differential expression analysis.

Output
The script will generate the following outputs:

All_candidates_annotated.csv: A CSV file with all genes and their differential expression statistics, ordered by adjusted p-value.
Plots: MA plot, gene count plot, heatmap, and volcano plot visualizations of the differential expression results.

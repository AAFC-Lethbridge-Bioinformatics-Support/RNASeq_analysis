##########################################################################################
# AAFC Bioinformatics conference and workshop / September 19-20, 2019
# Lethbridge Research Development Centre, AB, Canada									
# Day 2: Workshop / Friday, 20th September 2019		
# 
# Session on RNA-Seq and differential gene expression analyses
# 
# Arun Kommadath
# Biology Study Leader - Bioinformatics
# arun.kommadath@canada.ca / Tel: 403-782-8584
##########################################################################################


##########################################################################################
# Part B: Differential expression (DE) Analysis
##########################################################################################

# Reference: edgeR User's manual:
# http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

##########################################################################################


##########################################################################################
# Useful info on RStudio:
##########################################################################################

# Every line that begins with a # is a comment and is not executed and those that do 
# not begin with a # are code that is run.

# To run a specific line of code, point cursor to that line and click on the ->Run button 
# on the top right of this pane. Alternatively, hit Ctrl-Enter (or Cmd-Enter on Mac OS).

# To run the complete code in this script, click the ->Source button on the top right of 
# this pane.

# On the bottom right pane:
# - Any plot generated in the code will appear in the Plots tab 
# - Files and directory structure are visible and accessible from the Files tab
# - Help page of a command will appear on the Help tab. The command can be searched for 
#   in the search box. Alternatively, type ? before a command and hit Enter to invoke the
#   help page. For example, in the console pane below, the help page for the function 
#   "print" can be invoked by typing ?print 

# In the console below (as in the terminal), using tab to auto-complete commands and 
# variable names will help speed up typing and avoid errors.

##########################################################################################


##########################################################################################
### Define and set the working directory
##########################################################################################

# The code below will combine the path to your home directory to the path to the folder
# de_analysis and save it to a variable
# de_analysis_dir <- paste0( Sys.getenv("HOME"), "/rnaseq_analysis/de_analysis/")
de_analysis_dir <- paste0( Sys.getenv("HOME"), "/R/x86_64-pc-linux-gnu-library/RNA-Seq")

# Print that variable
print(de_analysis_dir)

# Print current working directory
getwd()

# Set your working directory to the path defined and stored in de_analysis_dir
setwd(de_analysis_dir)

##########################################################################################


##########################################################################################
### Set global options/variables
# Define global variables and custom functions, if any, at this point so that they are 
# accessible from anywhere further in the code. 
# Other global options may also be set, if applicable
##########################################################################################

# The following code prevents character vectors being converted to factors by default.
options(stringsAsFactors = FALSE);

# Set false discovery rate cutoff
FDR_cutoff <- 0.05 

##########################################################################################



##########################################################################################
### Load libraries/packages required for this analysis
#########################################################################################

# CRAN packages (https://cran.r-project.org/)
.libPaths("refGenome")
.libPaths()
.libPaths(c(.libPaths(), "/isilon/lethbridge-rdc/training/R_env"))
library(refGenome)
# For Gene and Splice Site Annotation Using Annotation Data from 'Ensembl' and 'UCSC'
# Genome Browsers [https://CRAN.R-project.org/package=refGenome]

# Bioconductor packages (https://www.bioconductor.org/)

library(edgeR)
# For differential expression (DE) analysis [DOI: 10.18129/B9.bioc.edgeR]
# [http://bioconductor.org/packages/release/bioc/html/edgeR.html]

library(RNASeqPower)
# For power analysis [DOI: 10.18129/B9.bioc.RNASeqPower]
# http://bioconductor.org/packages/release/bioc/html/RNASeqPower.html

library(dplyr)

library(calibrate)
##########################################################################################


##########################################################################################
### Sample information
##########################################################################################

# Read sample information
sample_info <- read.table(file="sample_info.txt", header = T) 
# header=T above indicates that the first row of the sample_info.txt file is a header
# Type ?read.table for more info on that function

print(sample_info)

# Print the dimensions of the object sample_info (number of rows and columns) 
dim(sample_info)      
##########################################################################################


##########################################################################################
### Read gene annotation
##########################################################################################

# function from refGenome library used to create an empty object of class 'ensemblGenome'.
gtf <- ensemblGenome()

# Reading and parsing GTF files into refGenome objects. NOTE: GTF file should be unzipped.
# read.gtf(gtf, useBasedir = F, filename = "../reference_genome/Listeria_monocytogenes_egd_e.ASM19603v1.44.gtf")
read.gtf(gtf, filename = "GCF_000240135.3_ASM24013v3_genomic.gtf")

# Copy content of gtf table to a data.frame object called lmono_gene_anno
# copying content of gtf table to a data.frame object.
fusarium_gene_anno <- gtf@ev$gtf
#filtering data.frame object to only keep reads that include these biotypes.
biotypes <- c("protein_coding","tRNA","pseudogene","rRNA", "tRNA_pseudogene")
fusarium_gene_anno_filtered <- fusarium_gene_anno %>% dplyr::filter(gene_biotype %in% biotypes)

# fusarium_gene_anno = fusarium_gene_anno[seq(1, nrow(fusarium_gene_anno_one), 6), ]

# View first 3 rows of data.frame object lmono_gene_anno (if n is not specified, it defaults to 6)
head(fusarium_gene_anno_filtered, n=3)
# Checking to see if there are dups
length(unique(fusarium_gene_anno_filtered$gene_id))
length(unique(fusarium_gene_anno_filtered$gene_id)) == nrow(fusarium_gene_anno_filtered)
# Returns TRUE= no duplicates


# Retrieve the number of rows and columns of an object
dim(fusarium_gene_anno_filtered)

# View column names (headers) of the data.frame object lmono_gene_anno
print( colnames(fusarium_gene_anno_filtered) )

# Keep only required columns (reordered) and reset row.names
fusarium_gene_anno <- fusarium_gene_anno_filtered[, c("gene_id", "seqid","start","end",
                                                      "strand","gene_biotype")]
row.names(fusarium_gene_anno) <- NULL

# View first 6 (when n is not specified) rows of data.frame object lmono_gene_anno following 
# the editing performed above 
head(fusarium_gene_anno)
dim(fusarium_gene_anno)
# Checking for dups
length(unique(fusarium_gene_anno$gene_id))
length(unique(fusarium_gene_anno$gene_id)) == nrow(fusarium_gene_anno)

# Build a contingency table of the counts at each combination of factor levels witin column
# gene_biotype
table(fusarium_gene_anno$gene_biotype, useNA="always")

##########################################################################################


##########################################################################################
## Process read counts
##########################################################################################

# read_counts_file <- "combined_counts_unstranded.txt"
read_counts_file <- "/home/AAFC-AAC/chengal/R/x86_64-pc-linux-gnu-library/RNA-Seq/featureCounts/merged_gene_counts.txt"
# The file combined_counts_unstranded.txt in the current working directory is the same file
# that would have been created after successfully running the bacterial_rnaseq_pipeline.sh 
# script earlier. It has been copied here to enable running the DE abalysis script even if 
# you haven't completed running the rnaseq_pipeline.sh script successfully.

# Read read_counts_file into a data.frame called read_counts. Print ?read.table for info on
# this command and what the arguments mean
read_counts <- read.table(file = read_counts_file, header = T, check.names = F, row.names = 1)
head(read_counts)
dim(read_counts) # Note the total number of genes 

# Note that the first 5 columns of read_counts (gene annotation) are not required in the 
# read_counts file for further processing. So, remove those columns
# in this case, only the first column of read_counts is not requred.
read_counts <- read_counts[ , -c(1:1) ]
head(read_counts)
dim(read_counts)

# Also note that the column names are based on the full path to the alignment files (BAM files). 
# We could remove the path and the .bam suffix to limit the column names to sample names.
# Print ?gsub and ?basename for the respective help pages of those commands.
#colnames(read_counts) <- gsub(".sorted.bam","", basename(colnames(read_counts)) )
colnames(read_counts) <- gsub(".out.bam","", basename(colnames(read_counts)) )
colnames(read_counts) <- gsub("Aligned","", basename(colnames(read_counts)) )
colnames(read_counts) <- gsub(".sortedByCoord","", basename(colnames(read_counts)) )
head(read_counts)
dim(read_counts)

##########################################################################################
## Code retrieved from: https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping/18540144#18540144
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

##########################################################################################
### Create DGEList object 
# The edgeR software requires data to be stored in a simple list-based data object called a 
# DGEList (DGE is an acronym for digital gene expression). Use the function DGEList to create 
# the DGEList object containing the read counts data from the read_counts object, a grouping 
# factor based on treatment from the sample_info object, and the gene annotation data from 
# lmono_gene_anno object. 
# Type ?DGEList for details on this function
##########################################################################################
head(read_counts)
read_counts_temp <- read_counts[moveme(names(read_counts), "192.CA-MKK2_OX10 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "191.CA-MKK2_OX10 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "190.CA-MKK2_OX10 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "189.CA-MKK2_OX9 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "188.CA-MKK2_OX9 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "187.CA-MKK2_OX9 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "186.MKK2_OX5 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "185.MKK2_OX5 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "104.MKK2_OX5 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "103.MKK2_OX4 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "102.MKK2_OX4 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "101.MKK2_OX4 first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "C5.WT first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "B5.WT first")]
read_counts_temp <- read_counts_temp[moveme(names(read_counts_temp), "A5.WT first")]
head(read_counts_temp)
head(read_counts)
read_counts <- read_counts_temp
dim(read_counts)

# Checking for dups, no dups here.
read_counts_temp <- tibble::rownames_to_column(read_counts_temp, "VALUE")
head(read_counts_temp)
length(unique(read_counts_temp$VALUE)) == nrow(read_counts_temp)


dge <- DGEList(counts = read_counts,
               group = sample_info$treatment,
               genes = fusarium_gene_anno )

read_counts
sample_info
dge$samples

# Checking for dups in read counts
# No dups here
read_counts_temp <- read_counts
read_counts_temp <- tibble::rownames_to_column(read_counts_temp, "VALUE")
head(read_counts_temp)
length(unique(read_counts_temp$VALUE)) == nrow(read_counts_temp)

print(dge) # The list has 3 elements: counts, samples and genes

# View summary of the library sizes and note the range of the library sizes
summary(dge$samples$lib.size)
# View barplot of the library sizes (in millions) and note the sample with the lowest or 
# highest library size
png('sample_sizes.png', width = 1000,height = 800)
par(mar = c(10, 10, 10, 5), mgp = c(7, 1, 0))
barplot(dge$samples$lib.size/1e6, ylab="Number of sequenced reads (in millions)", xlab = "Sample", 
        main = "Sample Sizes", las=2, names=row.names( dge$samples), cex.names = 0.9)
dev.off()

# Creating sample size plot with replications combined
head(dge$samples$lib.size)
(dge$samples$lib.size)
(dge$samples)

# creating treatment vector
treatment <- c("WT", "MKK2_OX4", "MKK2_OX5", "CA-MKK2_OX9", "CA-MKK2_OX10")  
# creating lib.size vector
lib.size <- dge$samples$lib.size
wt_lib.size <- lib.size[1:3]
head (wt_lib.size)
wt_lib.size <- sum (wt_lib.size)
wt_lib.size == (58303726 + 53826621 + 67167679)

lib.size <- dge$samples$lib.size
ox4_lib.size <- lib.size[4:6]
head (ox4_lib.size)
ox4_lib.size <- sum (ox4_lib.size)
ox4_lib.size == (67433556 + 61891582 + 55023876)

lib.size <- dge$samples$lib.size
ox5_lib.size <- lib.size[7:9]
head (ox5_lib.size)
ox5_lib.size <- sum (ox5_lib.size)
ox5_lib.size == (63322569 + 60794130 + 56386070)

lib.size <- dge$samples$lib.size
ox9_lib.size <- lib.size[10:12]
head (ox9_lib.size)
ox9_lib.size <- sum (ox9_lib.size)
ox9_lib.size == (66555318 + 64064753 + 60686661)

lib.size <- dge$samples$lib.size
ox10_lib.size <- lib.size[13:15]
head (ox10_lib.size)
ox10_lib.size <- sum (ox10_lib.size)
ox10_lib.size == (63770523 + 54494637 + 55907455)

lib.size <- c(wt_lib.size, ox4_lib.size, ox5_lib.size, ox9_lib.size, ox10_lib.size)

#Creating Data Frame using above vectors    
sample_sizes_combined_treatments <- data.frame(treatment, lib.size)

#Printing the above Data Frame    
sample_sizes_combined_treatments
sample_sizes_combined_treatments$treatment

#creating plot
png('sample_sizes_combined_treatments.png', width = 700,height = 700)
par(mar = c(10, 10, 10, 5), mgp = c(7, 1, 0))
barplot(sample_sizes_combined_treatments$lib.size/1e6, ylab="Number of sequenced reads (in millions)", xlab = "Treatment", 
        main = "Treatment Sizes", las = 2, names=row.names(sample_sizes_combined_treatments$treatment), cex.names = 0.5)
dev.off()

##########################################################################################



##########################################################################################
### Filtering low expressed genes
# First, we will identify and remove genes that are expressed at very low levels across the
# majority of the samples. These genes are likely to be of little relevance to the treatment
# under investigation and will interfere with some of the statistical approximations made 
# further on if included in the analysis. As a rule of thumb, genes are kept if they are 
# expressed above a certain level in all samples representing the condition of interest with
# the smallest sample size. The expression level filter should be based on counts-per-million
# (CPM) and not the read counts themselves as differences in library sizes between samples 
# need to be accounted for while filtering. CPM is a standardized measure of read counts 
# obtained by dividing each read count by the corresponding library size (in millions) thus
# enabling comparison of read abundance across libraries of different sizes. The function 
# cpm does the conversion from read counts to CPM. Usually a gene is considered expressed 
# if it has a read count of 5-10. In the sample with the smallest library
# size, we look at the read counts and corresponding CPM values to decide on the CPM-based 
# expression level filter to choose.
##########################################################################################
# Determine lowest library size sample (llss)
llss <- row.names(dge$samples)[ which(dge$samples$lib.size == min(dge$samples$lib.size)) ]
print(llss)

# Determine cpm corresponding to minimum reads required for a gene to be considered expressed
read_counts_llss <- dge$counts[ , llss ]
cpm_llss <- cpm( read_counts_llss )
read_count_expr_threshold <- 10
# find the CPM value of genes which have a read count equal to the read_count_expr_threshold
cpm_expr_threshold <- unique( cpm_llss[ which( read_counts_llss == read_count_expr_threshold ), ] )
print(cpm_expr_threshold)

# Filter out very lowly expressed genes, keeping genes that are expressed at a reasonable level.
# We will keep only those genes that are expressed above cpm_expr_threshold level in at least 
# 3 samples, which is the sample size of each group of samples in our experiment. This approach
# will allow genes expressed in one condition but not the other to also be considered for 
# further analysis.
group_size <- 3
dge_expressed <- dge[ rowSums(cpm(dge) > cpm_expr_threshold) >= group_size, ]
dim(dge_expressed$counts)

# Testing for duplicates
length(unique(dge_expressed$genes$gene_id)) == nrow(dge_expressed$genes)


##########################################################################################


##########################################################################################
### Re-compute the library sizes after filtering
# The library size will decrease slightly after removing these genes, so recalculating the 
# library sizes of the DGEList object after the filtering is recommended.
##########################################################################################
dge_expressed$samples$lib.size <- colSums(dge_expressed$counts)
#debugging - checking to see if samples and groups are matched.
dge_expressed$samples

##########################################################################################


##########################################################################################
### Data normalisation
# This step is required for the following reasons (quote from the edgeR user guide): 
# "RNA-seq provides a measure of the relative abundance of each gene in each RNA sample, but
# does not provide any measure of the total RNA output on a per-cell basis. This commonly 
# becomes important when a small number of genes are very highly expressed in one sample, 
# but not in another. The highly expressed genes can consume a substantial proportion of the 
# total library size, causing the remaining genes to be under-sampled in that sample. Unless 
# this RNA composition effect is adjusted for, the remaining genes may falsely appear to be
# down-regulated in that sample. The calcNormFactors function normalizes for RNA composition 
# by finding a set of scaling factors for the library sizes that minimize the log-fold changes
# between the samples for most genes. The default method for computing these scale factors uses 
# a trimmed mean of M-values (TMM) between each pair of samples (Robinson and Oshlack, 2010). 
# We call the product of the original library size and the scaling factor the effective library 
# size. The effective library size replaces the original library size in all downsteam analyses."
##########################################################################################

# Normalizing - Compute effective library sizes using the default TMM normalization:
norm_dge_expressed <- calcNormFactors(dge_expressed, method="TMM")

##########################################################################################


##########################################################################################
### Data exploration with hierarchical clustering and MDS plots
##########################################################################################
# Convert expressed read counts to counts per million
cpm_expressed <- cpm(norm_dge_expressed, log=FALSE)
# Convert to log scale (for plots like heatmaps)
log2cpm_expressed <- cpm(norm_dge_expressed, log=TRUE)

# Note the difference in range of expression on the log scale compared to original
summary(log2cpm_expressed)
summary(cpm_expressed)

dissimilarity_matrix <- dist(t( cpm( norm_dge_expressed$counts ) ), method = "euclidean")
hclust_object <- hclust( dissimilarity_matrix, method = "ward.D" )
plot( hclust_object, las=1 )
# Note distinct clustering of treatment samples and control samples indicating substantial DE

# Beautify figure
hclust_dendrogram <- as.dendrogram( hclust_object )
nodePar <- list( cex = 0.7, lab.cex = 0.8, pch = c(NA, 19), col = "blue")
edgePar <- list(col = 2:3, lwd = 2:1)
png('cluster_dendogram.png', width = 800,height = 800)
par(mar = c(10,5,5,2), mgp=c(6,1,0))  
plot( hclust_dendrogram, nodePar = nodePar, edgePar = edgePar,
      las=1, cex.main=1, cex.sub=0.9, cex.axis=0.9,
      main="Hierarchical clustering dendrogram",
      sub=paste("distance measure=euclidean", "agglomeration method=ward.D2", sep = "\n"))
dev.off()

# MDS plots
# Quote from the edgeR user guide: "The function plotMDS draws a multi-dimensional scaling plot of 
# the RNA samples in which distances correspond to leading log-fold-changes between each pair 
# of RNA samples. The leading log-fold-change is the average (root-meansquare) of the largest 
# absolute log-fold-changes between each pair of samples. This plot can be viewed as a type of
# unsupervised clustering”.
# The number of top genes used to calculate pairwise distances can be specified, which by default
# is 500. Also, by default, the first and second principal components (PC) are plotted. You can 
# set it to any pair of PCs but this function is limited to plotting on a 2 dimensional scatterplot.
par(cex=0.7, mar = c(5,5,5,5))
png('MDS_treatment.png', width = 800,height = 600)
plotMDS(norm_dge_expressed, main="Coloured by TREATMENT", las=1,
        col=rainbow( length(unique(sample_info$treatment)))[as.factor(sample_info$treatment) ])
dev.off()
png('MDS_batch.png', width = 800,height = 600)
plotMDS(norm_dge_expressed, main="Coloured by BATCH", las=1,
        col=rainbow( length(unique(sample_info$batch)))[as.factor(sample_info$batch) ])
dev.off()
##########################################################################################







### COMPARISON OF CA-MKK2OX to wildtype
##########################################################################################
### Statistical modelling of data and tests for DE
# Quote from the edgeR user guide: “Linear models are associated with normally distributed 
# data and have long been used to describe multifactor microarray experiments. Read counts
# however follow a negative binomial distribution. Generalized linear models (GLMs) are an 
# extension of classical linear models to non-normally distributed response data. GLMs specify
# probability distributions according to their mean-variance relationship and require variance
# parameters (dispersions) to be estimated. For general experiments (with multiple factors), 
# edgeR uses the Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions. 
# The CR method can be used to calculate a common dispersion for all the genes, trended 
# dispersion depending on the gene abundance, or separate dispersions for individual genes. 
# It takes care of multiple factors by fitting GLMs with a design matrix”.

##########################################################################################
### Define design matrix
##########################################################################################
sample_info_camkk2ox_vs_wt <- filter(sample_info, (treatment == "CA-MKK2_OX" | treatment == "wildtype"))
print(sample_info_camkk2ox_vs_wt)
TREATMENT <- as.factor(sample_info_camkk2ox_vs_wt$treatment )

# Define reference as control so that all DE results will be interpreted as changes in
# treatment with respect to control ( treament vs. control) samples.
TREATMENT <- relevel(TREATMENT, ref = "wildtype")
print(TREATMENT)

# To test for differences between treatment and control while adjusting for any differences
# between the batches, we add a term BATCH to the model defined as follows:
# convert batch information to factors
BATCH <- as.factor( sample_info_camkk2ox_vs_wt$batch )

# create design matrix without interaction term
design1 <- model.matrix( ~ TREATMENT + BATCH )
print(design1)

colnames(design1)[2] <- "CA-MKK2_OX"
print(design1)

# subsetting norm_dge_expressed so that its treatments match the treatments in sample_info.
norm_dge_expressed
norm_dge_expressed_camkk2ox_vs_wt <- norm_dge_expressed

# filtering counts - delete columns corresponding to WT
head(norm_dge_expressed_camkk2ox_vs_wt$counts)
norm_dge_expressed_camkk2ox_vs_wt$counts<-
  norm_dge_expressed_camkk2ox_vs_wt$counts[,-c(4,5,6,7,8,9)]
head(norm_dge_expressed_camkk2ox_vs_wt$counts)

#filtering samples
norm_dge_expressed_camkk2ox_vs_wt$samples
row.names.remove_camkk2_ox_vs_wt <- c("101.MKK2_OX4", "102.MKK2_OX4", "103.MKK2_OX4",
                                      "104.MKK2_OX5", "185.MKK2_OX5", "186.MKK2_OX5")
norm_dge_expressed_camkk2ox_vs_wt$samples <- norm_dge_expressed_camkk2ox_vs_wt$samples[
  !(row.names(norm_dge_expressed_camkk2ox_vs_wt$samples) %in% row.names.remove_camkk2_ox_vs_wt), ]
norm_dge_expressed_camkk2ox_vs_wt$samples

norm_dge_expressed_camkk2ox_vs_wt
dispersions_design1 <- estimateDisp(norm_dge_expressed_camkk2ox_vs_wt, design1, robust=TRUE)
dispersions_design1$common.dispersion


# Note from above the value of the common negative binomial dispersion. The square-root of this
# value is called the biological coefficient of variation (BCV) and its magnitude is indicative 
# of the biological variation between the replicate samples used in the experiment. From edgeR
# user guide: “Typical values for the common BCV for datasets arising from well-controlled 
# experiments are 0.4 for human data, 0.1 data on genetically identical model organisms or 0.01
# for technical replicates”.
bcv_design1 <- signif(sqrt(dispersions_design1$common.dispersion), 5)
bcv_design1

##########################################################################################
### DE contrasts: treatment versus control
##########################################################################################
# fit genewise negative binomial GLMs for design1
glmfit_design1 <- glmQLFit(dispersions_design1, design1)

# genewise statistical tests can be performed for a given coefficient or coefficient contrast
qlf_design1 <- glmQLFTest( glmfit_design1, coef = "CA-MKK2_OX" )

# The function topTags can be called to extract the top n DE genes ranked by p-value or 
# absolute log-fold change. Save all results by setting n to NULL. Use Benjamini-Hochberg (BH)
# method to adjust p-values for multiple testing. The output is set to be sorted by ordinary
# p-values. We do not provide a specified adjusted p-value cutoff now.
lrt_top_design1 <- topTags(qlf_design1, adjust.method = "BH", sort.by = "PValue",
                           p.value = 0.05, n=NULL)
row.names(lrt_top_design1$table) <- NULL
head(lrt_top_design1$table)
dim(lrt_top_design1$table)

# summary of DE genes based only on adjusted p-value threshold
summary(decideTests(qlf_design1, p.value = 0.05, lfc = 0))

# Summary of DE genes based on adjusted p-value threshold and a 2 fold change threshold.
# Note the numbers for up and down regulated genes have reduced.
summary(decideTests(qlf_design1, p.value = 0.05, lfc = log2(2)))

# Instead of defining a hard cut-off for fold change, it would be better to select a value
# based on the power of our experiment to detect changes in expression, which would be 
# proportional to the sample size and affected by factors like biological variablity among
# samples within groups. So, performing a power analylsis is useful.

##########################################################################################
### Power analysis
# Use function rnapower from package RNASeqPower to perform a basic power analysis
##########################################################################################

# View help function by entering ?rnapower
# The average depth of coverage for the transcript is one of the parameters required for 
# rnapower function. This can be estimated as follows:
# (# of sequenced bases) / (transcriptome size)
# Sequenced bases is the number of reads x read length (twice that if paired-end).
# An estimate of the transciptome size can be obtained from the GTF file by summing the
# difference between the start and end of each annotated gene. So, depth of sequencing in the 
# sample with the lowest library size will be:
depth <- min( dge$samples$lib.size )*100*2 / sum( fusarium_gene_anno$end - fusarium_gene_anno$start )

# Estimate power for a range of fold change values between 1 and 5 with other factors fixed
power_analysis <- rnapower(depth, n=group_size, cv=bcv_design1, effect=seq(1, 5, by=0.1),
                           alpha=0.05)

png('CAMKK2OX_vs_WT_power.png', width = 800,height = 600)
par(mar=c(8,8,8,4))
plot(names(power_analysis), power_analysis, las=1,
     xlab="Fold change", ylab="Power", main="Power analysis CAMKK2OX vs WT",
     sub="[parameters: n=3; alpha=0.05; coverage=5X]")
abline(v=1.84, h=0.9, col="red")
abline(v=1.71, h=0.8, col="blue")
dev.off()

# Note from the plot that fold changes of 1.6 can be detected at a power of 0.8
# Note from the plot that fold changes of 1.75 can be detected at a power of 0.9
fold_change_cutoff_p_8_camkk2_vs_wt=1.71
fold_change_cutoff_p_9_camkk2_vs_wt=1.84
fold_change_cutoff_3.5=3.5
fold_change_cutoff_2=2

# Define DE genes with these cut-offs
# Power = .8
summary(decideTests(qlf_design1, lfc = log2(fold_change_cutoff_p_8_camkk2_vs_wt), p.value = 0.05))
de_camkk2ox_vs_wt_.8 <- lrt_top_design1$table[ lrt_top_design1$table$FDR < FDR_cutoff &
                                         ( abs(as.numeric(lrt_top_design1$table$logFC)) >
                                             log2(fold_change_cutoff_p_8_camkk2_vs_wt) ), ]
# Power = .9
summary(decideTests(qlf_design1, lfc = log2(fold_change_cutoff_p_9_camkk2_vs_wt), p.value = 0.05))
de_camkk2ox_vs_wt_.9 <- lrt_top_design1$table[ lrt_top_design1$table$FDR < FDR_cutoff &
                                              ( abs(as.numeric(lrt_top_design1$table$logFC)) >
                                                  log2(fold_change_cutoff_p_9_camkk2_vs_wt) ), ]
# LFC = 3.5
summary(decideTests(qlf_design1, lfc = log2(fold_change_cutoff_3.5), p.value = 0.05))
de_camkk2ox_vs_wt_3.5 <- lrt_top_design1$table[ lrt_top_design1$table$FDR < FDR_cutoff &
                                                ( abs(as.numeric(lrt_top_design1$table$logFC)) >
                                                    log2(fold_change_cutoff_3.5) ), ]

head(de_camkk2ox_vs_wt_.8)
print(nrow(de_camkk2ox_vs_wt_.8))
head(de_camkk2ox_vs_wt_.9)
print(nrow(de_camkk2ox_vs_wt_.9))
head(de_camkk2ox_vs_wt_3.5)
print(nrow(de_camkk2ox_vs_wt_3.5))

##########################################################################################
### Heat map
# To visualize the expression differences among samples for the DE genes identified.
##########################################################################################
png('CAMKK2OX_vs_WT_heatmap_.8.png', width = 1000,height = 10000)
par(mar= c(10,10,10,10))
coolmap((log2cpm_expressed[ match(de_camkk2ox_vs_wt_.8$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map CAMKK2OX vs WT (Power = 0.8)")
dev.off()

png('CAMKK2OX_vs_WT_heatmap_.9.png', width = 1000,height = 10000)
par(mar= c(10,10,10,10))
coolmap((log2cpm_expressed[ match(de_camkk2ox_vs_wt_.9$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map CAMKK2OX vs WT (Power = 0.9)")
dev.off()

png('CAMKK2OX_vs_WT_heatmap_3.5.png', width = 1000,height = 10000)
par(mar= c(10,10,10,10))
coolmap((log2cpm_expressed[ match(de_camkk2ox_vs_wt_3.5$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map CAMKK2OX vs WT (FC Cutoff = 3.5)")
dev.off()


##########################################################################################
### Volcano plots
# To visualize log fold changes for different genes.
##########################################################################################
# Code based off https://www.r-bloggers.com/2014/05/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/
head(de_camkk2ox_vs_wt_.8)
png("CAMKK2_OX_vs_WT_volcano_.8.png", width = 700, height = 1500)
with(de_camkk2ox_vs_wt_.8, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot CAMKK2OX vs WT (Power = 0.8)", xlim=c(-20,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_camkk2ox_vs_wt_.8, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_camkk2ox_vs_wt_.8, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_camkk2ox_vs_wt_.8, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_camkk2ox_vs_wt_.8<- de_camkk2ox_vs_wt_.8 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_camkk2ox_vs_wt_.8, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.65))
dev.off()


head(de_camkk2ox_vs_wt_.9)
png("CAMKK2_OX_vs_WT_volcano_.9.png", width = 700, height = 1500)
with(de_camkk2ox_vs_wt_.9, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot CAMKK2OX vs WT (Power = 0.9)", xlim=c(-20,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_camkk2ox_vs_wt_.9, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_camkk2ox_vs_wt_.9, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_camkk2ox_vs_wt_.9, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_camkk2ox_vs_wt_.9<- de_camkk2ox_vs_wt_.9 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_camkk2ox_vs_wt_.9, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.65))
dev.off()


head(de_camkk2ox_vs_wt_3.5)
png("CAMKK2_OX_vs_WT_volcano_3.5.png", width = 700, height = 1500)
with(de_camkk2ox_vs_wt_3.5, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot CAMKK2OX vs WT (FC Cutoff = 3.5)", xlim=c(-20,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_camkk2ox_vs_wt_3.5, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_camkk2ox_vs_wt_3.5, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_camkk2ox_vs_wt_3.5, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_camkk2ox_vs_wt_3.5<- de_camkk2ox_vs_wt_3.5 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_camkk2ox_vs_wt_3.5, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.65))
dev.off()

##########################################################################################
### Tests for over-represented GO and KEGG pathway terms
# Performing gene ontology (GO) and KEGG pathway term enrichment tests on our lists of DE
# genes will give us an idea of the kind of biological processes that the genes are involved 
# in. Basically the enrichment test is a modified Fisher’s exact test (hypergeometric test)
# to determine GO or KEGG terms that are statistically over-represented in our set of DE 
# genes compared to particular background set of genes (universe). The choice of the universe
# is crucial and should include all genes that have a chance to be expressed. For our purpose, 
# the universe is selected to be the set of all genes identified as expressed. There are 
# several packages that perform these tests available on Bioconductor like GOstats or goana
# and kegga functions available from package “limma” that is already loaded as a dependency
# of “edgeR” or as web based tools like DAVID (https://david.ncifcrf.gov/) or Panther 
# (http://pantherdb.org/). We test the enrichments seperately for up- and down-regulated 
# genes.
##########################################################################################

# Upregulated gene_ids
upregulated_camkk2ox_vs_wt_ids_.8 <- de_camkk2ox_vs_wt_.8$gene_id[ de_camkk2ox_vs_wt_.8$logFC>0 ]
write.table( upregulated_camkk2ox_vs_wt_ids_.8, file="CAMKK2OX_vs_WT_upregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
# Upregulated genes_ids w/ info on PValues and logFC etc.
upregulated_camkk2ox_vs_wt_all_info_.8 <- filter(de_camkk2ox_vs_wt_.8, de_camkk2ox_vs_wt_.8$logFC>0)
write.table( upregulated_camkk2ox_vs_wt_all_info_.8, file="CAMKK2OX_vs_WT_upregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)

downregulated_camkk2ox_vs_wt_ids_.8 <- de_camkk2ox_vs_wt_.8$gene_id[ de_camkk2ox_vs_wt_.8$logFC<0 ]
write.table( downregulated_camkk2ox_vs_wt_ids_.8, file="CAMKK2OX_vs_WT_downregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
downregulated_camkk2ox_vs_wt_all_info_.8 <- filter(de_camkk2ox_vs_wt_.8, de_camkk2ox_vs_wt_.8$logFC<0)
write.table( downregulated_camkk2ox_vs_wt_all_info_.8, file="CAMKK2OX_vs_WT_downregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)



upregulated_camkk2ox_vs_wt_ids_.9 <- de_camkk2ox_vs_wt_.9$gene_id[ de_camkk2ox_vs_wt_.9$logFC>0 ]
write.table( upregulated_camkk2ox_vs_wt_ids_.9, file="CAMKK2OX_vs_WT_upregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
# Upregulated genes_ids w/ info on PValues and logFC etc.
upregulated_camkk2ox_vs_wt_all_info_.9 <- filter(de_camkk2ox_vs_wt_.9, de_camkk2ox_vs_wt_.9$logFC>0)
write.table( upregulated_camkk2ox_vs_wt_all_info_.9, file="CAMKK2OX_vs_WT_upregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)

downregulated_camkk2ox_vs_wt_ids_.9 <- de_camkk2ox_vs_wt_.9$gene_id[ de_camkk2ox_vs_wt_.9$logFC<0 ]
write.table( downregulated_camkk2ox_vs_wt_ids_.9, file="CAMKK2OX_vs_WT_downregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
downregulated_camkk2ox_vs_wt_all_info_.9 <- filter(de_camkk2ox_vs_wt_.9, de_camkk2ox_vs_wt_.9$logFC<0)
write.table( downregulated_camkk2ox_vs_wt_all_info_.9, file="CAMKK2OX_vs_WT_downregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)



upregulated_camkk2ox_vs_wt_ids_3.5 <- de_camkk2ox_vs_wt_3.5$gene_id[ de_camkk2ox_vs_wt_3.5$logFC>0 ]
write.table( upregulated_camkk2ox_vs_wt_ids_3.5, file="CAMKK2OX_vs_WT_upregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
# Upregulated genes_ids w/ info on PValues and logFC etc.
upregulated_camkk2ox_vs_wt_all_info_3.5 <- filter(de_camkk2ox_vs_wt_3.5, de_camkk2ox_vs_wt_3.5$logFC>0)
write.table( upregulated_camkk2ox_vs_wt_all_info_3.5, file="CAMKK2OX_vs_WT_upregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)

downregulated_camkk2ox_vs_wt_ids_3.5 <- de_camkk2ox_vs_wt_3.5$gene_id[ de_camkk2ox_vs_wt_3.5$logFC<0 ]
write.table( downregulated_camkk2ox_vs_wt_ids_3.5, file="CAMKK2OX_vs_WT_downregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
downregulated_camkk2ox_vs_wt_all_info_3.5 <- filter(de_camkk2ox_vs_wt_3.5, de_camkk2ox_vs_wt_3.5$logFC<0)
write.table( downregulated_camkk2ox_vs_wt_all_info_3.5, file="CAMKK2OX_vs_WT_downregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)


# The above files are saved to rnaseq_analysis/de_analysis. Click Refresh in the Files tab on 
# the bottom right panel if you do not see the files there). Export those files in zipped format
# to your local machine by selecting them from the Files tab on the bottom right panel and 
# clicking export (under More). Choose "Save as" option to save the exported files to a defined 
# folder on your local machine.

# Use those files to continue with overrepresentation analysis on http://pantherdb.org/
# We test the enrichments seperately for up- and down-regulated genes and using the list of 
# expressed genes as the background.





### COMPARISON OF MKK2OX to wildtype
##########################################################################################
### Statistical modelling of data and tests for DE
# Quote from the edgeR user guide: “Linear models are associated with normally distributed 
# data and have long been used to describe multifactor microarray experiments. Read counts
# however follow a negative binomial distribution. Generalized linear models (GLMs) are an 
# extension of classical linear models to non-normally distributed response data. GLMs specify
# probability distributions according to their mean-variance relationship and require variance
# parameters (dispersions) to be estimated. For general experiments (with multiple factors), 
# edgeR uses the Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions. 
# The CR method can be used to calculate a common dispersion for all the genes, trended 
# dispersion depending on the gene abundance, or separate dispersions for individual genes. 
# It takes care of multiple factors by fitting GLMs with a design matrix”.

##########################################################################################
### Define design matrix
##########################################################################################
# convert treatment information to factors
# filtering out sample info -- keeping only treatments wildtype and mkk2ox.
sample_info_mkk2ox_vs_wt <- filter(sample_info, (treatment == "wildtype" | treatment == "MKK2_OX"))
print(sample_info_mkk2ox_vs_wt)
TREATMENT <- as.factor( sample_info_mkk2ox_vs_wt$treatment )

# Define reference as control so that all DE results will be interpreted as changes in
# treatment with respect to control ( treament vs. control) samples.
TREATMENT <- relevel(TREATMENT, ref = "wildtype")
print(TREATMENT)

# To test for differences between treatment and control while adjusting for any differences
# between the batches, we add a term BATCH to the model defined as follows:
# convert batch information to factors
BATCH <- as.factor( sample_info_mkk2ox_vs_wt$batch )

# create design matrix without interaction term
design2 <- model.matrix( ~ TREATMENT + BATCH )
print(design2)

colnames(design2)[2] <- "MKK2_OX" 
print(design2)

# subsetting norm_dge_expressed so that its treatments match the treatments in sample_info.
norm_dge_expressed
norm_dge_expressed_mkk2ox_vs_wt <- norm_dge_expressed

# filtering counts
head(norm_dge_expressed_mkk2ox_vs_wt$counts)
norm_dge_expressed_mkk2ox_vs_wt$counts<-
  norm_dge_expressed_mkk2ox_vs_wt$counts[,-c(10,11,12,13,14,15)] #delete columns 1 and 4
head(norm_dge_expressed_mkk2ox_vs_wt$counts)

#filtering samples
norm_dge_expressed_mkk2ox_vs_wt$samples
row.names.remove_mkk2ox_vs_wt <- c("187.CA-MKK2_OX9", "188.CA-MKK2_OX9", "189.CA-MKK2_OX9", "190.CA-MKK2_OX10",
                                   "191.CA-MKK2_OX10", "192.CA-MKK2_OX10")
norm_dge_expressed_mkk2ox_vs_wt$samples <- norm_dge_expressed_mkk2ox_vs_wt$samples[
  !(row.names(norm_dge_expressed_mkk2ox_vs_wt$samples) %in% row.names.remove_mkk2ox_vs_wt), ]


norm_dge_expressed_mkk2ox_vs_wt
dispersions_design2 <- estimateDisp(norm_dge_expressed_mkk2ox_vs_wt, design2, robust=TRUE)
dispersions_design2$common.dispersion

# Note from above the value of the common negative binomial dispersion. The square-root of this
# value is called the biological coefficient of variation (BCV) and its magnitude is indicative 
# of the biological variation between the replicate samples used in the experiment. From edgeR
# user guide: “Typical values for the common BCV for datasets arising from well-controlled 
# experiments are 0.4 for human data, 0.1 data on genetically identical model organisms or 0.01
# for technical replicates”.
bcv_design2 <- signif(sqrt(dispersions_design2$common.dispersion), 5)
bcv_design2

##########################################################################################
### DE contrasts: treatment versus control
##########################################################################################
# fit genewise negative binomial GLMs for design1
glmfit_design2 <- glmQLFit(dispersions_design2, design2)

# genewise statistical tests can be performed for a given coefficient or coefficient contrast
qlf_design2 <- glmQLFTest( glmfit_design2, coef = "MKK2_OX" )

# lrt_top_design is a lot smaller than for the previous comparison!
# The function topTags can be called to extract the top n DE genes ranked by p-value or 
# absolute log-fold change. Save all results by setting n to NULL. Use Benjamini-Hochberg (BH)
# method to adjust p-values for multiple testing. The output is set to be sorted by ordinary
# p-values. We do not provide a specified adjusted p-value cutoff now.
lrt_top_design2 <- topTags(qlf_design2, adjust.method = "BH", sort.by = "PValue",
                           p.value = 0.05, n=NULL)
row.names(lrt_top_design2$table) <- NULL
head(lrt_top_design2$table)
dim(lrt_top_design2$table)

# summary of DE genes based only on adjusted p-value threshold
summary(decideTests(qlf_design2, p.value = 0.05, lfc = 0))

# Summary of DE genes based on adjusted p-value threshold and a 2 fold change threshold.
# Note the numbers for up and down regulated genes have reduced.
summary(decideTests(qlf_design2, p.value = 0.05, lfc = log2(2)))

# Instead of defining a hard cut-off for fold change, it would be better to select a value
# based on the power of our experiment to detect changes in expression, which would be 
# proportional to the sample size and affected by factors like biological variablity among
# samples within groups. So, performing a power analylsis is useful.

##########################################################################################
### Power analysis
# Use function rnapower from package RNASeqPower to perform a basic power analysis
##########################################################################################

# View help function by entering ?rnapower
# The average depth of coverage for the transcript is one of the parameters required for 
# rnapower function. This can be estimated as follows:
# (# of sequenced bases) / (transcriptome size)
# Sequenced bases is the number of reads x read length (twice that if paired-end).
# An estimate of the transciptome size can be obtained from the GTF file by summing the
# difference between the start and end of each annotated gene. So, depth of sequencing in the 
# sample with the lowest library size will be:
depth <- min( dge$samples$lib.size )*100*2 / sum( fusarium_gene_anno$end - fusarium_gene_anno$start )

# Estimate power for a range of fold change values between 1 and 5 with other factors fixed
power_analysis <- rnapower(depth, n=group_size, cv=bcv_design2, effect=seq(1, 5, by=0.1),
                           alpha=0.05)

png('MKK2OX_vs_WT_power.png', width = 800,height = 600)
par(mar=c(8,8,8,4))
plot(names(power_analysis), power_analysis, las=1,
     xlab="Fold change", ylab="Power", main="Power analysis MKK2OX vs WT",
     sub="[parameters: n=3; alpha=0.05; coverage=5X]")
abline(v=1.68, h=0.9, col="red")
abline(v=1.56, h=0.8, col="blue")
dev.off()


# Note from the plot that fold changes of 1.6 can be detected at a power of 0.8
# Note from the plot that fold changes of 1.75 can be detected at a power of 0.9
fold_change_cutoff_p_8_mkk_vs_wt=1.56
fold_change_cutoff_p_9_mkk_vs_wt=1.68

# Define DE genes with these cut-offs
summary(decideTests(qlf_design2, lfc = log2(fold_change_cutoff_p_8_mkk_vs_wt), p.value = 0.05))
de_mkk2ox_vs_wt_.8 <- lrt_top_design2$table[ lrt_top_design2$table$FDR < FDR_cutoff &
                                         ( abs(as.numeric(lrt_top_design2$table$logFC)) >
                                             log2(fold_change_cutoff_p_8_mkk_vs_wt) ), ]
summary(decideTests(qlf_design2, lfc = log2(fold_change_cutoff_p_9_mkk_vs_wt), p.value = 0.05))
de_mkk2ox_vs_wt_.9 <- lrt_top_design2$table[ lrt_top_design2$table$FDR < FDR_cutoff &
                                            ( abs(as.numeric(lrt_top_design2$table$logFC)) >
                                                log2(fold_change_cutoff_p_9_mkk_vs_wt) ), ]
summary(decideTests(qlf_design2, lfc = log2(fold_change_cutoff_3.5), p.value = 0.05))
de_mkk2ox_vs_wt_3.5 <- lrt_top_design2$table[ lrt_top_design2$table$FDR < FDR_cutoff &
                                            ( abs(as.numeric(lrt_top_design2$table$logFC)) >
                                                log2(fold_change_cutoff_3.5) ), ]
head(de_mkk2ox_vs_wt_.8)
print(nrow(de_mkk2ox_vs_wt_.8))
head(de_mkk2ox_vs_wt_.9)
print(nrow(de_mkk2ox_vs_wt_.9))
head(de_mkk2ox_vs_wt_3.5)
print(nrow(de_mkk2ox_vs_wt_3.5))

##########################################################################################
### Heat map
# To visualize the expression differences among samples for the DE genes identified.
##########################################################################################
png('MKK2OX_vs_WT_heatmap_power_.8.png', width = 700,height = 8000)
par(mar= c(15,15,15,15))
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_wt_.8$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs WT (Power = 0.8)")
dev.off()
png('MKK2OX_vs_WT_heatmap_power_.9.png', width = 700,height = 6000)
par(mar= c(15,15,15,15))
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_wt_.9$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs WT (Power = 0.9)")
dev.off()
png('MKK2OX_vs_WT_heatmap_LFC_3.5.png', width = 700,height = 1500)
par(mar= c(15,15,15,15))
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_wt_3.5$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs WT (FC Cutoff = 3.5)")
dev.off()



##########################################################################################
### Volcano plots
# To visualize log fold changes for different genes.
##########################################################################################
# Code based off https://www.r-bloggers.com/2014/05/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/

head(de_mkk2ox_vs_wt_.8)
png("MKK2_OX_vs_WT_volcano_.8.png", width = 700, height = 900)
with(de_mkk2ox_vs_wt_.8, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs WT (Power = 0.8)", xlim=c(-10,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_wt_.8, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_wt_.8, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_wt_.8, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_wt_.8<- de_mkk2ox_vs_wt_.8 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_wt_.8, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()



png("MKK2_OX_vs_WT_volcano_.9.png", width = 700, height = 900)
with(de_mkk2ox_vs_wt_.9, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs WT (Power = 0.9)", xlim=c(-10,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_wt_.9, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_wt_.9, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_wt_.9, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_wt_.9<- de_mkk2ox_vs_wt_.9 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_wt_.9, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()


png("MKK2_OX_vs_WT_volcano_.3.5.png", width = 700, height = 900)
with(de_mkk2ox_vs_wt_3.5, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs WT (FC Cutoff = 3.5)", xlim=c(-10,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_wt_3.5, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_wt_3.5, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_wt_3.5, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_wt_3.5<- de_mkk2ox_vs_wt_3.5 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_wt_3.5, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()

##########################################################################################
### Tests for over-represented GO and KEGG pathway terms
# Performing gene ontology (GO) and KEGG pathway term enrichment tests on our lists of DE
# genes will give us an idea of the kind of biological processes that the genes are involved 
# in. Basically the enrichment test is a modified Fisher’s exact test (hypergeometric test)
# to determine GO or KEGG terms that are statistically over-represented in our set of DE 
# genes compared to particular background set of genes (universe). The choice of the universe
# is crucial and should include all genes that have a chance to be expressed. For our purpose, 
# the universe is selected to be the set of all genes identified as expressed. There are 
# several packages that perform these tests available on Bioconductor like GOstats or goana
# and kegga functions available from package “limma” that is already loaded as a dependency
# of “edgeR” or as web based tools like DAVID (https://david.ncifcrf.gov/) or Panther 
# (http://pantherdb.org/). We test the enrichments seperately for up- and down-regulated 
# genes.
##########################################################################################
upregulated_mkk2ox_vs_wt_ids_.8 <- de_mkk2ox_vs_wt_.8$gene_id[ de_mkk2ox_vs_wt_.8$logFC>0 ]
write.table( upregulated_mkk2ox_vs_wt_ids_.8, file="MKK2OX_vs_WT_upregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_wt_all_info_.8 <- filter(de_mkk2ox_vs_wt_.8, de_mkk2ox_vs_wt_.8$logFC>0)
write.table( upregulated_mkk2ox_vs_wt_all_info_.8, file="MKK2OX_vs_WT_upregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_wt_ids_.8 <- de_mkk2ox_vs_wt_.8$gene_id[ de_mkk2ox_vs_wt_.8$logFC<0 ]
write.table( downregulated_mkk2ox_vs_wt_ids_.8, file="MKK2OX_vs_WT_downregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_wt_all_info_.8 <- filter(de_mkk2ox_vs_wt_.8, de_mkk2ox_vs_wt_.8$logFC<0)
write.table( downregulated_mkk2ox_vs_wt_all_info_.8, file="MKK2OX_vs_WT_downregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)



upregulated_mkk2ox_vs_wt_ids_.9 <- de_mkk2ox_vs_wt_.9$gene_id[ de_mkk2ox_vs_wt_.9$logFC>0 ]
write.table( upregulated_mkk2ox_vs_wt_ids_.9, file="MKK2OX_vs_WT_upregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_wt_all_info_.9 <- filter(de_mkk2ox_vs_wt_.9, de_mkk2ox_vs_wt_.9$logFC>0)
write.table( upregulated_mkk2ox_vs_wt_all_info_.9, file="MKK2OX_vs_WT_upregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_wt_ids_.9 <- de_mkk2ox_vs_wt_.9$gene_id[ de_mkk2ox_vs_wt_.9$logFC<0 ]
write.table( downregulated_mkk2ox_vs_wt_ids_.9, file="MKK2OX_vs_WT_downregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_wt_all_info_.9 <- filter(de_mkk2ox_vs_wt_.9, de_mkk2ox_vs_wt_.9$logFC<0)
write.table( downregulated_mkk2ox_vs_wt_all_info_.9, file="MKK2OX_vs_WT_downregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)



upregulated_mkk2ox_vs_wt_ids_3.5 <- de_mkk2ox_vs_wt_3.5$gene_id[ de_mkk2ox_vs_wt_3.5$logFC>0 ]
write.table( upregulated_mkk2ox_vs_wt_ids_3.5, file="MKK2OX_vs_WT_upregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_wt_all_info_3.5 <- filter(de_mkk2ox_vs_wt_3.5, de_mkk2ox_vs_wt_3.5$logFC>0)
write.table( upregulated_mkk2ox_vs_wt_all_info_3.5, file="MKK2OX_vs_WT_upregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_wt_ids_3.5 <- de_mkk2ox_vs_wt_3.5$gene_id[ de_mkk2ox_vs_wt_3.5$logFC<0 ]
write.table( downregulated_mkk2ox_vs_wt_ids_3.5, file="MKK2OX_vs_WT_downregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_wt_all_info_3.5 <- filter(de_mkk2ox_vs_wt_3.5, de_mkk2ox_vs_wt_3.5$logFC<0)
write.table( downregulated_mkk2ox_vs_wt_all_info_3.5, file="MKK2OX_vs_WT_downregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)


# expressed_genes <- dge_expressed$genes$gene_id # same as row.names(log2cpm_expressed)
# length(expressed_genes)
# write.table( expressed_genes, file="MKK2OX_vs_WT_expressed.txt", col.names = F, row.names = F, quote=F)


# The above files are saved to rnaseq_analysis/de_analysis. Click Refresh in the Files tab on 
# the bottom right panel if you do not see the files there). Export those files in zipped format
# to your local machine by selecting them from the Files tab on the bottom right panel and 
# clicking export (under More). Choose "Save as" option to save the exported files to a defined 
# folder on your local machine.

# Use those files to continue with overrepresentation analysis on http://pantherdb.org/
# We test the enrichments seperately for up- and down-regulated genes and using the list of 
# expressed genes as the background.





### COMPARISON OF MKK2OX to CA-MKK2OX
##########################################################################################
### Statistical modelling of data and tests for DE
# Quote from the edgeR user guide: “Linear models are associated with normally distributed 
# data and have long been used to describe multifactor microarray experiments. Read counts
# however follow a negative binomial distribution. Generalized linear models (GLMs) are an 
# extension of classical linear models to non-normally distributed response data. GLMs specify
# probability distributions according to their mean-variance relationship and require variance
# parameters (dispersions) to be estimated. For general experiments (with multiple factors), 
# edgeR uses the Cox-Reid profile-adjusted likelihood (CR) method in estimating dispersions. 
# The CR method can be used to calculate a common dispersion for all the genes, trended 
# dispersion depending on the gene abundance, or separate dispersions for individual genes. 
# It takes care of multiple factors by fitting GLMs with a design matrix”.

##########################################################################################
### Define design matrix
##########################################################################################
# convert treatment information to factors
sample_info_mkk2ox_vs_camkk2ox <- filter(sample_info, (treatment == "CA-MKK2_OX" | treatment == "MKK2_OX"))
print(sample_info_mkk2ox_vs_camkk2ox)
TREATMENT <- as.factor( sample_info_mkk2ox_vs_camkk2ox$treatment )

# Define reference as control so that all DE results will be interpreted as changes in
# treatment with respect to control ( treament vs. control) samples.
TREATMENT <- relevel(TREATMENT, ref = "CA-MKK2_OX")
print(TREATMENT)

# To test for differences between treatment and control while adjusting for any differences
# between the batches, we add a term BATCH to the model defined as follows:
# convert batch information to factors
BATCH <- as.factor( sample_info_mkk2ox_vs_camkk2ox$batch )

# create design matrix without interaction term
design3 <- model.matrix( ~ TREATMENT + BATCH )
print(design3)

colnames(design3)[2] <- "MKK2_OX"
print(design3)

# subsetting norm_dge_expressed so that its treatments match the treatments in sample_info.
norm_dge_expressed
norm_dge_expressed_mkk2ox_vs_camkk2ox <- norm_dge_expressed

# filtering counts - delete columns corresponding to WT
head(norm_dge_expressed_mkk2ox_vs_camkk2ox$counts)
norm_dge_expressed_mkk2ox_vs_camkk2ox$counts<-
  norm_dge_expressed_mkk2ox_vs_camkk2ox$counts[,-c(1,2,3)]
head(norm_dge_expressed_mkk2ox_vs_camkk2ox$counts)

#filtering samples
norm_dge_expressed_mkk2ox_vs_camkk2ox$samples
row.names.remove_mkk2ox_vs_camkk2ox <- c("A5.WT", "B5.WT", "C5.WT")
norm_dge_expressed_mkk2ox_vs_camkk2ox$samples <- norm_dge_expressed_mkk2ox_vs_camkk2ox$samples[
  !(row.names(norm_dge_expressed_mkk2ox_vs_camkk2ox$samples) %in% row.names.remove_mkk2ox_vs_camkk2ox), ]
norm_dge_expressed_mkk2ox_vs_camkk2ox$samples

norm_dge_expressed_mkk2ox_vs_camkk2ox
dispersions_design3 <- estimateDisp(norm_dge_expressed_mkk2ox_vs_camkk2ox, design3, robust=TRUE)
dispersions_design3$common.dispersion


# Note from above the value of the common negative binomial dispersion. The square-root of this
# value is called the biological coefficient of variation (BCV) and its magnitude is indicative 
# of the biological variation between the replicate samples used in the experiment. From edgeR
# user guide: “Typical values for the common BCV for datasets arising from well-controlled 
# experiments are 0.4 for human data, 0.1 data on genetically identical model organisms or 0.01
# for technical replicates”.
bcv_design3 <- signif(sqrt(dispersions_design3$common.dispersion), 5)
bcv_design3

##########################################################################################
### DE contrasts: treatment versus control
##########################################################################################
# fit genewise negative binomial GLMs for design1
glmfit_design3 <- glmQLFit(dispersions_design3, design3)

# genewise statistical tests can be performed for a given coefficient or coefficient contrast
qlf_design3 <- glmQLFTest( glmfit_design3, coef = "MKK2_OX" )

# The function topTags can be called to extract the top n DE genes ranked by p-value or 
# absolute log-fold change. Save all results by setting n to NULL. Use Benjamini-Hochberg (BH)
# method to adjust p-values for multiple testing. The output is set to be sorted by ordinary
# p-values. We do not provide a specified adjusted p-value cutoff now.
lrt_top_design3 <- topTags(qlf_design3, adjust.method = "BH", sort.by = "PValue",
                           p.value = 0.05, n=NULL)
row.names(lrt_top_design3$table) <- NULL
head(lrt_top_design3$table)
dim(lrt_top_design3$table)

# summary of DE genes based only on adjusted p-value threshold
summary(decideTests(qlf_design3, p.value = 0.05, lfc = 0))

# Summary of DE genes based on adjusted p-value threshold and a 2 fold change threshold.
# Note the numbers for up and down regulated genes have reduced.
summary(decideTests(qlf_design3, p.value = 0.05, lfc = log2(2)))

# Instead of defining a hard cut-off for fold change, it would be better to select a value
# based on the power of our experiment to detect changes in expression, which would be 
# proportional to the sample size and affected by factors like biological variablity among
# samples within groups. So, performing a power analylsis is useful.

##########################################################################################
### Power analysis
# Use function rnapower from package RNASeqPower to perform a basic power analysis
##########################################################################################

# View help function by entering ?rnapower
# The average depth of coverage for the transcript is one of the parameters required for 
# rnapower function. This can be estimated as follows:
# (# of sequenced bases) / (transcriptome size)
# Sequenced bases is the number of reads x read length (twice that if paired-end).
# An estimate of the transciptome size can be obtained from the GTF file by summing the
# difference between the start and end of each annotated gene. So, depth of sequencing in the 
# sample with the lowest library size will be:
depth <- min( dge$samples$lib.size )*100*2 / sum( fusarium_gene_anno$end - fusarium_gene_anno$start )

# Estimate power for a range of fold change values between 1 and 5 with other factors fixed
power_analysis <- rnapower(depth, n=group_size, cv=bcv_design3, effect=seq(1, 5, by=0.1),
                           alpha=0.05)

png('MKK2OX_vs_CAMKK2OX_power.png', width = 800,height = 600)
par(mar=c(8,8,8,4))
plot(names(power_analysis), power_analysis, las=1,
     xlab="Fold change", ylab="Power", main="Power analysis MKK2OX vs CAMKK2OX",
     sub="[parameters: n=3; alpha=0.05; coverage=5X]")
abline(v=1.61, h=0.9, col="red")
abline(v=1.52, h=0.8, col="blue")
dev.off()

# Note from the plot that fold changes of 1.6 can be detected at a power of 0.8
# Note from the plot that fold changes of 1.75 can be detected at a power of 0.9
fold_change_cutoff_p_8_mkk2_vs_camkk2=1.52
fold_change_cutoff_p_9_mkk2_vs_camkk2=1.61

# Define DE genes with these cut-offs
summary(decideTests(qlf_design3, lfc = log2(fold_change_cutoff_p_8_mkk2_vs_camkk2), p.value = 0.05))
de_mkk2ox_vs_camkk2ox_.8 <- lrt_top_design3$table[ lrt_top_design3$table$FDR < FDR_cutoff &
                                         ( abs(as.numeric(lrt_top_design3$table$logFC)) >
                                             log2(fold_change_cutoff_p_8_mkk2_vs_camkk2) ), ]
summary(decideTests(qlf_design3, lfc = log2(fold_change_cutoff_p_9_mkk2_vs_camkk2), p.value = 0.05))
de_mkk2ox_vs_camkk2ox_.9 <- lrt_top_design3$table[ lrt_top_design3$table$FDR < FDR_cutoff &
                                                  ( abs(as.numeric(lrt_top_design3$table$logFC)) >
                                                      log2(fold_change_cutoff_p_9_mkk2_vs_camkk2) ), ]
summary(decideTests(qlf_design3, lfc = log2(fold_change_cutoff_3.5), p.value = 0.05))
de_mkk2ox_vs_camkk2ox_3.5 <- lrt_top_design3$table[ lrt_top_design3$table$FDR < FDR_cutoff &
                                                  ( abs(as.numeric(lrt_top_design3$table$logFC)) >
                                                      log2(fold_change_cutoff_3.5) ), ]
head(de_mkk2ox_vs_camkk2ox_.8)
print(nrow(de_mkk2ox_vs_camkk2ox_.8))
head(de_mkk2ox_vs_camkk2ox_.9)
print(nrow(de_mkk2ox_vs_camkk2ox_.9))
head(de_mkk2ox_vs_camkk2ox_3.5)
print(nrow(de_mkk2ox_vs_camkk2ox_3.5))


##########################################################################################
### Heat map
# To visualize the expression differences among samples for the DE genes identified.
##########################################################################################
png('MKK2OX_vs_CAMKK2OX_heatmap_.8.png', width = 600,height = 10000)
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_camkk2ox_.8$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs CAMKK2OX (Power = 0.8)")
dev.off()

png('MKK2OX_vs_CAMKK2OX_heatmap_.9.png', width = 600,height = 7000)
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_camkk2ox_.9$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs CAMKK2OX (Power = 0.9)")
dev.off()

png('MKK2OX_vs_CAMKK2OX_heatmap_3.5.png', width = 600,height = 1500)
coolmap((log2cpm_expressed[ match(de_mkk2ox_vs_camkk2ox_3.5$gene_id, row.names(log2cpm_expressed)),]),
        keysize = 1, cexRow= .9, margins = c(10,10), main="Heat Map MKK2OX vs CAMKK2OX (FC Cutoff = 3.5)")
dev.off()


##########################################################################################
### Volcano plots
# To visualize log fold changes for different genes.
##########################################################################################
# Code based off https://www.r-bloggers.com/2014/05/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/
head(de_mkk2ox_vs_camkk2ox_.8)

png("MKK2_OX_vs_CAMKK2OX_volcano_.8.png", width = 700, height = 1000)
with(de_mkk2ox_vs_camkk2ox_.8, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs CAMKK2OX (Power = 0.8)", xlim=c(-8,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_camkk2ox_.8, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_camkk2ox_.8, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_camkk2ox_.8, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_camkk2ox_.8<- de_mkk2ox_vs_camkk2ox_.8 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_camkk2ox_.8, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()


png("MKK2_OX_vs_CAMKK2OX_volcano_.9.png", width = 700, height = 1000)
with(de_mkk2ox_vs_camkk2ox_.9, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs CAMKK2OX (Power = 0.9)", xlim=c(-8,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_camkk2ox_.9, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_camkk2ox_.9, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_camkk2ox_.9, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_camkk2ox_.9<- de_mkk2ox_vs_camkk2ox_.9 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_camkk2ox_.9, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()


head(de_mkk2ox_vs_camkk2ox_3.5)

png("MKK2_OX_vs_CAMKK2OX_volcano_.3.5.png", width = 700, height = 1000)
with(de_mkk2ox_vs_camkk2ox_3.5, plot(logFC, -log10(PValue), pch=20, main="Volcano Plot MKK2OX vs CAMKK2OX (FC Cutoff = 3.5)", xlim=c(-8,10)))
# Add colored points: red if -log10(PValue)>8, orange if log2FC>2, green if both)
with(subset(de_mkk2ox_vs_camkk2ox_3.5, -log10(PValue)>7 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(de_mkk2ox_vs_camkk2ox_3.5, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(de_mkk2ox_vs_camkk2ox_3.5, -log10(PValue)>7 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="green"))

library(calibrate)
# Trying to add labels
# Shortening gene_id
de_mkk2ox_vs_camkk2ox_3.5<- de_mkk2ox_vs_camkk2ox_3.5 %>% mutate_at(~gsub("gene-", "", .), .vars = 1)
with(subset(de_mkk2ox_vs_camkk2ox_3.5, -log10(PValue)>7 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_id, cex=.6))
dev.off()

##########################################################################################
### Tests for over-represented GO and KEGG pathway terms
# Performing gene ontology (GO) and KEGG pathway term enrichment tests on our lists of DE
# genes will give us an idea of the kind of biological processes that the genes are involved 
# in. Basically the enrichment test is a modified Fisher’s exact test (hypergeometric test)
# to determine GO or KEGG terms that are statistically over-represented in our set of DE 
# genes compared to particular background set of genes (universe). The choice of the universe
# is crucial and should include all genes that have a chance to be expressed. For our purpose, 
# the universe is selected to be the set of all genes identified as expressed. There are 
# several packages that perform these tests available on Bioconductor like GOstats or goana
# and kegga functions available from package “limma” that is already loaded as a dependency
# of “edgeR” or as web based tools like DAVID (https://david.ncifcrf.gov/) or Panther 
# (http://pantherdb.org/). We test the enrichments seperately for up- and down-regulated 
# genes.
##########################################################################################

upregulated_mkk2ox_vs_camkk2ox_ids_.8 <- de_mkk2ox_vs_camkk2ox_.8$gene_id[ de_mkk2ox_vs_camkk2ox_.8$logFC>0 ]
write.table( upregulated_mkk2ox_vs_camkk2ox_ids_.8, file="MKK2OX_vs_CAMKK2OX_upregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_camkk2ox_all_info_.8 <- filter(de_mkk2ox_vs_camkk2ox_.8, de_mkk2ox_vs_camkk2ox_.8$logFC>0)
write.table( upregulated_mkk2ox_vs_camkk2ox_all_info_.8, file="MKK2OX_vs_CAMKK2OX_upregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_camkk2ox_ids_.8 <- de_mkk2ox_vs_camkk2ox_.8$gene_id[ de_mkk2ox_vs_camkk2ox_.8$logFC<0 ]
write.table( downregulated_mkk2ox_vs_camkk2ox_ids_.8, file="MKK2OX_vs_CAMKK2OX_downregulated_ids_.8.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_camkk2ox_all_info_.8 <- filter(de_mkk2ox_vs_camkk2ox_.8, de_mkk2ox_vs_camkk2ox_.8$logFC<0)
write.table( downregulated_mkk2ox_vs_camkk2ox_all_info_.8, file="MKK2OX_vs_CAMKK2OX_downregulated_all_info_.8.txt", col.names = T, row.names = F, quote=F)




upregulated_mkk2ox_vs_camkk2ox_ids_.9 <- de_mkk2ox_vs_camkk2ox_.9$gene_id[ de_mkk2ox_vs_camkk2ox_.9$logFC>0 ]
write.table( upregulated_mkk2ox_vs_camkk2ox_ids_.9, file="MKK2OX_vs_CAMKK2OX_upregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_camkk2ox_all_info_.9 <- filter(de_mkk2ox_vs_camkk2ox_.9, de_mkk2ox_vs_camkk2ox_.9$logFC>0)
write.table( upregulated_mkk2ox_vs_camkk2ox_all_info_.9, file="MKK2OX_vs_CAMKK2OX_upregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_camkk2ox_ids_.9 <- de_mkk2ox_vs_camkk2ox_.9$gene_id[ de_mkk2ox_vs_camkk2ox_.9$logFC<0 ]
write.table( downregulated_mkk2ox_vs_camkk2ox_ids_.9, file="MKK2OX_vs_CAMKK2OX_downregulated_ids_.9.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_camkk2ox_all_info_.9 <- filter(de_mkk2ox_vs_camkk2ox_.9, de_mkk2ox_vs_camkk2ox_.9$logFC<0)
write.table( downregulated_mkk2ox_vs_camkk2ox_all_info_.9, file="MKK2OX_vs_CAMKK2OX_downregulated_all_info_.9.txt", col.names = T, row.names = F, quote=F)



upregulated_mkk2ox_vs_camkk2ox_ids_3.5 <- de_mkk2ox_vs_camkk2ox_3.5$gene_id[ de_mkk2ox_vs_camkk2ox_3.5$logFC>0 ]
write.table( upregulated_mkk2ox_vs_camkk2ox_ids_3.5, file="MKK2OX_vs_CAMKK2OX_upregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
upregulated_mkk2ox_vs_camkk2ox_all_info_3.5 <- filter(de_mkk2ox_vs_camkk2ox_3.5, de_mkk2ox_vs_camkk2ox_3.5$logFC>0)
write.table( upregulated_mkk2ox_vs_camkk2ox_all_info_3.5, file="MKK2OX_vs_CAMKK2OX_upregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)

downregulated_mkk2ox_vs_camkk2ox_ids_3.5 <- de_mkk2ox_vs_camkk2ox_3.5$gene_id[ de_mkk2ox_vs_camkk2ox_3.5$logFC<0 ]
write.table( downregulated_mkk2ox_vs_camkk2ox_ids_3.5, file="MKK2OX_vs_CAMKK2OX_downregulated_ids_3.5.txt", col.names = F, row.names = F, quote=F)
downregulated_mkk2ox_vs_camkk2ox_all_info_3.5 <- filter(de_mkk2ox_vs_camkk2ox_3.5, de_mkk2ox_vs_camkk2ox_3.5$logFC<0)
write.table( downregulated_mkk2ox_vs_camkk2ox_all_info_3.5, file="MKK2OX_vs_CAMKK2OX_downregulated_all_info_3.5.txt", col.names = T, row.names = F, quote=F)



# The above files are saved to rnaseq_analysis/de_analysis. Click Refresh in the Files tab on 
# the bottom right panel if you do not see the files there). Export those files in zipped format
# to your local machine by selecting them from the Files tab on the bottom right panel and 
# clicking export (under More). Choose "Save as" option to save the exported files to a defined 
# folder on your local machine.

# Use those files to continue with overrepresentation analysis on http://pantherdb.org/
# We test the enrichments seperately for up- and down-regulated genes and using the list of 
# expressed genes as the background.


##########################################################################################
### Creating excel sheets.
##########################################################################################
install.packages("writexl")
library("writexl")
write_xlsx(downregulated_camkk2ox_vs_wt_all_info_.8,"CAMKK2OX_vs_WT_dowregulated_all_info_.8.xlsx")
write_xlsx(downregulated_camkk2ox_vs_wt_all_info_.9,"CAMKK2OX_vs_WT_dowregulated_all_info_.9.xlsx")
write_xlsx(downregulated_camkk2ox_vs_wt_all_info_3.5,"CAMKK2OX_vs_WT_dowregulated_all_info_3.5.xlsx")

write_xlsx(downregulated_mkk2ox_vs_camkk2ox_all_info_.8,"MKK2OX_vs_CAMKK2OX_downregulated_all_info_.8.xlsx")
write_xlsx(downregulated_mkk2ox_vs_camkk2ox_all_info_.9,"MKK2OX_vs_CAMKK2OX_downregulated_all_info_.9.xlsx")
write_xlsx(downregulated_mkk2ox_vs_camkk2ox_all_info_3.5,"MKK2OX_vs_CAMKK2OX_downregulated_all_info_3.5.xlsx")

write_xlsx(downregulated_mkk2ox_vs_wt_all_info_.8,"MKK2OX_vs_WT_downregulated_all_info_.8.xlsx")
write_xlsx(downregulated_mkk2ox_vs_wt_all_info_.9,"MKK2OX_vs_WT_downregulated_all_info_.9.xlsx")
write_xlsx(downregulated_mkk2ox_vs_wt_all_info_3.5,"MKK2OX_vs_WT_downregulated_all_info_3.5.xlsx")



write_xlsx(upregulated_camkk2ox_vs_wt_all_info_.8,"CAMKK2OX_vs_WT_upregulated_all_info_.8.xlsx")
write_xlsx(upregulated_camkk2ox_vs_wt_all_info_.9,"CAMKK2OX_vs_WT_upregulated_all_info_.9.xlsx")
write_xlsx(upregulated_camkk2ox_vs_wt_all_info_3.5,"CAMKK2OX_vs_WT_upregulated_all_info_3.5.xlsx")

write_xlsx(upregulated_mkk2ox_vs_camkk2ox_all_info_.8,"MKK2OX_vs_CAMKK2OX_upregulated_all_info_.8.xlsx")
write_xlsx(upregulated_mkk2ox_vs_camkk2ox_all_info_.9,"MKK2OX_vs_CAMKK2OX_upregulated_all_info_.9.xlsx")
write_xlsx(upregulated_mkk2ox_vs_camkk2ox_all_info_3.5,"MKK2OX_vs_CAMKK2OX_upregulated_all_info_3.5.xlsx")

write_xlsx(upregulated_mkk2ox_vs_wt_all_info_.8,"MKK2OX_vs_WT_upregulated_all_info_.8.xlsx")
write_xlsx(upregulated_mkk2ox_vs_wt_all_info_.9,"MKK2OX_vs_WT_upregulated_all_info_.9.xlsx")
write_xlsx(upregulated_mkk2ox_vs_wt_all_info_3.5,"MKK2OX_vs_WT_upregulated_all_info_3.5.xlsx")



##########################################################################################
### Create a document for reference
##########################################################################################
# Go to File tab above (within RStudio) and click Knit document. Choose output format HTML.
# A report called de_analysis.html is created that is saved to the de_analysis folder and 
# also opens in a browser.
# It would be useful to track the session info in the document as a refernce to software 
# used and their versions.
sessionInfo()


##########################################################################################
### Copying your data to your local machine
##########################################################################################

# In the biocluster terminal, issue the following command to change to home directory and 
# create a zipped file of the rnaseq_analysis folder:
# (rnaseq) [bctraining145@biocomp-0-3 rnaseq_analysis]$ cd
# (rnaseq) [bctraining145@biocomp-0-3 ~]$ zip -r rnaseq_analysis.zip rnaseq_analysis

# Export the rnaseq_analysis.zip file (within home) to your local machine by selecting 
# it from the Files tab on the bottom right panel and clicking export (under More). Choose
# "Save as" option to save the exported file to a defined folder on your local machine.
# Unzip the file and explore the contents in detail.

##########################################################################################
### Clean exit from biocluster and RStudio server
##########################################################################################

# Log out from the biocluster session using the following commands:
# (rnaseq) [bctraining145@biocomp-0-3 rnaseq_analysis]$ conda deactivate
# (base) [bctraining145@biocomp-0-3 rnaseq_analysis]$ logout
# (base) [bctraining145@biocluster ~]$ logout


# Signout from RStudio: click top right button (hover the pointer over the top bar to find it).
##########################################################################################
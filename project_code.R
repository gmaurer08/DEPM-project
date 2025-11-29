
################################################################################
################ DIGITAL EPIDEMIOLOGY AND PRECISION MEDICINE ###################
################ ============ PROJECT - GROUP 06 =========== ###################
################################################################################


#### Libraries ####

library(BiocGenerics) 
library(DESeq2)
library(psych) 
library(NetworkToolbox)
library(ggplot2) 
library(GGally)
library(sna)
library(network)
library(TCGAbiolinks)
library(GenomicRanges)
library(SummarizedExperiment)
library(DT)

#### Downloading the data ####------------------------------------------------------------------------------

# File directory
# setwd("C:/Users/galax/OneDrive/Dokumente/University/Magistrale/DEPM/Project")
setwd("C:/Users/lamor/OneDrive/Documents/LAURITA/Sapienza/DE")
proj <- "TCGA-BLCA" # Bladder Urothelial Carcinoma
dir.create(file.path(proj))

# Primary tumor
rna.query.tumor <- TCGAbiolinks::GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Solid tissue normal
rna.query.normal <- TCGAbiolinks::GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

# Download the data
GDCdownload(rna.query.tumor, directory = "Tumor", method = "api")
GDCdownload(rna.query.normal, directory = "Normal", method = "api")

# Directory Primary Tumor: Normal/TCGA-BLCA/Transcriptiome_Profiling/Gene_Expression_Quantification
# Directory Solid Tissue Normal: Tumor/TCGA-BLCA/Transcriptiome_Profiling/Gene_Expression_Quantification

# Prepare data
rna.data.tumor <- GDCprepare(rna.query.tumor, directory = "Tumor")
rna.expr.data.tumor <- assay(rna.data.tumor)
rna.data.normal <- GDCprepare(rna.query.normal, directory = "Normal")
rna.expr.data.normal <- assay(rna.data.normal)

# Inspect the data
genes.info.tumor <- BiocGenerics::as.data.frame(rowRanges(rna.data.tumor))
genes.info.normal <- BiocGenerics::as.data.frame(rowRanges(rna.data.normal))
all(na.omit(genes.info.normal) == na.omit(genes.info.tumor))

# Loading and analyzing clinical data
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
table(clinical.query$ajcc_pathologic_stage)
boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )

# View all the data
View(rna.expr.data.normal)
View(rna.expr.data.tumor)
View(clinical.query)

# Dimensions
dim(rna.expr.data.tumor)
dim(rna.expr.data.normal)

# Saving the data
save(rna.query.normal, rna.query.tumor, rna.data.normal, rna.data.tumor, rna.expr.data.normal, rna.expr.data.tumor,
genes.info.tumor, genes.info.normal, genes.info.tumor, genes.info.normal, clinical.query, file="initial-project-data.RData")


#### 1. Pre-process the data  --------------------------------------------------------------------------------------

# 1.1 Data Pre-processing & Matching

# participant ids (first 12 characters of the barcode) to match patients
tumor.patients <- substr(colnames(rna.expr.data.tumor), 1, 12)
normal.patients <- substr(colnames(rna.expr.data.normal), 1, 12)

# the intersection (patients that have BOTH Tumor and Normal data) 
common.patients <- intersect(tumor.patients, normal.patients)
message("Number of matched patients: ", length(common.patients))
cat("Only 19 patiens matched") #lets validate this with the teacher

# Identify columns corresponding to these common patients
tumor.idx <- match(common.patients, tumor.patients)
normal.idx <- match(common.patients, normal.patients)

# Subset the expression matrices to keep only matched patients
# Note: We must ensure the order of columns matches for paired analysis later
matched.tumor.data <- rna.expr.data.tumor[, tumor.idx]
matched.normal.data <- rna.expr.data.normal[, normal.idx]

# Check if barcodes align 
print(head(colnames(matched.tumor.data))) #all good
print(head(colnames(matched.normal.data))) #all good

# Combine into a single matrix for processing
count.matrix <- cbind(matched.tumor.data, matched.normal.data) #Columns: First all Tumors, then all Normals

# 1.2 Filtering Zeros 

# erase the genes (rows) with at least one zero value 
genes.to.keep <- rowSums(count.matrix == 0) == 0 # check which rows have ANY zero in the matched dataset
filtered.counts <- count.matrix[genes.to.keep, ]

message("Original gene count: ", nrow(count.matrix))
message("Genes remaining after removing rows with zeros: ", nrow(filtered.counts))
cat("the original genes count was 60.660, and after doing the cleaning of data, we remained with 17.223, which is 28% of the original data")

# Update the separated matrices with the filtered genes
tumor.final <- matched.tumor.data[genes.to.keep, ] 
normal.final <- matched.normal.data[genes.to.keep, ]


#### 2. DIFFERENTIAL EXPRESSION ANALYSIS (DEGs) -----------------------------------------------------------------

# 2.1 Metadata (Experimental Design)
# 'n' Tumor samples followed by 'n' Normal samples in 'filtered.counts'
# n = length(common.patients) will be the 19 from the results above
n_samples <- length(common.patients)

col_data <- data.frame(condition = factor(rep(c("Tumor", "Normal"), each = n_samples)), row.names = colnames(filtered.counts)
)

all(rownames(col_data) == colnames(filtered.counts))# Check if columns match. PD: It's true, so we can continue

# 2.2 Create DESeqDataSet Object
# DESeq2 requires raw counts (integers), which we have from "STAR - Counts"
dds <- DESeqDataSetFromMatrix(
  countData = filtered.counts,
  colData = col_data,
  design = ~ condition
)

# Set "Normal" as the reference level so "Tumor" is compared TO "Normal"
# (Positive FC means Upregulated in Tumor)
dds$condition <- relevel(dds$condition, ref = "Normal")

# 2.3 Run DESeq Analysis
# This performs normalization, dispersion estimation, and the Wald test
dds <- DESeq(dds)

# 2.4 Extract Results
# We use alpha = 0.05 for p-value optimization
res <- results(dds, alpha = 0.05)
summary(res)

cat("LFC > 0 (up) shows that 3,775 genes have a positive Log Fold Change. So,these genes are more active (over-expressed) in the Tumor samples compared to normal tissue.")
cat("LFC < 0 (down) shows that 3,242 genes have a negative Log Fold Change. So,those genes are less active (silenced or suppressed) in the Tumor samples compared to normal tissue.")
cat("outliers [1] is 0, So, the data is clean and no single patient is messing up the results, ouliers[2] is also 0, which makes sense because we deleted the rows (genes) with 0")

cat("I'm still worried that all those genes are a lot and we need to filter them more - Laura")

# 2.3. Filtering DEGs (Applying Thresholds)

# Project Guideline
# 1. Adjusted P-value (padj) <= 0.05
# 2. |Fold Change| >= 1.2
#    Note: DESeq2 outputs log2(FoldChange).
#    log2(1.2) approx 0.263
#    log2(1/1.2) approx -0.263

fc_threshold_linear <- 5 #iterating it to get hundreds of genes and not thousands
log2_fc_threshold <- log2(fc_threshold_linear)
padj_threshold <- 0.05

# Create logical vector for significance
is_significant <- (res$padj <= padj_threshold & !is.na(res$padj)) & (abs(res$log2FoldChange) >= log2_fc_threshold)

# Extract the significant genes
degs_list <- res[is_significant, ]

# Order them by significance (lowest p-value at top)
degs_list <- degs_list[order(degs_list$padj), ]

cat("For thresholds: padj <= ", padj_threshold, " | |FC| >= ", fc_threshold_linear)
cat("we, got ",  nrow(degs_list), "number of DEGs identified.")

# Save the names of the DEGs for the next tasks (Co-expression networks)
deg_names <- rownames(degs_list)

#2.4 Visualization: Volcano Plot 


volcano_data <- as.data.frame(res) #the results to a data frame for plotting
volcano_data$diffexpressed <- "No"
# Upregulated (green)
volcano_data$diffexpressed[volcano_data$log2FoldChange > log2_fc_threshold & volcano_data$padj < padj_threshold] <- "UP"
# Downregulated (Blue)
volcano_data$diffexpressed[volcano_data$log2FoldChange < -log2_fc_threshold & volcano_data$padj < padj_threshold] <- "DOWN"

# plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "green")) +
  geom_vline(xintercept = c(-log2_fc_threshold, log2_fc_threshold), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), col = "black", linetype = "dashed") +
  labs(title = "Tumor vs Normal (TCGA-BLCA)",
       subtitle = paste("Thresholds: |FC| >", fc_threshold_linear, "& adj.p <", padj_threshold),
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.title = element_blank())

cat("Downregulated (blue) seems denser and has more points than the upregulated")
cat("it seems that loss of gene expression is a dominant feature of this bladder cancer dataset when looking at high-fold changes. Mainly, there are many genes being turned off.")

##Save DEG Results ####
write.csv(as.data.frame(degs_list), file = "DEGs_List.csv")


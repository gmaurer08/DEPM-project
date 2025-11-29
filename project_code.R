
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
library(biomaRt) #new from laura

#### Downloading the data ####------------------------------------------------------------------------------

# File directory
# setwd("C:/Users/galax/OneDrive/Dokumente/University/Magistrale/DEPM/Project")
#setwd("C:/Users/lamor/OneDrive/Documents/LAURITA/Sapienza/DE")-->LAURA
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
# 2. |Fold Change| >= 5
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

# Save the names of the DEGs for the next tasks 
deg_names <- rownames(degs_list) #this is for the next steps
message("Here are the DEG names:", deg_names)

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

##Save DEG Results
write.csv(as.data.frame(degs_list), file = "DEGs_List.csv")

# Subset the original data to keep the significant DEGs :D
network.tumor.data <- tumor.final[deg_names, ]
network.normal.data <- normal.final[deg_names, ]

### 3. Co-expression networks -----------------------------------------------------------------------

#3.1 Before staring, we transform the data using log2.
tumor.log <- log2(network.tumor.data + 1) # Adding +1 avoids log(0) errors
normal.log <- log2(network.normal.data + 1)

#3.2 Correlation Matrices (Pearson correlation)
# We use t() because cor() calculates correlation between columns, and our genes are rows.
tumor.cor <- cor(t(tumor.log), method = "pearson")
normal.cor <- cor(t(normal.log), method = "pearson")

# 3.3 Build Binary Adjacency Matrices
# "Binary adjacency matrix where a_ij=0 if |rho| < threshold" 

# Define Correlation Threshold
cor_thresh <- 0.75 #standar measure

# Create Adjacency Matrices (1 if connected, 0 if not)
# We use abs() because strong negative correlation is also a valid biological link.
tumor.adj <- ifelse(abs(tumor.cor) > cor_thresh, 1, 0)
normal.adj <- ifelse(abs(normal.cor) > cor_thresh, 1, 0)

# Remove self-loops (a gene is always correlated 1.0 with itself, but that's not a network link)
diag(tumor.adj) <- 0
diag(normal.adj) <- 0


cat("Tumor Network Density: ", round(mean(tumor.adj), 3))
cat("Normal Network Density: ", round(mean(normal.adj), 3))



# 3.4 Scale-Free Property and Degree ------------------------------------------

# "Compute the degree index and check if the network is a scale free network" [cite: 32]

# Calculate Degree (k) = Number of connections per gene
k.tumor <- rowSums(tumor.adj)
k.normal <- rowSums(normal.adj)

# Function to Plot Degree Distribution (Log-Log Plot)
# A straight line descending implies a Scale-Free Network (Power Law)
plot_scale_free <- function(degrees, title) {
  # Frequency table of degrees
  deg_counts <- table(degrees)
  
  # X axis: Degree (k), Y axis: Probability P(k)
  x <- as.numeric(names(deg_counts))
  y <- as.numeric(deg_counts) / sum(deg_counts)
  
  # Remove degree 0 to avoid log(0) errors
  valid <- x > 0
  x <- x[valid]
  y <- y[valid]
  
  # Plot
  plot(log10(x), log10(y), pch = 19, col = "darkblue",
       main = title,
       xlab = "Log10 (Degree k)", 
       ylab = "Log10 (P(k))")
  
  # Fit a linear model to check R-squared (goodness of fit)
  fit <- lm(log10(y) ~ log10(x))
  abline(fit, col = "red", lwd = 2)
  
  # Display R-squared on plot
  r2 <- round(summary(fit)$r.squared, 3)
  legend("topright", legend = paste("R2 =", r2), bty = "n")
}

# Generate Plots (Required for Report "Expected figures: Degree distribution" )
par(mfrow = c(1, 2)) # Side by side
plot_scale_free(k.tumor, "Tumor Network Degree Dist.")
plot_scale_free(k.normal, "Normal Network Degree Dist.")
par(mfrow = c(1, 1)) # Reset


cat("The Tumor Network is a scale-free Network. We can see that the tumor network is organized hierarchically. Most genes have very few connections, but a few Hubs control everything. This structure is efficient but fragile if we target the hubs (with drugs).")
cat("On the other hand, the Normal Network shows a cluster of points on the far right (high degree), and that breaks the pattern. So, the normal network is not scale-free under these parameters. It is much denser. In healthy tissue, gene regulation is often tighter and more redundant, leading to many genes having high connectivity, rather than just a few select hubs.")


# 3.5 Hub Identification & Comparison

# Calculate cutoff for top 5% (5% of the nodes with highest degree values)
top_5_pct <- ceiling(length(deg_names) * 0.05)

# sort genes by degree and pick top N
hubs.tumor <- names(sort(k.tumor, decreasing = TRUE)[1:top_5_pct]) #descending order
hubs.normal <- names(sort(k.normal, decreasing = TRUE)[1:top_5_pct])

# Compare hubs sets related to the two condition.
common_hubs <- intersect(hubs.tumor, hubs.normal)
unique_tumor_hubs <- setdiff(hubs.tumor, hubs.normal)
unique_normal_hubs <- setdiff(hubs.normal, hubs.tumor)

cat("Total Hubs per network (top 5%): ", top_5_pct)
cat("Shared Hubs: ", length(common_hubs))
cat("Tumor-Specific Hubs: ", length(unique_tumor_hubs))
cat("Normal-Specific Hubs: ", length(unique_normal_hubs))
 
cat("There is a huge 'rewiring' of the network. Only 4 out of 36 'Hub' genes remained central in both healthy and cancer tissues.")
cat("In the normal tissue, 32 specific hubs were responsible for maintaining healthy function (homeostasis). In the cancer tissue, these genes lost their central connectivity, suggesting a collapse of the normal regulatory systems.")
cat("On the contrary, 32 new hubs emerged exclusively in the tumor network. These genes, which were peripheral or dormant in healthy tissue, appear to have 'hijacked' the network to drive the pathological state of the cancer.")
cat("Conclusion: Bladder Cancer is not just caused by individual genes changing, but by a complete restructuring of the gene communication network.")


#write.csv(data.frame(Gene=unique_tumor_hubs), "Hubs_Specific_to_Tumor.csv") 
#write.csv(data.frame(Gene=unique_normal_hubs), "Hubs_Specific_to_Normal.csv")

# Print first few tumor-specific hubs to check
print(unique_tumor_hubs) #these are the IDS but we want the Gene Symbols

#3.6. IDs to gene symbols 

clean_ids <- gsub("\\..*", "", unique_tumor_hubs)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # connect to the Ensembl Database

#Fetch the Gene Symbols
#ask for the ID and the "hgnc_symbol" (the human readable name)
gene_names <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = clean_ids,
  mart = mart
)

gene_names <- gene_names[gene_names$hgnc_symbol != "", ] #filter empty names

print(gene_names)
cat("The names of the tumor genes are: ")
print(gene_names$hgnc_symbol)

#write.csv(gene_names, "Tumor_Hub_Symbols.csv", row.names = FALSE) #save

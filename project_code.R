
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
library(network) #new
library(mclust) #5C
library(SNFtool) #5C
library(igraph) #5A


#### Downloading the data ####------------------------------------------------------------------------------

# File directory
#setwd("C:/Users/galax/OneDrive/Dokumente/University/Magistrale/DEPM/Project")
#setwd("C:/Users/lamor/OneDrive/Documents/LAURITA/Sapienza/DE")
setwd("/Users/home/Desktop")
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
#save(rna.query.normal, rna.query.tumor, rna.data.normal, rna.data.tumor, rna.expr.data.normal, rna.expr.data.tumor,
#genes.info.tumor, genes.info.normal, genes.info.tumor, genes.info.normal, clinical.query, file="initial-project-data.RData")


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


# 2.3. Filtering DEGs (Applying Thresholds)

# Project Guideline
# 1. Adjusted P-value (padj) <= 0.05
# 2. |Fold Change| >= 5
#    Note: DESeq2 outputs log2(FoldChange).
#    log2(1.2) approx 0.263
#    log2(1/1.2) approx -0.263

fc_threshold_linear <- 1.4 #iterating it to get hundreds of genes and not thousands
log2_fc_threshold <- log2(fc_threshold_linear)
padj_threshold <- 1e-8

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


cat("Downregulated seems denser and has more points than the upregulated")
cat("A loss of gene expression (genes being silenced or suppressed) is a dominant and characteristic feature of this specific cancer subtype compared to normal tissue, since the tails are not symetric")

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


# 3.5 Hub Identification & Comparison------------------------

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
 
cat("There is a huge 'rewiring' of the network. Only 4 out of 34 'Hub' genes remained central in both healthy and cancer tissues.")
cat("In the normal tissue, 30 specific hubs were responsible for maintaining healthy function (homeostasis). In the cancer tissue, these genes lost their central connectivity, suggesting a collapse of the normal regulatory systems.")
cat("On the contrary, 30 new hubs emerged exclusively in the tumor network. These genes, which were peripheral or dormant in healthy tissue, appear to have 'hijacked' the network to drive the pathological state of the cancer.")
cat("Conclusion: Bladder Cancer is not just caused by individual genes changing, but by a complete restructuring of the gene communication network.")


#write.csv(data.frame(Gene=unique_tumor_hubs), "Hubs_Specific_to_Tumor.csv") 
#write.csv(data.frame(Gene=unique_normal_hubs), "Hubs_Specific_to_Normal.csv")

# Print first few tumor-specific hubs to check
print(unique_tumor_hubs) #these are the IDS but we want the Gene Symbols

#3.6. IDs to gene symbols ---------------------------

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

# 3.7 Tumor Network view ----------------------------

# Create the 'network' object 'tumor.adj' is the binary adjacency matrix from Task 3.3
net.tumor <- network(tumor.adj, 
                     matrix.type = "adjacency", 
                     directed = FALSE)

# 2. attributes for plotting 
net.tumor %v% "type" = ifelse(network.vertex.names(net.tumor) %in% hubs.tumor, "Hub", "Non-Hub") # type (Hub/Non-Hub) using the hubs list from point 3.5
net.tumor %v% "color" = ifelse(net.tumor %v% "type" == "Hub", "tomato", "blue") #color based on type

# Tumor network plot
tumor_network_plot <- ggnet2(net.tumor, 
                             color = "type",                
                             palette = c("Hub" = "tomato", "Non-Hub" = "blue"),
                             size = "degree",                 
                             size.cut = 5,                    
                             node.alpha = 0.8,
                             edge.size = 0.1,                 
                             edge.alpha = 0.4,
                             mode = "fruchtermanreingold",    
                             legend.size = 10) +
  guides(size = "none") +
  labs(title = "Tumor co-expression Network (DEGs)",
       subtitle = paste("Hubs (5%) highlighted. Correlation threshold: |rho| >", cor_thresh)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Display and Save the plot
print(tumor_network_plot)
ggsave("Tumor_Coexpression_Network.png", plot = tumor_network_plot, width = 8, height = 8, dpi = 300)


# 3.8 Normal network view --------------------

# network' object from the binary adjacency matrix from Task 3.3
net.normal <- network(normal.adj, 
                      matrix.type = "adjacency", 
                      directed = FALSE)

#assign attributes
net.normal %v% "type" = ifelse(network.vertex.names(net.normal) %in% hubs.normal, "Hub", "Non-Hub") # 'normal.adj' is the binary adjacency matrix from point 3.3
net.normal %v% "color" = ifelse(net.normal %v% "type" == "Hub", "tomato", "green") #  color based on type

# plot the normal network
normal_network_plot <- ggnet2(net.tumor, 
                             color = "type",                 
                             palette = c("Hub" = "tomato", "Non-Hub" = "green"), 
                             size = "degree",                 
                             size.cut = 5,                    
                             node.alpha = 0.8,
                             edge.size = 0.1,                 
                             edge.alpha = 0.4,
                             mode = "fruchtermanreingold",    
                             legend.size = 10) +
  guides(size = "none") +
  labs(title = "Normal Co-expression Network (DEGs)",
       subtitle = paste("Hubs (5%) highlighted. Correlation threshold: |rho| >", cor_thresh)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Display and Save the plot
print(normal_network_plot)
ggsave("normal_Coexpression_Network.png", plot = normal_network_plot, width = 8, height = 8, dpi = 300)


# 3.9 Plotting the Hub Subnetwork (Tumor) -------------------------------------

# We used these variables: tumor.adj (binary matrix), hubs.tumor (hub names), deg_names (all gene names)

# rows/columns corresponding to Hubs
hubs.tumor.indices <- which(deg_names %in% hubs.tumor)
hubs.tumor.names <- deg_names[hubs.tumor.indices]

# neighbors (nodes connected to any hub in the binary matrix)
# rowSums checks which columns in the tumor.adj matrix have a connection (a '1') 
# when restricted to the hub rows.
is_neighbor <- colSums(tumor.adj[hubs.tumor.indices, ]) > 0 
neighbor.names <- deg_names[is_neighbor]

# final Subnetwork definition: Hubs + Neighbors
subnet.nodes <- unique(c(hubs.tumor.names, neighbor.names))

# Sub-Adjacency Matrix
hub.tumor.adj <- tumor.adj[subnet.nodes, subnet.nodes]

#Sub-Network Object
net.hub.tumor <- network(hub.tumor.adj, 
                         matrix.type="adjacency", 
                         directed = FALSE)

#assign attributes
net.hub.tumor %v% "type" = ifelse(network.vertex.names(net.hub.tumor) %in% hubs.tumor.names, "Hub", "Neighbor")
net.hub.tumor %v% "color" = ifelse(net.hub.tumor %v% "type" == "Hub", "tomato", "blue")

#Subnetwork plot
hub_subnetwork_plot <- ggnet2(net.hub.tumor,  
                              color = "type",                 
                              palette = c("Hub" = "tomato", "Neighbor" = "blue"),
                              size = "degree",                 
                              size.cut = 5,
                              alpha = 0.9,
                              edge.color = "grey", 
                              edge.alpha = 0.6,
                              edge.size = 0.15, 
                              label.color = "black", 
                              label.size = 2.5) +
  guides(size = "none") +
  labs(title = "Tumor Hub Subnetwork and First Neighbors", 
       subtitle = paste("Showing connectivity of", length(hubs.tumor.names), "Tumor-specific Hubs"))

print(hub_subnetwork_plot)
ggsave("Tumor_Hub_Subnetwork.png", plot = hub_subnetwork_plot, width = 8, height = 8, dpi = 300)




#### 4. DIFFERENTIAL CO-EXPRESSION NETWORK -----------------------------------------------


# 4.1 Computing Fisher's Z-transformation and Z-scores

# The goal is to find gene pairs whose correlation changes significantly 
# between tumor and normal conditions

# Fisher's Z-transformation: 0.5 * ln((1+rho)/(1-rho))   (rho: Pearson correlation)
# This make correlations approximately normally distributed

fisher_z <- function(rho) {
  return(0.5*log((1+rho)/(1-rho)))
}

# Transform both correlation matrices
z.tumor <- fisher_z(tumor.cor)
z.normal <- fisher_z(normal.cor)

# Handle extreme cases (when correlation is 1 or -1)
z.tumor[is.infinite(z.tumor)] <- NA
z.normal[is.infinite(z.normal)] <- NA

# For each gene pair compute how much the correlation changed
# Z = (Z_tumor - Z_normal) / sqrt(1/(n_tumor-3) + 1/(n_normal-3))

n_tumor <- ncol(tumor.log)  # number of tumor samples
n_normal <- ncol(normal.log)  # number of normal samples

# Calculate differential Z-score matrix
z_diff <- (z.tumor - z.normal) / sqrt(1/(n_tumor-3) + 1/(n_normal-3))
z_diff[is.na(z_diff)] <- 0 # Replace NAs with 0 (no significant change)
cat("Range of Z_diff scores:", round(min(z_diff), 2), "to", round(max(z_diff), 2), "\n")

# 4.2 Build binary differential adjacency matrix

z_threshold <- 3 # from the assignment, |Z_diff|>3 means that correlation change is statistically significant

# Create adjacency matrix: 1 if connection changed significantly, 0 otherwise
diff.adj <- ifelse(abs(z_diff) > z_threshold, 1, 0)
diag(diff.adj) <- 0 # Remove self-loops

cat("Threshold: |Z_diff|>", z_threshold, "\n")
cat("Total possible edges:", length(deg_names) * (length(deg_names)-1)/2, "\n")
cat("Significant differential edges:", sum(diff.adj)/2, "\n")  # divide by 2 because matrix is symmetric
cat("Network Density:", round(mean(diff.adj), 4), "\n\n")

# 4.3 Degree distribution and scale-free analysis

k.diff <- rowSums(diff.adj) # Calculate degree for each gene in the differential network

# Printing the statistics
cat("Degree Distribution Statistics\n")
cat("Min degree:", min(k.diff),"    Max degree:", max(k.diff), "\n")
cat("Mean degree:", round(mean(k.diff),2), "     Median degree:", median(k.diff), "\n")
cat("Genes with degree 0:", sum(k.diff==0), "\n\n")

# Plot degree distribution (log-log plot for scale-free assessment)
plot_scale_free_diff <- function(degrees, title) {
  
  # Create frequency table
  deg_counts <- table(degrees)
  
  # Prepare data
  x <- as.numeric(names(deg_counts))
  y <- as.numeric(deg_counts)/sum(deg_counts)
  
  # remove degree 0 to avoid log(0)
  valid <- x>0
  x<-x[valid]
  y<-y[valid]
  
  # check if we have enough data points
  if(length(x)<3) {
    plot(1, 1, type="n", main=title, xlab="Log10 (Degree k)", ylab="Log10 (P(k))")
    text(1, 1, "Insufficient data for scale-free analysis", col="red")
    return(NULL)
  }
  
  # log-log plot
  plot(log10(x), log10(y), pch=19, col="purple", main=title,
       xlab="Log10 (Degree k)", ylab="Log10 (P(k))", cex=1.5)
  
  # Fit linear model
  fit <- lm(log10(y) ~ log10(x))
  abline(fit, col="red", lwd=2)
  
  # Display statistics
  r2 <- round(summary(fit)$r.squared, 3)
  slope <- round(coef(fit)[2], 3)
  
  legend("topright", legend=c(paste("R^2 =", r2), paste("Slope =", slope)), bty="n", cex=1.1)
  
  return(list(r2=r2, slope=slope))
}

# Generate the plot
par(mfrow=c(1, 1))
diff_stats <- plot_scale_free_diff(k.diff, "Differential Network Degree Distribution")

# Interpretation
if(!is.null(diff_stats)) {
  if(diff_stats$r2 > 0.8) {
    cat("The differential network appears to be scale-free (R^2>0.8)\n")
    # this indicates a hierarchical organization with hub genes
  } else {
    cat("The differential network does NOT follow a clear scale-free topology (R^2<0.8)\n")
    # this may indicate a more distributed network structure
  }
}

# 4.4 Hub Identification in Differential Network

# Calculate top 5% threshold
top_5_pct_diff <- ceiling(length(deg_names)*0.05)
# Identify hubs (top 5% by degree)
hubs.diff <- names(sort(k.diff, decreasing=TRUE)[1:top_5_pct_diff])

# Print top hubs
cat("Top 5% threshold:", top_5_pct_diff, "genes\n")

# Display top 10 hubs with their degrees
cat("Top 10 Differential Hubs:\n")
top_diff_degrees <- sort(k.diff, decreasing = TRUE)[1:10]
print(data.frame(Gene=names(top_diff_degrees), Degree=as.numeric(top_diff_degrees)))

# 4.5 Comparison with Task 3 Hubs

# Compare with tumor hubs from Task 3
overlap_tumor <- intersect(hubs.diff, hubs.tumor)
unique_to_diff <- setdiff(hubs.diff, hubs.tumor)
unique_to_tumor <- setdiff(hubs.tumor, hubs.diff)

cat("Shared hubs (differential and tumor-specific):", length(overlap_tumor), "\n")
cat("Unique to differential network:", length(unique_to_diff), "\n")
cat("Unique to tumor network (Task 3):", length(unique_to_tumor), "\n\n")

# Compare with normal hubs from Task 3
overlap_normal <- intersect(hubs.diff, hubs.normal)
unique_to_normal <- setdiff(hubs.normal, hubs.diff)

cat("Shared hubs (differential and normal-specific):", length(overlap_normal), "\n")
cat("Unique to normal network (Task 3):", length(unique_to_normal), "\n\n")

# Compare with common hubs from Task 3
overlap_common <- intersect(hubs.diff, common_hubs)

cat("Shared hubs (differential and common hubs of both networks):", length(overlap_common), "\n\n")

# Create summary visualization
cat("SUMMARY TABLE\n")
comparison_summary <- data.frame(
  Hub_Set = c("Tumor-Specific (Task 3)", "Normal-Specific (Task 3)", "Common (Task 3)", "Differential (Task 4)"),
  Count = c(length(hubs.tumor), length(hubs.normal), length(common_hubs), length(hubs.diff)),
  Overlap_with_Diff = c(length(overlap_tumor), length(overlap_normal), length(overlap_common), NA))
print(comparison_summary)

# 4.6 Convert Differential Hub IDs to Gene Symbols

# Clean Ensembl IDs (remove version numbers)
clean_diff_ids <- gsub("\\..*", "", hubs.diff)

# Fetch gene symbols using biomaRt
gene_names_diff <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = clean_diff_ids,
  mart = mart
)

# Filter out empty symbols
gene_names_diff <- gene_names_diff[gene_names_diff$hgnc_symbol != "", ]

cat("Differential Hub Gene Symbols\n")
print(gene_names_diff[, c("ensembl_gene_id", "hgnc_symbol")])
cat("\nGene symbols:", paste(gene_names_diff$hgnc_symbol, collapse = ", "), "\n")

# Save results
write.csv(gene_names_diff, "Differential_Hub_Symbols.csv", row.names=FALSE)
write.csv(data.frame(Gene=hubs.diff, Degree=k.diff[hubs.diff]), 
          "Hubs_Differential_Network.csv", row.names=FALSE)

# Separate: Enrichr analysis


#### 5.A PATIENT SIMILARITY NETWORK (PSN)  -----------------------------------------------

## From section 3 we already have the log-transformed genes x samples matrix
## tumor data: tumor.log - tumor.log <- log2(network.tumor.data + 1); network.tumor.data <- tumor.final[deg_names, ]; tumor.final <- matched.tumor.data[genes.to.keep, ] 
## normal data: normal.log - similarly to the above 

## Build Patient–Patient similarity matrix (PSN)
## We compute sample-sample similarity 
psn.cor <- cor(tumor.log, method = "pearson")

# Pre-Louvain preprocessing: Louvain expects non-negative weights
psn.cor[psn.cor < 0] <- 0   # ensure that possible negative similarities would be set to 0

## Sparsification of the network: via thresholding weak similarities 
sim_thresh <- 0.86    # 0.86 seems to be the most reasonable threshold, allowing for 2 well-seprated communities and 3 disconnected nodes
psn.cor.thresh <- psn.cor
psn.cor.thresh[psn.cor.thresh < sim_thresh] <- 0

## Build igraph PSN object
g.psn <- graph_from_adjacency_matrix(
  psn.cor.thresh,
  mode    = "undirected",
  weighted = TRUE,
  diag    = FALSE
)

# Quick inspection of the graph
cat("Patient Similarity Network (Tumor):")
g.psn
cat("Summary of Patient Similarity Network (Tumor):")
summary(g.psn)

V(g.psn)$name <- substr(V(g.psn)$name, 1, 12) # shorten the patient names

# Vizualizations 
# A simple one of the graph after thresholding
set.seed(42) # reproducibility for the layout

# Fruchterman-Reingold force-directed layout: It tries to place similar/connected nodes close together and reduce edge crossings.
coords <- layout_with_fr(g.psn, weights = E(g.psn)$weight)

png("psn_tumor_simple.png", width = 900, height = 900, res = 100)

plot(
  g.psn,
  layout = coords,
  vertex.size = 8,
  vertex.label = NA,  # hide labels for clarity
  vertex.color = "skyblue",
  edge.width = E(g.psn)$weight * 2, # scale by weight
  edge.color = rgb(0, 0, 0, 0.2),   # semi-transparent edges
  main = "Patient Similarity Network \n (Tumor)",
  frame        = FALSE,   # no box
  margin       = 0        # igraph's own margin around the layout
)
dev.off()

#### 5.B LOUVAIN COMMUNITY DETECTION  -----------------------------------------------
set.seed(42) # reproducibility, to fix Louvain randomness
louvain.res <- cluster_louvain(g.psn, weights = E(g.psn)$weight)

# Basic info
cat("Louvain communities:")
print(louvain.res)

membership_vec <- membership(louvain.res)
table(membership_vec)   # community sizes
cat("As we can observe there are two prominent communities and 3 singletons.")

modularity(louvain.res) # modularity score
cat("Given our network size, modularity of 0.269 indicate moderate community structure")

# Set community membership as node attributes
V(g.psn)$community <- membership_vec

# Community membership
comm <- V(g.psn)$community

# Colorblind-safe pallette - cyan and yellow used for the biggest communities 
paul_tol_bright <- c(
  "#66CCEE", # cyan
  "#CCBB44", # yellow
  "#228833", # green
  "#EE6677", # red/pink
  "#AA3377", # purple
  "#BBBBBB"  # grey
)

# Assign 1 color per community (in numeric order)
comm_ids <- sort(unique(comm))
pal <- setNames(paul_tol_bright[seq_along(comm_ids)], comm_ids)

V(g.psn)$color <- pal[as.character(comm)]


# Simple plot of the communities
png("psn_tumor_communities.png", width = 1200, height = 900, res = 150)
layout_psn <- layout_with_fr(g.psn)

plot(
  g.psn,
  layout = layout_psn,
  vertex.size  = 10,
  vertex.label = NA,
  vertex.color = V(g.psn)$color,
  edge.width = E(g.psn)$weight * 2,
  main = "Patient Similarity Network (Tumor, DEGs)\nLouvain Community Structure"
)
dev.off()

# Enhanced vizualization
tg <- as_tbl_graph(g.psn)

set.seed(42)
louvain_plot <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.2) +
  geom_node_point(aes(color = factor(community)), size = 4) +
  scale_color_manual(values = pal, name = "Community") +
  scale_edge_width(range = c(0.1, 2)) +
  theme_void() +
  ggtitle("Patient Similarity Network (Tumor, DEGs) \n Louvain Communities")

louvain_plot <- louvain_plot +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

ggsave("psn_tumor_louvain.png", louvain_plot,
       width = 6, height = 4, dpi = 300,
       bg = "white")

# Barplot to highlight the community size
comm_sizes <- sizes(louvain.res)

png("psn_tumor_barplot.png", width = 1200, height = 900, res = 150)
par(mar = c(5, 6, 4, 2) + 0.1)  
barplot(
  comm_sizes,
  col = pal[names(comm_sizes)],
  main = "Community Sizes (Louvain on PSN)",
  xlab = "Community",
  ylab = "Number of Patients",
  ylim = c(0, max(comm_sizes) * 1.15),
  cex.main = 1.7,
  cex.lab = 1.5, 
  cex.axis = 1.1   
)

dev.off()

# Vizualization with the patient codes
tg <- as_tbl_graph(g.psn)

tg <- tg %>%
  activate(nodes) %>%
  mutate(label = name)

set.seed(42)
louvain_patients_plot <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.2) +
  geom_node_point(aes(color = factor(community)), size = 4) +
  geom_node_text(
    aes(label = label, color = "darkgray"),
    repel = FALSE,
    size = 3,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal, name = "Community") +
  scale_edge_width(range = c(0.1, 2)) +
  theme_void() +
  ggtitle("Patient Similarity Network \n (Louvain Communities)") +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    plot.title  = element_text(hjust = 0.5, margin = margin(b = 10)),
    legend.margin = margin(10, 10, 10, 10)
  )
ggsave("psn_tumor_louvain_patients.png", louvain_patients_plot,
       width = 6, height = 4, dpi = 300,
       bg = "white")

#### 5.C SIMILARITY NETWORK FUSION -----------------------------------------------

# Prepare input matrices
# Expression matrix for SNF (samples x features)
# tumor.log: genes x samples 
tumor.log.short <- tumor.log 
colnames(tumor.log.short) <- substr(colnames(tumor.log), 1, 12)  
pat_ids <- colnames(tumor.log.short)
expr.mat <- t(tumor.log.short)  # now: rows = patients, columns = DEG expression features

#length(intersect(substr(colnames(tumor.log.short),1,12), pat_ids)) # verify if the ids match
#rownames(expr.mat) <- pat_ids  # shortens the row names into just the patient ids

# Mutation matrix for the same patients
# Built a mutation matrix 'mut.data' with: rows = mutation features, columns = tumor patients
query.mut <- GDCquery(
  project = proj, 
  data.category = "Simple Nucleotide Variation", 
  data.type = "Masked Somatic Mutation"
)

GDCdownload(query.mut)
mut.raw <- GDCprepare(query.mut)   # long MAF-like table; rows: genes, cols: features; genes for one patient are grouped and patients are concatanated one by one
mut.maf <- mut.raw   

mut.maf$Tumor_Sample_Barcode <- substr(mut.maf$Tumor_Sample_Barcode, 1, 12) # modify to the patient-id like before

# Verify patients from our matched list that are not present in the mut.maf - no mutational data avail.
missing_in_maf <- setdiff(pat_ids, unique(mut.maf$Tumor_Sample_Barcode))
missing_in_maf # "TCGA-CU-A0YN" has no data avail.

# The patient has no mutational data but still has sequencing data so we will keep him 
#exclude_id <- "TCGA-CU-A0YN"
#mut_pat_ids <- pat_ids[pat_ids != exclude_id]
mut.maf_sub <- mut.maf[mut.maf$Tumor_Sample_Barcode %in% pat_ids, ]

# Use gene symbol as mutation feature
# Each row in mut.maf is one mutation event -> turn into patient x gene table; we force a row creation even for a patient missing mut data
mut_table <- with(
  mut.maf_sub,
  table(
    factor(Tumor_Sample_Barcode, levels = pat_ids),  # ensures one row per pat_id
    Hugo_Symbol
  )
)

dim(mut_table)

# Binary mutation matrix: 1 if gene mutated in that patient
mut.mat <- (mut_table > 0) * 1
mut.mat <- as.matrix(mut.mat)

# Reorder to match expression matrix
mut.mat <- mut.mat[pat_ids, , drop = FALSE]

# Check alignment
stopifnot(rownames(expr.mat) == rownames(mut.mat))

# Normalization
expr.norm <- standardNormalization(as.matrix(expr.mat))
mut.norm <- standardNormalization(as.matrix(mut.mat))

# Build the affinity/PSN matrices for each layer
# Parameters for SNF (you can tune K, alpha, T)
K <- 4    # number of neighbors
alpha <- 0.5   # scaling factor (as in SNFtool examples)
T <- 20    # number of iterations

# Distance matrices (samples x samples) - the diagonal are non-zero numbers not equal to 1
dist.expr <- dist2(expr.norm, expr.norm)  # Euclidean distance
dist.mut <- dist2(mut.norm, mut.norm)

# Affinity / similarity matrices - these ARE the PSNs for each layer
W.expr <- affinityMatrix(dist.expr, K, alpha)  # expression-based PSN
W.mut <- affinityMatrix(dist.mut,  K, alpha)  # mutation-based PSN

#rownames(W.expr) <- pat_ids
#colnames(W.expr) <- pat_ids
#rownames(W.mut)  <- pat_ids
#colnames(W.mut)  <- pat_ids

# Build SNF
W.list <- list(expr = W.expr, mut = W.mut) 
W.fused <- SNF(W.list, K, T)
# W.fused is the fused patient similarity network (samples x samples)

# Build igraph object from fused matrix
g.fused <- graph_from_adjacency_matrix(
  W.fused,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Apply Louvain on fused PSN
set.seed(42)
louvain.fused <- cluster_louvain(g.fused, weights = E(g.fused)$weight)
membership_fused <- membership(louvain.fused)

# Give names to membership_fused
names(membership_fused) <- pat_ids

# Community sizes
table(membership_fused)
modularity(louvain.fused)

# Visualization of community structure
set.seed(42)

png("snf_communities_simple.png", width = 1200, height = 900, res = 150)
layout_fused <- layout_with_fr(g.fused)

plot(
  g.fused,
  layout = layout_fused,
  vertex.size = 6,
  vertex.label = NA,
  vertex.color = as.factor(membership_fused),
  edge.width = E(g.fused)$weight * 5,
  edge.color = "black",
  main = "Fused Patient Similarity Network (Expression + Mutation)\nLouvain Communities"
)
dev.off()

# Plot with patient codes
V(g.fused)$name <- rownames(expr.mat)              # same patients
V(g.fused)$label <- substr(V(g.fused)$name, 1, 12)  # patient IDs
V(g.fused)$cluster <- membership_fused

png("snf_communities_patients.png", width = 1200, height = 900, res = 150)
plot(
  g.fused,
  layout = layout_fused,
  vertex.size = 6,
  vertex.label = V(g.fused)$label,
  vertex.label.cex = 0.9,
  vertex.color = as.factor(V(g.fused)$cluster),
  edge.width = E(g.fused)$weight * 5,
  edge.color = "black",
  main = "Fused PSN (Expression + Mutation)\n Louvain Communities"
)
dev.off()

# Verify that the order of patients is identical:
names(membership_vec) <- pat_ids      # tumor PSN communities
names(membership_fused) <- pat_ids    # fused PSN communities
stopifnot(names(membership_vec) == names(membership_fused))

# Cross-tabulation of communities
table(Expression_Only = membership_vec,
      Fused = membership_fused)

ari <- adjustedRandIndex(membership_vec, membership_fused)
cat("Adjusted Rand Index between expression-only and fused communities:", ari, "\n")

comm_sizes_expr <- table(membership_vec)
comm_sizes_fused <- table(membership_fused)

png("psn_barplots_expr_vs_fused.png", width = 1800, height = 900, res = 150, bg = "white")

par(mfrow = c(1, 2),
    mar = c(5, 6, 4, 2) + 0.1) 
n_slots <- 5

# Expression-only PSN
barplot(
  comm_sizes_expr,
  col  = pal[names(comm_sizes_expr)],
  main = "Expression-Only PSN",
  xlab = "Community",
  ylab = "Number of Patients",
  ylim = c(0, max(comm_sizes_expr) * 1.2),
  xlim = c(0, n_slots + 1),
  cex.main = 1.5,
  cex.lab  = 1.4,
  cex.axis = 1.2
)

# Fused PSN
barplot(
  comm_sizes_fused,
  col  = pal[names(comm_sizes_fused)],
  main = "Fused PSN: Expression + Mutation",
  xlab = "Community",
  ylab = "Number of Patients",
  ylim = c(0, max(comm_sizes_fused) * 1.3),
  xlim = c(0, n_slots + 1),
  cex.main = 1.5,
  cex.lab  = 1.4,
  cex.axis = 1.2
)

# Reset layout
par(mfrow = c(1, 1))

dev.off()


#Optional tasks ---------------------------

#----------------------------------------------
#Centrality Comparison (Degree vs. Betweenness)
#----------------------------------------------

# Betweenness Centrality  
b.tumor <- sna::betweenness(net.tumor, gmode = "graph")  # how often a node lies on the shortest path between two other nodes.
print(b.tumor)
names(b.tumor) <- network::network.vertex.names(net.tumor) 

# Betweenness (Top 5%)
top_5_pct <- ceiling(network::network.size(net.tumor) * 0.05) 
hubs.betweenness <- names(sort(b.tumor, decreasing = TRUE)[1:top_5_pct])
print(hubs.betweenness)

# Overlap with Degree Hubs (hubs.tumor)
common_hubs_ci <- intersect(hubs.tumor, hubs.betweenness)

message(paste("Total nodes:", network::network.size(net.tumor)))
message(paste("Total Hubs per method:", top_5_pct))
message(paste("Overlap between Degree Hubs and Betweenness Hubs:", length(common_hubs_ci)))
message(paste("Percentage Overlap:", round(length(common_hubs_ci) / top_5_pct * 100, 1), "%"))

# Print the names of the most critical overlapping hubs
print(paste("Overlapping Hubs:", common_hubs_ci))

cat("The 12 genes identified in the overlap are the most critical nodes in the Tumor Network. 
They simultaneously possess: high Connectivity (degree) where are the center of their local gene
expression clusters, and global influence (betweenness) where lie on the shortest communication
paths across the entire network.")
cat("The other 22 hubs, mean that they do not overlap and have different focus. However, to build a drug
for the cancer, there should be focus on the overlap genes (12) + the other 11 betweenness-only to block
communication and coordination between different cancer pathways.")

#get the symbols
clean_ids <- gsub("\\..*", "", common_hubs_ci)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # connect to the Ensembl Database

#Fetch the Gene Symbols
#ask for the ID and the "hgnc_symbol" (the human readable name)
gene_names_comparisson <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = clean_ids,
  mart = mart
)

gene_names_comparisson <- gene_names_comparisson[gene_names_comparisson$hgnc_symbol != "", ] #filter empty names

print(gene_names_comparisson)
cat("The names of the critical genes are: ")
print(gene_names_comparisson$hgnc_symbol)


#---------------------------------------------------
# DIFFERENTIAL CO-EXPRESSED NETWORK - ANALYSIS
# BINARY + SIGNED DIFFERENTIAL SUBNETWORKS (POS/NEG)
#---------------------------------------------------

# Build signed adjacency matrices
# positive differential edges: z_diff > z_threshold  (gains of correlation in tumor vs normal)
pos_diff.adj <- ifelse(z_diff>z_threshold, 1, 0)
diag(pos_diff.adj) <- 0

# negative differential edges: z_diff < -z_threshold (loss or reversed correlation in tumor vs normal)
neg_diff.adj <- ifelse(z_diff < -z_threshold, 1, 0)
diag(neg_diff.adj) <- 0

cat("Positive differential edges (count, directed pairwise):", sum(pos_diff.adj)/2, "\n")
cat("Negative differential edges (count, directed pairwise):", sum(neg_diff.adj)/2, "\n\n")

# Positive and negative degrees per gene
pos_degree <- rowSums(pos_diff.adj)
neg_degree <- rowSums(neg_diff.adj)

deg_df <- data.frame(
  Ensembl_with_version = deg_names,
  Ensembl = gsub("\\..*","",deg_names),
  PosDegree = pos_degree,
  NegDegree = neg_degree,
  TotalDiffDegree = pos_degree+neg_degree,
  stringsAsFactors = FALSE
)

# save degree table
write.csv(deg_df, "Differential_PosNeg_Degrees_raw.csv", row.names = FALSE)

# Identify hubs separately in positive vs negative networks
topN <- top_5_pct_diff # top 5% found earlier

# handle cases with fewer unique non-zero degrees than topN
nonzero_pos <- sum(pos_degree>0)
nonzero_neg <- sum(neg_degree>0)

topN_pos <- min(topN, nonzero_pos)
topN_neg <- min(topN, nonzero_neg)

top_pos_hubs <- names(sort(pos_degree, decreasing=TRUE)[1:topN_pos])
top_neg_hubs <- names(sort(neg_degree, decreasing=TRUE)[1:topN_neg])

cat("Top positive hubs (#):", length(top_pos_hubs), "\n")
cat("Top negative hubs (#):", length(top_neg_hubs), "\n\n")

# Save lists
write.csv(data.frame(Ensembl_with_version=top_pos_hubs, PosDegree=pos_degree[top_pos_hubs]), "Top_Positive_Hubs.csv", row.names=FALSE)
write.csv(data.frame(Ensembl_with_version=top_neg_hubs, NegDegree=neg_degree[top_neg_hubs]), "Top_Negative_Hubs.csv", row.names = FALSE)

# Merge and annotate hubs with gene symbols (biomaRt)
# Prepare list of unique Ensembl ids to query
hubs_all <- unique(c(top_pos_hubs, top_neg_hubs))
hubs_all_clean <- gsub("\\..*", "", hubs_all)

# Query biomaRt (re-use existing 'mart' connection)
if(!exists("mart")) {
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
}

hub_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = hubs_all_clean,
  mart = mart
)

# Merge back to produce annotated hub tables
# create mapping (may be fewer than hubs if some have no HGNC symbol)
map_df <- data.frame(Ensembl = hub_symbols$ensembl_gene_id,
                     HGNC = hub_symbols$hgnc_symbol,
                     Description = hub_symbols$description,
                     stringsAsFactors = FALSE)

# attach to positive hubs
pos_hubs_df <- data.frame(
  Ensembl_with_version = top_pos_hubs,
  Ensembl=gsub("\\..*", "", top_pos_hubs),
  PosDegree=pos_degree[top_pos_hubs],
  stringsAsFactors=FALSE
)
pos_hubs_df <- merge(pos_hubs_df, map_df, by="Ensembl", all.x=TRUE)

# attach to negative hubs
neg_hubs_df <- data.frame(
  Ensembl_with_version=top_neg_hubs,
  Ensembl=gsub("\\..*", "", top_neg_hubs),
  NegDegree=neg_degree[top_neg_hubs],
  stringsAsFactors=FALSE
)
neg_hubs_df <- merge(neg_hubs_df, map_df, by="Ensembl", all.x=TRUE)

# Save annotation files
write.csv(pos_hubs_df, "Top_Positive_Hubs_annotated.csv", row.names=FALSE)
write.csv(neg_hubs_df, "Top_Negative_Hubs_annotated.csv", row.names=FALSE)

# Also save combined hub summary
hub_summary <- merge(deg_df, map_df, by.x = "Ensembl", by.y = "Ensembl", all.x=TRUE)
hub_summary <- hub_summary[order(-hub_summary$TotalDiffDegree), ]
write.csv(hub_summary, "Differential_Hub_Summary_annotated.csv", row.names=FALSE)

# Build and plot hub subnetworks (separate for positive and negative)
# Helper function: create subnetwork for hubs + first neighbours and plot
plot_hub_subnet <- function(adj_mat, hub_list, deg_names, filename_prefix,
                            color_palette = c("Hub"="tomato","Neighbor"="grey")) {
  
  # restrict to nodes that are either hubs or neighbors of hubs (first neighbors)
  hub_indices <- which(deg_names %in% hub_list)
  if (length(hub_indices) == 0) {
    message("No hubs found for ", filename_prefix)
    return(NULL)
  }
  
  # neighbors of these hubs
  is_neighbor <- colSums(adj_mat[hub_indices, , drop = FALSE]) > 0
  nodes_keep <- unique(c(which(deg_names %in% hub_list), which(is_neighbor)))
  nodes_names <- deg_names[nodes_keep]
  
  # build adjacency
  sub_adj <- adj_mat[nodes_names, nodes_names, drop = FALSE]
  
  # build network object
  net_sub <- network(sub_adj, matrix.type = "adjacency", directed = FALSE)
  
  # assign attributes
  net_sub %v% "type" <- ifelse(
    network.vertex.names(net_sub) %in% hub_list, 
    "Hub", 
    "Neighbor"
  )
  net_sub %v% "display_name" <- network.vertex.names(net_sub)
  
  sub_deg <- degree(net_sub)
  
  # Bigger nodes for hubs, tiny nodes for neighbors
  node_sizes <- ifelse(net_sub %v% "type" == "Hub", sub_deg, 0.05)
  
  net_sub %v% "size" <- node_sizes
  
  # plotting
  p <- ggnet2(net_sub,
              color = "type",
              color.palette = color_palette,
              size = "size",
              size.cut = 4,
              #label = TRUE,
              #label.size = 3,
              edge.size = 0.2,
              edge.alpha = 0.4,
              mode = "fruchtermanreingold") +
    labs(
      title = paste0(filename_prefix, " (Hubs + 1st neighbors)"),
      subtitle = paste0("Nodes: ", network::network.size(net_sub),
                        " | Hubs: ", length(hub_list))
    )
  
  # print + save
  print(p)
  ggsave(paste0(filename_prefix, ".png"), plot=p, width=8, height=8, dpi=300)
  
  return(list(net=net_sub, plot=p))
}

# Call plotting function
pos_plot_info <- plot_hub_subnet(pos_diff.adj,top_pos_hubs,deg_names,"Positive_Diff_Hub_Subnetwork")
neg_plot_info <- plot_hub_subnet(neg_diff.adj,top_neg_hubs,deg_names,"Negative_Diff_Hub_Subnetwork",color_palette = c("Hub"="steelblue","Neighbor"="grey"))

# Simple overlap and combined-ranking
# Genes that are hubs for positive vs negative (disjoint or overlapping)
overlap_pos_neg <- intersect(top_pos_hubs, top_neg_hubs)

cat("Overlap between positive-hubs and negative-hubs:", length(overlap_pos_neg), "\n")
if(length(overlap_pos_neg)>0) {
  print(overlap_pos_neg)
}

# Combined ranking by difference of pos vs neg degree (useful to find polarizing genes)
deg_df$PosRank <- rank(-deg_df$PosDegree, ties.method = "min")
deg_df$NegRank <- rank(-deg_df$NegDegree, ties.method = "min")
deg_df$PolarScore <- deg_df$PosDegree - deg_df$NegDegree  # positive: more gained pos links; negative: more lost/neg links

combined_ranking <- deg_df[order(-abs(deg_df$PolarScore), -deg_df$TotalDiffDegree), ]
write.csv(combined_ranking, "Differential_Hub_CombinedRanking.csv", row.names=FALSE)

#---------------------------------------------------
# Computation of the Patient Similarity Network 
# using normal gene expression profile
#---------------------------------------------------

# Patient–Patient similarity (correlation between normal samples)
# normal.log from section 3
psn.cor.normal <- cor(normal.log, method = "pearson")  # patients x patients

# Clean up: remove self-similarity and negative edges
diag(psn.cor.normal) <- 0
psn.cor.normal[psn.cor.normal < 0] <- 0

# Threshold to sparsify: we apply the same applied to the tumor PSN previouisly
sim_thresh_normal <- 0.78
psn.cor.normal.thresh <- psn.cor.normal
psn.cor.normal.thresh[psn.cor.normal.thresh < sim_thresh_normal] <- 0

# Build PSN of normal samples
g.psn.normal <- graph_from_adjacency_matrix(
  psn.cor.normal.thresh,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

summary(g.psn.normal)

# Check: tumor and normal columns correspond to the same patients, same order
stopifnot(
  substr(colnames(network.tumor.data), 1, 12) ==
    substr(colnames(network.normal.data), 1, 12)
)

set.seed(123)
louvain.normal <- cluster_louvain(g.psn.normal, weights = E(g.psn.normal)$weight)
membership_normal <- membership(louvain.normal)

# Community sizes & modularity
table(membership_normal)
modularity(louvain.normal)

# Visualize
set.seed(123)
layout_normal <- layout_with_fr(g.psn.normal)

plot(
  g.psn.normal,
  layout = layout_normal,
  vertex.size = 10,
  vertex.label = NA,
  vertex.color = as.factor(membership_normal),
  edge.width = E(g.psn.normal)$weight,
  main = "Normal PSN (DEG Expression)\nLouvain Communities"
)

V(g.psn.normal)$name <- colnames(network.normal.data)
V(g.psn.normal)$label <- substr(V(g.psn.normal)$name, 1, 12)

plot(
  g.psn.normal,
  layout = layout_normal,
  vertex.size = 10,
  vertex.label = V(g.psn.normal)$label,
  vertex.label.cex = 0.5,
  vertex.color = as.factor(membership_normal),
  edge.width = E(g.psn.normal)$weight,
  main = "Normal PSN with Louvain Communities"
)

# Enhanced vizualization
tg <- as_tbl_graph(g.psn.normal)

set.seed(123)
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.2) +
  geom_node_point(aes(color = factor(membership_normal)), size = 4) +
  scale_color_manual(values = pal, name = "Community") +  # <-- matching colors
  scale_edge_width(range = c(0.1, 2)) +
  theme_void() +
  ggtitle("Normal PSN with Louvain Communities")



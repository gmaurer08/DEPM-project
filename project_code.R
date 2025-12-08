
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


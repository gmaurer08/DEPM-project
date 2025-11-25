
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

#### Downloading the data ####

# File directory
# setwd("C:/Users/galax/OneDrive/Dokumente/University/Magistrale/DEPM/Project")
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
# save(rna.query.normal, rna.query.tumor, rna.data.normal, rna.data.tumor, rna.expr.data.normal, rna.expr.data.tumor,
#     genes.info.tumor, genes.info.normal, genes.info.tumor, genes.info.normal, clinical.query, file="initial-project-data.RData")



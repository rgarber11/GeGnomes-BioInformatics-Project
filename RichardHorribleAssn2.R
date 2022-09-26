#Dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("biomaRt", "DESeq2", "EDASeq"))
if (!require("data.table")) install.packages("data.table")
library(tidyverse)
library(data.table)
library(biomaRt)
library("DESeq2")
library("EDASeq")
library("ggplot2")
#Load Necessary Data
counts_data <- fread("GSE207751_PBMC_counts.csv")
metadata <- fread("GSE207751_PBMC_metadata.csv")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#Use biomaRt to convert to from Ensembl to hgnc
mapped_data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = counts_data$gene_id, mart=mart)
#Remove ones without map
mapped_data <-  data.table(mapped_data[!(is.na(mapped_data$hgnc_symbol) | mapped_data$hgnc_symbol==""), ])
#combine tables
setkey(mapped_data, "ensembl_gene_id")
setkey(counts_data, "gene_id")

counts_data_filtered_hgnc <- (mapped_data[counts_data, nomatch=0])
#Remove data where min samples is 0, since it messes with log
counts_data_filtered_hgnc$Min <- apply(counts_data_filtered_hgnc[, 3:58], MARGIN = 1, FUN = min, na.rm = TRUE)
counts_data_filtered_hgnc <- counts_data_filtered_hgnc[Min != 0]
counts_data_filtered_hgnc <- counts_data_filtered_hgnc[, !"Min", with=FALSE]
#If gene_lengths have been found, no need to refind them
if(file.exists("gene_lengths.csv")) {
  gene_lengths_table <- fread("gene_lengths.csv")
} else {
  gene_lengths <- getGeneLengthAndGCContent(rangedcounts_data_filtered_hgnc$ensembl_gene_id, org="hsa")
  gene_lengths_table <- setDT(as.data.frame(gene_lengths), keep.rownames = "ensembl_gene_id")
  fwrite(gene_lengths_table, "gene_lengths.csv")
}
#Normalize Counts
gene_lengths_table$length <- gene_lengths_table$length / 1000;
gene_lengths_table <- gene_lengths_table[, c("ensembl_gene_id", "length")]
setkey(gene_lengths_table, "ensembl_gene_id")
tpm_counts_data <- gene_lengths_table[counts_data_filtered_hgnc, nomatch=0]
tpm_counts_data <- tpm_counts_data %>% mutate(across(c(4:59), .fns= ~./length))
tpm_counts_data <- tpm_counts_data %>% mutate(across(c(4:59), .fns=~./(sum(.)/1000000)))
tpm_counts_data <- tpm_counts_data[, !"length", with=FALSE]

#Log specify
log_tpm_counts <- tpm_counts_data
log_tpm_counts[, 3:58] <- log(log_tpm_counts[, 3:58], 2)
log_tpm_counts$Min <- apply(log_tpm_counts[, 3:58], MARGIN=1, FUN=min, na.rm=TRUE)
log_tpm_counts$Max <- apply(log_tpm_counts[, 3:58], MARGIN=1, FUN=max, na.rm=TRUE)
log_tpm_counts$Range <- log_tpm_counts[, Max-Min]
# Create density plot
density_plot_table <- log_tpm_counts[, c("ensembl_gene_id", "hgnc_symbol", "Range")]
density_plot <- ggplot(density_plot_table, aes(x=Range)) + geom_density()
density_plot
#PCA Plot
counts_df <- (data.frame(counts_data_filtered_hgnc[, !c("hgnc_symbol"), with=FALSE])) %>% tibble::column_to_rownames("ensembl_gene_id")
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                               colData = metadata,
                               design = ~ condition)
vsd <- vst(dds, blind=FALSE)
pca_plot <- plotPCA(vsd, intgroup="condition")
pca_plot
dups <- counts_data_filtered_hgnc[duplicated(counts_data_filtered_hgnc$hgnc_symbol), ]

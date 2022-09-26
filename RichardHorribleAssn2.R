if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("biomaRt","org.Hs.eg.db", "DESeq2", "EDASeq"))
if (!require("data.table")) install.packages("data.table")
library(tidyverse)
library(data.table)
library(magrittr)
library(org.Hs.eg.db)
library(biomaRt)
library("DESeq2")
library("EDASeq")
library("ggplot2")
countsdata <- fread("GSE207751_PBMC_counts.csv")
normbatchcorrecteddata <- fread("GSE207751_IMSA_PBMC_normandbatchcorrected.csv")
metadata <- fread("GSE207751_PBMC_metadata.csv")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapped_data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = countsdata$gene_id, mart=mart)
map2_data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = normbatchcorrecteddata$V1, mart=mart)
mappeddataAnnotation <- mapIds(
  org.Hs.eg.db,
  keys = normbatchcorrecteddata$V1,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list"
)
mapped_data_org_df <- mappeddataAnnotation %>% tibble::enframe(name = "Ensembl", value = "Symbol") %>% tidyr::unnest(cols = Symbol)
mapped_data_org_df_filtered <- mapped_data_org_df[complete.cases(mapped_data_org_df),]

mapped_data <-  data.table(mapped_data[!(is.na(mapped_data$hgnc_symbol) | mapped_data$hgnc_symbol==""), ])
map2_data <-  data.table(map2_data[!(is.na(map2_data$hgnc_symbol) | map2_data$hgnc_symbol==""), ])
setkey(mapped_data, "ensembl_gene_id")
setkey(map2_data, "ensembl_gene_id")
setkey(countsdata, "gene_id")
setkey(normbatchcorrecteddata, "V1")
normdatafilteredhgnc <- (map2_data[normbatchcorrecteddata, nomatch=0])
rangednormdatafilteredhgnc <- normdatafilteredhgnc
rangednormdatafilteredhgnc[, 3:58] <- log(rangednormdatafilteredhgnc[, 3:58], 2)
rangednormdatafilteredhgnc$Max <- apply(rangednormdatafilteredhgnc[, 3:58], MARGIN = 1, FUN = max, na.rm = TRUE)
rangednormdatafilteredhgnc$Min <- apply(rangednormdatafilteredhgnc[, 3:58], MARGIN = 1, FUN = min, na.rm = TRUE)
rangednormdatafilteredhgnc$Range <- rangednormdatafilteredhgnc[, Max-Min]
rangednormdatafilteredhgnc <- rangednormdatafilteredhgnc[, c("ensembl_gene_id", "hgnc_symbol", "Range")]
p <- ggplot(rangednormdatafilteredhgnc, aes(x=Range)) + geom_density()
unloggedhgnc <- normdatafilteredhgnc
unloggedhgnc$Max <- apply(unloggedhgnc[, 3:58], MARGIN = 1, FUN = max, na.rm = TRUE)
unloggedhgnc$Min <- apply(unloggedhgnc[, 3:58], MARGIN = 1, FUN = min, na.rm = TRUE)
unloggedhgnc$Range <- unloggedhgnc[, Max-Min]
unloggedhgnc <- unloggedhgnc[, c("ensembl_gene_id", "hgnc_symbol", "Range")]
w <- ggplot(unloggedhgnc, aes(x=Range, ..scaled..)) + geom_density()
countsdatafilteredhgnc <- (mapped_data[countsdata, nomatch=0])
rangedcountsdatafilteredhgnc <- countsdatafilteredhgnc
rangedcountsdatafilteredhgnc$Max <- apply(rangedcountsdatafilteredhgnc[, 3:58], MARGIN = 1, FUN = max, na.rm = TRUE)
rangedcountsdatafilteredhgnc <- rangedcountsdatafilteredhgnc[Max != 0]
rangedcountsdatafilteredhgnc <- rangedcountsdatafilteredhgnc[, !"Max", with=FALSE]
rangedcountsdatafilteredhgnc[, 3:58] <- log(rangedcountsdatafilteredhgnc[, 3:58], 2)
rangedcountsdatafilteredhgnc$Max <- apply(rangedcountsdatafilteredhgnc[, 3:58], MARGIN = 1, FUN = max, na.rm = TRUE)
rangedcountsdatafilteredhgnc$Min <- apply(rangedcountsdatafilteredhgnc[, 3:58], MARGIN = 1, FUN = min, na.rm = TRUE)
rangedcountsdatafilteredhgnc$Range <- rangedcountsdatafilteredhgnc[, Max-Min]
rangedcountsdatafilteredhgnc <- rangedcountsdatafilteredhgnc[, c("ensembl_gene_id", "hgnc_symbol", "Range")]
s <- ggplot(rangedcountsdatafilteredhgnc, aes(x=Range)) + geom_density()
genelengths <- getGeneLengthAndGCContent(countsdatafilteredhgnc$ensembl_gene_id, "hsa")
# countsdf <- (data.frame(countsdata)) %>% tibble::column_to_rownames("gene_id")
# dds <- DESeqDataSetFromMatrix(countData = countsdf,
#                               colData = metadata,
#                               design = ~ condition)

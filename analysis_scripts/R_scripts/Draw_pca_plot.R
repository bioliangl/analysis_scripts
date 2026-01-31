PCA_analysis_TPM <- function(rawdata, sampleinfo, ntop = 2000) {
  library(ggplot2)
  library(ggsci)
  ## load data
  tpm_data <- data.matrix(read.csv(rawdata, 
                                   header = TRUE, 
                                   row.names = 1))
  tpm_group <- read.csv(sampleinfo, header = T)
  ## filter data
  keep <- rowSums(tpm_data > 1) >= ncol(tpm_data) * 0.2
  keep_data <- tpm_data[keep, ]
  log_tpm <- log2(keep_data + 1)
  ## get ntop variant gene
  gene_variances <- apply(log_tpm, 1, var)
  n_top_genes <- min(ntop, nrow(log_tpm))
  top_var_indices <- order(gene_variances, decreasing = TRUE)[1:n_top_genes]
  hvg_matrix <- log_tpm[top_var_indices, ]
  ## PCA
  pca_result <- prcomp(t(hvg_matrix), center = TRUE, scale. = TRUE)
  pcaData <- data.frame(
    sample = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    group = tpm_group$group
  )
  percentVar <- round(100 * (pca_result$sdev ^ 2 / sum(pca_result$sdev ^ 2))[1:2])
  ## Plot
  ggplot(pcaData, aes(
    x = PC1,
    y = PC2,
    color = group,
    text = paste0("Sample: ", sample, "<br>Group: ", group)
  )) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      x = paste0("PC1 (", percentVar[1], "%)"),
      y = paste0("PC2 (", percentVar[2], "%)"),
      color = "Group"
    ) +
    scale_color_npg() +
    coord_fixed(ratio = 1.2) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 1, colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black", face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(size = 12, face = "bold")
    )
}

PCA_analysis_Count <- function(rawdata, sampleinfo, ntop = 2000) {
  library(ggplot2)
  library(ggsci)
  library(DESeq2)
  
  ## load data
  count_data <- data.matrix(read.csv(rawdata, header = TRUE, row.names = 1))
  count_group <- read.csv(sampleinfo, header = T)
  ## PCA
  dds <- DESeqDataSetFromMatrix(round(count_data), colData = count_group, design = ~
                                  group)
  vsd <- vst(dds)
  
  pcaData <- plotPCA(
    vsd,
    intgroup = c('group'),
    ntop = 2000,
    returnData = TRUE
  )
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ## Plot
  ggplot(pcaData, aes(
    x = PC1,
    y = PC2,
    color = group,
    text = paste0("Sample: ", name, "<br>Group: ", group)
  )) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      x = paste0("PC1 (", percentVar[1], "%)"),
      y = paste0("PC2 (", percentVar[2], "%)"),
      color = "Group"
    ) +
    scale_color_npg() +
    coord_fixed(ratio = 1.2) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 1, colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black", face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(size = 12, face = "bold")
    )
}
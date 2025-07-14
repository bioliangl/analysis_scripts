library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(ggsci)

gene_count <- read.csv("gene_count_matrix.csv", header = T, row.names = 1)
gene_group <- read.table("sample.info", header = T)
gene_group$treatment <- factor(gene_group$treatment, levels = c("R", "T", "H"))
dds <- DESeqDataSetFromMatrix(round(gene_count), colData = gene_group, design= ~treatment)
vsd <- vst(dds)

pcaData <- plotPCA(vsd, intgroup = c('sample', 'treatment'),ntop = 2000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = treatment, fill = treatment,
                    group = treatment)) +
  geom_point(size = 1.5) +
  geom_mark_ellipse(aes(fill = treatment),
                    expand = unit(2, "mm"),
                    alpha  = 0.2) +
  geom_text_repel(aes(label = sample),
                  size = 3,
                  segment.size = 0.5,
                  show.legend = FALSE) +
  labs(x = paste0("PC1 (", percentVar[1], "%)"),
       y = paste0("PC2 (", percentVar[2], "%)")) +
  scale_color_npg() +
  scale_fill_npg() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(linewidth = 1.2),
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  coord_fixed()
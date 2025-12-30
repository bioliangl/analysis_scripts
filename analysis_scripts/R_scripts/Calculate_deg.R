Calculate_DEG <- function(gene_count,group){
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = gene_count,colData = group,design = ~ treatment)
  dds<-DESeq(dds)
  res_p<-results(dds)
  res_p <- res_p[order(res_p$padj, res_p$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  res_p[which(res_p$log2FoldChange >= 1 & res_p$padj < 0.05),'sig'] <- 'up'
  res_p[which(res_p$log2FoldChange <= -1 & res_p$padj < 0.05),'sig'] <- 'down'
  res_p[which(abs(res_p$log2FoldChange) <= 1 | res_p$padj >= 0.05),'sig'] <- 'none'
  res_p$sig[is.na(res_p$sig)] <- 'none'
  up_gene <- row.names(res_p[res_p$sig=="up",])
  down_gene <- row.names(res_p[res_p$sig=="down",])
  result <- list(res=res_p,up_gene=up_gene,down_gene=down_gene)
  return(result)
}

GetedgeRresult <- function(count,bcv){
  library(edgeR)
  DGElist <- DGEList(counts = count,group = 1:2)
  DGElist <- DGElist[rowSums(cpm(DGElist)>0.1) >= 1,,keep.lib.sizes=FALSE]
  DGElist <- calcNormFactors(DGElist)
  et <- exactTest(DGElist,dispersion = bcv^2)
  result <- cbind(et$table,DGElist$count)
  result[which(result$logFC > 1 & result$PValue < 0.05),'sig'] <- 'up'
  result[which(result$logFC < -1 & result$PValue < 0.05),'sig'] <- 'down'
  result$sig[is.na(result$sig)] <- 'none'
  return(result)
}

Calculate_go_kegg <- function(org,gene_list){
  library(clusterProfiler)
  require(org, character.only = TRUE)
  org_db <- get(org)
  data("pathway2gene")
  data("pathway2name")
  go_result <-  enrichGO(gene = gene_list, OrgDb = org_db,
                         keyType = "GID",ont = "ALL",
                         pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = F)
  kegg_result <- enricher(gene_list, 
                          TERM2GENE = pathway2gene, 
                          TERM2NAME = pathway2name, 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 1)
  
  result <- list(go_result=go_result,kegg_result=kegg_result)
  return(result)
}
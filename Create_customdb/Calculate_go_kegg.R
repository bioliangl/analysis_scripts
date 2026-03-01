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

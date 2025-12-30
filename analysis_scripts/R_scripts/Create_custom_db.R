make_customdb <- function(eggnog_anno, genus, species){
  library(devtools)
  library(tidyverse)
  library(AnnotationForge)
  if(!file.exists('makeorg_info.RData')){
    library(jsonlite)
    library(purrr)
    library(RCurl)
    ## get KEGG info about plants
    if(T){
      kegg_info <- read.table("https://rest.kegg.jp/list/organism", header = F, sep = "\t", fill = T, quote = "")
      kegg_info_plant <- kegg_info[grep("Plants",kegg_info$V4),]
      plant_pathway_info <- vector()
      
      for (i in 1:length(kegg_info_plant$V2)){
        file_url <- paste0("https://rest.kegg.jp/list/pathway/", kegg_info_plant[i,2])
        tmp_info <- try(read.table(file_url, header = F, sep = "\t", fill = T, quote = ""), silent = TRUE)
        if(!inherits(tmp_info, "try-error") && nrow(tmp_info) > 0) {
          tmp_K <- sub(kegg_info_plant[i,2], "ko", tmp_info$V1) %>% unique()
          plant_pathway_info <- append(plant_pathway_info, tmp_K) 
        }
        if(i %% 10 == 0) {
          Sys.sleep(1)
        }
      }
      plant_kegg_list <- plant_pathway_info %>% unique()
    }
    ## get KEGG related information
    if(T){
      kegg_josn_file <- read_lines("https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=")
      pathway2name <- tibble(Pathway = character(), Name = character())
      ko2pathway <- tibble(Ko = character(), Pathway = character())
      kegg <- fromJSON(kegg_josn_file)
      for (a in seq_along(kegg[["children"]][["children"]])) {
        A <- kegg[["children"]][["name"]][[a]]
        
        for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
          B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
          
          for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
            pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
            
            pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
            pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
            pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
            
            kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
            
            kos <- str_match(kos_info, "K[0-9]*")[, 1]
            
            ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
          }
        }
      }
    }
    ## save RData
    save(ko2pathway, pathway2name, plant_kegg_list, file = "makeorg_info.RData")
  }else{
    load("makeorg_info.RData")
  }
  ## load file
  if(T){
    emapper_ann <- read_delim(eggnog_anno,"\t", comment = "##",
                              escape_double = FALSE, trim_ws = TRUE)
  }
  ## get Go annotation
  if(T){
    emapper_useful <-  emapper_ann
    emapper_useful[emapper_useful == "-"] <-  NA
    colnames(emapper_useful)[1] <- "query"
    eggNOG_line_with_go <- emapper_useful$GOs != ""
    eggNOG_annotations_go <- strsplit(emapper_useful$GOs[eggNOG_line_with_go],",")
    emapper_gene2go <-  data.frame(GID = rep(emapper_useful$query[eggNOG_line_with_go],
                                             times = sapply(eggNOG_annotations_go, length)),
                                   GO = unlist(eggNOG_annotations_go)) %>% na.omit() %>% dplyr::distinct()
    emapper_gene2go <- data.frame(emapper_gene2go, EVIDENCE = "IEA")
  }
  ## get KEGG annotation
  if(T){
    eggNOG_line_with_ko <-  emapper_useful$KEGG_ko != ""
    eggNOG_annotations_ko <-  strsplit(emapper_useful$KEGG_ko[eggNOG_line_with_ko],",")
    emapper_gene2ko <-  data.frame(GID = rep(emapper_useful$query[eggNOG_line_with_ko],
                                             times = sapply(eggNOG_annotations_ko, length)),
                                   Ko = unlist(eggNOG_annotations_ko)) %>% na.omit() %>% dplyr::distinct()
    emapper_gene2ko$Ko <-  str_sub(emapper_gene2ko$Ko,4,-1)
    
    emapper_gene2pathway_byko <-  left_join(emapper_gene2ko,ko2pathway, by = "Ko", relationship = "many-to-many") %>%
      dplyr::select(GID, Pathway) %>% na.omit() %>% dplyr::distinct()
    
    eggNOG_line_with_pathway <-  emapper_useful$KEGG_Pathway != ""
    eggNOG_annotations_pathway <-  strsplit(emapper_useful$KEGG_Pathway[eggNOG_line_with_pathway],",")
    gene2pathway_eggNOG_raw <-  data.frame(GID = rep(emapper_useful$query[eggNOG_line_with_pathway],
                                                     times = sapply(eggNOG_annotations_pathway,length)),
                                           Pathway = unlist(eggNOG_annotations_pathway)) %>% na.omit() %>% dplyr::distinct()
    gene2pathway_eggNOG <-  gene2pathway_eggNOG_raw[grep("ko",gene2pathway_eggNOG_raw$Pathway),]  
    emapper_gene2pathway <-  dplyr::bind_rows(emapper_gene2pathway_byko,gene2pathway_eggNOG) %>% dplyr::distinct()
    
    pathway2gene <- emapper_gene2pathway[emapper_gene2pathway$Pathway %in% plant_kegg_list,] %>% dplyr::select(Pathway, GID)
  }
  ## makedb
  if(T){
    makeOrgPackage(go=emapper_gene2go,
                   version="0.1",
                   maintainer="xx<xx@xx>",
                   author="xx",
                   outputDir = ".",
                   tax_id = 666,
                   genus = genus,
                   species = species,
                   goTable="go")
  }
  ## add pathway2gene and pathway2name
  if(T){
    R_pkg_name <- paste0("org", ".", paste0(substr(genus, 1, 1),species, sep = ""), ".", "eg", ".", "db")
    setwd(R_pkg_name)
    devtools::load_all(".")
    usethis::use_data(pathway2gene, pathway2name, overwrite = T)
    devtools::unload(R_pkg_name)
    remotes::install_local(".", force = TRUE)
  }
}
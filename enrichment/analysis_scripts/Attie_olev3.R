#' Attie Ole wrapper for Allez
#' Ole! Ole! Ole!
#' 
#' @usage
#' Run this function for the Allez analysis of the WGCNA modules
#' click the "source on save" button at the top of the script window and
#' then click the save button. The function should appear in the functions
#' window of the enviornment
#' 
#' @details 
#' requires the output for the WGCNA and the Gene info file
#' requires the Allez_attie_v2 function
#' 
#' @param verbose TRUE/FALSE
#' determines the level of detail to describe regarding the function
#' default is TRUE
#' 
#' @param flnms_choose TRUE/FALSE
#'setting whether the use wishes to
#' specify output file names. if not, default names will be used:
#' 
#' @param setenv TRUE/FALSE
#' it chooses whether you want to set the working directory
#'  for writing the files 
#'  default is TRUE
#' 
#' @param gene_set is the variable for the list of genes to analyze
#' can be passed from other functions. If NULL, it will prompt
#' the user to load it. 
#' 
#' @param gene_universe is the variable that defines the set of genes
#' against which to test enrichment. If NULL it will prompt the user
#' to select the corresponding files. 
#' 
#' @param spec_type is the variable setting which species to look in
#' default is "mouse" but can be "human"
#' 
#' @param tst_tail is a numeric setting the tails of the test in allez we 
#' use. The 1 tailed is default and is "T". 2 tailed is "F"
#' 
#' @param flnms_choose is a TRUE/FALSE setting whether the use wishes to
#' specify output file names. if not, default names will be used:
#' Attie_Allez_output_significant.csv
#' 
#' @param pcutoff is a numeric value
#' it sets the cutoff value for significance. default is 0.05
#' 
#' @param Lowersetsize is a numeric value 
#' sets the minimum number of genes considered for
#' each GO term to be included
#' 
#' @param Uppersetsize is a numeric value
#' sets the maximum number of genes considered for
#' each GO term to be included
#' 
#' @src_path_choose TRUE/FALSE
#' defines whether to choose the path to the Allez_attie_v3 script for 'source on save'
#' if FALSE, it will call src_path as the path to the script
#' if TRUE, it will prompt the user to select the Allez_attie_v2 script 
#' default is TRUE
#' 
#' @src_path is the path to the Allez_attie_v3 script. it is called only
#' if src_path_choose is TRUE. Default is NULL
#' 
#' A note on normalcy:
#' I've found that Allez will return modules with 'significant' enrichment that have
#' few genes overlap, which seems weird to me. I allow the user to set a threshold
#' for number of overlapping genes as well as the fraction of the whole set that
#' are enriched. For the latter, I would consider a very small overlap (say 15 in 450 genes)
#' to be "significantly enriched" even if Allez returns a significant p value. Consequently
#' I allow the user to set the thresholds for these things. 
#' 
#' @param numb.overlap sets the minimum number of overlapping genes to be considered 'normal'
#' default is 3
#' 
#' @param frxn_overlap sets the minimum fraction of overlapping genes for the enrichment
#' to be considered 'normal'. Default is 0.05
#' 
#' @author Chris Emfinger, PhD. Reanimated by Herbert West using enriched elixirs.
#' 

#start the function========================================================================================================
Attie_Ole_v3 <- function(
    gene_set = NULL, 
    gene_universe= NULL,
    setenv = TRUE,
    spec_type = "mouse",
    tst_tail = "T",
    flnms_choose = FALSE,
    pcutoff = 0.05,
    Lowersetsize = 5,
    Uppersetsize = 500,
    src_path_choose = TRUE,
    src_path = NULL,
    verbose = TRUE,
    frxn_overlap = 0.05,
    numb.overlap = 3
    ){

#set env=====================================================================
  d<-.GlobalEnv
  
#load libraries============================================================================================================
  if (verbose == TRUE){
    print("...loading libraries...", quote=FALSE)
  }
  
  pckg_missing <- c("You are missing the following packages: ")
  numb_missing <- 0
  if (!require("BiocManager", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " BiocManager")
    numb_missing <- numb_missing + 1
  }
  if (!require("rstudioapi", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " rstudioapi")
    numb_missing <- numb_missing + 1
  }
  if (!require("WGCNA", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " WGCNA")
    numb_missing <- numb_missing + 1
  }
  if (!require("dplyr", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " dplyr")
    numb_missing <- numb_missing + 1
  }
  if (!require("tidyverse", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " tidyverse")
    numb_missing <- numb_missing + 1
  }
  if (!require("shiny", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " shiny")
    numb_missing <- numb_missing + 1
  }
  if (!require("DT", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " DT")
    numb_missing <- numb_missing + 1
  }
  if (!require("shinyFiles", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " shinyFiles")
    numb_missing <- numb_missing + 1
  }
  if (!require("allez", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " allez")
    numb_missing <- numb_missing + 1
  }
  if (!require("org.Mm.eg.db", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " org.Mm.eg.db")
    numb_missing <- numb_missing + 1
  }
  if (!require("allez", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " allez")
    numb_missing <- numb_missing + 1
  }
  if (numb_missing == 0){
    require("rstudioapi")
    require("WGCNA")
    require("dplyr")
    require("tidyverse")
    require("shiny")
    require("shinyFiles")
    require("DT")
    require("allez")
  }
  if (numb_missing >= 1){
    print(paste0("ERROR: You need ", numb_missing, " packages."), quote = FALSE)
    print(pckg_missing, quote = FALSE)
    break
  }
  if (verbose == TRUE){
    print("...done...", quote=FALSE)
  }

  
#set env ==================================================================================================================
  if (setenv == TRUE){
    print("...now you'll select the working directory...", quote=FALSE)
    setwd(selectDirectory())
  }
d <- .GlobalEnv

#src_path ==================================================================================================================

if (verbose == TRUE){
  print("...source saving the Allez_attie_v3 script...", quote=FALSE)
}
if (src_path_choose == TRUE){
  print("...now you'll select the Allez_attie_v3 script...", quote=FALSE)
  src_path <- file.choose()
  source(src_path)
}
if (src_path_choose == FALSE){
  print("...sourcing the Allez attie function...", quote=FALSE)
  source(src_path)
}
if (verbose == TRUE){
  print("...done...", quote=FALSE)
}


#load the data=============================================================================================================
print("...now you'll load the data...", quote=FALSE)
#load the modules file
if (is.null(gene_set)==TRUE){
  print("...select the multi-tissue modules file...", quote=FALSE)
  gene_set <-  data.frame(read.csv(file.choose()),
                          stringsAsFactors = FALSE)
  rownames(gene_set)<-gene_set$Accession
}
if (is.null(gene_set)!=TRUE){
  print("...select the multi-tissue modules file...", quote=FALSE)
  gene_set <-  data.frame(read.csv(gene_set),
                          stringsAsFactors = FALSE)
  rownames(gene_set)<-gene_set$Accession
}

if (is.null(gene_universe) == TRUE){
  print("...select the 'universe' file...", quote=FALSE)
  gene_universe <-  data.frame(read.csv(file.choose()),
                               stringsAsFactors = FALSE)
}  
if (is.null(gene_universe) != TRUE){
  print("...select the 'universe' file...", quote=FALSE)
  gene_universe <-  data.frame(read.csv(gene_universe),
                               stringsAsFactors = FALSE)
}  

if (spec_type == "mouse"){
  lib.v <- "org.Mm.eg"
  lib.n <- AnnotationDbi::as.list(org.Mm.egSYMBOL)
}
if (spec_type == "human"){
  lib.v <- "org.Hs.eg"
  lib.n <- AnnotationDbi::as.list(org.Hs.egSYMBOL)
}
if (verbose == TRUE){
  print("...data loaded...", quote=FALSE)
}
  
#determine the enrichments for each set====================================================================================

#set counters:
cntr <- 0
cntr2 <- 0
cntr3 <- 0

#set dummy variables
gene_set2 <- gene_set
gene_set3 <- gene_set
gene_set4 <- gene_set
gene_set5 <- gene_set[1,]
nriched <- data.frame(matrix(nrow=1,ncol=2),stringsAsFactors = FALSE)
colnames(nriched)<-c("module","enriched?")
nriched2 <- nriched
enrich.strt <- data.frame(matrix(nrow=1, ncol=12), stringsAsFactors = FALSE)
colnames(enrich.strt) <- c("GO_ID", 
                           "Term", 
                           "p.value", 
                           "p.adj", 
                           "z.score",
                           "set.size",
                           "set.mean",
                           "set.sd",
                           "num.overlap",
                           "Genes",
                           "module",
                           "tissue")
#assign enrich.strt
#assign("enrich.strt", enrich.strt, envir = d)

#start repeat loops 
repeat{
  cntr <- cntr+1
  cntr2 <- 0
  print(paste0("...now on ", unique(gene_set$module)[cntr],", ",
               cntr, " of ",length(unique(gene_set$module)),"..."), quote = FALSE)
  gene_set2 <- subset(gene_set, module == unique(gene_set$module)[cntr])
  
  #start inner loop 1
  uni <- data.frame(gene_universe$symbol, stringsAsFactors = FALSE)
  colnames(uni)<-"universe"
  gset <- data.frame(gene_set2$gene_symbol, stringsAsFactors = FALSE)
  colnames(gset) <- "set"
  enrich.temp <- NULL
  
  #note- 'tryCatch' is to deal with the problem of errors int he allez
  #functions themselves. When they occur in normal runs, it breaks the
  #functions linked to them. This way, it will run them and if there
  #are errors it will instead return NULL but not break the whole function

  enrich.temp <- tryCatch({Allez_attie_v3(gene_set = gset, 
                 gene_universe = uni,
                 spec_type = spec_type,
                 tst_tail = tst_tail,
                 pcutoff = pcutoff,
                 Lowersetsize = Lowersetsize,
                 Uppersetsize = Uppersetsize,
                 verbose = FALSE)}, 
                 error = function(e){
                   print("no enrichment in tissue", quote=FALSE)
                   return(NULL)})

  
  if ((is.null(enrich.temp)==TRUE)){
    nriched$module <- gene_set2$module[1]
    nriched$`enriched?`<-"NO" 
  }
  if ((is.null(enrich.temp)==FALSE) && nrow(enrich.temp)==0){
    nriched$module <- gene_set2$module[1]
    nriched$`enriched?`<-"NO" 
  }
  if ((is.null(enrich.temp)==FALSE) && nrow(enrich.temp)!=0){
  
  nriched$module <- gene_set2$module[1]
  nriched$`enriched?`<-"YES"
  
  mods <- data.frame(matrix(nrow=nrow(enrich.temp),ncol=1),stringsAsFactors = FALSE)
  colnames(mods) <- "module"
  tiss <- data.frame(matrix(nrow=nrow(enrich.temp),ncol=1),stringsAsFactors = FALSE)
  colnames(tiss) <- "tissue"
  mods$module <- gene_set2$module[1]
  if (length(unique(gene_set2$tissue)) > 1){
  tiss$tissue <- "all"
  }
  if (length(unique(gene_set2$tissue)) <= 1){
    tiss$tissue <- gene_set2$tissue[1]
  }
  enrich.temp <- cbind(enrich.temp,mods, tiss)
  rownames(enrich.temp)<-paste0(enrich.temp$GO_ID,"_",
                                enrich.temp$module,"_", 
                                enrich.temp$tissue)
  
  enrich.strt <- rbind(enrich.strt, enrich.temp)
  #assign("enrich.strt", enrich.strt, envir = d)
  
  if (length(unique(gene_set2$tissue_origin))<=1){
    cat("...your module has only one tissue or class...\n")
  }
  
  #start inner repeat loop
  if (length(unique(gene_set2$tissue_origin))>1){
  repeat{
    cntr2 <- cntr2+1
    print(paste0("...now on the tissue ", unique(gene_set2$tissue_origin)[cntr2],", ",
                 cntr2, " of ",length(unique(gene_set2$tissue_origin))," tissues in module ",
                 gene_set2$module[1]), quote = FALSE)
    
    gene_set3 <- subset(gene_set2, tissue_origin == unique(gene_set$tissue_origin)[cntr2])
    
    uni <- data.frame(gene_universe$symbol, stringsAsFactors = FALSE)
    colnames(uni)<-"universe"
    gset <- data.frame(gene_set3$gene_symbol, stringsAsFactors = FALSE)
    colnames(gset) <- "set"
    enrich.temp <- NULL
    
    enrich.temp <- tryCatch({Allez_attie_v3(gene_set = gset, 
                   gene_universe = uni,
                   spec_type = spec_type,
                   tst_tail = tst_tail,
                   pcutoff = pcutoff,
                   Lowersetsize = Lowersetsize,
                   Uppersetsize = Uppersetsize,
                   verbose = FALSE)},
                   error=function(e){
                     print("no enrichment in tissue", quote=FALSE)
                     return(NULL)})

    
    if ((is.null(enrich.temp)==FALSE) && nrow(enrich.temp)!=0){
      
      mods <- data.frame(matrix(nrow=nrow(enrich.temp),ncol=1),stringsAsFactors = FALSE)
      colnames(mods) <- "module"
      tiss <- data.frame(matrix(nrow=nrow(enrich.temp),ncol=1),stringsAsFactors = FALSE)
      colnames(tiss) <- "tissue"
      mods$module <- gene_set2$module[1]
      tiss$tissue <- gene_set3$tissue_origin[1]
      
      enrich.temp <- cbind(enrich.temp,mods, tiss)
      rownames(enrich.temp)<-paste0(enrich.temp$GO_ID,"_",
                                    enrich.temp$module,"_", 
                                    enrich.temp$tissue)
      
      enrich.strt <- rbind(enrich.strt, enrich.temp)
      #assign("enrich.strt", enrich.strt, envir = d)
    }
    
    nriched2 <- rbind(nriched2, nriched)
    #end of 2nd loop
    if (cntr2 == length(unique(gene_set2$tissue_origin))){break}
  }
  }
  
  }
  

  
  #end loop
  if (cntr==length(unique(gene_set$module))){break}
}

cat("...writing the summary files...\n")

#write the all containing file
all_enriched <- enrich.strt[-1,]
write.csv(all_enriched, file="All_enrichments.csv", row.names = FALSE)
cat("...writing all enrichments...\n")

#find the weird ones
with_overlap <- subset(all_enriched, num.overlap > 0)
with_overlap <- subset(with_overlap, is.na(with_overlap$num.overlap)!=TRUE)

overlap_fraction <- data.frame(with_overlap[1])
rownames(overlap_fraction)<-rownames(with_overlap)
colnames(overlap_fraction)<-"fraction_overlap"
overlap_fraction$fraction_overlap <- with_overlap$num.overlap/with_overlap$set.size 
print(head(overlap_fraction))
with_overlap <- cbind (with_overlap, overlap_fraction)
print(head(with_overlap))

with_overlap <- subset(with_overlap, fraction_overlap > frxn_overlap)
print(head(with_overlap))
with_overlap <- subset(with_overlap, num.overlap > numb.overlap)
print(head(with_overlap))

#write the ones that seem "normal" 
cat("...writing 'normal' enrichments...\n")
write.csv(with_overlap, file="enriched_looks_normal.csv", row.names = FALSE)

#write the ones that don't seem "normal"
cat("...writing 'abnormal' enrichments...\n")
not_normal <- subset(all_enriched, !(rownames(all_enriched) %in% rownames(with_overlap)))
write.csv(not_normal, file="enriched_seems_abnormal.csv", row.names = FALSE)

if(any(with_overlap$tissue == "all")){
#write the ones that have tissue subset
  cat("...writing specific tissue module enrichments...\n")
with_overlap2 <- subset(with_overlap, tissue != "all")
write.csv(with_overlap2, file="enriched_looks_normal_no_all.csv", row.names = FALSE)
}

#end of function===========================================================================================================  
  }


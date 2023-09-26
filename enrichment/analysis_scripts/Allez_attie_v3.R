#' Allez! Allez! Allez! *waves the cloth at the bull*
#' 
#' Allez_attie()
#' This is for running Allez enrichments using the allez packages
#' 
#' @usage 
#' click "Source on save" at the top of the window then save the file
#' you should see it appear in the "functions" portion of the environment
#' window
#' 
#' type Allez_attie() into the command line. 
#' If you want to turn off messages, set verbose to FALSE
#' if you want to alter the p value cutoff, adjust the pcutoff to
#' the desired value
#' numerics and TRUE/FALSE are entered without quotations
#' 
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
#' @param verbose is a TRUE/FALSE determining whether a lot of the 
#' steps are indicated. 
#' 
#' @author Chris Emfinger, PhD. Award winning author of 'Handbook for
#' the Recently Diseased'
#' 
#' ==========================================================================================

Allez_attie_v3 <- function(gene_set = NULL, 
                        gene_universe = NULL,
                        spec_type = "mouse",
                        tst_tail = "T",
                        pcutoff = 0.05,
                        Lowersetsize = 5,
                        Uppersetsize = 500,
                        verbose = TRUE){
  
# load libraries
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
  
  


#load data ================================================================================================================
  print("...now you'll load the data...", quote=FALSE)
  
  if (is.null(gene_set)==TRUE){
    print("...select the list of genes of interest...", quote=FALSE)
    gene_set <-  data.frame(read.csv(file.choose(),
                                     header=FALSE),
                            stringsAsFactors = FALSE)
    colnames(gene_set)<-"set"
  }
  
  if (is.null(gene_universe) == TRUE){
    print("...select the universe...", quote=FALSE)
    gene_universe <-  data.frame(read.csv(file.choose(),
                                          header=FALSE),
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

  
# process the data=========================================================================================================
  if (verbose == TRUE){
    print("...processing data...", quote=FALSE)
  }
  
  #enrich_allez <- function(File, FileSep, alter, namev, alterTotal, nameTotalv, Local, LocalSep, outprefix, lib.v, side,  Lowersetsize, Uppersetsize, pcut ){
    
    #formatting files
  if (verbose == TRUE){
    print("...formatting files...", quote=FALSE)
  }
    
    #load all of the genes
    allgnames <- data.frame(unique(t(data.frame(lib.n))))
    colnames(allgnames)<-"all_genes"
    
    #id genes in the 'universe' segment
    namev <- gene_universe
    colnames(namev)<-"Universe"
    namev <- as.list(namev$Universe)
    names(namev)<-namev
    rownames(allgnames)<-allgnames$all_genes
    
    #set the gene list
    gene_set<-as.list(gene_set$set)
    names(gene_set)<-gene_set
    
    #define the "Score" input for the Allez feature
    Score <- data.frame(rep(0, length(as.matrix(unique(namev)))))
    colnames(Score)<-"Score"
    rownames(Score)<-unique(names(namev))
    Score[names(gene_set),1] <- 1
    rnms <- rownames(Score)
    Score <- as.list(Score$Score)
    names(Score)<-rnms
    Score<-Score[unique(names(Score))]
    Score<-unlist(Score, use.names = TRUE)
    rnms2 <- names(Score)
    Score<-as.numeric(Score)
    names(Score)<-rnms2
    #names(Score)<-rnms
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
#run the Allez function=========================================================================
    if (verbose == TRUE){
      print("...running the Allez function...", quote=FALSE)
    }
    
    Out<-allez(score=Score,lib=lib.v,idtype="SYMBOL",locallist=namev)
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #extract the p values & scores
    if (verbose == TRUE){
      print("...extracting the scores & p values...", quote=FALSE)
    }
    Mat <- Out$setscores[,c("Term","set.mean","set.sd","set.size","z.score")]
    message(c("one tailed p value? ", tst_tail))
    
    # two tailed
    if(tst_tail=="F"){
      print("...2 test p values calculating...", quote=FALSE)
      Mat$p.value <- pnorm(-abs(Mat$z.score))}
    
    # one tailed
    if (verbose == TRUE){
      print("...calculating one-tailed p values...", quote=FALSE)
    }
    
    if(tst_tail=="T"){
      prb <- pnorm(Mat$z.score)# one tailed
      Mat$p.value <- ifelse(1-prb>prb, prb, 1-prb)*2
    }
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #correcting p value
    if (verbose == TRUE){
      print("...correcting p values...", quote=FALSE)
    }
   
    Mat$p.adj <- p.adjust(Mat$p.value, method="BH")
    
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #filtering by set counts
    if (verbose == TRUE){
      print("...filtering by gene set counts...", quote=FALSE)
    }
    
    Mat <- Mat[which(Mat$set.size>Lowersetsize),]
    Mat <- Mat[which(Mat$set.size<Uppersetsize),]
    
    MatOut <- Mat[order(Mat$p.value),c("Term","p.value","p.adj","z.score","set.size","set.mean","set.sd")]
    
    message("sets with size < ",Lowersetsize, " or > ", Uppersetsize, " are not considered" )
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #add GO IDs
    if (verbose == TRUE){
      print("...adding GO term IDs...", quote=FALSE)
    }
    
    LocalOut <- MatOut[which(is.na(MatOut[,"Term"])),]
    MatOut2 <-  cbind(rownames(MatOut), MatOut)
    LocalOut2 <- cbind(rownames(LocalOut), LocalOut)
    colnames(MatOut2)[1] = colnames(LocalOut2)[1] = "GO_ID"
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #defining the GO terms for which the genes appear
    if (verbose == TRUE){
      print("...defining the GO terms for which the genes appear...", quote=FALSE)
      
    }
    Mat.Cats <- rownames(MatOut2)
    Mat.Aux <- Out$aux$set.data
    gInData <- names(which(Score==1))
    Mat.Aux <- Mat.Aux[which(Mat.Aux$symbol %in% gInData),]
      
    MaxCat <- BiocGenerics::pmin(length(Mat.Cats), 1000)
    if (verbose == TRUE){
      print("...done...", quote=FALSE)
    }
    
    #determining the number of genes in each category
    if (verbose == TRUE){
      print("...determining the number of overlapping genes in each category...", quote=FALSE)
    }
    
    GenesInCats <- lapply(1:MaxCat, function(x) {
        useg <- Mat.Aux[which(Mat.Aux$go_id == Mat.Cats[x]),"symbol"]
        numOL <- length(useg)
        gcats <- paste0(useg, collapse=", ")
        return(list(gcats, numOL))
      })
    
    NumInCats <- unlist(sapply(GenesInCats, function(x) x[2]))
      if(length(Mat.Cats) > 1000) {NumInCats <- c(NumInCats, rep(0, length(Mat.Cats) - MaxCat))}
      names(NumInCats) <- Mat.Cats
      GenesInCats <- unlist(sapply(GenesInCats, function(x) x[1]))
      if(length(Mat.Cats) > 1000) {GenesInCats <- c(GenesInCats, rep("", length(Mat.Cats) - MaxCat))}
      names(GenesInCats) <- Mat.Cats
      if (verbose == TRUE){
        print("...done...", quote=FALSE)
      }
      
    #integrating the overlapping genes
      MatOut2 <- data.frame(MatOut2[Mat.Cats,], num.overlap = NumInCats[Mat.Cats], Genes = GenesInCats[Mat.Cats])
    
    Mat.p <- MatOut2[which(MatOut2$p.adj <= pcutoff),]
    Local.p <- LocalOut2[which(LocalOut2$p.adj <= pcutoff),]
    
    #assign the output 
    #assign("enrich.temp",Mat.p, envir = d)
    enrich.temp <- Mat.p
    return(enrich.temp)
    print("...Allez_attie_v2 function complete.", quote=FALSE)
   
    
    #Out <- list(Allres=MatOut2, Localres=LocalOut2, SigAllres=Mat.p, SigLocalres=Local.p)
  #}  
  
#end of function===========================================================================================================
}

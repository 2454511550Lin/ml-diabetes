#testing v3 of Allez Attie & Attie Ole v3
#with reworked file selection options

gene_set1 <- file.choose() # your gene set file showing the module gene assingments
gene_universe1 <- file.choose() # your gene universe file showing symbols of all genes measured in your experiment
setenv1 <- FALSE # if TRUE it will ask you to set the directory it will write files to
spec_type1 <- "mouse" #options are mouse or human
tst_tail1 <- "T" #default here uses a 1-tailed P value which is the Allez default option
flnms_choose1 <- FALSE #uses the default filenames. If TRUE it will let you save with specific file names. I don't recommend changing that
pcutoff1 <- 0.05 #the adjusted p cutoff
Lowersetsize1 <- 5 #the minimum number of genes in an enrichment term
Uppersetsize1 <- 500 #the maximum number of genes in an enrichment term
src_path_choose1 <- FALSE #if TRUE, it will prompt you to graphically select the Allez_attive_v3.R script
src_path1 <- file.choose() #the path to the Allez_attie_v3.R script
verbose1 <- TRUE #do you want to spit out a lot of messages to let you know how the thing's running?
fraction_overlap <- 0.05 #the minimum fraction of a term's genes found in the module
num.overlap <- 3 #minimum number of genes in the module that map to that enrichment term

#note that some of the options I mentioned
Attie_Ole_v3 (
  gene_set = gene_set1, 
  gene_universe= gene_universe1,
  setenv = setenv1,
  spec_type = spec_type1,
  tst_tail = tst_tail1,
  flnms_choose = flnms_choose1,
  pcutoff = pcutoff1,
  Lowersetsize = Lowersetsize1,
  Uppersetsize = Uppersetsize1,
  src_path_choose = src_path_choose1,
  src_path = src_path1,
  verbose = verbose1,
  frxn_overlap = fraction_overlap,
  numb.overlap = num.overlap
)

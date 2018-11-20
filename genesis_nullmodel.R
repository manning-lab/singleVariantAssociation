# genesis_nullmodel.R
# Description: Generate a null model using the Genesis package; this code is mainly derived from code written by Jen Brody for the DNANexus platform.
# Inputs:
# genotype.file : a genotype file, only sample ids will be used (.gds)
# phenotype.file : delimited text file with columns for each of id.col, outcome, and covariates (.csv or .tsv) 
# outcome.name : column name of the outcome to be tested (string)
# outcome.type : type of outcome (string = continuous, dichotomous)
# covariate.string : comma separated list of covariates to include (string, ex: age,sex)
# sample.file : a text file with one sample id per line (matching both the genotype file and phenotype file) to be included in analysis (.txt)
# label : prefix for output file (string)
# kinship.matrix : a file containing a matrix of relatedness between samples, this should be numeric with entries i,j = relatedness measure for sample i and sample j (.RDa)
# id.col : sample id column in phenotype file (string)
# 
# Ouputs:
# <label>_null.RDa : a null model generated from the desired samples (.RDa)
# 
# 
# Author : Tim Majarian tmajaria@broadinstitute.org The Broad Institute
# Inspiration from Jen Brody, Stephanie M. Gogarten

input_args <- commandArgs(trailingOnly=T)
genotype.file <- input_args[1]
phenotype.file <- input_args[2]
outcome.name <- input_args[3]
outcome.type <-  input_args[4]
covariate.string <- input_args[5]
conditional.string <- input_args[6]
ivars.string <- input_args[7]
group.var <- input_args[8]
sample.file <- input_args[9]
label <- input_args[10]
kinship.matrix <- input_args[11]
id.col <- input_args[12]

## Test inputs ##
setwd("/Users/tmajaria/Documents/projects/public_workflows/singleVariantAssociation/test_inputs/")
genotype.file <- "1KG_phase3_subset_chrX.gds"
phenotype.file <- "1KG_phase3_subset.phenotypes.csv"
outcome.name <- "t2d"
outcome.type <-  "Dichotomous"
covariate.string <- "ancestry,sex,last_exam_age,bmi"
conditional.string <- "NA"
ivars.string <- "NA"
group.var <- "NA"
sample.file <- "1KG_phase3_subset.sampleids.txt"
label <- "1KG_phase3_subset_chrX"
kinship.matrix <- "1KG_phase3_subset.grm.RData"
id.col <- "sample_id"

# BiocManager::install(c("SeqArray","SeqVarTools", "GENESIS", "GWASTools"))
#################

# Load required libraries
suppressMessages(library(data.table))
suppressMessages(library(SeqArray))
suppressMessages(library(GWASTools))
suppressMessages(library(GENESIS))

# MKL function from topmed DCC pipeline, credit Stephanie M. Gogarten (https://github.com/UW-GAC/analysis_pipeline/blob/a530949a49ed99c4bf24e0531ca8f75ddfbf0ca9/TopmedPipeline/R/countThreads.R)
countThreads <- function() {
  nSlots <- Sys.getenv("NSLOTS")
  nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
  if (nThreads == 0) nThreads <- 1
  message(paste("Running with", nThreads,"thread(s)."))
  Sys.setenv(MKL_NUM_THREADS=nThreads)
  nThreads
}

# used for MKL
countThreads()

# Parse the covariate string
if (covariate.string == "NA"){
  covariates = c()
} else {
  covariates <- unlist(strsplit(covariate.string,","))
}

# Check for group_var, heterogenous variances category to use
if (!(group.var == "NA")){
  if (!(group.var %in% covariates)){
    covariates <- c(covariates,group.var)
  }
}
# Check if we have ivars input, add to covars
if (!ivars.string == "NA") {
  ivars <- unlist(strsplit(ivars.string, ","))
  covariates = c(covariates,ivars.string)
}

# If this is conditional, combine with covariates
if (!(conditional.string == "NA")) {
  conditional = unlist(strsplit(conditional.string,","))
  
  # check that the conditional covars start with a letter, if not, add chr
  conditional.edited <- c()
  for (c in conditional){
    if(grepl("[[:digit:]]", substr(c, 1, 1))){
      conditional.edited <- c(conditional.edited, sub(":","\\.",paste("chr",c,sep="")))
    } else {
      conditional.edited <- c(conditional.edited, sub(":","\\.",c))
    }
  }
  
  # combine with covariates
  covariates = c(covariates,conditional.edited)
}

## Load phenotype data
phenotype.data <- fread(
  phenotype.file,
  header = T,
  stringsAsFactors = F, 
  data.table = F)

# show the covariants
print(covariates)
print(colnames(phenotype.data))

# Correct the outcome column if we have a continuous variable
if (tolower(outcome.type) == "continuous"){
  phenotype.data[,outcome.name] <- as.numeric(phenotype.data[,outcome.name])
}

# If we had to change the conditional names, change the fields of the phenotype file
if (!(conditional.string == "NA")) {
  colnames(phenotype.data)[match(conditional,colnames(phenotype.data))] <- conditional.edited
}

# Make sure other continuous variables have numeric columns
for (cur_cov in covariates){
  if (length(unique(phenotype.data[,cur_cov])) >= length(phenotype.data[,1])/2){
    phenotype.data[,cur_cov] <- as.numeric(phenotype.data[,cur_cov])
  }
}

# Remove any duplicate sample ids from the phenotype file
phenotype.data <- phenotype.data[!duplicated(phenotype.data[,id.col]),]

# Gather all the columns that we will need (sample ids, outcome, covariates)
all_vals <- c(id.col,outcome.name,covariates)

# Remove any samples that do not a value for any covariate or outcome
phenotype.slim <- na.omit(as.data.frame(phenotype.data[,all_vals,drop=F]))

# Read in the list of sample ids to be used in analysis
if (!(sample.file == "NA")){
  sample.ids <- unique(readLines(sample.file))
} else {
  sample.ids = phenotype.slim[,id.col]
}

# Only keep phenotypes from these samples
phenotype.slim <- phenotype.slim[phenotype.slim[,id.col] %in% sample.ids,na.omit(all_vals,drop=F)]

# Make row names the sample ids
row.names(phenotype.slim) <- phenotype.slim[,id.col]

# Correct the column name of the sample ids
names(phenotype.slim)[names(phenotype.slim) == id.col] <- "scanID"

# Phenotypes must be in same order in both phenotype file and the genotype file
# Open the genotype file
genotype.data <- seqOpen(genotype.file)

# Extract the sample ids in the correct order
genotype.ids <- seqGetData(genotype.data, "sample.id")

# Remove any sample phenotypes that are not in the genotype file
phenotype.slim <- phenotype.slim[phenotype.slim$scanID %in% genotype.ids,]

# Filter again over the samples in the genotype file
seqSetFilter(genotype.data,sample.id=phenotype.slim$scanID)

# Ensure that we have the right order
genotype.ids <- seqGetData(genotype.data, "sample.id")

# Reorder the phenotype file to match the genotype sample id order
phenotype.slim <- phenotype.slim[match(genotype.ids,phenotype.slim$scanID),,drop=F]

# Close the genotype file
seqClose(genotype.data)

# load the kinship matrix (this should be in .RDa format)
kinship <- get(load(kinship.matrix))

# Remove any samples from the phenotype file that are not in the relatedness matrix
phenotype.slim <- phenotype.slim[phenotype.slim$scanID %in% row.names(kinship),,drop=F]

# Remove any samples from the relatedness matrix that are not in the phenotype file
kinship = kinship[row.names(kinship) %in% phenotype.slim$scanID,colnames(kinship) %in% phenotype.slim$scanID]

# Reorder the relatedness matrix to match the phenotype/genotype sample id order
kinship = kinship[match(phenotype.slim$scanID,row.names(kinship)),match(phenotype.slim$scanID,colnames(kinship))]
kinship = as.matrix(kinship)

# Sample data to the correct format for Genesis
scan.frame <- ScanAnnotationDataFrame(phenotype.slim)

# Fit the null model
if (group.var == "NA"){
  nullmod <- fitNullModel(scan.frame,
                       outcome = outcome.name,
                       covars = covariates,
                       family = ifelse(tolower(outcome.type) == "dichotomous", "binomial", "gaussian"),
                       cov.mat = kinship)
} else {
  nullmod <- fitNullModel(scan.frame,
                       outcome = outcome.name,
                       covars = covariates,
                       family = "gaussian",
                       group.var = group.var,
                       cov.mat = kinship)
}

# Save the null model to the output file
save(nullmod,file=paste(label,"_null.RDa",sep=""))

# show memory stats
ms <- gc()
print(paste0("memory: ", ms[1,6]+ms[2,6], " mb"))

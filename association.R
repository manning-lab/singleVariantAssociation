### association.R
# Description: This function performs an association test to generate p-values for each variant included.
# Inputs:
# gds.file : a genotype file containing data for all samples are variants to be tested (.gds)
# null.file : output of the function *fitNull* or a pregenerated null model (.RDa)
# label : prefix for output filename (string)
# test : statistical test (Score or Wald, default = Score)
# mac : minimum minor allele count for variants to be included in analysis (int, default = 5)
# ivars.string : comma separated list of interaction variables
# variant.range : comma separated list of variant ranges in format <chr#>:pos_start-pos_end
# Outputs:
# assoc : an RData file of associations results (.RData)

# Load packages
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)
library(dplyr)
library(tidyr)

# Parse input arguments
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
mac <- as.numeric(input_args[5])
ivars.string <- input_args[6]
variant.range.file <- input_args[7]

# these are from the DCC pipeline, credit -> S. Gogarten
.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    rename_(allele="alt") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}

# Print input arguments
print_ <- function(to_print) {
	print(paste(to_print, collapse=" "))
}

print_(c("gds.file", gds.file))
print_(c("null.file", null.file))
print_(c("test", test))
print_(c("mac", mac))
print_(c("ivars.string", ivars.string))
print_(c("variant.range.file", variant.range.file))


# Load nullfile
load(null.file)

# Open gds file
gds.data <- seqOpen(gds.file)

# Filter by desired MAC
seqSetFilter(gds.data,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)
gds.freq <- seqAlleleFreq(gds.data, .progress=TRUE)
gds.maf <- pmin(gds.freq, 1-gds.freq)
gds.mac.filt <- 2 * gds.maf * (1-gds.maf) * length(nullmod$scanID) >= mac

# If no snps remain, return empty
if(sum(gds.mac.filt, na.rm = TRUE)==0) {
	print("No SNPs pass MAC filter. Finished Association Step")
	assoc <- NA

# Else move on to association testing
} else {

	# Filter to snps with mac greater than threshold
	seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)
  
  # Filter to only passing variants
  var.ids <- seqGetData(gds.data,"variant.id")
  filt <- seqGetData(gds.data, "annotation/filter")
  var.ids <- var.ids[filt == "PASS"]
  seqSetFilter(gds.data, variant.id = var.ids, action="intersect", verbose=TRUE)

  # stop again if we have no passing variants
  if (length(var.ids) == 0){
  	print("No SNPs pass MAC filter. Finished Association Step")
	  assoc <- NA
	} else {

	  # Organize data for output
	  snps.pos <- .expandAlleles(gds.data)[,c(1,2,3,4,5)]
	  names(snps.pos) <- c("id","chr","pos","ref","alt")
	  
		# Filter by variant range if it exists
		if (!(variant.range.file == "NA")){
			variant.range <- fread(variant.range.file, data.table = F, stringsAsFactors = F)
			var.tokeep.id <- c()
			for (rng.ind in seq(1,nrow(variant.range))){
				rng <- variant.range[rng.ind,]
			  cur.var <- subset(snps.pos, chr == rng[,1] & pos >= as.numeric(rng[,2]) & pos <= as.numeric(rng[,3]))
			  var.tokeep.id <- c(var.tokeep.id, cur.var$id)
			}
			var.tokeep.id <- unique(var.tokeep.id)
			
			seqSetFilter(gds.data, variant.id=var.tokeep.id, action="intersect", verbose=TRUE)
			snps.pos <- snps.pos[snps.pos$id %in% var.tokeep.id,]
			snps.pos <- snps.pos[,!(names(snps.pos) == "chr")]
		}
		
		# Print the number of snps were working with
		print("Filtered SNPs")
		print(nrow(snps.pos))

		# Genotype data to the correct format
		gds.geno.data <- SeqVarData(gds.data)

		# Run association test
		if (ivars.string == "NA"){
			assoc <- assocTestMM(genoData = gds.geno.data, nullMMobj = nullmod, test = test)
		} else {
			assoc <- assocTestMM(genoData = gds.geno.data, nullMMobj = nullmod, test = test, ivars=unlist(strsplit(ivars.string, ",")))
		}
		print("Finished Association Step")
		print(dim(assoc))
		snps.pos <- snps.pos[,!(names(snps.pos) == "chr")]
		assoc <- merge(snps.pos, assoc, by.x = "id", by.y = "snpID")
		
		# get genotype counts for variants w/ MAF < 5% and P < 0.01 if we have a dichotomous trait
		if (!is.null(nullmod$outcome) & length(unique(nullmod$outcome)) == 2){
		
  		# get case/control
  		pheno <- data.frame(sample = as.character(nullmod$scanID), outcome = nullmod$outcome, stringsAsFactors = F)
  		
  		# get the variants that pass both maf and pval threshold 
  		assoc.top_var <- assoc[(assoc$MAF < 0.05 & assoc[,paste0(test,".pval")] < 0.01), "id"]
  		
  		# set filter
  		seqSetFilter(gds.data, variant.id = assoc.top_var, sample.id = pheno$sample)
  		
  		# get genotypes
  		geno <- altDosage(gds.data)
  		geno.ctrl <- geno[row.names(geno) %in% pheno[pheno$outcome == 0, "sample"],]
  		geno.case <- geno[row.names(geno) %in% pheno[pheno$outcome == 1, "sample"],]
  		rm(geno)
  		
  		# get counts per geno
  		geno.ctrl.counts <- apply(geno.ctrl, 2, function(x) sum(x == 0, na.rm = T))
  		geno.ctrl.counts <- data.frame(variant.id = names(geno.ctrl.counts), homref = as.character(geno.ctrl.counts), stringsAsFactors = F)
  		geno.ctrl.counts$het <- as.character(apply(geno.ctrl, 2, function(x) sum(x == 1, na.rm = T)))
  		geno.ctrl.counts$homalt <- as.character(apply(geno.ctrl, 2, function(x) sum(x == 2, na.rm = T)))
  		
  		geno.case.counts <- apply(geno.case, 2, function(x) sum(x == 0, na.rm = T))
  		geno.case.counts <- data.frame(variant.id = names(geno.case.counts), homref = as.character(geno.case.counts), stringsAsFactors = F)
  		geno.case.counts$het <- as.character(apply(geno.case, 2, function(x) sum(x == 1, na.rm = T)))
  		geno.case.counts$homalt <- as.character(apply(geno.case, 2, function(x) sum(x == 2, na.rm = T)))
  		
  		# get to right format
  		geno.counts <- data.frame(
  		  variant.id = geno.ctrl.counts$variant.id, 
  		  homref = paste0(geno.case.counts$homref, "/", geno.ctrl.counts$homref),
  		  het = paste0(geno.case.counts$het, "/", geno.ctrl.counts$het),
  		  homalt = paste0(geno.case.counts$homalt, "/", geno.ctrl.counts$homalt),
  		  stringsAsFactors = F
  		)
  		
  		rm(geno.ctrl.counts)
  		rm(geno.case.counts)
  		rm(geno.ctrl)
  		rm(geno.case)
  		
  		# merge with assoc results
  		assoc <- merge(assoc, geno.counts, by.x = "id", by.y = "variant.id", all.x = T)
  		assoc[is.na(assoc)] <- ""
		}
	}
}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))


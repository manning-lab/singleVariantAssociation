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

# Parse input arguments
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
mac <- as.numeric(input_args[5])
ivars.string <- input_args[6]
variant.range.file <- input_args[7]

## Test inputs ##
# setwd("/Users/tmajaria/Documents/projects/public_workflows/singleVariantAssociation/test_inputs/")
# gds.file <- "1KG_phase3_subset.gds"
# null.file <- "1KG_phase3_subset_null.RDa"
# label <- "1KG_phase3_subset"
# test <- "Score"
# mac <- 20
# ivars.string <- "NA"
# variant.range.file <- "NA"
#################

# Load packages
suppressMessages(library(GWASTools))
suppressMessages(library(GENESIS))
suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))



# these are from the DCC pipeline, credit -> S. Gogarten
.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             filter=seqGetData(gds,"annotation/filter"),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
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
gds.mac.filt <- 2 * gds.maf * (1-gds.maf) * length(nullmod$sample.id) >= mac

# If no snps remain, return empty
if(sum(gds.mac.filt, na.rm = TRUE)==0) {
	print("No SNPs pass MAC filter. Finished Association Step")
	assoc <- NA

# Else move on to association testing
} else {

	# Filter to snps with mac greater than threshold
	seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)
  
  # Filter to only passing variants
  var.data <- .expandAlleles(gds.data)
  var.data <- var.data[var.data$filter == "PASS", ]
  var.data <- var.data[,c("variant.id", "chromosome","position","ref","alt")]
  seqSetFilter(gds.data, variant.id = var.data$variant.id, action="intersect")

  # stop again if we have no passing variants
  if (nrow(var.data) == 0){
  	print("No SNPs pass MAC filter. Finished Association Step")
	  assoc <- NA
	} else {

		# Filter by variant range if it exists
		if (!(variant.range.file == "NA")){
			print("Applying variant filter to this set of variants.")
			variant.range <- fread(variant.range.file, data.table = F, stringsAsFactors = F)
			var.tokeep.id <- c()
			for (rng.ind in seq(1,nrow(variant.range))){
				rng <- variant.range[rng.ind,]
			  cur.var <- subset(var.data, chromosome == rng[,1] & position >= as.numeric(rng[,2]) & position <= as.numeric(rng[,3]))
			  var.tokeep.id <- c(var.tokeep.id, cur.var$id)
			}
			var.tokeep.id <- unique(var.tokeep.id)
			
			seqSetFilter(gds.data, variant.id = var.tokeep.id, action="intersect")
			var.data <- var.data[var.data$variant.id %in% var.tokeep.id,]
		}
		
		# Print the number of snps were working with
		print("Filtered SNPs")
		print(nrow(var.data))
		if(nrow(var.data)==0) {
			print("Zero variants for analysis")
	  		assoc <- NA
		} else {
		  # get to right format
			gds.geno.data <- SeqVarData(gds.data)
			
			# make iterator for association testing
			iterator <- SeqVarBlockIterator(gds.geno.data)
	
			# Run association test
			if (ivars.string == "NA"){
				assoc <- assocTestSingle(iterator, null.model = nullmod, test = test, sparse = FALSE)
			} else {
				assoc <- assocTestSingle(iterator, null.model = nullmod, test = test, ivars=unlist(strsplit(ivars.string, ",")), sparse = FALSE)
			}
			print("Finished Association Step")
			print(dim(assoc))
			
			# merge ref/alt with assoc statistics
			assoc <- merge(var.data, assoc, by.x = c("variant.id","chromosome","position"), by.y = c("variant.id","chr","pos"), all.y = T)
			
			# add mac
			all.mac <- list()
			for(n in unique(assoc$allele.index)){
			  cur.ids <- assoc[assoc$allele.index == n, c("variant.id", "ref","alt")]
			  seqSetFilter(gds.data, variant.id = cur.ids$variant.id, sample.id = nullmod$sample.id)
			  cur.dat <- alleleCount(gds.data, n=n, use.names = T)
			  cur.dat <- data.frame(variant.id = names(cur.dat), mac = cur.dat, stringsAsFactors = F)
			  all.mac[[n]] <- merge(cur.dat, cur.ids, by.x = "variant.id", by.y = "variant.id")
			}
			all.mac <- do.call(rbind, all.mac)
			assoc <- merge(assoc, all.mac, by.x = c("variant.id", "ref", "alt"), by.y = c("variant.id", "ref", "alt"), all.x = T)
			assoc <- assoc[!(is.na(assoc[,paste0(test,".pval")])),]
			
			# get genotype counts for variants w/ MAF < 5% and P < 0.01 if we have a dichotomous trait
			if (!is.null(nullmod$outcome) & length(unique(nullmod$outcome)) == 2) {
			
		  		# get case/control
  				pheno <- data.frame(sample = as.character(nullmod$sample.id), outcome = nullmod$outcome, stringsAsFactors = F)
  				names(pheno) <- c("sample","outcome")
  		
  				# get the variants that pass both maf and pval threshold 
  				assoc.top_var <- assoc[(assoc$freq < 0.05 & assoc[,paste0(test,".pval")] < 0.01), "variant.id"]
  		    
  				# stop if we dont have any variants
  				if (length(assoc.top_var) < 2 ){
  				  assoc$homref.case <- ""
  				  assoc$homref.control <- ""
  				  assoc$het.case <- ""
  				  assoc$het.control <- ""
  				  assoc$homalt.case <- ""
  				  assoc$homalt.control <- ""
  				} else {
  				  
    				# set filter
    				seqSetFilter(gds.data, variant.id = assoc.top_var, sample.id = pheno$sample)
    		
    				# get genotypes
    				geno <- altDosage(gds.data)
    				geno.ctrl <- geno[row.names(geno) %in% pheno[pheno$outcome == 0, "sample"],]
    				geno.case <- geno[row.names(geno) %in% pheno[pheno$outcome == 1, "sample"],]
    				rm(geno)
    		
    				# get counts per geno
    				geno.ctrl.counts <- apply(geno.ctrl, 2, function(x) sum(x == 0, na.rm = T))
    				geno.ctrl.counts <- data.frame(variant.id = names(geno.ctrl.counts), homref = as.numeric(as.character(geno.ctrl.counts)), stringsAsFactors = F)
    				geno.ctrl.counts$het <- as.numeric(as.character(apply(geno.ctrl, 2, function(x) sum(x == 1, na.rm = T))))
    				geno.ctrl.counts$homalt <- as.numeric(as.character(apply(geno.ctrl, 2, function(x) sum(x == 2, na.rm = T))))
    		
    				geno.case.counts <- apply(geno.case, 2, function(x) sum(x == 0, na.rm = T))
    				geno.case.counts <- data.frame(variant.id = names(geno.case.counts), homref = as.numeric(as.character(geno.case.counts)), stringsAsFactors = F)
    				geno.case.counts$het <- as.numeric(as.character(apply(geno.case, 2, function(x) sum(x == 1, na.rm = T))))
    				geno.case.counts$homalt <- as.numeric(as.character(apply(geno.case, 2, function(x) sum(x == 2, na.rm = T))))
    		
    				# get to right format
    				geno.counts <- data.frame(
    				  variant.id = geno.ctrl.counts$variant.id, 
    		  		  homref = paste0(geno.case.counts$homref, "/", geno.ctrl.counts$homref),
    	  			  het = paste0(geno.case.counts$het, "/", geno.ctrl.counts$het),
    				  homalt = paste0(geno.case.counts$homalt, "/", geno.ctrl.counts$homalt),
    				  stringsAsFactors = F
    				)
  
		      	geno.counts <- data.frame(
 		          variant.id = geno.ctrl.counts$variant.id, 
		          homref.case = geno.case.counts$homref,
		          homref.control = geno.ctrl.counts$homref,
		          het.case = geno.case.counts$het,
		          het.control = geno.ctrl.counts$het,
		          homalt.case = geno.case.counts$homalt,
 		          homalt.control = geno.ctrl.counts$homalt,
 		          stringsAsFactors = F
    				)
    		
    				rm(geno.ctrl.counts)
    				rm(geno.case.counts)
    				rm(geno.ctrl)
    				rm(geno.case)
    		
    				# merge with assoc results
    				assoc <- merge(assoc, geno.counts, by.x = "variant.id", by.y = "variant.id", all.x = T)
    				assoc[is.na(assoc)] <- ""
  				}
			}
		}
	}
}

# add marker id column
assoc$MarkerName <- paste(assoc$chromosome, assoc$position, assoc$ref, assoc$alt, sep = "-")
assoc$minor.allele <- "alt"
# assoc$minor.allele <- ifelse(assoc$freq < 0.5, "alt", "ref")

# fix the names of assoc
assoc <- assoc[,c("MarkerName","chromosome","position","ref","alt","minor.allele","freq",paste0(test,".pval"),"n.obs",paste0(test,".Stat"),ifelse(test == "Score", c("Score","Score.SE"), c("Est","Est.SE")),"mac","homref.case","homref.control","het.case","het.control","homalt.case","homalt.control")]  
names(assoc) <- c("MarkerName","chr","pos","ref","alt","minor.allele","maf","pvalue","n",paste0(test,".Stat"),ifelse(test == "Score", c("Score","Score.SE"), c("Est","Est.SE")),"mac","homref.case","homref.control","het.case","het.control","homalt.case","homalt.control")


## save assoc object
save(assoc, file=paste0(label, ".assoc.RData"))


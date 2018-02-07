### association.R
# Description: This function performs an association test to generate p-values for each variant included.
# Inputs:
# gds.file : a genotype file containing data for all samples are variants to be tested (.gds)
# null.file : output of the function *fitNull* or a pregenerated null model (.RDa)
# label : prefix for output filename (string)
# test : statistical test (Score or Wald, default = Score)
# mac : minimum minor allele count for variants to be included in analysis (int, default = 5)
# Outputs:
# assoc : an RData file of associations results (.RData)

# Load packages
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(data.table)

# Parse input arguments
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
label <- input_args[3]
test <- input_args[4]
mac <- as.numeric(input_args[5])
ivars.string <- input_args[6]

# Print input arguments
print_ <- function(to_print) {
	print(paste(to_print, collapse=" "))
}

print_(c("gds.file", gds.file))
print_(c("null.file", null.file))
print_(c("test", test))
print_(c("mac", mac))
print_(c("ivars.string", ivars.string))


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

	# Organize data for output
	pos <- seqGetData(gds.data,"position")
	ref <- as.character(ref(gds.data))
	alt <- as.character(unlist(alt(gds.data)))
	snps.pos <- cbind(pos,ref,alt)
	
	# Print the number of snps were working with
	print("Filtered SNPs")
	print(dim(snps.pos))

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
	assoc <- cbind(snps.pos, assoc)

}
## save assoc object
save(assoc, file=paste(label, ".assoc.RData", sep=""))


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

# Load packages
suppressMessages(library(GWASTools))
suppressMessages(library(GENESIS))
suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Rcpp))



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
seqSetFilter(gds.data,sample.id=nullmod$sample.id, action="intersect", verbose=TRUE)

# old, slow way of doing this
gds.freq <- seqAlleleFreq(gds.data, .progress=TRUE)

# new, faster way of doing this
# cppFunction("
#     double calc_freq(IntegerVector x)
#             {
#             int len=x.size(), n=0, n0=0;
#             for (int i=0; i < len; i++)
#             {
#             int g = x[i];
#             if (g != NA_INTEGER)
#             {
#             n++;
#             if (g == 0) n0++;
#             }
#             }
#             return double(n0) / n;
#             }")
# 
# # seqParallelSetup()
# gds.freq <- seqParallel(TRUE, gds.data, FUN = function(f) {
#   seqApply(f, "genotype", FUN=calc_freq, as.is="double", margin="by.variant")
# }, split = "by.variant")
# seqParallelSetup(FALSE,verbose=F)
# 
# gds.freq <- as.numeric(as.character(gds.freq[1:length(gds.freq)-1]))

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
  var.data <- var.data[,c("variant.id", "chromosome","position","ref","alt","allele.index")]
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
      seqParallelSetup(3)
      
      # Run association test
      if (ivars.string == "NA"){
        assoc <- assocTestSingle(iterator, null.model = nullmod, test = test, sparse = FALSE)
      } else {
        assoc <- assocTestSingle(iterator, null.model = nullmod, test = test, ivars=unlist(strsplit(ivars.string, ",")), sparse = FALSE)
      }
      seqParallelSetup(FALSE)
      
      print("Finished Association Step")
      print(dim(assoc))
      
      # merge ref/alt with assoc statistics
      assoc <- merge(var.data, assoc, by.x = c("variant.id","chromosome","position","allele.index"), by.y = c("variant.id","chr","pos","allele.index"), all.y = T)
      
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
        assoc.top_var <- assoc[(assoc$freq < 0.05 & assoc[,paste0(test,".pval")] < 0.01), c("variant.id","allele.index")]
        # assoc.top_var <- assoc[(assoc$freq < 0.05 & assoc[,paste0(test,".pval")] < 0.0005), c("variant.id","allele.index")]
        
        # stop if we dont have any variants
        if (length(assoc.top_var) == 0){
          assoc$homref.case <- ""
          assoc$homref.control <- ""
          assoc$het.case <- ""
          assoc$het.control <- ""
          assoc$homalt.case <- ""
          assoc$homalt.control <- ""
        } else if (length(assoc.top_var) == 1){
          # set filter
          seqSetFilter(gds.data, variant.id = assoc.top_var, sample.id = pheno$sample)
          
          # get genotypes
          geno <- alleleDosage(gds.data, n=assoc.top_var$allele.index)
          geno.ctrl <- geno[row.names(geno) %in% pheno[pheno$outcome == 0, "sample"],]
          geno.case <- geno[row.names(geno) %in% pheno[pheno$outcome == 1, "sample"],]
          
          # get counts per geno
          geno.ctrl.counts <- length(geno.ctrl[geno.ctrl == 0])
          geno.ctrl.counts <- data.frame(variant.id = colnames(geno), homref = as.numeric(as.character(geno.ctrl.counts)), stringsAsFactors = F)
          geno.ctrl.counts$het <- length(geno.ctrl[geno.ctrl == 1])
          geno.ctrl.counts$homalt <- length(geno.ctrl[geno.ctrl == 2])
          
          geno.case.counts <- length(geno.case[geno.case == 0])
          geno.case.counts <- data.frame(variant.id = colnames(geno), homref = as.numeric(as.character(geno.case.counts)), stringsAsFactors = F)
          geno.case.counts$het <- length(geno.case[geno.case == 1])
          geno.case.counts$homalt <- length(geno.case[geno.case == 2])
          
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
          rm(geno)
          
          # merge with assoc results
          assoc <- merge(assoc, geno.counts, by.x = "variant.id", by.y = "variant.id", all.x = T)
          assoc[is.na(assoc)] <- ""
          
        } else {
          # loop through alt allele indices
          all.counts <- list()
          for(n in unique(assoc.top_var$allele.index)){
            cur.ids <- assoc.top_var[assoc.top_var$allele.index == n, "variant.id"]
            
            # set filter
            seqSetFilter(gds.data, variant.id = cur.ids, sample.id = pheno$sample)
            
            if (length(cur.ids) == 1){
              # get genotypes
              geno <- alleleDosage(gds.data, n=n)
              geno.ctrl <- geno[row.names(geno) %in% pheno[pheno$outcome == 0, "sample"],]
              geno.case <- geno[row.names(geno) %in% pheno[pheno$outcome == 1, "sample"],]
              
              # get counts per geno
              geno.ctrl.counts <- length(geno.ctrl[geno.ctrl == 0])
              geno.ctrl.counts <- data.frame(variant.id = colnames(geno), homref = as.numeric(as.character(geno.ctrl.counts)), stringsAsFactors = F)
              geno.ctrl.counts$het <- length(geno.ctrl[geno.ctrl == 1])
              geno.ctrl.counts$homalt <- length(geno.ctrl[geno.ctrl == 2])
              
              geno.case.counts <- length(geno.case[geno.case == 0])
              geno.case.counts <- data.frame(variant.id = colnames(geno), homref = as.numeric(as.character(geno.case.counts)), stringsAsFactors = F)
              geno.case.counts$het <- length(geno.case[geno.case == 1])
              geno.case.counts$homalt <- length(geno.case[geno.case == 2])
              
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
              rm(geno)
              all.counts[[n]] <- geno.counts
            } else {
              
              # get genotypes
              geno <- alleleDosage(gds.data, n=n)
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
              all.counts[[n]] <- geno.counts
            }
          }
          geno.counts <- do.call(rbind, all.counts)
          
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
assoc$minor.allele <- ifelse(assoc$freq < 0.5, "alt", "ref")

# fix the names of assoc
assoc <- assoc[,c("MarkerName","chromosome","position","ref","alt","minor.allele","freq",paste0(test,".pval"),"n.obs",paste0(test,".Stat"),ifelse(test == "Score", c("Score","Score.SE"), c("Est","Est.SE")),"mac","homref.case","homref.control","het.case","het.control","homalt.case","homalt.control")]  
names(assoc) <- c("MarkerName","chr","pos","ref","alt","minor.allele","maf","pvalue","n",paste0(test,".Stat"),ifelse(test == "Score", c("Score","Score.SE"), c("Est","Est.SE")),"mac","homref.case","homref.control","het.case","het.control","homalt.case","homalt.control")


## save assoc object
save(assoc, file=paste0(label, ".assoc.RData"))


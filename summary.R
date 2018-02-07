# summary.R
# Description: Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all and the top associated variants.
# Inputs:
# pval : the p-value column in the output of assocTest, this should be the statistical test with ".pval" appended (string, Score -> Score.pval, Wald -> Wald.pval)
# pval.threshold : p-value threshold for the returning top associations, top association output will include only variants with a p-value less than the threshold (float, default = 0.0001)
# label : prefix for output filename (string)
# assoc.files : comma separated list of association results, output of assocTest (string)

# Check if required packages are installed (sourced from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them)
packages <- c("qqman","data.table","stringr")
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install,repos='http://cran.us.r-project.org')

# Load packages
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
pval <- input_args[1]
pval.threshold <- as.numeric(input_args[2])
label <- input_args[3]
assoc.files <- unlist(strsplit(input_args[4],","))

# Stop if no assoc files
if (length(assoc.files) == 0){
  fwrite(list(),paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
  fwrite(list(), paste(label, ".topassoc.csv", sep=""),row.names=F)
  pdf(paste(label,"_association_plots.pdf",sep=""),width=8,height=8)
  dev.off()
  
} else {

  # Prep for association files
  assoc.compilation <- c() 
  
  # Loop through association files
  for (i in seq(1,length(assoc.files))) {
    load(assoc.files[i])
    
    # Check that the file is not empty
    if (!is.na(assoc)[1]){
      assoc <- assoc[!is.na(assoc[,pval]),]
      print(dim(assoc))
      
      # Write the results out to a master file
      if (i == 1) {
        write.table(assoc,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
      } else {
        write.table(assoc,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F, append=TRUE)
      }	
    }
  }
  
  # Read master file back in
  assoc.compilation <- fread(paste(label, ".assoc.csv", sep=""),sep=",",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
  
  # Make sure the columns are in the right format
  assoc.compilation$chr <- as.numeric(as.character(assoc.compilation$chr))
  assoc.compilation$pos <- as.numeric(as.character(assoc.compilation$pos))
  assoc.compilation$P <- as.numeric(as.character(assoc.compilation[,pval]))
  
  # Write out the top results
  fwrite(assoc.compilation[assoc.compilation[,pval] < pval.threshold, ], paste(label, ".topassoc.csv", sep=""), sep=",", row.names = F)
  
  # QQ plots by maf
  png(filename = paste(label,"_association_plots.png",sep=""),width = 11, height = 11, units = "in", res=400, type = "cairo")
  par(mfrow=c(3,3))
  
  # All variants
  qq(assoc.compilation$P,main="All variants")
  
  # Common variants
  qq(assoc.compilation$P[assoc.compilation$MAF>0.05],main="Variants with MAF>0.05")
  
  # Rare/Low frequency variants
  qq(assoc.compilation$P[assoc.compilation$MAF<=0.05],main="Variants with MAF<=0.05")
  
  # Manhattan plots by maf
  # All variants
  manhattan(assoc.compilation,chr="chr",bp="pos",p="P", main="All variants")
  
  # Common variants
  manhattan(assoc.compilation[assoc.compilation$MAF>0.05,],chr="chr",bp="pos",p="P", main="Variants with MAF>0.05")
  
  # Rare/Low frequency variants
  manhattan(assoc.compilation[assoc.compilation$MAF<=0.05,],chr="chr",bp="pos",p="P", main="Variants with MAF<=0.05")
  dev.off()
}

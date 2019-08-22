# summary.R
# Description: Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all and the top associated variants.
# Inputs:
# pval : the p-value column in the output of assocTest, this should be the statistical test with ".pval" appended (string, Score -> Score.pval, Wald -> Wald.pval)
# pval.threshold : p-value threshold for the returning top associations, top association output will include only variants with a p-value less than the threshold (float, default = 0.0001)
# label : prefix for output filename (string)
# assoc.files : comma separated list of association results, output of assocTest (string)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
pval.threshold <- as.numeric(input_args[1])
label <- input_args[2]
assoc.files <- unlist(strsplit(input_args[3],","))

## Test inputs ##
# pval.threshold <- "0.05"
# label <- "testing"
# assoc.files <- c("1KG_phase3_subset_chrX.assoc.RData")
#################

# Load packages
suppressMessages(library(qqman))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))

# generate the QQ plot (from AManning, JWessel)
lam = function(x,p=.5){
  x = x[!is.na(x)]
  x.quantile <- quantile(x,p)
  round((qchisq(1-x.quantile,1)/qchisq(p,1)),2)
}

qqpval = function(x, main="", col="black"){
  x<-sort(-log10(x))
  n<-length(x)
  for.xaxis <- qexp(ppoints(n))/log(10)
  plot(x=for.xaxis[seq(round(n*.9),n)], y=x[seq(round(n*.9),n)], 
       xlab="Expected", ylab="Observed", main=main ,col=col ,cex=.8, bg= col, 
       pch = 19,
       xlim=c(1,max(for.xaxis)),
       ylim=c(1,max(x)),bty="n")
  abline(0,1, lty=2)
}

plotQQMH <- function(data, pval, maf, chr, pos, filename) {
  # need to check if we have variants to make all of the plots
  which.plot <- "all"
  
  if (length(data[data$maf < 0.01,pval]) == 0){
    which.plot <- "single"
  } else if (length(data[data$maf < 0.01 & data$mac > 20,pval]) == 0 || length(data[data$mac < 20,pval]) == 0) {
    which.plot <- "double"
  } else {
    which.plot <- "all"
  }
  
  # make sure data types are right
  data$pos <- as.numeric(as.character(data[,pos]))
  data$pval <- as.numeric(data[,pval])
  cols <- brewer.pal(8,"Dark2")
  
  if (which.plot == "single"){
    png(filename = filename, width = 12, height = 4, units = "in", res=400)
    layout(matrix(c(1,2,2),nrow=1,byrow = T))
    
    qqpval(data[,pval], main = "All variants", col=cols[8])
    legend('topleft',c(paste0('GC = ',lam(data[,pval]))),col=c(cols[8]),pch=c(21))
    manhattan(data,chr=chr, bp=pos, p=pval, main="All variants")
    
    dev.off()
    
  } else if (which.plot == "double"){
    png(filename = filename, width = 8, height = 12, units = "in", res=400)
    layout(matrix(c(1,2,3,3,4,4),nrow=3,byrow = T))
    
    qqpval(data[data$maf >= 0.01,pval], main = "MAF>=1%", col=cols[8])
    legend('topleft',c(paste0('GC = ',lam(data[data$maf >= 0.01,pval]))),col=c(cols[8]),pch=c(21))
    qqpval(data[data$maf < 0.01,pval], main = "MAF<1%", col=cols[1])
    legend('topleft',c(paste0('GC = ',lam(data[data$maf < 0.01,pval]))),col=c(cols[8]),pch=c(21)) 
    manhattan(data[data$maf >= 0.01,],chr=chr, bp=pos, p=pval, main="MAF>=1%")
    manhattan(data[data$maf < 0.01,],chr=chr, bp=pos, p=pval, main="MAF<1% & MAC>20")
    
    dev.off()
    
  } else if (which.plot == "all"){
    png(filename = filename, width = 12, height = 16, units = "in", res=400)
    layout(matrix(c(1,2,3,4,4,4,5,5,5,6,6,6),nrow=4,byrow = T))
    cols <- brewer.pal(8,"Dark2")
    
    qqpval(data[data$maf >= 0.01,pval], main = "MAF>=1%", col=cols[8])
    legend('topleft',c(paste0('GC = ',lam(data[data$maf >= 0.01,pval]))),col=c(cols[8]),pch=c(21))
    qqpval(data[data$maf < 0.01 & data$mac > 20,pval], main = "MAF<1% & MAC>20", col=cols[1])
    legend('topleft',c(paste0('GC = ',lam(data[data$maf < 0.01 & data$mac > 20,pval]))),col=c(cols[8]),pch=c(21)) 
    qqpval(data[data$mac <= 20,pval], main = "MAC<=20", col=cols[7])
    legend('topleft',c(paste0('GC = ',lam(data[data$mac <= 20,pval]))),col=c(cols[7]),pch=c(21))
    manhattan(data[data$maf >= 0.01,],chr=chr, bp=pos, p=pval, main="MAF>=1%")
    manhattan(data[data$maf < 0.01 & data$mac > 20,],chr=chr, bp=pos, p=pval, main="MAF<1% & MAC>20")
    manhattan(data[data$mac <= 20,],chr=chr, bp=pos, p=pval, main="MAC<=20")
    
    dev.off()
    
  }
}

top.file <- paste0(label, ".topassoc.csv")
all.file <- paste0(label, ".assoc.csv")
png.file <- paste0(label,"_association_plots.png")

# Stop if no assoc files
if (length(assoc.files) == 0){
  fwrite(list(), file = top.file, sep=",", row.names=F)
  fwrite(list(), file = all.file, sep=",", row.names=F)
  png(filename = png.file, width = 1, height = 1, units = "in", res=1)
  dev.off()
  
} else {

  # Prep for association files
  assoc.compilation <- c() 
  j <- 1
  
  # Loop through association files
  for (i in seq(1,length(assoc.files))) {
    load(assoc.files[i])
    
    # Check that the file is not empty
    if (!is.na(assoc)[1]){
      assoc <- assoc[!is.na(assoc$pvalue),]
      print(dim(assoc))
      
      # Write the results out to a master file
      if (j == 1) {
        write.table(assoc, all.file, sep=",", row.names=F, quote = FALSE)
      } else {
        write.table(assoc, all.file, col.names=FALSE, sep=",", row.names=F, quote = FALSE, append=TRUE)
      }	
      j <- j + 1
    }
  }
  
  # Read master file back in
  assoc.compilation <- fread(all.file, sep=",", header=T, stringsAsFactors=F, data.table=F)
  
  # Make sure the columns are in the right format
  assoc.compilation[assoc.compilation$chr == "X", "chr"] <- 23
  assoc.compilation[assoc.compilation$chr == "Y", "chr"] <- 24
  assoc.compilation[assoc.compilation$chr == "MT", "chr"] <- 25
  assoc.compilation$chr <- as.numeric(as.character(assoc.compilation$chr))
  assoc.compilation$pos <- as.numeric(as.character(assoc.compilation$pos))
  assoc.compilation$pvalue <- as.numeric(as.character(assoc.compilation$pvalue))

  # QQ plot
  plotQQMH(assoc.compilation, "pvalue", "maf", "chr", "pos", png.file)
}

# Write out the top results
fwrite(assoc.compilation[assoc.compilation[,"pvalue"] < pval.threshold, ], top.file, sep=",", row.names = F, quote = FALSE, na = "")

# write out all results
fwrite(assoc.compilation, all.file, sep=",", row.names = F, quote = FALSE, na = "")

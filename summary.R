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
pval.threshold <- "0.05"
label <- "1KG_phase3_subset"
assoc.files <- c("1KG_phase3_subset_chrX.assoc.RData")
#################

# Load packages
suppressMessages(library(qqman))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

# Stop if no assoc files
if (length(assoc.files) == 0){
  fwrite(list(),paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
  fwrite(list(), paste(label, ".topassoc.csv", sep=""),row.names=F)
  pdf(paste(label,"_association_plots.pdf",sep=""),width=8,height=8)
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
        write.table(assoc,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F,quote = FALSE)
      } else {
        write.table(assoc,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F,quote = FALSE, append=TRUE)
      }	
      j <- j + 1
    }
  }
  
  # Read master file back in
  assoc.compilation <- fread(paste(label, ".assoc.csv", sep=""),sep=",",header=T,stringsAsFactors=F,data.table=F)
  
  # Make sure the columns are in the right format
  assoc.compilation[assoc.compilation$chr == "X", "chr"] <- 23
  assoc.compilation[assoc.compilation$chr == "Y", "chr"] <- 24
  assoc.compilation[assoc.compilation$chr == "MT", "chr"] <- 25
  assoc.compilation$chr <- as.numeric(as.character(assoc.compilation$chr))
  assoc.compilation$pos <- as.numeric(as.character(assoc.compilation$pos))
  assoc.compilation$pvalue <- as.numeric(as.character(assoc.compilation$pvalue))

  # generate the QQ plot (from J Wessel)
  qqpval2 = function(x, ymin, main="", col="black"){
    x<-sort(-log(x[x>0],10))
    n<-length(x)
    ymax <- -log(ymin,10)

    plot(x=qexp(ppoints(n))/log(10), y=x, xlab="Expected", ylab="Observed", main=main ,col=col ,cex=.8, bg= col, pch = 21, ylim=c(0,ymax))
    abline(0,1, lty=2)
  }
  
  # QQ without identity line
  qqpvalOL = function(x, col="blue"){
    x<-sort(-log(x[x>0],10))
    n<-length(x)
    points(x=qexp(ppoints(n))/log(10), y=x, col=col, cex=.8, bg = col, pch = 21)
  }
  
  # get the right colors
  library(RColorBrewer)
  cols <- brewer.pal(8,"Dark2")
  
  # calculate control
  lam = function(x,p=.5){
    x = x[!is.na(x)]
    chisq <- qchisq(1-x,1)
    round((quantile(chisq,p)/qchisq(p,1)),2)
  }
  
  # QQ plot
  png(filename = paste(label,"_association_plots.png",sep=""),width = 11, height = 11, units = "in", res=400, type = "cairo")
  # par(mfrow=c(1,2))
  layout(matrix(c(1,2,3,3),nrow=2,byrow = T))

  min.p <- min(assoc.compilation$pvalue)
  qqpval2(assoc.compilation$pvalue,col=cols[8],ymin=min.p)
  legend('topleft',c(paste0('ALL ',lam(assoc.compilation$pvalue))),col=c(cols[8]),pch=c(21))
  
  qqpval2(assoc.compilation[assoc.compilation$maf>=0.05,"pvalue"],col=cols[1],ymin=min.p)

  qqpvalOL(assoc.compilation[assoc.compilation$maf < 0.05,"pvalue"],col=cols[2])
  legend('topleft',c(paste0('MAF >= 5%  ',lam(assoc.compilation[assoc.compilation$maf>=0.05,"pvalue"])),
                     paste0('MAF < 5%  ',lam(assoc.compilation[assoc.compilation$maf < 0.05,"pvalue"]))
  ),
  col=c(cols[1],cols[2]),pch=c(21,21))

  manhattan(assoc.compilation,chr="chr",bp="pos",p="pvalue", main="All variants")
  dev.off()
}

# Write out the top results
fwrite(assoc.compilation[assoc.compilation[,"pvalue"] < pval.threshold, ], paste(label, ".topassoc.csv", sep=""), sep=",", row.names = F, quote = FALSE, na = "")

# write out all results
fwrite(assoc.compilation,paste(label, ".assoc.csv", sep=""), sep=",", row.names = F, quote = FALSE, na = "")

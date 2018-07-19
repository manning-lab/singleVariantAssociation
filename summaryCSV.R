# summaryCSV.R
# Description: Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all and the top associated variants.
# Inputs:
# pval : the p-value column in the output of assocTest, this should be the statistical test with ".pval" appended (string, Score -> Score.pval, Wald -> Wald.pval)
# pval.threshold : p-value threshold for the returning top associations, top association output will include only variants with a p-value less than the threshold (float, default = 0.0001)
# label : prefix for output filename (string)
# assoc.file : csv

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
assoc.file <- input_args[4]

assoc <- fread(assoc.file, data.table = F, stringsAsFactors = F)
assoc$MarkerName <- apply(assoc,1,function(x){paste("chr",sub(" +","",x["chr"]),"-",sub(" +","",x["pos"]),"-",x["ref"],"-",x["alt"],sep="")})

# Make sure the columns are in the right format
assoc$chr <- as.numeric(as.character(assoc$chr))
assoc$pos <- as.numeric(as.character(assoc$pos))
assoc$P <- as.numeric(as.character(assoc[,pval]))

# Write out the top results
fwrite(assoc[assoc[,pval] < pval.threshold, ], paste(label, ".topassoc.csv", sep=""), sep=",", row.names = F, quote = FALSE)

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

min.p <- min(assoc[,pval])
qqpval2(assoc[,pval],col=cols[8],ymin=min.p)
legend('topleft',c(paste0('ALL ',lam(assoc[,pval]))),col=c(cols[8]),pch=c(21))

qqpval2(assoc[assoc$MAF>=0.05,pval],col=cols[1],ymin=min.p)

qqpvalOL(assoc[assoc$MAF < 0.05,pval],col=cols[2])
legend('topleft',c(paste0('MAF >= 5%  ',lam(assoc[assoc$MAF>=0.05,pval])),
                   paste0('MAF < 5%  ',lam(assoc[assoc$MAF < 0.05,pval]))
),
col=c(cols[1],cols[2]),pch=c(21,21))

manhattan(assoc,chr="chr",bp="pos",p="P", main="All variants")
dev.off()


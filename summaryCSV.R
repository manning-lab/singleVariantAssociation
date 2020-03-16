# Check if required packages are installed (sourced from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them)
packages <- c("qqman","data.table","stringr","RColorBrewer")
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install,repos='http://cran.us.r-project.org')

# Load packages
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
pval <- input_args[1]
pval.threshold <- as.numeric(input_args[2])
label <- input_args[3]
assoc.files <- unlist(strsplit(input_args[4], ","))

# load results
assoc <- do.call(rbind, lapply(assoc.files, fread, data.table = F, stringsAsFactors = F))
assoc$MarkerName <- apply(assoc,1,function(x){paste("chr",sub(" +","",x["chr"]),"-",sub(" +","",x["pos"]),"-",x["ref"],"-",x["alt"],sep="")})

# get right pval column
if (pval == "NA"){
  pval <- names(assoc)[grep("pval",names(assoc))][1]
}

assoc$chr <- sub("chr", "", assoc$chr)

# Make sure the columns are in the right format
if (any(assoc$chr == "X")) assoc[assoc$chr == "X", "chr"] <- 23
if (any(assoc$chr == "Y")) assoc[assoc$chr == "Y", "chr"] <- 24
if (any(assoc$chr == "M")) assoc[assoc$chr == "M", "chr"] <- 25
assoc$chr <- as.numeric(as.character(assoc$chr))
assoc$pos <- as.numeric(as.character(assoc$pos))
assoc$P <- as.numeric(as.character(assoc[,pval]))

# remove 0, NA, or Inf pvalues
assoc <- assoc[!is.na(assoc[[pval]]),]
assoc <- assoc[assoc[[pval]] > 0,]
assoc <- assoc[!is.infinite(assoc[[pval]]),]

# write out all results
fwrite(assoc, paste0(label, ".all.assoc.csv"), sep=",", row.names = F)

# Write out the top results
top.assoc <- assoc[assoc[,pval] < pval.threshold, ]
if (nrow(top.assoc) == 0){
  fwrite(list(), paste0(label, ".top.assoc.csv"), sep=",", row.names = F)
} else {
  fwrite(top.assoc, paste0(label, ".top.assoc.csv"), sep=",", row.names = F)  
}

# generate the QQ plot (from J Wessel)
qqpval2 = function(x, ymin, main="", col="black"){
  x <- x[!is.na(x)]
  x<-sort(-log(x[x>0],10))
  n<-length(x)
  ymax <- -log(ymin,10)
  
  plot(x=qexp(ppoints(n))/log(10), y=x, xlab="Expected", ylab="Observed", main=main ,col=col ,cex=.8, bg= col, pch = 21, ylim=c(0,ymax))
  abline(0,1, lty=2)
}

# QQ without identity line
qqpvalOL = function(x, col="blue"){
  x <- x[!is.na(x)]
  x<-sort(-log(x[x>0],10))
  n<-length(x)
  points(x=qexp(ppoints(n))/log(10), y=x, col=col, cex=.8, bg = col, pch = 21)
}

# get the right colors
cols <- brewer.pal(8,"Dark2")

# calculate control
lam = function(x,p=.5){
  x = x[!is.na(x)]
  chisq <- qchisq(1-x,1)
  round((quantile(chisq,p)/qchisq(p,1)),2)
}

# QQ plot

png(filename = paste0(label,".association.plots.png"),width = 11, height = 11, units = "in", res=400, type = "cairo")
layout(matrix(c(1,2,3,3),nrow=2,byrow = T))

min.p <- min(assoc[,pval], na.rm = T)
qqpval2(assoc[,pval],col=cols[8],ymin=min.p)
legend('topleft',c(paste0('ALL ',lam(assoc[,pval]))),col=c(cols[8]),pch=c(21))

qqpval2(assoc[assoc$freq>=0.05,pval],col=cols[1],ymin=min.p)

qqpvalOL(assoc[assoc$freq < 0.05,pval],col=cols[2])
legend('topleft',c(paste0('MAF >= 5%  ',lam(assoc[assoc$freq>=0.05,pval])),
                   paste0('MAF < 5%  ',lam(assoc[assoc$freq < 0.05,pval]))
),
col=c(cols[1],cols[2]),pch=c(21,21))

manhattan(assoc[!is.na(assoc[[pval]]),],chr="chr",bp="pos",p="P", main="All variants")
dev.off()


FROM robbyjo/r-mkl-bioconductor:3.4.3-16.04-2018.1

MAINTAINER tmajaria@broadinstitute.org

RUN apt-get update && apt-get -y install git dstat

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN echo "install.packages(c('data.table','qqman','stringr','RColorBrewer','dplyr','tidyr'), repos='http://cran.us.r-project.org')" > install.R && \
	echo "source('https://bioconductor.org/biocLite.R')" >> install.R && \
	echo "biocLite(c('SeqArray','SeqVarTools','GENESIS','GWASTools'))" >> install.R && \
	R --vanilla < install.R && \
	rm install.R

RUN git clone https://github.com/manning-lab/singleVariantAssociation.git && \
	cd ./singleVariantAssociation && \
	git pull origin master

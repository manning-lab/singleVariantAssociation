FROM robbyjo/r-mkl-bioconductor:3.4.3-16.04-2018.1
RUN apt-get update
RUN apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"

RUN git clone https://github.com/manning-lab/singleVariantAssociation.git
RUN cd ./singleVariantAssociation && git pull origin master
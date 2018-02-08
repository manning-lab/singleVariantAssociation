FROM uwgac/r343-topmed:master
RUN apt-get update
RUN apt-get -y install git

RUN git clone https://github.com/manning-lab/singleVariantAssociation.git
RUN cd ./singleVariantAssociation && git pull origin master
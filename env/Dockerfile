FROM continuumio/miniconda
MAINTAINER Katie Evans <katiesevans9@gmail.com>

#COPY conda.yml .
#RUN \
#   conda env update -n root -f conda.yml \
#&& conda clean -a

RUN conda install -c bioconda bcftools=1.9
RUN conda install -c conda-forge pandoc=2.12
RUN conda install -c bioconda plink=1.90b6.12
RUN conda install -c bioconda bedtools=2.29.2
RUN conda install -c conda-forge r=3.6.0



# install other tools not avalible on conda cloud
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

RUN conda install -c conda-forge mscorefonts
FROM continuumio/miniconda
MAINTAINER Katie Evans <katiesevans9@gmail.com>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a


# install other tools not avalible on conda cloud
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

RUN conda install -c conda-forge mscorefonts
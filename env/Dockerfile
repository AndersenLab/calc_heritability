FROM rocker/tidyverse
MAINTAINER Katie Evans <katiesevans9@gmail.com>

RUN Rscript -e "install.packages(c('boot', 'sommer', 'lme4', 'futile.logger'), dependencies=TRUE, repos='http://cran.us.r-project.org')"




# FROM continuumio/miniconda
# MAINTAINER Katie Evans <katiesevans9@gmail.com>

# RUN conda install -c bioconda bcftools=1.9 plink=1.90b6.12 bedtools=2.29.2
# RUN conda install -c conda-forge pandoc=2.12 r mscorefonts
# RUN conda install -c conda-forge r-tidyverse

# #r-data.table r-lme4 r-knitr r-markdown r-futile.logger

# # RUN Rscript -e "install.packages('boot', dependencies=TRUE, repos='http://cran.us.r-project.org')"
# # RUN Rscript -e "install.packages('sommer', dependencies=TRUE, repos='http://cran.us.r-project.org')"


# # install other tools not avalible on conda cloud
# RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

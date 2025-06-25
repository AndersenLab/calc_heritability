FROM mambaorg/micromamba:2.3.0

USER root
# install other tools not avalible on conda cloud
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps libxrender1 libcairo2-dev 

COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml .

RUN micromamba install -y -n base -f conda.yml \
&& micromamba clean -a

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN Rscript -e "install.packages('valr', dependencies=TRUE, repos='http://cran.us.r-project.org')"

RUN micromamba install -c conda-forge r-fuzzyjoin
RUN micromamba install -c bioconda bioconductor-iranges

RUN Rscript -e "remotes::install_version('dplyr', version = '1.0.4', repos = 'http://cran.us.r-project.org')"

RUN micromamba install -c conda-forge mscorefonts
RUN micromamba clean --all --yes

USER $MAMBA_USER
ENV PATH "/opt/conda/bin:$PATH"
ENTRYPOINT source /usr/local/bin/_entrypoint.sh
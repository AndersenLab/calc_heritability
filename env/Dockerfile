FROM r-base:4.5.1

RUN apt update \
    && apt install -y r-cran-lme4 \
                      r-cran-tidyverse \
                      r-cran-dt \
                      r-cran-genetics \
                      r-cran-ggbeeswarm \
                      r-cran-ggrepel \
                      r-cran-markdown \
                      r-cran-plotly \
                      r-cran-rspectra \
                      r-cran-cowplot \
                      r-cran-futile.logger \
                      procps \
    && apt clean
RUN Rscript -e "install.packages('sommer')"

RUN apt-get install -y autoconf \
                       automake \
                       make \
                       gcc \
                       perl \
                       zlib1g-dev \
                       libbz2-dev \
                       liblzma-dev \
                       libcurl4-gnutls-dev \
                       libssl-dev \
                       libncurses5-dev \
    && apt clean

RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 -O bcftools.tar.bz2 \
    && tar -xjvf bcftools.tar.bz2 \
    && cd bcftools-1.22 && make && make prefix=/usr/local/bin install \
    && ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools \
    && cd ../ && rm -r bcftools-1.22

RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 -O htslib.tar.bz2 \
    && tar -xjvf htslib.tar.bz2 \
    && cd htslib-1.22 \
    && make && make install \
    && cd ../ && rm -r htslib-1.22

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250615.zip -O plink.zip \
    && unzip plink.zip \
    && mv plink /usr/bin/ \
    && rm plink.zip
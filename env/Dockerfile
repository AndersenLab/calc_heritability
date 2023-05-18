FROM rocker/tidyverse:3.6.0
MAINTAINER Katie Evans <katiesevans9@gmail.com>

RUN Rscript -e "install.packages(c('futile.logger', 'remotes', 'devtools'), dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_github('lme4/lme4', build_vignettes=FALSE)"
RUN Rscript -e "remotes::install_version('sommer', version = '4.0.4', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_version('boot', version = '1.3-28', repos = 'http://cran.us.r-project.org')"

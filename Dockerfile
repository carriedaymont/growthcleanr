FROM rocker/tidyverse:latest

MAINTAINER dlchudnov@mitre.org

WORKDIR /usr/src/app

COPY LICENSE /LICENSE
COPY README.md /README.md

RUN R -e "install.packages(c( \
    'argparser', \
    'bit64', \
    'data.table', \
    'doParallel', \
    'dplyr', \
    'foreach', \
    'Hmisc', \
    'plyr' \
    ))"

RUN R -e "library(devtools); \
    devtools::install_github('mitre/growthcleanr')"

ADD exec/gcdriver.R /usr/local/bin/
RUN chmod +x /usr/local/bin/gcdriver.R

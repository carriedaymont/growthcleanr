FROM rocker/tidyverse:latest

MAINTAINER dlchudnov@mitre.org

WORKDIR /usr/src/app

COPY LICENSE /LICENSE
COPY README.md /README.md

RUN R -e "remotes::install_github('mitre/growthcleanr', dependencies = TRUE)"

ADD exec/gcdriver.R /usr/local/bin/
RUN chmod +x /usr/local/bin/gcdriver.R

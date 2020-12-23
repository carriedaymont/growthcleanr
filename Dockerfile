FROM rocker/tidyverse:latest

LABEL maintainer="Daniel Chudnov <dlchudnov@mitre.org>"

WORKDIR /app

COPY LICENSE /LICENSE
COPY README.md /README.md

RUN R -e "devtools::install_github('mitre/growthcleanr', dependencies = TRUE, ref = 'main')"

ADD exec/gcdriver.R /usr/local/bin/
RUN chmod +x /usr/local/bin/gcdriver.R

RUN addgroup --gid 1202 gcuser && adduser --system --uid 1002 --ingroup gcuser gcuser
RUN chown -R gcuser:gcuser /app

USER gcuser

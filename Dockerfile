FROM ghcr.io/rocker-org/tidyverse:latest

LABEL maintainer="Daniel Chudnov <dlchudnov@mitre.org>"

WORKDIR /app

COPY LICENSE /LICENSE
COPY README.md /README.md

RUN addgroup --gid 1202 gcuser && adduser --system --uid 1002 --ingroup gcuser gcuser

RUN mkdir /app/R_libs
RUN echo "R_LIBS_USER=/app/R_libs" > /home/gcuser/.Renviron
RUN echo ".libPaths(c('/app/R_libs', .libPaths()))" > /home/gcuser/.Rprofile

RUN R -e "install.packages('growthcleanr', dependencies = TRUE, lib='/app/R_libs')"

ADD exec/gcdriver.R /usr/local/bin/
RUN chmod ugo+rx /usr/local/bin/gcdriver.R
RUN chown -R gcuser:gcuser /app

USER gcuser

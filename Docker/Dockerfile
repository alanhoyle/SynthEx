FROM rocker/r-ver:3.5.2
# FROM rocker/r-base:latest

LABEL maintainer="Alan Hoyle <alanh@unc.edu>"
LABEL description="The SynthEx Copy Number Variant caller"

RUN apt-get update
RUN apt-get install -y libcurl4-openssl-dev  libssl-dev python zlib1g-dev

ENV bedtools_version 2.28.0

ADD https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools /usr/local/bin
RUN chmod 755 /usr/local/bin/bedtools

RUN R -q -e "install.packages(c('devtools', 'covr', 'roxygen2', 'testthat','optparse','data.table'), repos = 'https://cloud.r-project.org/')" &&\
rm -rf /tmp/*

RUN Rscript -e ' install.packages("BiocManager", repos = "https://cloud.r-project.org/") ; BiocManager::install() ; BiocManager::install("DNAcopy")'

# Docker isn't smart enough to know if something has changed in the github repository at build time.
# To force Docker to build everything after the following line (ARG CACHEBUST=1),
# build using the following command:
#     docker build -t synthex --build-arg CACHEBUST=$(date +%s)

ARG CACHEBUST=1
RUN Rscript -e 'devtools::install_github("alanhoyle/SynthEx")'

COPY Rprofile.site /etc/R/
COPY SynthEx.R /SynthEx.R

CMD ["Rscript", "/SynthEx.R", "--help"]


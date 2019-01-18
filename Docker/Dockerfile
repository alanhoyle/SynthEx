from rocker/r-base

run apt-get update
run apt-get install -t unstable -y libcurl4-openssl-dev  libssl-dev python
# run R -e "install.packages('devtools')"
MAINTAINER Alan Hoyle <alanh@unc.edu


#RUN git clone https://github.com/ropensci/git2r.git &&\
#R CMD INSTALL --configure-args="--with-libssl-include=/usr/lib/" git2r &&\
#rm -rf git2r /tmp/*

ENV bedtools_version 2.27.1

ADD https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools-${bedtools_version}.tar.gz /tmp/

RUN cd /tmp/ \
    && tar xzvf bedtools-${bedtools_version}.tar.gz \
    && cd bedtools2 \
    && make \
    && make install \
    && rm /tmp/bedtools-${bedtools_version}.tar.gz \
    && rm -rf /tmp/bedtools2/test \
    && rm -rf /tmp/bedtools2/src \
    && rm -rf /tmp/bedtools2/docs

RUN R -q -e "install.packages(c('devtools', 'covr', 'roxygen2', 'testthat','optparse'), repos = 'https://cloud.r-project.org/')" &&\
rm -rf /tmp/*

RUN Rscript -e ' install.packages("BiocManager", repos = "https://cloud.r-project.org/") ; BiocManager::install() ; BiocManager::install("DNAcopy")'
RUN Rscript -e 'devtools::install_github("alanhoyle/SynthEx")'

COPY Rprofile.site /etc/R/
COPY SynthEx.R /SynthEx.R

ENTRYPOINT ["Rscript", "/SynthEx.R"]
CMD ["--help"]

# CMD ["R"]

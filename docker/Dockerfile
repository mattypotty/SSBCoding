# BASE IMAGE
FROM debian:latest
ENV HOME="/"

# BASE SOFTWARE
WORKDIR "/"
RUN apt-get update && \
    apt-get upgrade && \
    apt-get install -y git wget zip unzip vim build-essential 

# R
WORKDIR "/"
RUN apt-get install -y gfortran libreadline-dev libbz2-dev liblzma-dev libpcre2-dev libcurl4-openssl-dev
RUN wget https://cran.r-project.org/src/base/R-4/R-4.5.1.tar.gz && \
    tar xfz R-4.5.1.tar.gz; rm R-4.5.1.tar.gz
WORKDIR "/R-4.5.1"
RUN ./configure --with-x=no --enable-R-shlib --enable-memory-profiling && \
    make -s ; make install
ENV PATH /R-4.5.1/bin:$PATH
ENV R_HOME /R-4.5.1
ENV R_LIBS_USER /R-4.5.1/library
ENV LD_LIBRARY_PATH /R-4.5.1/lib:$LD_LIBRARY_PATH

# R PACKAGES
# diversitree: libgsl-dev
RUN apt-get install -y libgsl-dev
# "this.path", "ape", "phangorn", "castor", "caper", "phytools", "diversitree", "MCMCglmm", "coda", "MCMCvis", "mcmcr", "styler"
WORKDIR "/"
RUN R -e "install.packages('this.path',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ape',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('phangorn',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('castor',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('caper',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('phytools',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('diversitree',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('MCMCglmm',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('coda',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('MCMCvis',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('mcmcr',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('styler',dependencies=TRUE,repos='http://cran.rstudio.com/')"

# FINISHING
WORKDIR "/"
COPY dockerfile_SSB.txt /dockerfile_SSB.txt
COPY test.sh /test.sh
COPY start.sh /start.sh
ENV HOME="/"
CMD sleep infinity

# Command for building Docker image:
# cp Dockerfile dockerfile_ssb.txt; docker build -t sswiston/ssb:1 . --progress=plain > ssb.log 2>&1
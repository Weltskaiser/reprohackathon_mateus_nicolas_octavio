##############################################################
# STAGE 1 — compile featureCounts 1.4.6-p3
##############################################################
FROM debian:stretch-slim AS fcbuild
ENV DEBIAN_FRONTEND=noninteractive

# install from  legacy repository (Stretch EOL)
RUN printf 'Acquire::Check-Valid-Until "false";\nAcquire::AllowInsecureRepositories "true";\n' \
      > /etc/apt/apt.conf.d/99archive && \
    printf 'deb http://archive.debian.org/debian stretch main contrib non-free\n' \
      > /etc/apt/sources.list && \
    printf 'deb http://archive.debian.org/debian-security stretch/updates main contrib non-free\n' \
      >> /etc/apt/sources.list && \
    apt-get -o Acquire::Check-Valid-Until=false -o Acquire::AllowInsecureRepositories=true update && \
    apt-get install -y --no-install-recommends \
      build-essential wget ca-certificates tar make zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Compile only featureCounts
RUN wget -O /tmp/subread-1.4.6-p3-source.tar.gz \
      "https://downloads.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-source.tar.gz" && \
    tar -xzf /tmp/subread-1.4.6-p3-source.tar.gz -C /tmp && \
    cd /tmp/subread-1.4.6-p3-source/src && \
    make -f Makefile.Linux featureCounts && \
    install -m 0755 featureCounts /featureCounts && \
    strip /featureCounts || true


##############################################################
# STAGE 2 — final image with R 3.4.1, DESeq2 1.16.1,
#            Python 3.5, Cutadapt 1.11, Bowtie 0.12.7 e featureCounts
##############################################################
FROM rocker/r-ver:3.4.1
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8
ENV R_DEFAULT_INTERNET_TIMEOUT=60

# Acquire files from legacy  Debian Stretch repo
RUN printf 'Acquire::Check-Valid-Until "false";\nAcquire::AllowInsecureRepositories "true";\n' \
      > /etc/apt/apt.conf.d/99archive && \
    printf 'deb http://archive.debian.org/debian stretch main contrib non-free\n' \
      > /etc/apt/sources.list && \
    printf 'deb http://archive.debian.org/debian-security stretch/updates main contrib non-free\n' \
      >> /etc/apt/sources.list && \
    apt-get -o Acquire::Check-Valid-Until=false -o Acquire::AllowInsecureRepositories=true update

# system dependences to compile R + toolchain Python 3.5
RUN apt-get install -y --no-install-recommends \
      build-essential gfortran pkg-config \
      libxml2-dev libssl-dev libcurl4-openssl-dev \
      zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libreadline-dev libicu-dev \
      libblas-dev liblapack-dev \
      wget ca-certificates tar unzip \
      python3.5 python3.5-dev python3-pip && \
    rm -rf /var/lib/apt/lists/*

# pip/setuptools compatible and Cutadapt 1.11
RUN python3.5 -m pip install --upgrade "pip==20.3.4" "setuptools<50" wheel "cython<0.29" && \
    python3.5 -m pip install "cutadapt==1.11"

# Bowtie 0.12.7 (download + symlinks robust)
RUN wget -O /tmp/bowtie.zip \
      "https://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip" && \
    unzip /tmp/bowtie.zip -d /opt && \
    rm /tmp/bowtie.zip && \
    bash -lc 'set -euo pipefail; d=$(find /opt -maxdepth 1 -type d -name "bowtie*"); \
              install -d /usr/local/bin; \
              for f in "$d"/bowtie* ; do [ -x "$f" ] && ln -sf "$f" /usr/local/bin/; done'

# R user library
RUN mkdir -p /usr/local/lib/R/site-library && chmod 777 /usr/local/lib/R/site-library && \
    bash -lc "echo \".libPaths(c('/usr/local/lib/R/site-library', .libPaths()))\" >> /usr/local/lib/R/etc/Rprofile.site"

# Repos (CRAN snapshot legacy + Bioconductor 3.5 oficial)
# CRAN via PPM (snapshot de 2017-10-15 for compatibility with R 3.4.1)
ENV CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2017-10-15
# Bioconductor 3.5
ENV BIOC_SOFT=https://bioconductor.org/packages/3.5/bioc
ENV BIOC_ANN=https://bioconductor.org/packages/3.5/data/annotation
ENV BIOC_EXP=https://bioconductor.org/packages/3.5/data/experiment

# Install dependences for CRAN + BiocInstaller e DESeq2 1.16.1
RUN R -e "options(repos=c(CRAN='${CRAN_SNAPSHOT}', BioCsoft='${BIOC_SOFT}', BioCann='${BIOC_ANN}', BioCexp='${BIOC_EXP}'), \
                 download.file.method='libcurl'); \
          message('Repos usados:'); print(getOption('repos')); \
          invisible(available.packages()); \
          install.packages(c('Rcpp','locfit','BH','RcppArmadillo','Matrix'), type='source'); \
          install.packages('BiocInstaller'); \
          BiocInstaller::biocLite('DESeq2', suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE); \
          if (!'DESeq2' %in% rownames(installed.packages())) stop('DESeq2 não foi instalado')"

# Copy featureCounts for  build
COPY --from=fcbuild /featureCounts /usr/local/bin/featureCounts

# Smoke tests for build
RUN cutadapt --version && \
    bowtie --version && \
    featureCounts -v && \
    R --version && \
    Rscript -e "library(DESeq2); cat('DESeq2 version:', as.character(packageVersion('DESeq2')), '\n')"

WORKDIR /data
CMD ["/bin/bash"]

# FROM python:3.9
FROM continuumio/miniconda3:4.12.0

LABEL \
    author="Shenghui huang" \
    maintainer="Shenghui huang" \
    email="hsh-me@outlook.com" \
    description="Docker Image for ASTK"

# Install system dependencies
RUN apt update -y && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev \
     libreadline-dev build-essential gosu libexpat-dev &&  \
     apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/* && \
     mkdir -p /home/software/astk /project


# install R packages
RUN conda install -y -c https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/ --override-channels \
    r-base=4.1.3 r-upsetr=1.4.0  r-argparse=2.1.6 r-ggnewscale=0.4.7 r-tidyverse=1.3.2 \
    r-ggplot2=3.3.6 r-biocmanager=1.30.18 r-usethis=2.1.6 r-ggraph=2.0.6 && \
    conda install -y -c https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/ --override-channels \
    bedtools=2.30.0 && \
    Rscript -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');\
    BiocManager::install(c('clusterProfiler', 'org.Mm.eg.db', 'tximport', 'ComplexHeatmap', \
        'org.Hs.eg.db', 'simplifyEnrichment','universalmotif', 'memes', 'DESeq2'))" 

# install meme-suite
ADD meme-5.4.1.tar.gz /home/software/
RUN cd /home/software/meme-5.4.1 && eval "$(conda shell.bash hook)" && conda activate base && \
    ./configure --prefix=/home/software/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && make install 

ENV PATH="${PATH}:/home/software/meme/bin:/home/software/meme/libexec/meme-5.4.1"

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
COPY setup.py /home/software/astk/setup.py
COPY astk /home/software/astk/astk

# Install astk
RUN python -m pip install --upgrade pip -i https://mirrors.bfsu.edu.cn/pypi/web/simple && \
    pip install -e /home/software/astk/ -i https://mirrors.bfsu.edu.cn/pypi/web/simple && \
    chmod a+x /usr/local/bin/docker-entrypoint.sh


WORKDIR /project

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

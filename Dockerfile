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


COPY setup.py /home/software/astk/setup.py
COPY astk /home/software/astk/astk

# Install astk
RUN python -m pip install --upgrade pip -i https://mirrors.bfsu.edu.cn/pypi/web/simple && \
    pip install -e /home/software/astk/ -i https://mirrors.bfsu.edu.cn/pypi/web/simple

# install R packages
RUN conda install -y -c https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/ --override-channels \
    r-base=4.1.3 r-upsetr=1.4.0  r-argparse=2.1.6 r-ggnewscale=0.4.7 r-tidyverse=1.3.2 \
    r-ggplot2=3.3.6 r-biocmanager=1.30.18 r-usethis=2.1.6

RUN Rscript -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');\
    BiocManager::install(c('clusterProfiler', 'org.Mm.eg.db', 'tximport', 'ComplexHeatmap', \
        'org.Hs.eg.db', 'simplifyEnrichment','universalmotif', 'memes', 'DESeq2'))" 

# install meme-suite        
ADD meme-5.4.1.tar.gz /home/software/
RUN cd /home/software/meme-5.4.1 && eval "$(conda shell.bash hook)" && conda activate base && \
    ./configure --prefix=/home/software/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && make install 

ENV PATH="${PATH}:/home/software/meme/bin:/home/software/meme/libexec/meme-5.4.1"

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod a+x /usr/local/bin/docker-entrypoint.sh

WORKDIR /project

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]



# bzip2=1.0.8 ncurses=6.3    
# # install R packages
# RUN conda install -y -c bioconda --override-channels  bioconductor-annotationdbi=1.56.1 \
#     bioconductor-complexheatmap=2.10.0 bioconductor-clusterprofiler=4.2.0 \
#     bioconductor-org.mm.eg.db=3.14.0  bioconductor-org.hs.eg.db=3.14.0 \
#     bioconductor-simplifyenrichment=1.4.0 bioconductor-universalmotif=1.12.3 \
#     bioconductor-tximport=1.22.0 meme=5.4.1 bioconductor-memes=1.2.0 


# docker run --rm astk meta -o metadata/fb_e11_based -repN 2 \
#     -p1 data/quant/fb_e11.5_rep*/quant.sf \
#     -p2 data/quant/fb_e1[2-6].5_rep*/quant.sf  data/quant/fb_p0_rep*/quant.sf \
#     -gn fb_e11_12 fb_e11_13 fb_e11_14 fb_e11_15 fb_e11_16 fb_e11_p0

FROM mambaorg/micromamba:1.5.1

RUN micromamba install -y -n base -c conda-forge -c bioconda\
       crispresso2 \
       benchling-sdk \
       psycopg2~=2.9.3 \
       python-dotenv~=0.20.0 \
       pandas~=1.4.3 \
       numpy~=1.23.0 \
       seaborn~=0.11.2 \
       matplotlib~=3.5.2 \
       requests \
       quilt3 \
       altair && \
    micromamba clean --all --yes

COPY / /home/ubuntu/bin/tbOnT/

WORKDIR /home/ubuntu/bin/tbOnT/

CMD python automated_run.py
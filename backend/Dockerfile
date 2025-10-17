FROM continuumio/miniconda3
MAINTAINER Andrey Alekseenko <andrey.alekseenko@stonybrook.edu>

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y ncbi-blast+ gcc g++ cmake gperf check pkg-config && apt-get clean

COPY . /code
RUN /bin/bash ./code/bootstrap.sh && conda clean --all

# set work directory
WORKDIR /code
ENTRYPOINT ["./entrypoint.sh"]

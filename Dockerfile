FROM python:3.12-slim

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential libxml2-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    python -m pip install --upgrade pip && \
    pip install cellxgene

ENTRYPOINT ["cellxgene"]

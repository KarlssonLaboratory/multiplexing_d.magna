FROM rocker/tidyverse:4.5.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    cmake \
    less \
    git \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN which cmake && cmake --version

# Install R-packages
COPY docker/install.R /tmp/install.R
RUN Rscript /tmp/install.R

WORKDIR /data
CMD ["/bin/bash"]

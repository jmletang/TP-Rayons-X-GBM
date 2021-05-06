FROM tbaudier/cours_jm:v1.6

# Adapted from https://github.com/ContinuumIO/docker-images/blob/master/miniconda3/debian/Dockerfile
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# hadolint ignore=DL3008
RUN apt-get update -q && \
    apt-get install -q -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mercurial \
        subversion \
        wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

USER root

COPY . ${HOME}

RUN cp /root/.profile /home/$NB_USER/.profile
RUN cat /home/$NB_USER/.profile
RUN chown --recursive $NB_UID:users /home/$NB_USER/

USER $NB_UID

ENV PATH ${HOME}/miniconda/bin:$PATH

CMD [ "/bin/bash" ]

# Leave these args here to better use the Docker build cache
ARG CONDA_VERSION=py38_4.9.2
ARG CONDA_MD5=122c8c9beb51e124ab32a0fa6426c656

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
    echo "${CONDA_MD5}  miniconda.sh" > miniconda.md5 && \
    if ! md5sum --status -c miniconda.md5; then exit 1; fi && \
    mkdir -p /opt && \
    sh miniconda.sh -b -p ${HOME}/miniconda && \
    rm miniconda.sh miniconda.md5 && \
    find ${HOME}/miniconda/ -follow -type f -name '*.a' -delete && \
    find ${HOME}/miniconda/ -follow -type f -name '*.js.map' -delete && \
    ${HOME}/miniconda/bin/conda clean -afy

# Adapted from https://pythonspeed.com/articles/activate-conda-dockerfile/
RUN conda env create -f ${HOME}/environment.yml
SHELL ["conda", "run", "-n", "example-environment", "/bin/bash", "-c"]

ENV JUPYTER_ENABLE_LAB=yes

USER ${USER}
RUN cat /home/$NB_USER/.profile
SHELL ["conda", "run", "-n", "example-environment", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "example-environment"]

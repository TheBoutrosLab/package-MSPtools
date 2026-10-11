ARG MINIFORGE_VERSION=26.1.1-2
ARG UBUNTU_VERSION=24.04
ARG CONDA_ENV_PATH=/opt/conda/envs/msptools

FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS builder

ARG CONDA_ENV_PATH

RUN mamba create -qy -p ${CONDA_ENV_PATH} \
    -c bioconda \
    -c conda-forge \
    python=3.10 \
    pip \
    pyopenms==3.2.0 && \
    mamba clean -afy

ENV PATH="${CONDA_ENV_PATH}/bin:${PATH}"

WORKDIR /opt/msptools
COPY . .
RUN python -m pip install --no-cache-dir . && \
    python -c "import pyopenms"

FROM ubuntu:${UBUNTU_VERSION} AS final

ARG CONDA_ENV_PATH

RUN apt-get update && \
    apt-get install -y --no-install-recommends libglib2.0-0t64 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder ${CONDA_ENV_PATH} ${CONDA_ENV_PATH}

ENV CONDA_ENV_PATH="${CONDA_ENV_PATH}" \
    OPENMS_DATA_PATH="${CONDA_ENV_PATH}/share/OpenMS" \
    PATH="${CONDA_ENV_PATH}/bin:${PATH}"

RUN python -c "import pyopenms" && \
    msptools --version

RUN groupadd -g 500001 bldocker && \
    useradd -r -u 500001 -g bldocker bldocker
USER bldocker

LABEL org.opencontainers.image.source="https://github.com/TheBoutrosLab/package-MSPtools"
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later"
ENTRYPOINT ["msptools"]

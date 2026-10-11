FROM python:3.10-slim

WORKDIR /opt/msptools
COPY . .
RUN python -m pip install --no-cache-dir .

RUN groupadd -g 500001 bldocker && useradd -r -u 500001 -g bldocker bldocker
USER bldocker

LABEL org.opencontainers.image.source="https://github.com/TheBoutrosLab/package-MSPtools"
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later"
ENTRYPOINT ["msptools"]

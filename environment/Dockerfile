# hash:sha256:b23c353c207e6bdd02f78d047d0831dabb710c6ad3d8f5c64f5dab80d09bcb73
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ENV R_LIBS_USER=/opt/renv/library
ENV RENV_PATHS_LIBRARY=/opt/renv/library
ENV RENV_PATHS_CACHE=/opt/renv/cache
ENV RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall

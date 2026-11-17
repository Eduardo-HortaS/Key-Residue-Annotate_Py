# =========================================================================
# nf-core/kra Dockerfile
# =========================================================================
# Dockerfile for Key Residue Annotate (KRA) pipeline
# Based on micromamba for fast Conda environment resolution
# =========================================================================

FROM mambaorg/micromamba:1.5.8

LABEL authors="Eduardo Horta S."
LABEL description="Docker image containing all software requirements for nf-core/kra pipeline"

# Install environment from environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy Python scripts to project directory
COPY --chown=$MAMBA_USER:$MAMBA_USER *.py /opt/kra/

# Add to PYTHONPATH so modules can import utilities
ENV PYTHONPATH="/opt/kra:$PYTHONPATH"

# Set working directory
WORKDIR /data

# Activate environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import pyhmmer; import Bio; import pandas" || exit 1

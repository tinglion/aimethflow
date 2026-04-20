FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

# ── System packages ───────────────────────────────────────────────────────────
# deadsnakes PPA 提供 Ubuntu 24.04 上的 Python 3.10
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common \
        samtools parallel \
        curl wget ca-certificates \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update && apt-get install -y --no-install-recommends \
        python3.10 python3.10-venv python3.10-dev \
    && rm -rf /var/lib/apt/lists/*

# ── Python venv (3.10.14) ─────────────────────────────────────────────────────
RUN python3.10 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# ── fastp (static binary) ─────────────────────────────────────────────────────
ARG FASTP_VERSION=0.23.4
RUN wget -qO /usr/local/bin/fastp \
        "https://github.com/OpenGene/fastp/releases/download/v${FASTP_VERSION}/fastp" \
    && chmod +x /usr/local/bin/fastp

# ── bwa-mem2 (pre-built bundle) ───────────────────────────────────────────────
ARG BWAMEM2_VERSION=2.2.1
RUN wget -qO- \
        "https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2" \
    | tar -xj -C /usr/local/bin/ --strip-components=1 --wildcards "*/bwa-mem2*"

# ── MethylDackel (static binary) ──────────────────────────────────────────────
# 如果版本有变，到 https://github.com/dpryan79/MethylDackel/releases 确认文件名
ARG METHYLDACKEL_VERSION=0.6.0
RUN wget -qO /usr/local/bin/MethylDackel \
        "https://github.com/dpryan79/MethylDackel/releases/download/${METHYLDACKEL_VERSION}/MethylDackel-${METHYLDACKEL_VERSION}_linux_x86_64" \
    && chmod +x /usr/local/bin/MethylDackel

# ── Python packages ───────────────────────────────────────────────────────────
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir \
        umi_tools \
        bwameth \
        jupyterlab ipykernel \
    && pip install --no-cache-dir -r /tmp/requirements.txt

# ── Application code ──────────────────────────────────────────────────────────
WORKDIR /app
COPY src/ /app/src/

ENV PYTHONPATH=/app/src

# Default: bash; overridden per service in docker-compose.yml
CMD ["bash"]

FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

# ── mamba（conda 加速求解器）─────────────────────────────────────────────────
RUN conda install -n base -c conda-forge mamba -y && conda clean -afy

# ── conda 环境 Python 3.10.14 ─────────────────────────────────────────────────
RUN mamba create -n meth python=3.10.14 -y && conda clean -afy

# ── 生信工具（bioconda）──────────────────────────────────────────────────────
RUN mamba install -n meth -y -c bioconda -c conda-forge \
        fastp=0.23.4 \
        bwameth \
        samtools \
        bedtools \
        methyldackel \
        umi_tools \
        pysam=0.22.0 \
    && conda clean -afy

# ── Jupyter（conda-forge，避免 pip SSL 问题）──────────────────────────────────
RUN mamba install -n meth -y -c conda-forge \
        jupyterlab ipykernel \
    && conda clean -afy

# ── 应用代码 + Python 依赖（pyproject.toml 单一来源）────────────────────────
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src/ ./src/
RUN /opt/conda/envs/meth/bin/pip install --no-cache-dir ".[ml]" -i https://pypi.tuna.tsinghua.edu.cn/simple

ENV PATH="/opt/conda/envs/meth/bin:$PATH"
ENV PYTHONPATH=/app/src

CMD ["bash"]

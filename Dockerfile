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

# ── Python 依赖（pip）────────────────────────────────────────────────────────
SHELL ["conda", "run", "-n", "meth", "/bin/bash", "-c"]
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir \
        jupyterlab ipykernel \
    && pip install --no-cache-dir -r /tmp/requirements.txt

# ── 应用代码 ──────────────────────────────────────────────────────────────────
WORKDIR /app
COPY src/ /app/src/

ENV PATH="/opt/conda/envs/meth/bin:$PATH"
ENV PYTHONPATH=/app/src

CMD ["bash"]

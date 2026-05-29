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

# ── Python 依赖（pyproject.toml 单一来源）───────────────────────────────────
# 拆成两步以充分利用 Docker 层缓存：
#   Step 1 仅复制 pyproject.toml + 最小占位 src/，安装第三方依赖
#          → 只要 pyproject.toml 不变，这一层就命中缓存（即使 src/ 改动也不会重装依赖）
#   Step 2 复制真实 src/，因为我们用 editable 模式安装，无需再次 pip install
#          → src/ 改动只是普通文件复制，重建在秒级完成
WORKDIR /app
COPY pyproject.toml ./
RUN mkdir -p src/scripts src/utils src/methylation \
    && touch src/__init__.py src/scripts/__init__.py src/utils/__init__.py src/methylation/__init__.py \
    && /opt/conda/envs/meth/bin/pip install --no-cache-dir -e ".[ml]" \
         -i https://pypi.tuna.tsinghua.edu.cn/simple \
    && rm -rf src

# ── 应用代码（editable 安装已通过 .pth 文件指向 /app/src，复制真实源码即可）──
COPY src/ ./src/

ENV PATH="/opt/conda/envs/meth/bin:$PATH"
ENV PYTHONPATH=/app/src

CMD ["bash"]

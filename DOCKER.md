# Docker 部署说明

## 设计决策

| 项目 | 决策 |
|------|------|
| 基础镜像 | `continuumio/miniconda3` |
| 生信工具 | conda（bioconda + conda-forge），与服务器环境一致 |
| Python 依赖 | conda env `meth` (3.10.14) + pip |
| GPU | 不启用 |
| 参考基因组 | 外部挂载，不打包 |
| data / temp / result | 外部挂载，不打包 |
| Jupyter | docker-compose 中独立 service，与分析容器共享卷 |

---

## 镜像层结构

```
continuumio/miniconda3
  ├─ [conda] mamba（求解加速）
  ├─ [mamba] create -n meth python=3.10.14    ~400 MB
  ├─ [mamba] bioconda 生信工具                 ~800 MB
  │          fastp / bwameth / samtools /
  │          bedtools / methyldackel / umi_tools / pysam
  ├─ [pip]   requirements.txt + jupyterlab    ~500 MB
  └─ [COPY]  src/                              ~10 MB
             （data / temp / result 运行时挂载）
```

---

## 工具安装方式

### conda（bioconda + conda-forge）

| 工具 | 版本 | 用途 |
|------|------|------|
| fastp | 0.23.4 | 质控 trimming |
| bwameth | latest | 亚硫酸盐比对 |
| samtools | latest | BAM 处理 |
| bedtools | latest | 区域过滤 |
| methyldackel | latest | 甲基化提取 |
| umi_tools | latest | UMI 去重 |
| pysam | 0.22.0 | Python BAM 读取 |

### pip（conda env 内）

```
jupyterlab ipykernel   # Jupyter 环境
-r requirements.txt    # 项目 Python 依赖
```

---

## 挂载卷规划

| 宿主机路径（示例）       | 容器内路径          | 只读 | 说明                      |
| ----------------------- | ------------------- | ---- | ------------------------- |
| `./data`                | `/app/data`         | 是   | panel BED、样本 Excel、候选区域 |
| `./config.toml`         | `/app/config.toml`  | 是   | 项目配置                  |
| `/data/meth/OC_tissue`  | `/mnt/OC_tissue`    | 是   | 组织 bedGraph             |
| `/data/meth/OC_plasma`  | `/mnt/OC_plasma`    | 是   | 血浆 BAM / bedGraph       |
| `/data/meth/ref`        | `/mnt/ref`          | 是   | hg38 + bwa-meth 索引      |
| `/data/meth/temp`       | `/mnt/temp`         |      | 中间文件                  |
| `/data/meth/result`     | `/mnt/result`       |      | 分析结果                  |

`config.toml` 中路径需改为容器内路径，例如：

```toml
> tissue_dir = "/mnt/OC_tissue"
> plasma_dir = "/mnt/OC_plasma"
> reference_genome = "/mnt/ref/hg38.fa"
> [directories]
> temp_dir   = "/mnt/temp"
> result_dir = "/mnt/result"
> ```

---

## 端口

| 端口 | 服务 |
|------|------|
| 8888 | JupyterLab |

---

## docker-compose 服务规划

```yaml
services:
  analysis:          # 全流程分析（python run_analysis.py）
    build: .
    volumes: [data, config, OC_tissue, OC_plasma, ref, temp, result]
    command: python /app/src/scripts/run_analysis.py --steps 1-19

  jupyter:           # 分步查看结果
    build: .
    ports: ["8888:8888"]
    volumes: [data, config, OC_tissue, OC_plasma, ref, temp, result, notebooks]
    command: jupyter lab --ip=0.0.0.0 --no-browser --allow-root
```

---

## 环境变量

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `NCBI_API_KEY` | *(空)* | PubMed 速率提升（可选） |

---

## 文件清单（待创建）

```
meth/
├── Dockerfile            # 镜像构建
├── docker-compose.yml    # 一键启动（analysis + jupyter）
├── .dockerignore         # 排除 data/ temp/ result/ __pycache__
├── DOCKER.md             # 本文件
└── notebooks/            # Jupyter 分步分析 notebook
    ├── 01_qc.ipynb
    ├── ...
    └── 19_report.ipynb
```

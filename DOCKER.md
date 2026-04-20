# Docker 部署说明

## 设计决策

| 项目 | 决策 |
|------|------|
| 基础镜像 | `ubuntu:24.04` |
| Python 依赖 | pip only，无 conda |
| GPU | 不启用 |
| 参考基因组 | 外部挂载，不打包 |
| data / temp / result | 外部挂载，不打包 |
| Jupyter | docker-compose 中独立 service，与分析容器共享卷 |

---

## 镜像层结构

```
ubuntu:24.04
  ├─ [apt]  python3.12 + 编译依赖 + 系统库           ~200 MB
  ├─ [apt/binary]  fastp / bwa-mem2 / samtools       ~150 MB
  │              / MethylDackel / parallel
  ├─ [pip]  requirements.txt + umi_tools + bwameth   ~600 MB
  │         + jupyterlab
  └─ [COPY] src/                                      ~10 MB
            （data / temp / result 运行时挂载）
```

---

## 工具安装方式

### apt（Ubuntu 24.04 官方仓库）

```
samtools          # ≥1.19，Ubuntu 24.04 仓库自带
parallel          # gnu-parallel
zlib1g-dev libcurl4-openssl-dev libssl-dev  # 编译依赖
```

### 预编译二进制（GitHub Release 直接下载）

| 工具 | 来源 |
|------|------|
| fastp | OpenGene/fastp releases |
| bwa-mem2 | bwa-meth/bwa-mem2 releases |
| MethylDackel | dpryan79/MethylDackel releases |

### pip install

```
umi_tools bwameth          # 生信 Python 工具
jupyterlab ipykernel       # Jupyter 环境
-r requirements.txt        # 项目依赖
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

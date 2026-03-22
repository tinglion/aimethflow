# 生信项目 - 甲基化分析

## 项目简介

本项目用于分析甲基化相关的生信数据，包含一系列脚本工具用于处理和分析FastQ格式的测序数据，以及进行甲基化分析和可视化。

## 数据目录结构

数据根目录: `F:\meth`

包含以下路径:

- `OC_tissue`: 组织样本数据
- `OC_plasma`: 血浆样本数据

## 环境配置

### 重要说明
本项目**仅支持在Linux环境或Windows Subsystem for Linux (WSL)**中运行，因为生物信息学工具（如samtools、Bismark等）在Windows原生环境中可能不可用。

### 1. 在Linux环境中设置

#### 1.1 创建conda环境

```bash
# 创建名为meth的conda环境（Python 3.10）
conda create -n meth -y python=3.10

# 激活环境
conda activate meth

# 安装生物信息学工具和Python依赖
conda install -y -c bioconda -c conda-forge \
    fastp=0.23.4 \
    bwa-mem2=2.2.1 \
    samtools=1.19 \
    pysam=0.22.0 \
    libdeflate=1.20 \
    htslib=1.20 \
    pandas=2.2.0 \
    numpy=1.26.2 \
    matplotlib=3.8.2 \
    seaborn=0.13.0 \
    biopython=1.84
   #  methyldackel=0.6.2 \
   #  gnu-parallel=20230422 \

# 安装项目及其依赖
pip install -e .
```

### 2. 在Windows电脑上使用WSL

### 3. 安装Bismark（用于真实甲基化分析）

在Linux或WSL环境中安装：

```bash
# 安装Bismark及其依赖 bioconda
conda install  bismark bowtie2 samtools
conda install -y -c bioconda bwa=0.7.17 methyldackel=1.0.0
```

### 4. 参考基因组准备（仅需 1 次）

```bash
# 创建参考基因组目录
mkdir -p /mnt/f/meth/hg38
cd /mnt/f/meth/hg38

# 下载hg38参考基因组（WSL2中直接下载，比Windows快）
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# 构建bwa-mem2索引（32线程，约30分钟）
export BWA_MEM2_MAX_MEM=40G
bwa-mem2 index -p hg38 hg38.fa

# bwa index默认用IS算法，内存占用≤8GB，WSL2无压力
bwa index -p hg38 hg38.fa

# 验证索引（生成5个索引文件）
ls hg38.*.bt2*  # 能看到hg38.0.bt2、hg38.1.bt2等文件
```

## 开发指南

### 项目结构

```txt
meth/
├── src/
│   ├── methylation/       # 甲基化分析核心模块
│   ├── scripts/           # 命令行脚本工具
│   └── utils/             # 工具模块
├── config.toml            # 配置文件
├── pyproject.toml         # Python 包配置
├── README.md              # 项目文档
└── requirements.txt       # 依赖列表
```

## 生信分析工具和功能

### 1. 数据质量控制

#### fq-inspect

快速查看.fq.gz文件的统计信息，包括序列数量、碱基数、长度分布和质量分数等。

**使用方法:**

```bash
fq-inspect <file_path>
```

**示例:**

```bash
fq-inspect "F:\meth\OC_plasma\B10001.R1.fq.gz"
```

#### 质量控制工具

- **FastQC**: 用于对测序数据进行质量评估
- **MultiQC**: 整合多个样本的质量控制报告

### 2. 数据预处理

#### fq-dump

用于从.fq.gz文件中提取序列数据的脚本，支持指定开始位置和读取数量。

**使用方法:**

```bash
fq-dump <file_path> [--start <position>] [--count <number>]
```

**参数:**

- `--start`: 读取开始位置，默认为0
- `--count`: 读取数量，默认为10

**示例:**

```bash
# 从位置5开始，读取3条序列
fq-dump "F:\meth\OC_plasma\B10001.R1.fq.gz" --start 5 --count 3
```

#### 预处理工具

- **Trimmomatic**: 用于修剪低质量的 reads 和接头序列
- **Cutadapt**: 用于移除接头序列和低质量碱基

### 3. 序列比对

#### 常用比对工具

- **Bismark**: 专门用于甲基化测序数据的比对工具
- **Bowtie2**: 快速的通用序列比对工具
- **BSMAP**: 用于 bisulfite 处理后测序数据的比对

### 4. 甲基化分析

#### 分析工具

- **Bismark methylation extractor**: 从比对结果中提取甲基化信息
- **MethylKit**: R包，用于甲基化数据分析和差异甲基化检测
- **BSseq**: Bioconductor包，用于处理和分析 bisulfite 测序数据

### 5. 差异甲基化分析

#### 分析工具

- **DMRfinder**: 用于鉴定差异甲基化区域
- **methylSig**: 用于差异甲基化位点和区域的统计分析
- **edgeR**: 用于差异表达分析，也可用于差异甲基化分析

### 6. 结果可视化

#### 可视化工具

- **IGV**: 基因组浏览器，用于查看甲基化数据
- **methylKit**: 提供甲基化数据的可视化功能
- **ggplot2**: R包，用于创建自定义的可视化图表

### 7. 下游分析

#### 功能分析

- **GO分析**: 基因本体论分析
- **KEGG通路分析**: 分析差异甲基化基因富集的通路
- **基因注释**: 将甲基化位点与基因关联起来

### 8. 脚本工具开发

本项目已开发和集成了以下脚本工具，以支持甲基化分析的各个环节。

#### 已开发的工具

##### methylation-extractor

从Bismark等工具的比对结果中提取甲基化信息。

**使用方法:**

```bash
methylation-extractor <input_file> <output_file>
```

**参数:**

- `input_file`: 比对结果文件路径（.bam格式）
- `output_file`: 输出文件路径（.bedmethyl或.txt格式）

**示例:**

```bash
methylation-extractor sample.bam methylation_results.txt
```

#### 真实甲基化提取流程

要进行真实的甲基化提取，需要先将FastQ文件比对到参考基因组上生成BAM文件。以下是完整的流程：

1. **安装Bismark**

使用conda安装Bismark及其依赖：

```bash
# 安装Bismark及其依赖
conda install -c bioconda bismark bowtie2 samtools
```

2. **下载参考基因组**

```bash
# 下载人类参考基因组（hg38）
bismark_genome_preparation --path_to_genome /path/to/genome
```

3. **比对FastQ文件**

```bash
# 单端测序数据
bismark --genome /path/to/genome <read1.fq.gz>

# 双端测序数据
bismark --genome /path/to/genome -1 <read1.fq.gz> -2 <read2.fq.gz>
```

4. **提取甲基化信息**

```bash
# 使用Bismark的甲基化提取工具
bismark_methylation_extractor --bedGraph --counts --gzip sample.bam

# 或使用本项目的methylation-extractor工具
methylation-extractor sample.bam methylation_results.txt
```

5. **运行完整分析**

```bash
# 使用run-analysis命令运行完整分析流程
run-analysis F:\meth\OC_tissue F:\meth\OC_plasma
```

#### 注意事项

- 真实的甲基化分析需要大量的计算资源和时间，特别是比对步骤
- 对于大型FastQ文件，建议使用高性能计算集群或云服务
- 本项目的run-analysis脚本会自动检查是否存在BAM文件，如果不存在则使用模拟数据
- 要获得准确的甲基化分析结果，建议使用真实的BAM文件

##### diff-methylation

比较两个样本的甲基化数据，找出差异甲基化位点。

**使用方法:**

```bash
diff-methylation <sample1_file> <sample2_file> <output_file> [--threshold <value>]
```

**参数:**

- `sample1_file`: 第一个样本的甲基化数据文件
- `sample2_file`: 第二个样本的甲基化数据文件
- `output_file`: 差异甲基化分析结果输出文件
- `--threshold`: 甲基化差异阈值，默认为0.2

**示例:**

```bash
diff-methylation sample1_methylation.txt sample2_methylation.txt diff_result.txt --threshold 0.3
```

##### methylation-visualizer

可视化甲基化数据，生成甲基化水平分布图。

**使用方法:**

```bash
methylation-visualizer <input_file> <output_file> [--type <plot_type>]
```

**参数:**

- `input_file`: 甲基化数据文件路径
- `output_file`: 输出图片文件路径（支持.png, .pdf, .svg格式）
- `--type`: 图表类型，可选值：density, scatter, heatmap，默认为density

**示例:**

```bash
# 生成密度图
methylation-visualizer methylation_results.txt density_plot.png --type density

# 生成散点图
methylation-visualizer methylation_results.txt scatter_plot.png --type scatter

# 生成热图
methylation-visualizer methylation_results.txt heatmap_plot.png --type heatmap
```

## 配置文件

项目使用 `config.toml` 文件进行配置，存储项目相关的配置信息。以下是配置文件的主要部分：

### 配置文件结构

```toml
# 项目基本信息
name = "甲基化分析项目"
version = "1.0.0"

# 数据目录配置
[data]
root_dir = "F:\\meth"
tissue_dir = "F:\\meth\\OC_tissue"
plasma_dir = "F:\\meth\\OC_plasma"

# 目录配置
[directories]
temp_dir = "temp"
result_dir = "result"
log_dir = "logs"

# FastQ分析配置
[fastq]
progress_interval = 1000000  # 每处理多少条序列显示进度
sample_size = 1000          # 随机抽样大小
random_seed = 42            # 随机种子，确保结果可重复

# 甲基化提取配置
[extraction]
min_quality = 20            # 最小甲基化质量分数
min_depth = 10              # 最小覆盖深度

# 差异甲基化分析配置
[diff_methylation]
threshold = 0.2             # 甲基化差异阈值
p_value = 0.05              # 统计显著性水平
min_depth = 20              # 最小覆盖深度

# 可视化配置
[visualization]
plot_types = ["density", "scatter", "heatmap"]  # 图表类型
fig_width = 12              # 图表宽度
fig_height = 8              # 图表高度
dpi = 300                   # DPI
font_size = 12              # 字体大小

# 日志配置
[logging]
level = "INFO"               # 日志级别
format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"  # 日志格式
save_to_file = true         # 是否保存到文件
console_output = true       # 是否输出到控制台
```

### 配置说明

1. **数据目录配置**：设置数据的根目录和各个样本类型的目录路径。
2. **目录配置**：设置临时文件、结果文件和日志文件的存储目录。
3. **FastQ分析配置**：
   - `sample_size`：设置FastQ文件分析时的随机抽样大小，默认1000条序列
   - `random_seed`：设置随机种子，确保抽样结果可重复
4. **甲基化提取配置**：设置提取甲基化信息时的质量控制参数。
5. **差异甲基化分析配置**：设置差异分析的阈值和统计参数。
6. **可视化配置**：设置结果可视化的图表参数。
7. **日志配置**：设置日志的级别和格式。

### 随机抽样说明

为了高效处理大型FastQ文件，项目采用了随机抽样的方法进行质量分析：

1. **抽样原理**：从FastQ文件中随机抽取指定数量的序列进行分析，而不是分析整个文件。
2. **抽样参数**：
   - `sample_size`：控制抽样的序列数量，默认为1000条
   - `random_seed`：确保抽样结果可重复，默认为42
3. **优势**：
   - 显著减少分析时间，特别是对于大型FastQ文件
   - 提供代表性的质量统计信息
   - 支持对多个文件进行快速分析

通过修改 `config.toml` 文件中的 `[fastq]` 部分，可以调整抽样参数以适应不同的分析需求。

## 完整分析流程

使用 `run-analysis` 命令可以执行完整的甲基化分析流程：

```bash
run-analysis
```

该命令会按照配置文件中的设置，执行数据质量控制、预处理、甲基化提取、差异分析和可视化等步骤。

## 注意事项

- 本项目使用 Python 3.10 版本
- 所有脚本都需要在激活虚拟环境后运行
- 数据文件较大时，建议使用分块加载的方式处理，避免内存不足
- 执行分析时，临时数据文件会统一放到 `temp` 目录下，重要结果文件放到 `result` 目录下，log 放到 `logs` 目录下，每次执行都单独创建一个带时间戳的子目录

# 卵巢癌血浆 cfDNA 靶向甲基化分析

卵巢癌早期诊断项目——基于血浆 cfDNA 的靶向甲基化测序分析流程，涵盖从原始 FastQ 到机器学习诊断分类的全流程（17 个步骤）。

## 分析目标

- 检测血浆游离 DNA（cfDNA）中与卵巢癌相关的 CpG 位点差异甲基化
- 利用 88 个候选区域（来自组织样本分析）构建诊断特征 panel
- 训练良恶性二分类模型（LR / SVM / RF / XGB / LGBM，LOO-CV 评估）
- 对恶性 + 交界性样本进行亚型多分类（高级别浆液性癌 / 子宫内膜样癌 / 透明细胞癌 / 低级别浆液性癌 / 黏液性）

---

## 环境要求

> **仅支持 Linux / WSL2**，生信工具（bismark、bwameth、MethylDackel 等）不支持 Windows 原生环境。

### 1. 创建 conda 环境

```bash
conda create -n meth python=3.10 -y
conda activate meth

# 生信工具（bioconda）
conda install -y -c bioconda -c conda-forge \
    fastp=0.23.4 \
    bwameth \
    samtools \
    bedtools \
    methyldackel \
    bismark \
    bowtie2 \
    umi_tools \
    pysam=0.22.0

# Python 依赖
pip install -e ".[ml]" -i https://pypi.tuna.tsinghua.edu.cn/simple
```

`[ml]` 包含 `xgboost` 和 `lightgbm`；不需要 ML 步骤时可省略：`pip install -e .`

### 2. 参考基因组索引（首次运行）

```bash
# bwameth 索引（血浆流程，推荐）
bwameth.py index /data/sting/meth/hg38.fa

# bismark 索引（组织流程）
bismark_genome_preparation /data/sting/meth/
```

### 3. 辅助文件

| 文件 | 说明 |
| --- | --- | | 探针捕获区域 BED |
| `data/OC_tissue_88_regions.hg38.txt` | 候选甲基化区域（88 个） |
| `data/卵巢癌fq编号与样本信息对应表v3.xlsx` | 样本分组 Excel（含分类/亚型列） |
| `/mnt/f/meth/hg38.repeat.bed` | hg38 重复序列 BED |
| `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` | 启动子区域 BED |

---

## 项目结构

```text
├── config.toml                    # 全局配置
├── requirements.txt               # pip 依赖
├── pyproject.toml                 # 包配置
├── data/                          # 参考数据（BED、Excel）
├── src/
│   ├── methylation/               # 核心分析模块
│   │   ├── diff_methylation.py    # 差异甲基化分析
│   │   ├── sample_grouping.py     # Excel 样本分组
│   │   └── methylation_visualizer.py
│   ├── scripts/
│   │   ├── run_analysis.py        # 主入口（CLI）
│   │   └── modules/
│   │       ├── alignment_utils.py     # 比对 / 去重 / 靶区过滤 / 提取
│   │       ├── fastq_utils.py         # fastp 预处理
│   │       ├── feature_utils.py       # 特征选择（Step 14）
│   │       ├── classifier_utils.py    # 诊断分类（Step 15）
│   │       ├── subtype_utils.py       # 亚型分类（Step 16）
│   │       └── report_utils.py        # 报告生成（Step 17）
│   ├── test/
│   │   └── test_parse_sample_info.py
│   └── utils/
│       ├── config.py
│       └── logger.py
└── docs/
```

---

## 分析流程（17 步）

### 血浆模式（`--plasma`，推荐）

| 步骤 | Key | 说明 | 工具 |
| --- | --- | --- | --- |
| 4 | `preprocess` | fastp QC + UMI 提取 | fastp |
| 5 | `align` | 双端比对 | bwameth（默认）/ bismark |
| 6 | `dedup` | UMI 去重 | umi_tools dedup --paired |
| 7 | `target` | On-target 过滤 | bedtools intersect |
| 7.5 | `mbias` | M-bias 检查（可选） | MethylDackel mbias |
| 8 | `extract` | CpG 甲基化提取 | MethylDackel extract |
| 9 | `filter` | 过滤（重复序列/深度/候选区域） | 自研脚本 + bedtools |
| 10 | `diff` | 差异甲基化分析（Wilcoxon） | 自研脚本 |
| 11 | `visualize` | 密度图 / 散点图 / 热图 | matplotlib |
| 14 | `feature_select` | 特征选择（LASSO + RF） | scikit-learn |
| 15 | `classify` | 诊断分类器 LOO-CV | LR / SVM / RF / XGB / LGBM |
| 16 | `subtype_classify` | 亚型多分类 LOO-CV | RF + SVM OvR |
| 17 | `report` | 汇总 HTML + JSON 报告 | 自研脚本 |

### 组织模式（无 `--plasma`）

步骤 4–11 可用（bismark + deduplicate_bismark）；步骤 14–17 需要 Excel 含良/恶性分组。

---

## 快速开始

### 全流程一键运行

```bash
run-analysis --plasma /mnt/f/meth/OC_plasma/
```

### 分步执行

```bash
# Step 4：预处理（生成 workdir，后续步骤复用）
run-analysis --plasma --step preprocess /mnt/f/meth/OC_plasma/

# Step 5–17：指定 workdir 依次执行
run-analysis --plasma --step align          --workdir temp/20260406
run-analysis --plasma --step dedup         --workdir temp/20260406
run-analysis --plasma --step target        --workdir temp/20260406
run-analysis --plasma --step mbias         --workdir temp/20260406
run-analysis --plasma --step extract       --workdir temp/20260406
run-analysis --plasma --step filter        --workdir temp/20260406
run-analysis --plasma --step diff          --workdir temp/20260406
run-analysis --plasma --step visualize     --workdir temp/20260406
run-analysis --plasma --step feature_select --workdir temp/20260406
run-analysis --plasma --step classify      --workdir temp/20260406
run-analysis --plasma --step subtype_classify --workdir temp/20260406
run-analysis --plasma --step report        --workdir temp/20260406

# 查看所有步骤及当前进度
run-analysis --list-steps --workdir temp/20260406
```

### 常用选项

| 选项 | 说明 |
| --- | --- |
| `--plasma` | 血浆模式（bwameth + UMI dedup） |
| `--step <name>` | 只运行指定步骤（分步模式） |
| `--workdir <path>` | 指定已有临时目录（分步模式必须） |
| `--skip-report` | 跳过报告生成 |
| `--clean-temp` | 完成后清理临时目录 |
| `--list-steps` | 打印步骤列表及完成状态 |
| `--config <path>` | 使用指定配置文件 |

---

## 配置文件 `config.toml`

### 关键路径

```toml
[data]
plasma_dir         = "/mnt/f/meth/OC_plasma"
sample_info_excel  = "data/卵巢癌fq编号与样本信息对应表v3.xlsx"
probe_panel_bed    = "data/panel_B10308.bed"
candidate_regions  = "data/OC_tissue_88_regions.hg38.txt"

[extraction]
reference_genome   = "/data/sting/meth/hg38.fa"

[tools]
conda_env_bin      = "/path/to/miniconda3/envs/meth/bin"
```

### 血浆参数 `[plasma]`

```toml
aligner       = "bwameth"   # 或 "bismark"
umi_length    = 6
trim_front1   = 12
trim_front2   = 12
trim_tail2    = 5           # R2 3' poly-C 裁剪，根据 M-bias SVG 调整
nOT           = "0,0,0,20"  # MethylDackel OT 链位置偏差校正
nOB           = "0,0,17,0"  # MethylDackel OB 链位置偏差校正
diff_delta_beta = 0.10      # 差异分析 |Δβ| 阈值
```

> **M-bias 调整建议**：先运行 `--step mbias` 查看 SVG，若 R2 3' 端甲基化率骤降，相应增大 `trim_tail2` / `nOT` / `nOB`，然后从 `extract` 步骤重跑（不需要重新比对）。

### 中间文件清理 `[cleanup]`

```toml
delete_after_dedup  = false   # dedup 后是否删除 sorted.bam
delete_after_target = false   # target 后是否删除 dedup.bam
```

### 线程配置 `[threads]`

```toml
fastp           = 8
bwameth         = 8
samtools_sort   = 4
methyldackel    = 4
```

### 机器学习 `[ml]`

```toml
# 特征选择
feature_selection_mode = "union"   # "union" | "lasso" | "rf"
min_feature_coverage   = 0.5       # CpG 位点最低样本覆盖率
lasso_cv_folds         = 5

# 分类器参数
rf_n_estimators   = 500
lr_C              = 1.0
svm_C             = 1.0
xgb_n_estimators  = 200
xgb_max_depth     = 4
xgb_learning_rate = 0.1
lgbm_n_estimators = 200
lgbm_max_depth    = 4
lgbm_learning_rate = 0.1
```

---

## 输出文件

| 目录 | 文件 | 说明 |
| --- | --- | --- |
| `result/` | `analysis_report.html` | 全流程汇总报告（含可选 ML 节） |
| `result/` | `analysis_report.json` | 结构化报告数据 |
| `result/` | `diff_methylation.txt` | 差异甲基化位点列表 |
| `result/` | `feature_importance.csv` | LASSO + RF 特征得分 |
| `result/` | `panel_cpg.csv` | 最终入选 CpG 位点 |
| `result/` | `classifier_summary.csv` | 5 个分类器 AUC/Acc/Sens/Spec |
| `result/` | `classifier_predictions.csv` | 逐样本预测结果 |
| `result/` | `roc_curve.png` | ROC 曲线对比图 |
| `result/` | `confusion_matrix.png` | 混淆矩阵（2 行 × 3 列） |
| `result/` | `subtype_summary.csv` | 亚型分类 RF + SVM 汇总 |
| `result/` | `subtype_predictions.csv` | 亚型逐样本预测 |
| `result/` | `subtype_roc.png` | 亚型 OvR ROC 曲线 |
| `temp/` | `*.bwameth_pe.sorted.bam` | 比对 BAM |
| `temp/` | `*.dedup.bam` | 去重 BAM |
| `temp/` | `*_CpG.bedGraph` | 单样本甲基化 bedGraph |
| `temp/` | `*.mbias_OT.svg` | M-bias 图（OT 链） |

---

## 样本 Excel 格式

`data/卵巢癌fq编号与样本信息对应表v3.xlsx` 须包含以下列（支持多 sheet）：

| 列名 | 内容 |
| --- | --- |
| `fq_R1` | R1 文件名，如 `B10001.R1.fq.gz` |
| `fq_R2` | R2 文件名 |
| `分类` | `良性` / `恶性` / `交界性` |
| `亚型` | `高级别浆液性癌` 等（恶性/交界性填写） |

- `交界性` 不参与二分类训练，但参与亚型分类
- 同一样本在多个 sheet 出现时自动去重

---

## 辅助工具

```bash
# 快速查看 FastQ 统计
fq-inspect /mnt/f/meth/OC_plasma/B10001.R1.fq.gz

# 从 BAM 文件生成 panel.bed（若无捕获 BED）
bedtools genomecov -ibam sample.dedup.bam -bg \
  | awk '$4>10' | bedtools merge -d 100 \
  | awk 'length($1)<6' > panel.bed

# 单元测试
pytest src/test/ -v
```

---

## 注意事项

- `pysam` 仅支持 Linux/WSL，Windows 本机无法安装
- XGBoost / LightGBM 为可选依赖，未安装时自动跳过，其余分类器正常运行
- 步骤间状态通过 `temp/<workdir>/pipeline_state.json` 传递，重跑单步时不需要重新指定输入文件
- dedup / target 步骤均支持幂等重跑：输出文件已存在时自动跳过
- 修改 `nOT`/`nOB`/`trim_tail2` 后只需从 `--step extract` 重跑，不需要重新比对

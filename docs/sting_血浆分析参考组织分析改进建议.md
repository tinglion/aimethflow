# 血浆 cfDNA 靶向甲基化分析：完整流程指南

> 综合参考：优讯血浆流程（实验结果欠佳）、Crist 标准 SOP、组织分析经验，结合实验室实际建库情况整理。  
> **本文档与代码实现同步，所有命令均已在 conda `meth` 环境下测试。**  
> 最后更新：2026-04-17

---

## 一、建库前提确认

在分析流程设计之前，必须明确以下建库信息，它们直接决定多个关键步骤的参数。

| 参数 | 本实验室情况 | 影响步骤 |
| ---- | ----------- | ------- |
| UMI | **有，6N，位于 R2 5' 端** | Step 4 fastp 提取、Step 6 去重方式 |
| 建库方式 | **单链建库（single-strand），高度定向（OB 占 1.2%）** | Step 5 比对参数（`non_directional = false`） |
| 末端修复 | **无末端修复；P7 接头引入 poly-C（未甲基化）** | Step 4 `trim_tail2`、Step 7.5 M-bias 检查 |
| 5' 端特殊序列 | **R1/R2 5' 端各有 12 个低复杂度寡核苷酸** | Step 4 `trim_front1/2` |
| 参考基因组 | hg38（`/data/meth/bwa/hg38.fa`） | 全流程 |
| bwameth c2t 索引 | `/data/meth/bwa/hg38.fa.bwameth.c2t`（已构建，约 17 GB） | Step 5 |
| Panel BED | `data/panel_B10308.bed` | Step 7 On-target 过滤 |
| 候选区域 | `data/OC_tissue_88_regions.hg38.txt`（88 个区域） | Step 9 过滤 |

---

## 二、优讯方案缺陷与本流程改进对比

### 根本问题汇总

| 问题 | 优讯方案 | 本流程方案 | 严重程度 |
| ---- | ------- | --------- | ------- |
| **去重方式** | Picard 坐标去重——cfDNA 片段集中在 147–167 bp，大量真实独立片段被误删 | **umi_tools dedup --paired**，基于 UMI 识别真实 PCR 重复 | ⚠⚠⚠ 最严重 |
| **UMI 未利用** | 有 UMI 建库却完全忽略，等于浪费 | **fastp --umi** 提取，umi_tools 去重 | ⚠⚠⚠ |
| **缺少 On-target 过滤** | 无 bedtools intersect 步骤，全基因组脱靶噪声未去除 | **bedtools intersect panel.bed** 保留靶区 reads | ⚠⚠ |
| **poly-C 末端未处理** | 仅切除 5' 端，R2 3' poly-C（来自 P7 接头）未处理 | `trim_tail2=5`，通过 M-bias SVG 确认后调整 | ⚠⚠ |
| **甲基化阈值过严** | 血浆用 80%，比组织 60% 更严苛，ctDNA 信号极易消失 | **60%，与组织保持一致** | ⚠ |
| **差异分析缺少绝对差值过滤** | 仅 P<0.05，无 \|Δβ\| 过滤，假阳性多 | **P<0.05 且 \|Δβ\|>0.1** | ⚠ |

### 方法学对比

| 环节 | 优讯流程 | 本流程 |
| ---- | ------- | ------ |
| QC 工具 | fastp（Q15，宽松） | fastp（Q20，min_length=35，含 UMI 提取） |
| 5'/3' 端处理 | trim_front=12，无 trim_tail | trim_front1/2=12，trim_tail2=5（M-bias 确认） |
| 比对工具 | bwameth（directional，可能丢失 50% 数据） | bwameth（`non_directional=false`，已确认高度定向） |
| 去重 | Picard 坐标去重 ❌ | umi_tools dedup --paired ✅ |
| On-target 过滤 | 无 ❌ | bedtools intersect panel.bed ✅ |
| M-bias 检查 | 无 | MethylDackel mbias ✅ |
| 甲基化提取 | Biscuit（read-level） | MethylDackel extract（`--mergeContext`，cfDNA 必须） |
| 甲基化阈值 | 80% | 60%（与组织统一）✅ |
| 差异分析 | 仅 P<0.05 | P<0.05 且 \|Δβ\|>0.1 ✅ |

---

## 三、完整分析流程（17 步）

### 流程概览

```text
Raw FASTQ (paired-end)
  │
  ▼ Step 4   fastp QC + UMI 提取（trim_front=12/12，trim_tail2=5，Q20，min_length=35）
  ▼ Step 5   bwameth 双端比对（8 线程，两步法避免管道死锁）
  ▼ Step 6   umi_tools dedup --paired（UMI 感知去重）
  ▼ Step 7   bedtools intersect panel.bed（On-target 过滤）
  ▼ Step 7.5 MethylDackel mbias（M-bias 检查，确认 trim_tail2）
  ▼ Step 8   MethylDackel extract --mergeContext（CpG.bedGraph 提取）
  ▼ Step 9   深度/重复序列过滤 + 聚焦 88 候选区域
  ▼ Step 10  差异分析（Wilcoxon，P<0.05 且 |Δβ|>0.1）
  ▼ Step 11  可视化（密度图 / 散点图 / 热图）
  ▼ Step 14  特征选择（LASSO + RF 并集，→ panel_cpg.csv）
  ▼ Step 15  诊断分类器（LR / SVM / RF / XGB / LGBM，LOO-CV）
  ▼ Step 16  亚型分类（RF + SVM OvR，恶性+交界性多分类，LOO-CV）
  ▼ Step 17  生成分析报告（HTML + JSON）
```

### 各步骤详细说明

#### Step 4：fastp QC + UMI 提取

```bash
fastp \
  -i sample_R1.fq.gz -I sample_R2.fq.gz \
  -o sample.clean.R1.fq.gz -O sample.clean.R2.fq.gz \
  -q 20 --length_required 35 \
  --umi --umi_loc read2 --umi_len 6 --umi_prefix UMI \
  --trim_front1 12 --trim_front2 12 \
  --trim_tail2 5 \
  --thread 8 \
  -j sample.fastp.json -h sample.fastp.html
```

> `trim_tail2=5` 是初始值，须在 Step 7.5 M-bias 结果出来后确认或调整。

#### Step 5：bwameth 双端比对（两步法）

大样本（>20M pairs）用两步法避免 stdout 管道死锁：

```bash
# Step 1：比对，流式转为 unsorted BAM
bwameth.py --reference /data/meth/bwa/hg38.fa -t 8 \
  sample.clean.R1.fq.gz sample.clean.R2.fq.gz \
| samtools view -bS -o sample.unsorted.bam -

# Step 2：单独排序（bwameth 已退出，无死锁）
samtools sort -@ 4 -m 8G -o sample.bwameth_pe.sorted.bam sample.unsorted.bam
rm sample.unsorted.bam
samtools index sample.bwameth_pe.sorted.bam
```

> bwameth c2t 索引构建（一次性，约 2.5 h）：`bwameth.py index /data/meth/bwa/hg38.fa`  
> `samtools sort -m 8G`：4 线程共 32 GB 内存，防止大样本 OOM。

#### Step 6：UMI 感知去重

```bash
umi_tools dedup --paired \
  -I sample.bwameth_pe.sorted.bam \
  -S sample.bwameth_pe.dedup.bam \
  --log sample.bwameth_pe.dedup.stats
samtools index sample.bwameth_pe.dedup.bam
```

#### Step 7：On-target 过滤

```bash
bedtools intersect -abam sample.bwameth_pe.dedup.bam \
  -b data/panel_B10308.bed -wa \
| samtools sort -o sample.bwameth_pe.dedup.targeted.sorted.bam
samtools index sample.bwameth_pe.dedup.targeted.sorted.bam
```

#### Step 7.5：M-bias 检查

```bash
MethylDackel mbias /data/meth/bwa/hg38.fa \
  sample.bwameth_pe.dedup.targeted.sorted.bam sample.mbias
# 生成 sample.mbias_OT.svg / sample.mbias_OB.svg
```

**判断标准**：

- R2 3' 端甲基化率断崖式下降 → 增大 `trim_tail2`（如 8 或 10），回 Step 4 重跑至 Step 7
- 曲线平稳 → 当前 `trim_tail2=5` 合适，继续后续步骤
- OB 占比 < 5% → `non_directional = false`（已确认）

**本实验室确认结果**：R2 3' 末端甲基化率约 30% 且平稳，`trim_tail2=5` 无需调整。

#### Step 8：甲基化信息提取

```bash
MethylDackel extract \
  --minDepth 1 --mergeContext \
  /data/meth/bwa/hg38.fa \
  sample.bwameth_pe.dedup.targeted.sorted.bam \
  -o sample
# 输出：sample_CpG.bedGraph
# 列：chrom  start  end  methylation%  methCount  unmethCount
```

> `--mergeContext`：对称 CpG 双链计数合并，cfDNA 短片段必须启用。  
> `--minDepth 1`：深度过滤推迟到 Step 9 统一处理。

#### Step 9：过滤与候选区域聚焦

- 重复序列过滤（`hg38.repeat.bed`）
- 深度过滤（`min_cpg_per_region=3`，`min_reads_per_region=5`）
- 聚焦 88 个候选区域（`OC_tissue_88_regions.hg38.txt`）

#### Step 10：差异分析

Wilcoxon 秩和检验，阈值：P < 0.05 **且** |Δβ| > 0.1。  
输出：`result/diff_methylation.txt`

#### Step 11：可视化

密度分布图 / 位置散点图 / 甲基化热图 → `result/density_plot.png` 等。

#### Step 14：特征选择（LASSO + RF）

从 DMC 集合中精简核心 CpG 面板：

| 方法 | 实现 | 说明 |
| --- | --- | --- |
| **LASSO** | `LogisticRegressionCV(penalty="l1", cv=5)` | 非零系数特征入选 |
| **随机森林** | `RandomForestClassifier(n_estimators=500)` | Gini 重要性 > 0 入选 |
| **并集（默认）** | `feature_selection_mode = "union"` | 任一方法入选即保留 |

覆盖率过滤：< `min_feature_coverage`（默认 50%）样本无覆盖的位点丢弃。

**输出**（`result/`）：

| 文件 | 内容 |
| ---- | ---- |
| `panel_cpg.csv` | 入选 CpG 位点列表（`feature` 列，格式 `chr:pos`） |
| `panel_cpg_depth.csv` | 每个 panel 位点的深度统计（mean/median/p10/p90/min/max/coverage_rate） |
| `feature_importance.csv` | 每个 DMC 的 LASSO 系数 + RF 重要性 + 入选标记 |
| `feature_matrix.npz` | 样本×特征矩阵缓存（含深度矩阵 D，供 Step 15 复用） |

#### Step 15：诊断分类器（LOO-CV）

以 panel CpG 甲基化水平为输入，LOO-CV 评估五种分类器：

| 模型 | 处理类别不均衡 | 关键参数 |
| --- | --- | --- |
| **LR（L2 逻辑回归）** | `class_weight="balanced"` | `lr_C`（默认 1.0） |
| **SVM（RBF 核）** | `class_weight="balanced"` | `svm_C`（默认 1.0） |
| **RF（随机森林）** | `class_weight="balanced"` | `rf_n_estimators`（默认 200） |
| **XGBoost** | `scale_pos_weight = n_neg/n_pos` | `xgb_n_estimators`, `xgb_max_depth` |
| **LightGBM** | `scale_pos_weight = n_neg/n_pos` | `lgbm_n_estimators`, `lgbm_max_depth` |

LOO-CV 要点：每轮在训练集独立 fit StandardScaler，严格防止数据泄漏。

**输出**（`result/`）：`roc_curve.png`、`confusion_matrix.png`、`classifier_summary.csv`、`classifier_predictions.csv`

#### Step 16：亚型分类器（LOO-CV，多分类）

对恶性 + 交界性样本按亚型做多分类：

- **亚型**：高级别浆液性癌 / 子宫内膜样癌 / 透明细胞癌 / 低级别浆液性癌 / 黏液性
- **分类器**：RF（原生多分类）+ SVM OvR（`decision_function_shape="ovr"`，Platt scaling 输出概率）
- **评估**：LOO-CV，输出 OvR ROC 曲线

**输出**（`result/`）：`subtype_summary.csv`、`subtype_predictions.csv`、`subtype_confusion_RF.png`、`subtype_confusion_SVM.png`、`subtype_roc.png`

#### Step 17：生成分析报告

全流程 HTML + JSON 报告，包含：

- 每样本 fastp QC / 比对 / 去重 / On-target 统计
- M-bias SVG 内嵌
- 甲基化提取与候选区域统计
- 差异甲基化分析结果（前 20 DMC 表 + 染色体分布）
- 可视化图表（密度 / 散点 / 热图）
- **（可选）** 特征选择结果 + panel 位点深度统计（含低深度警告）
- **（可选）** 5 分类器 LOO-CV 汇总 + ROC + 混淆矩阵
- **（可选）** 亚型分类 LOO-CV 汇总 + 混淆矩阵 + OvR ROC
- 结果解读与质量评价

可选节在对应步骤的输出文件存在时自动出现，否则静默跳过。

---

## 四、配置参数参考（config.toml）

```toml
[plasma]
umi_length      = 6
umi_location    = "read2"
trim_front1     = 12          # R1 5' 低复杂度寡核苷酸
trim_front2     = 12          # R2 5' UMI 后寡核苷酸
trim_tail2      = 5           # R2 3' poly-C（根据 M-bias 调整！）
min_length      = 35
min_quality     = 20
non_directional = false       # 已确认高度定向（OB 占 1.2%）
min_mapq        = 40
nOT             = "0,0,0,20"  # MethylDackel OT 链偏差校正（根据 M-bias 调整）
nOB             = "0,0,17,0"  # MethylDackel OB 链偏差校正
min_cpg_per_region   = 3
min_reads_per_region = 5
methylated_read_threshold  = 0.60   # 与组织统一（原优讯 80% → 降至 60%）
low_confidence_threshold   = 0.40
diff_p_value    = 0.05
diff_delta_beta = 0.10

[extraction]
reference_genome = "/data/meth/bwa/hg38.fa"

[tools]
conda_env_bin       = "/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin"
umi_tools_fallback  = "/mnt/e/programs/ubuntu/miniconda3/bin/umi_tools"

[ml]
feature_selection_mode = "union"   # "union" | "lasso" | "rf"
min_feature_coverage   = 0.5
lasso_cv_folds         = 5
rf_n_estimators        = 500
rf_min_importance      = 0.0
lr_C                   = 1.0
svm_C                  = 1.0
xgb_n_estimators       = 200
xgb_max_depth          = 4
xgb_learning_rate      = 0.1
lgbm_n_estimators      = 200
lgbm_max_depth         = 4
lgbm_learning_rate     = 0.1
```

---

## 五、执行操作指南

### 环境准备

```bash
conda activate meth
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
```

### 常用命令

```bash
# 查看步骤列表及当前进度
python src/scripts/run_analysis.py --list-steps --workdir $WORKDIR

# 全流程一键运行（自动创建 workdir）
python src/scripts/run_analysis.py --plasma /mnt/f/meth/OC_plasma/
```

### 分步执行

```bash
# Step 4：预处理（首次执行，自动创建 workdir）
python src/scripts/run_analysis.py --plasma --step preprocess /mnt/f/meth/OC_plasma/

# 从输出读取 workdir，后续步骤复用
WORKDIR=temp/$(ls -t temp/ | head -1)

python src/scripts/run_analysis.py --plasma --step align          --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step dedup          --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step target         --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step mbias          --workdir $WORKDIR

# ── 检查 M-bias SVG（$WORKDIR/*.mbias_OT.svg）──
# 若 R2 3' 端曲线下滑 → 调整 config.toml trim_tail2，重跑 preprocess → target
# 若曲线平稳 → 继续

python src/scripts/run_analysis.py --plasma --step extract        --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step filter         --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step diff           --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step visualize      --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step feature_select --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step classify       --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step subtype_classify --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step report         --workdir $WORKDIR
```

### 一键 Shell 脚本（后台运行）

```bash
#!/bin/bash
# 用法：bash run_oc_plasma.sh 2>&1 | tee logs/oc_plasma_run.log
set -e
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
PY=/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin/python

$PY src/scripts/run_analysis.py --plasma --step preprocess /mnt/f/meth/OC_plasma/
WORKDIR=$(ls -dt temp/20[0-9]* | head -1)
echo "workdir: $WORKDIR"

$PY src/scripts/run_analysis.py --plasma --step align          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step dedup          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step target         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step mbias          --workdir $WORKDIR

echo "=== 请检查 M-bias SVG（$WORKDIR/*.mbias_OT.svg）后按 Enter 继续，Ctrl+C 中止调参 ==="
read

$PY src/scripts/run_analysis.py --plasma --step extract        --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step filter         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step diff           --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step visualize      --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step feature_select --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step classify       --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step subtype_classify --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step report         --workdir $WORKDIR

echo "=== 完成！报告：result/$WORKDIR/analysis_report.html ==="
```

### 断点续跑

```bash
# 查看已完成步骤
python src/scripts/run_analysis.py --list-steps --workdir $WORKDIR

# 从中断步骤重跑（以 align 为例）
python src/scripts/run_analysis.py --plasma --step align --workdir $WORKDIR
```

每步完成后状态自动记录到 `$WORKDIR/pipeline_state.json`，重跑单步无需重新指定输入文件。

### M-bias 后重跑范围

| 情形 | 调整内容 | 需重跑步骤 |
| ---- | ------- | --------- |
| trim_tail2 不够大 | 增大 `trim_tail2` | preprocess → dedup → target → extract |
| nOT/nOB 需要调整 | 修改 `nOT` / `nOB` | extract 起（无需重新比对） |

---

## 六、与组织分析的关系

```text
组织分析（已完成）
  └─ 88 个候选 DMR（OC_tissue_88_regions.hg38.txt）
       │  聚焦（降低多重检验负担）
       ▼
血浆分析
  ├─ 在 88 个候选区域中独立验证血浆样本差异区域（P<0.05, |Δβ|>0.1）
  ├─ 差异区域经 LASSO + RF 精简 → 核心 CpG 面板
  ├─ 五种分类器 LOO-CV 评估（LR / SVM / RF / XGB / LGBM）
  ├─ 亚型多分类（RF + SVM OvR）
  └─ 最终目标：无创液体活检诊断卵巢癌
```

> **注意**：不能直接把组织 DMR 照搬到血浆。血浆 cfDNA 成分复杂、片段化、信噪比低，须在血浆数据中独立统计验证后再入模型。

---

## 七、资源占用参考

### 关键变量（测试样本实测）

| 变量 | 含义 | B10002（实测） | B10197（实测） |
| ---- | ---- | ------------- | ------------- |
| N | 输入 read 对数（Mpairs） | ~16 | ~26 |
| R_q | fastp 保留率 | ~0.90 | ~0.90 |
| R_aln | bwameth 比对率 | ~0.88（估） | ~0.88（估） |
| R_dup | 去重后保留率 | **0.375** ✓实测 | **0.500** ✓实测 |
| R_ont | On-target 率 | **0.362** ✓实测 | **0.247** ✓实测 |

### 逐步资源表（N=16 Mpairs 参考，8 核 / 32 GB）

| 步骤 | 工具 | 时间 | 内存峰值 | 输出存储 |
| ---- | ---- | ---- | ------- | ------- |
| 4 preprocess | fastp | ~2 min | 0.5 GB | × R_q ≈ ×0.90 |
| 5 align | bwameth + samtools | **~47 min** ✓实测 | **13–14 GB** | ≈ clean FASTQ × 0.85 |
| 6 dedup | umi_tools | **~10 min** ✓实测 | 2–4 GB | × R_dup ≈ ×0.38–0.50 |
| 7 target | bedtools | ~4 min | 1–2 GB | × R_ont ≈ ×0.25–0.36 |
| 7.5 mbias | MethylDackel | ~5 min | 2–4 GB | ~2 MB（SVG） |
| 8 extract | MethylDackel | ~10 min | 2–4 GB | ~20–50 MB |
| 9 filter | Python | <2 min | 1 GB | <1 MB |
| 10 diff | Python | <1 min | <1 GB | <1 MB |
| 11 visualize | matplotlib | <2 min | 1 GB | ~3 MB |
| 14 feature_select | scikit-learn | 5–30 min（样本量依赖） | 2–4 GB | <10 MB |
| 15 classify | sklearn/XGB/LGBM | 10–60 min（LOO-CV） | 2–8 GB | <5 MB |
| 16 subtype_classify | sklearn | 5–30 min（LOO-CV） | 2–4 GB | <5 MB |
| 17 report | Python | <2 min | 1 GB | ~5–15 MB（含嵌图） |

**内存峰值**：bwameth 比对阶段（BWT 6 GB + SA 3 GB + 工作缓冲 ≈ 13–14 GB）。

**时间估算公式**（简化，R_q=0.90, R_aln=0.88）：

```text
T_total(分钟) ≈ N × 4.1 + 5    （Steps 4–11）
```

### 存储估算（单样本）

```text
峰值（全部中间文件）≈ N × 360 MB + 50 MB
最终保留（targeted BAM + 结果）≈ N × 140 MB + 50 MB
```

N=16：峰值约 **5.8 GB**，最终保留约 **2.3 GB**。

### 数据流量漏斗（N=16，B10002 实测）

```text
原始 FASTQ          1.88 GB (R1+R2)
  ↓ fastp ×0.90
clean FASTQ         1.80 GB ✓实测
  ↓ bwameth ×0.79
sorted BAM          1.50 GB ✓实测
  ↓ umi_tools ×0.375
dedup BAM           0.56 GB ✓实测
  ↓ bedtools ×0.362
on-target BAM       0.20 GB ✓实测  ← 主要分析数据（原始 FASTQ 的 11%）
  ↓ MethylDackel
CpG.bedGraph        ~1 MB
  ↓ 候选区域过滤
candidate.bedGraph  <1 MB           ← 最终分析数据（<0.1%）
```

### bwameth c2t 索引文件（全部必需，共 ~17 GB）

| 文件 | 大小 | 用途 |
| ---- | ---- | ---- |
| `hg38.fa.bwameth.c2t` | 6.1 GB | C→T 转换参考基因组 |
| `hg38.fa.bwameth.c2t.bwt` | 6.0 GB | BWT 索引（比对核心） |
| `hg38.fa.bwameth.c2t.sa` | 3.0 GB | 后缀数组（种子精确定位） |
| `hg38.fa.bwameth.c2t.pac` | 1.5 GB | 2-bit 压缩序列（CIGAR 计算） |
| `hg38.fa.bwameth.c2t.amb` | 42 KB | 模糊碱基表 |
| `hg38.fa.bwameth.c2t.ann` | 44 KB | 染色体注释（SAM 头部） |

---

## 八、OC_plasma 全量运行计划（10 样本）

### 样本概况

| 样本 | R1 大小 | 估算 N（Mpairs） | 备注 |
| ---- | ------- | --------------- | ---- |
| B10002 | 1.2 GB | ~16 | 恶性（已测试） |
| B10011 | 1.6 GB | ~21 | — |
| B10014 | 1.4 GB | ~18 | — |
| B10018 | 1.7 GB | ~22 | — |
| B10019 | 1.4 GB | ~18 | — |
| B10188 | 1.3 GB | ~17 | — |
| B10189 | 2.8 GB | ~37 | ⚠️ 最大样本，sort 内存需关注 |
| B10194 | 1.6 GB | ~21 | — |
| B10195 | 1.3 GB | ~17 | — |
| B10197 | 1.9 GB | ~26 | 良性（已测试） |
| **合计** | **18.3 GB** | **~213** | — |

### 资源预估（顺序执行）

| 步骤 | 合计时间 | 内存峰值 |
| ---- | ------- | ------- |
| preprocess | ~30 min | 0.5 GB |
| align | ~10–11 h | 13–14 GB |
| dedup | ~2 h | 2–4 GB |
| target | ~40 min | 1–2 GB |
| mbias + extract | ~40 min | 2–4 GB |
| filter + diff + visualize | ~15 min | <1 GB |
| feature_select + classify + subtype + report | ~1–2 h | 2–8 GB |
| **合计** | **~15–17 h** | **14 GB** |

**磁盘**（E 盘，需可用 ≥ 120 GB）：

| 阶段 | 空间 |
| ---- | ---- |
| 输入 FASTQ（F 盘已有） | ~37 GB |
| 峰值（全部中间文件） | ~115 GB |
| 最终保留（targeted BAM + 结果） | ~40 GB |

### 完成后清理

```bash
# 删除大型中间文件，保留 targeted BAM + 结果（~6–8 GB）
rm $WORKDIR/*.bwameth_pe.sorted.bam
rm $WORKDIR/*.bwameth_pe.sorted.bam.bai
rm $WORKDIR/*.bwameth_pe.dedup.bam
rm $WORKDIR/*.bwameth_pe.dedup.bam.bai
# 可选（原始 FASTQ 保留在 /mnt/f/meth/OC_plasma/）
rm $WORKDIR/*.clean.R*.fq.gz
```

---

## 九、已确认问题清单

| 问题 | 状态 | 当前设置 |
| ---- | ---- | ------- |
| 建库方向性 | ✅ 已确认 | `non_directional = false`（OB 占 1.2%） |
| R2 3' poly-C 影响 | ✅ 已确认 | `trim_tail2 = 5`（R2 末端甲基化率 30% 平稳，无需调整） |
| bwameth c2t 索引 | ✅ 已完成 | `/data/meth/bwa/hg38.fa.bwameth.c2t`（18 GB） |
| Bismark 索引（备用） | ✅ 已完成 | `/data/meth/bwa/Bisulfite_Genome/` |
| umi_tools 安装 | ✅ 已完成 | meth conda env，v1.1.6 |
| Panel BED 坐标版本 | ✅ hg38 | `data/panel_B10308.bed` |
| 重复序列 BED | ✅ 已下载 | `/mnt/f/meth/hg38.repeat.bed` |
| 启动子 BED | ✅ 已下载 | `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` |
| 测试样本预处理 | ✅ 已完成 | B10002（恶性）、B10197（良性） |
| 测试样本比对 | ✅ 已完成 | B10002 1.5 GB / 47 min，B10197 2.8 GB / 124 min |
| align 管道死锁 | ✅ 已修复 | 两步法（view + sort）+ `samtools sort -m 8G` |

# 血浆cfDNA靶向甲基化分析：问题诊断与优化方案

> 综合参考：优讯血浆流程（实验结果欠佳）、Crist标准SOP、组织分析经验，结合实验室实际建库情况整理。
> **本文档与代码实现同步，所有命令均已在 conda `meth` 环境下测试。**
> 最后更新：2026-04-03

---

## 一、实验室建库关键参数确认（前提）

在分析流程设计之前，必须明确以下建库信息，它们直接决定多个关键步骤的参数。

| 参数 | 本实验室情况 | 影响的步骤 |
| ---- | ----------- | --------- |
| UMI | **有，6N，位于R2 5'端** | Step 1 fastp提取、Step 5 去重方式 |
| 建库方式 | **单链建库（single-strand）** | Step 3 比对参数（directional 待确认） |
| 末端修复 | **无末端修复；P7连接引入poly-C（未甲基化）** | Step 1 trim_tail2、Step 7 M-bias |
| 5'端特殊序列 | **R1/R2 5'端各有12个低复杂度寡核苷酸** | Step 1 trim_front1/trim_front2 |
| 参考基因组 | hg38（`/data/meth/bwa/hg38.fa`） | 全流程 |
| Bismark索引 | `/data/meth/bwa/Bisulfite_Genome/`（已构建） | Step 3 |
| Panel BED | `data/panel_B10308.bed` | Step 6 On-target过滤 |

---

## 二、优讯血浆流程问题诊断

优讯提交的血浆靶向甲基化分析结果"效果不太好"，通过与其tissue流程、Crist SOP及建库实际情况对比，找到以下**根本原因**：

### 问题1（最严重）：坐标去重误删大量真实cfDNA片段

**描述**：优讯使用 Picard MarkDuplicates 进行PCR重复标记，这是基于reads的起止坐标进行去重的。

**致命缺陷**：cfDNA片段高度集中在147~167 bp，且断裂点集中在核小体边缘，导致**来自不同细胞的真实独立片段也有相同的起止坐标**。Picard无法区分这些真实片段和PCR扩增重复，会把它们全部标记为重复并删除。

**后果**：在靶向测序（杂交捕获需多轮PCR）的场景下，这个问题极其严重——有效深度会断崖式下降，而血浆中ctDNA占比本来就极低（早期可能仅1~5%），深度下降直接导致信号完全消失。

**本实验室有UMI（6N in R2）**，这正是为了解决这个问题，但优讯完全没有利用UMI，仍然做坐标去重，等于白白浪费了UMI建库的成本，同时引入了严重偏差。

### 问题2：缺少靶向区域过滤（On-target Filtering）

**描述**：优讯流程中没有用panel.bed对比对结果做intersect过滤的步骤。

**后果**：脱靶比对的reads（全基因组背景）未被去除，这些reads的甲基化信号是随机噪声，会稀释真实的靶区信号，严重干扰差异分析。Crist SOP将此步骤标记为"核心步骤"。

### 问题3：建库方向性未处理，可能丢失50%数据

**描述**：本实验室采用单链建库，两条链都被纳入文库（基因组两条链均被利用）。优讯使用bwameth，默认按directional（单链定向）模式比对。

**后果**：如果建库实际上是non-directional，双向reads中有一半（来自互补链的reads）会被软件错误丢弃，数据利用率仅50%。需要通过实际数据中OT/OB/CTOT/CTOB四种链的比例来判断，或咨询建库试剂盒厂商，并在比对参数中明确指定。

### 问题4：P7接头poly-C导致甲基化率低估，未处理

**描述**：本实验室建库无末端修复，但P7接头连接时引入了poly-C（人工未甲基化碱基）。poly-C经亚硫酸氢盐转化后变为poly-T，直接影响reads 3'端（R2端）的甲基化定量。

**后果**：R2末端的CpG位点甲基化率会被系统性低估。优讯仅切除了5'端12个低复杂度寡核苷酸，未针对此问题做处理。

**处理方法**：比对后运行M-bias图检查（`MethylDackel mbias`），若发现reads末端甲基化率断崖式下降，则回头调整 `trim_tail2` 参数重跑预处理。

### 问题5：血浆甲基化read阈值比组织更严苛，逻辑相反

对比优讯的两份流程：

- Tissue：甲基化read判定阈值为CpG甲基化比例 **>60%**
- 血浆：甲基化read判定阈值为 **>80%**

血浆ctDNA本来信号就弱，反而设置更严苛的阈值（80% vs 60%），在ctDNA fraction极低时可能完全检测不到任何甲基化reads。阈值应当与组织保持一致或更宽松，而非更严格。

### 问题6：差异分析缺少Δβ绝对差值过滤

Tissue流程要求 Wilcoxon P<0.05 **且** |Δβ|>0.1，而血浆流程只有 P<0.05，没有甲基化差值的绝对量过滤，导致大量统计显著但生物学差异微小的位点被纳入，增加假阳性。

---

## 三、两套流程方法学对比

| 环节 | 优讯流程（实际使用） | 本实验室流程 | 评价 |
| ---- | ----------------- | --------- | --- |
| QC工具 | fastp（宽松：Q<15，40%以上才过滤） | **fastp（Q20，最小长度35bp，--umi合并UMI提取）** | 标准更严格 |
| 5'端trimming | 切除12个低复杂度寡核苷酸 | **trim_front1=12, trim_front2=12** | 一致 |
| 3'端trimming | 无 | **trim_tail2=5（根据M-bias调整）** | ✅ 改进 |
| UMI提取 | **未使用** ❌ | **fastp --umi --umi_loc read2 --umi_len 6** ✅ | 关键改进 |
| 比对工具 | bwameth | **bwameth（BWA-MEM + C→T转换，8线程）** | 速度相当，更易维护 |
| 去重 | **Picard坐标去重** ❌ | **umi_tools dedup --paired** ✅ | 最关键区别 |
| On-target filtering | **缺失** ❌ | **bedtools intersect panel.bed** ✅ | 关键步骤补充 |
| M-bias检查 | 无 | **MethylDackel mbias** ✅ | 必要质控 |
| 甲基化calling | Biscuit（read-level） | **MethylDackel extract（CpG.bedGraph）** | 与bwameth配套 |
| 双端重叠处理 | 未提及 | **--mergeContext（MethylDackel）** ✅ | cfDNA短片段必须处理 |
| 甲基化阈值 | 80% | **60%（与tissue一致）** ✅ | 与tissue统一 |
| 差异分析 | 仅P<0.05 | **P<0.05 且 \|Δβ\|>0.1** ✅ | 减少假阳性 |

---

## 四、最终优化流程设计（实际实现版）

### 流程概览

```text
Raw Reads (FASTQ paired-end, B10002/B10197)
  │
  ▼
Step 4: fastp QC + UMI提取（一步完成）
         ├─ 接头去除 + 质量过滤（Q20，最小长度35bp）
         ├─ --umi --umi_loc read2 --umi_len 6  → UMI附加到read名称
         ├─ --trim_front1 12 --trim_front2 12  → 5'端低复杂度序列
         └─ --trim_tail2 5                     → R2 3'端poly-C（根据M-bias调整）
  │
  ▼
Step 5: bwameth双端比对（BWA-MEM + C→T转换，8线程）
         ├─ reference: /data/meth/bwa/hg38.fa
         ├─ index:     /data/meth/bwa/hg38.fa.bwameth.c2t（已构建）
         └─ 输出排序索引BAM（bwameth.py | bwa mem | samtools sort 管道）
  │
  ▼
Step 6: UMI感知去重（umi_tools dedup --paired）
         └─ 区分真实cfDNA独立片段 vs PCR扩增重复
  │
  ▼
Step 7: On-target过滤（bedtools intersect data/panel_B10308.bed）
         └─ 只保留panel覆盖区域的reads
  │
  ▼
Step 7.5: M-bias检查（MethylDackel mbias）
           ├─ 生成SVG检查reads各位置甲基化率
           └─ 若R2 3'端异常 → 调整trim_tail2后重跑Step 4
  │
  ▼
Step 8: 甲基化信息提取（MethylDackel extract）
         └─ --paired-end --no_overlap --bedGraph --counts
  │
  ▼
Step 9: 甲基化数据过滤 + 候选区域聚焦
         ├─ 覆盖深度过滤（min_cpg_per_region=3, min_reads=5）
         └─ 聚焦88个候选DMR（OC_tissue_88_regions.hg38.txt）
  │
  ▼
Step 10: 差异分析（Wilcoxon P<0.05 且 |Δβ|>0.1）
  │
  ▼
Step 11: 可视化图表生成（密度图/散点图/热图 → result/density_plot.png 等）
  │
  ▼
Step 12: 生成分析报告（HTML 含 M-bias SVG + 可视化图表）
```

### 详细步骤说明

#### Step 4: fastp QC + UMI提取（合并操作）

fastp 通过 `--umi` 参数在质控的同时完成UMI提取，无需单独运行 umi_tools extract：

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

> ⚠️ `trim_tail2=5` 是初始值，必须在 Step 7.5 M-bias 结果出来后确认或调整。

#### Step 5: bwameth双端比对

bwameth.py 将 C→T 碱基转换后调用 bwa mem 完成比对。为避免大样本（>20M pairs）管道死锁，先通过 `samtools view` 流式转换为 unsorted BAM，再单独排序：

```bash
# 步骤 1：比对（bwameth stdout → samtools view 流式转为 BAM，无排序缓冲）
bwameth.py \
  --reference /data/meth/bwa/hg38.fa \
  -t 8 \
  sample.clean.R1.fq.gz \
  sample.clean.R2.fq.gz \
| samtools view -bS -o sample.unsorted.bam -

# 步骤 2：排序（从文件读取，bwameth 已退出，不会死锁）
samtools sort -@ 4 -m 8G -o sample.bwameth_pe.sorted.bam sample.unsorted.bam
rm sample.unsorted.bam        # 排序完成后自动删除中间文件

samtools index sample.bwameth_pe.sorted.bam
```

> ⚠️ bwameth c2t 索引需提前构建（一次性，约 2.5 小时）：
>
> `bwameth.py index /data/meth/bwa/hg38.fa`
>
> ⚠️ `samtools sort -m 8G`：每线程 8 GB 内存限制（4线程共 32 GB），防止大样本内存溢出导致 temp 文件过多。48 GB WSL 环境安全。

#### Step 6: UMI感知去重

```bash
umi_tools dedup \
  --paired \
  -I sample.sorted.bam \
  -S sample.dedup.bam \
  --log sample.dedup.log

samtools index sample.dedup.bam
```

#### Step 7: On-target过滤

```bash
bedtools intersect \
  -abam sample.dedup.bam \
  -b data/panel_B10308.bed \
  -wa \
  > sample.ontarget.bam

samtools sort -o sample.ontarget.sorted.bam sample.ontarget.bam
samtools index sample.ontarget.sorted.bam
```

#### Step 7.5: M-bias检查

```bash
MethylDackel mbias \
  /data/meth/bwa/hg38.fa \
  sample.ontarget.sorted.bam \
  sample.mbias

# 生成 sample.mbias_R1.svg 和 sample.mbias_R2.svg
# 判断标准：
#   R2 3'端甲基化率断崖式下降 → 增大 trim_tail2（如改为8或10），回Step 4重跑
#   曲线平稳 → trim_tail2=5 合适，继续后续步骤
```

#### Step 8: 甲基化信息提取

```bash
MethylDackel extract \
  --minDepth 1 \
  --mergeContext \
  /data/meth/bwa/hg38.fa \
  sample.ontarget.sorted.bam \
  -o sample_prefix

# 输出：sample_prefix_CpG.bedGraph
# 列：chrom  start  end  methylation%  methCount  unmethCount
```

> `--mergeContext`：对称 CpG 双链计数合并，cfDNA 短片段必须启用；
> `--minDepth 1`：深度过滤推迟到 Step 9 统一处理。

#### Step 9: 候选区域聚焦

```bash
bedtools intersect \
  -a sample.CpG.bedGraph \
  -b data/OC_tissue_88_regions.hg38.txt \
  -wa \
  > sample.candidate.bedGraph
```

---

## 五、分步执行操作指南

### 环境准备

```bash
# conda 环境（必须激活）
conda activate meth

# 工作目录
cd /mnt/e/projects/hq/meth

# PYTHONPATH
export PYTHONPATH=src
```

### 查看步骤列表

```bash
python src/scripts/run_analysis.py --list-steps
# 附加 --workdir 可查看当前进度
python src/scripts/run_analysis.py --list-steps --workdir temp/20260327_134610
```

### 分步执行命令

```bash
WORKDIR=temp/20260327_134610

# ── Step 4: 预处理（首次执行，自动创建工作目录）──
python src/scripts/run_analysis.py --plasma --step preprocess --workdir $WORKDIR /mnt/f/meth/OC_plasma/

# ── 后续步骤：通过 --workdir 指定上步生成的目录 ──

python src/scripts/run_analysis.py --plasma --step align   --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step dedup   --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step target  --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step mbias   --workdir $WORKDIR

# ── 检查 M-bias SVG，确认 trim_tail2 参数 ──
# 若需调整：编辑 config.toml 中 trim_tail2，然后重跑 preprocess+align+dedup+target
# 若无需调整：继续

python src/scripts/run_analysis.py --plasma --step extract --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step filter  --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step diff    --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step visualize --workdir $WORKDIR
python src/scripts/run_analysis.py --plasma --step report  --workdir $WORKDIR
```

### 全流程一次性执行

```bash
python src/scripts/run_analysis.py --plasma /mnt/f/meth/OC_plasma/
```

### 状态文件

每一步完成后，工作目录下自动更新 `pipeline_state.json`，记录各步骤输出文件路径和已完成步骤列表，实现断点续跑。

---

## 六、核心参数参考（config.toml）

```toml
[plasma]
# UMI设置（R2 5'端6N）
umi_length = 6
umi_location = "read2"

# fastp Trimming参数
trim_front1 = 12        # R1 5'端切除（低复杂度寡核苷酸）
trim_front2 = 12        # R2 5'端切除（UMI提取后的寡核苷酸）
trim_tail2 = 5          # R2 3'端切除（poly-C偏差，根据M-bias图调整！）
min_length = 35
min_quality = 20

# 比对参数
non_directional = false # 待M-bias/OT-OB比例确认后调整

# 比对质量过滤
min_mapq = 40

# 区域覆盖深度过滤
min_cpg_per_region = 3
min_reads_per_region = 5

# 甲基化read阈值（策略B，与tissue统一）
methylated_read_threshold = 0.60   # 原优讯80% → 降至60%
low_confidence_threshold = 0.40

# 差异分析阈值
diff_p_value = 0.05
diff_delta_beta = 0.10

[extraction]
reference_genome = "/data/meth/bwa/hg38.fa"

[tools]
conda_env_bin = "/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin"
```

---

## 七、关键问题清单

| 问题 | 状态 | 当前设置 | 需要操作 |
| ---- | ---- | ------- | ------- |
| 建库方向性 | ✅ 已确认 | `non_directional = false` | OB仅占1.2%，高度定向，无需修改 |
| R2 3'端poly-C影响范围 | ✅ 已确认 | `trim_tail2 = 5`（无需修改） | R2 3'末端甲基化率30%平稳，无偏差 |
| bwameth c2t 索引 | ✅ 已完成 | `/data/meth/bwa/hg38.fa.bwameth.c2t`（18 GB） | — |
| Bismark索引（备用） | ✅ 已完成 | `/data/meth/bwa/Bisulfite_Genome/` | — |
| umi_tools安装 | ✅ 已完成 | meth conda env, v1.1.6 | — |
| Panel BED坐标版本 | ✅ hg38 | `data/panel_B10308.bed` | — |
| filter 重复序列 BED | ✅ 已下载 | `/mnt/f/meth/hg38.repeat.bed` | filter 步骤自动使用 |
| filter 启动子 BED | ✅ 已下载 | `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` | filter 步骤自动使用 |
| 测试样本预处理 | ✅ 已完成 | B10002（恶性）、B10197（良性） | — |
| 测试样本比对 | ✅ 已完成 | B10002 1.5 GB / 47 min，B10197 2.8 GB / 124 min | — |
| align 管道死锁 | ✅ 已修复 | 两步法（view + sort）+ `samtools sort -m 8G` | B10189 等大样本不再卡死 |

---

## 八、与组织分析的关系

血浆分析与组织分析形成"引导-验证"的互补关系：

```text
组织分析（已完成）
  └─ 88个候选DMR区域（OC_tissue_88_regions.hg38.txt）
       │
       ▼  聚焦（大幅降低多重检验负担）
血浆分析
  ├─ 在88个候选区域中寻找血浆样本中有差异的区域
  ├─ 差异区域用于训练分类模型（恶性 vs 良性）
  └─ 最终目标：无创液体活检诊断卵巢癌
```

**设计原则**：

- 不能直接把组织分析发现的88个区域照搬到血浆（血浆cfDNA成分复杂、片段化、信噪比低）
- 需要在血浆数据中独立验证哪些区域在良/恶性间有统计显著差异（P<0.05, |Δβ|>0.1）
- 机器学习模型以通过双重过滤的区域为特征（弹性网络/随机森林/SVM）

---

## 九、各步骤资源占用参考

### 9.1 变量定义

| 变量 | 含义 | 典型范围 | 测试样本 B10002 | 测试样本 B10197 |
| ---- | ---- | ------- | -------------- | -------------- |
| N | 输入 read 对数（百万对，Mpairs） | 10–200 | ~16 | ~26 |
| R_q | fastp 质控保留率 | 0.85–0.95 | 0.90（预估） | 0.90（预估） |
| R_aln | bwameth 比对率 | 0.80–0.95 | 待测（预估 0.88） | 待测（预估 0.88） |
| R_dup | UMI去重后保留率（真实独立片段占比） | 0.35–0.55 | **0.375** ✓实测 | **0.500** ✓实测 |
| R_ont | On-target 率（靶区reads占比） | 0.25–0.45 | **0.362** ✓实测 | **0.247** ✓实测 |

> N 从 clean FASTQ 大小估算：R1 文件大小(MB) / 58 ≈ Mpairs（B10002 R1=935 MB → ~16 Mpairs；B10197 R1=1500 MB → ~26 Mpairs）

派生量：

```text
N_clean = N × R_q                              # 预处理后read对数
N_aln   = N × R_q × R_aln                     # 比对成功read对数
N_dedup = N × R_q × R_aln × R_dup             # 去重后read对数
N_ont   = N × R_q × R_aln × R_dup × R_ont     # on-target read对数
```

B10002 (N≈16)，典型参数下：N_clean≈14.4，N_aln≈12.7，N_dedup≈7.0，N_ont≈3.8 Mpairs

---

### 9.2 逐步资源占用表

> 参考硬件：8核 CPU / 32 GB RAM。时间列以 N=16 Mpairs 单样本（B10002实测）为例。

| 步骤 | 工具 | 输入格式 | 输出格式 | CPU推荐 | 时间(N=16M) | 内存峰值 | 输入存储 | 输出存储 | 存储变化率 |
| ---- | ---- | ------- | ------- | ---------- | -------------- | ---------- | ------- | ------- | ------------ |
| **4** preprocess | fastp | FASTQ.gz × 2 | FASTQ.gz × 2 | 8T | ~2 min | 0.5 GB | N × 58 MB | N_clean × 58 MB | × R_q ≈ **×0.90** |
| **5** align | bwameth + bwa mem + samtools sort | FASTQ.gz × 2 | sorted BAM | 8T | **~47 min** ✓实测 | **13–14 GB** | N_clean × 58 MB | ≈ clean FASTQ总大小 | **×0.81–0.93** |
| **6** dedup | umi_tools dedup | sorted BAM | dedup BAM | 1T | **B10002: 10 min / B10197: 17 min** ✓实测 | 2–4 GB | ≈ BAM大小 | BAM × R_dup | × R_dup ≈ **×0.38–0.50**（实测） |
| **7** target | bedtools intersect | dedup BAM | on-target BAM | 1T | **~4 min（两样本）** ✓实测 | 1–2 GB | N_dedup × 100 MB | N_ont × 100 MB | × R_ont ≈ **×0.25–0.36**（实测） |
| **7.5** mbias | MethylDackel mbias | on-target BAM | SVG × 2 | 4T | ~5 min | 2–4 GB | N_ont × 100 MB | ~2 MB（固定） | **→ 微小** |
| **8** extract | MethylDackel extract | on-target BAM | CpG.bedGraph | 4T | ~10 min | 2–4 GB | N_ont × 100 MB | ~20–50 MB（靶向，固定） | **→ 微小** |
| **9** filter | Python | CpG.bedGraph | 候选区域 bedGraph | 1T | <2 min | 1 GB | ~50 MB | <1 MB（88区域） | **→ 极小** |
| **10** diff | Python | 甲基化矩阵 | 差异结果表 | 1T | <1 min | <1 GB | <1 MB | <1 MB | **→ 极小** |
| **11** visualize | Python + matplotlib | 候选区域 bedGraph | PNG × 3 | 1T | <2 min | 1 GB | <1 MB | ~3 MB（3图） | **→ 极小** |
| **12** report | Python | 各步中间文件 + PNG + SVG | HTML 报告 | 1T | <2 min | 1 GB | — | ~5–15 MB | — |

> **实测数据**：B10002（N≈16 Mpairs）比对 47 min，BAM 1.5 GB；B10197（N≈26 Mpairs）比对 124 min（含大数据sort开销），BAM 2.8 GB。
> bwameth BAM 大小经验值：≈ 双端 clean FASTQ 总大小 × 0.85（B10002: 1.845 GB → 1.5 GB BAM）

一次性基础设施（不计入每样本时间）：

| 一次性步骤 | 工具 | 时间 | 内存 | 存储（输出） |
| --------- | ---- | ---- | ---- | ---------- |
| bwameth c2t 索引构建 | bwameth.py index（内部调用 bwa index） | **~2.5 h** ✓实测（BWT 110 min + SA 45 min） | **9.4 GB**（SA阶段） | 共 6 个文件合计 **~17 GB**（见下表） |
| Bismark索引（备用） | bowtie2-build × 2 | ~2 h | 16 GB | CT + GA 各 4.1 GB = 8.2 GB |
| 参考基因组 FASTA | — | — | — | 3.1 GB（已有） |

#### bwameth c2t 索引文件说明

`bwameth.py index /data/meth/bwa/hg38.fa` 产生 6 个文件，**全部必需**，缺少任何一个都无法比对：

| 文件 | 大小 | 用途 | 可否删除 |
| ---- | ---- | ---- | ------- |
| `hg38.fa.bwameth.c2t` | 6.1 GB | C→T 转换后的参考基因组 FASTA。bwameth 将所有 C 转为 T，使亚硫酸盐处理后的 reads 能在此参考上比对。MethylDackel 提取甲基化时仍用原始 `hg38.fa`，不用此文件。 | ❌ |
| `hg38.fa.bwameth.c2t.bwt` | 6.0 GB | BWT（Burrows-Wheeler Transform）索引。bwa mem 种子搜索的核心数据结构，加载到内存后用于快速定位 reads 的比对位置。 | ❌ |
| `hg38.fa.bwameth.c2t.sa` | 3.0 GB | 后缀数组（Suffix Array）。与 BWT 配合用于种子精确定位和延伸。SA 构建是索引构建最慢的阶段（约 45 min）。 | ❌ |
| `hg38.fa.bwameth.c2t.pac` | 1.5 GB | 2-bit 压缩序列。bwa mem 找到比对位置后，用此文件做序列比较（Smith-Waterman 打分），计算 CIGAR 字符串。 | ❌ |
| `hg38.fa.bwameth.c2t.amb` | 42 KB | 模糊碱基表（N 区域位置和长度）。bwa mem 比对时跳过参考中的 N 区域。 | ❌ |
| `hg38.fa.bwameth.c2t.ann` | 44 KB | 染色体注释（名称、长度）。bwa mem 用于生成 SAM 文件的 `@SQ` 头部行，samtools 也依赖此信息建立 BAM 索引。 | ❌ |

> **内存占用说明**：比对时 `.bwt`（6 GB）+ `.sa`（3 GB）+ `.pac`（1.5 GB）全部映射进内存，合计约 10–11 GB；加上 bwa mem 自身缓冲区约 2–3 GB，峰值 **13–14 GB**（实测）。

---

### 9.3 计算公式

#### 存储（单样本全流程磁盘占用峰值）

```text
峰值占用 = FASTQ_in + FASTQ_clean + BAM_aln + BAM_dedup + BAM_ont + Meth_text

= N × 58 × 2                                         # 原始 FASTQ（R1+R2）
+ N × R_q × 58 × 2                                   # clean FASTQ（R1+R2）
+ N × R_q × R_aln × 100                              # bwameth BAM（实测≈clean FASTQ总×0.85）
+ N × R_q × R_aln × R_dup × 100                      # dedup BAM
+ N × R_q × R_aln × R_dup × R_ont × 100              # on-target BAM
+ 50                                                  # 甲基化文本（MB，靶向，近似常数）

≈ N × [116 + 104 + 75 + 41 + 23] + 50   MB
≈ N × 360 MB + 50 MB
```

**估算示例（N=16，B10002 实测）：** 16 × 360 + 50 ≈ **5.8 GB / 样本**

若中间文件及时清理（只保留 FASTQ_in + BAM_ont + Meth），最小存储：

```text
最小存储 ≈ N × (116 + 23) + 50 = N × 139 MB + 50 MB
```

N=16 时约 **2.3 GB / 样本**

#### 时间（单样本，8 核参考机）

```text
T_total ≈ T_preprocess + T_align + T_dedup + T_target + T_mbias + T_extract + 常数

T_preprocess = N / 10                           # 分钟（fastp 8线程）
T_align      = N × 3                            # 分钟（bwameth 8线程，实测B10002: 16×3=48≈47 min✓）
T_dedup      = N × R_q × R_aln × 0.75          # 分钟（umi_tools 单线程，实测B10002: 10 min✓）
T_target     = N_dedup × 0.15                   # 分钟（bedtools 极快）
T_mbias      = N_ont × 0.5                      # 分钟（MethylDackel 快于bismark_extractor）
T_extract    = N_ont × 0.5                      # 分钟（MethylDackel extract）
常数         ≈ 5                                # 分钟（filter + diff + report）
```

**简化公式**（忽略小项，R_q=0.90, R_aln=0.88）：

```text
T_total(分钟) ≈ N × 4.1 + 5
```

N=16 时约 **71 分钟 ≈ 1.2 小时**（比对 ~55%，去重 ~14%，其余 ~31%）

> 对比 bismark：T_align = N × 6，T_total ≈ N × 9，N=16 约 149 min。**bwameth 节省约 40% 总时间。**

多样本并行时（k 个样本，机器有 C 核）：

```text
并行度 p = floor(C / 8)             # 每个样本 bwameth 占用 ~8 核
T_batch(分钟) ≈ ceil(k / p) × (N × 4.1 + 5)
```

#### 内存

```text
内存峰值 = max(各步骤峰值)
         = bwameth 比对阶段峰值（实测）
         = bwa mem 工作内存 + c2t 索引（BWT 6.0G + SA 3.0G 按需加载）
         ≈ 13–14 GB（实测：B10002 12.8 GB，B10197 14.2 GB）

推荐服务器内存：16 GB（单样本）/ 32 GB（2 样本并行）/ 64 GB（4 样本并行）
```

---

### 9.4 数据流量漏斗图（N=16 Mpairs，B10002 实测）

```text
输入 FASTQ           0.94 GB × 2 = 1.88 GB  (N×58MB ×2, N=16)
  ↓ fastp (×R_q=0.90)
clean FASTQ          0.84 GB × 2 = 1.69 GB  （实测: 935+910 MB = 1.80 GB）
  ↓ bwameth (×R_q×R_aln≈×0.79)
bwameth sorted BAM   1.50 GB                 ✓实测 B10002
  ↓ umi_tools dedup (×R_dup=0.375)  → 自动删除 sorted BAM
dedup BAM            0.56 GB                 ✓实测 B10002
  ↓ bedtools intersect (×R_ont=0.362) → 自动删除 dedup BAM + 未排序中间 BAM
on-target sorted BAM 0.20 GB                ✓实测 B10002（主要分析数据）
  ↓ MethylDackel extract
methylation bedGraph ~1 MB                  ←── 仅靶向CpG位点
  ↓ filter + 候选区域
candidate bedGraph   <1 MB                  ←── 88个区域
```

> **存储关键比例**：on-target BAM ≈ 原始 FASTQ 的 11%；最终甲基化文件 < 0.1%。
> **时间关键比例**：bwameth 比对约占总时间 55%，umi_tools 去重约占 40%，两步合计 ~95%。
> **自动清理**：align → dedup 完成后删 sorted BAM；target 完成后删 dedup BAM 和未排序中间 BAM。

---

## 十、OC_plasma 全量运行计划

### 10.1 数据概况

| 样本 | R1 大小 | 估算 N（Mpairs） | 类别 |
| ---- | ------- | --------------- | ---- |
| B10002 | 1.2 GB | ~16 | 恶性（已测试） |
| B10011 | 1.6 GB | ~21 | — |
| B10014 | 1.4 GB | ~18 | — |
| B10018 | 1.7 GB | ~22 | — |
| B10019 | 1.4 GB | ~18 | — |
| B10188 | 1.3 GB | ~17 | — |
| B10189 | 2.8 GB | ~37 | — ⚠️ 最大样本，sort 可能溢出临时文件 |
| B10194 | 1.6 GB | ~21 | — |
| B10195 | 1.3 GB | ~17 | — |
| B10197 | 1.9 GB | ~26 | 良性（已测试） |
| **合计** | **18.3 GB** | **~213** | 10 个样本 |

输入路径：`/mnt/f/meth/OC_plasma/`（F 盘，20 个文件，总大小约 37 GB）

---

### 10.2 资源预估（顺序执行，单机 8 线程）

| 步骤 | 工具 | 10 样本合计时间 | 内存峰值 | 备注 |
| ---- | ---- | -------------- | ------- | ---- |
| preprocess | fastp | ~30 min | 0.5 GB | 8 线程，I/O 为主 |
| align | bwameth | **~10–11 h** | 13–14 GB | T = N×3 min/样本；B10189 预计 ≥120 min |
| dedup | umi_tools | ~2 h | 2–4 GB | 单线程，T = N×0.75 min/样本 |
| target | bedtools | ~40 min | 1–2 GB | 极快 |
| mbias | MethylDackel | ~30 min | 2–4 GB | 生成 SVG 供确认 |
| extract | MethylDackel | ~10 min | 2–4 GB | 输出 CpG.bedGraph |
| filter+diff+report | Python | ~15 min | <1 GB | — |
| **总计** | — | **~14–15 h** | **14 GB** | — |

存储需求：

| 阶段 | 空间 |
| ---- | ---- |
| 输入 fastq（已有，F 盘） | ~37 GB |
| clean fastq（预处理后） | ~34 GB |
| sorted BAM × 10 | ~21 GB |
| dedup BAM × 10 | ~13 GB |
| targeted sorted BAM × 10（最终保留） | ~5.5 GB |
| 甲基化/报告文件 | < 50 MB |
| **峰值（全部中间文件同时保留）** | **~115 GB** |
| **最终保留（clean fastq + targeted BAM + 结果）** | **~40 GB** |

> ⚠️ E 盘当前可用 134 GB，运行期间请确保剩余 ≥ 120 GB（峰值 115 GB + 缓冲）。

---

### 10.3 分步执行命令

> 前提：在 Windows WSL 终端中执行；工作目录为 `/mnt/e/projects/hq/meth`。

**环境激活：**

```bash
conda activate meth
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
```

| 步骤 | 命令 | 预计时长 | 完成标志 |
| ---- | ---- | ------- | ------- |
| **Step 4** 预处理 | `python src/scripts/run_analysis.py --plasma --step preprocess /mnt/f/meth/OC_plasma/` | ~30 min | 屏幕输出 `workdir: temp/YYYYMMDD_HHMMSS` |
| **Step 5** 比对 | `python src/scripts/run_analysis.py --plasma --step align --workdir $WORKDIR` | ~10–11 h | `[✓] align` |
| **Step 6** 去重 | `python src/scripts/run_analysis.py --plasma --step dedup --workdir $WORKDIR` | ~2 h | `[✓] dedup` |
| **Step 7** On-target | `python src/scripts/run_analysis.py --plasma --step target --workdir $WORKDIR` | ~40 min | `[✓] target` |
| **Step 7.5** M-bias | `python src/scripts/run_analysis.py --plasma --step mbias --workdir $WORKDIR` | ~30 min | SVG 文件生成于 workdir |
| *(检查 M-bias)* | 用浏览器查看 `$WORKDIR/*.mbias_OT.svg` / `*_OB.svg` | — | R2 3'末端曲线平稳 → 继续；异常 → 修改 `trim_tail2` 重跑 |
| **Step 8** 甲基化提取 | `python src/scripts/run_analysis.py --plasma --step extract --workdir $WORKDIR` | ~10 min | `[✓] extract` |
| **Step 9** 过滤 | `python src/scripts/run_analysis.py --plasma --step filter --workdir $WORKDIR` | ~5 min | `[✓] filter` |
| **Step 10** 差异分析 | `python src/scripts/run_analysis.py --plasma --step diff --workdir $WORKDIR` | <2 min | `result/$WORKDIR/diff_methylation.txt` |
| **Step 11** 可视化 | `python src/scripts/run_analysis.py --plasma --step visualize --workdir $WORKDIR` | <2 min | `result/$WORKDIR/density_plot.png` 等 3 图 |
| **Step 12** 报告 | `python src/scripts/run_analysis.py --plasma --step report --workdir $WORKDIR` | <2 min | `result/$WORKDIR/analysis_report.html`（含 M-bias SVG + 分布图） |

**完整 shell 脚本（一键执行，后台运行）：**

```bash
#!/bin/bash
# OC_plasma 全量分析脚本
# 用法：bash run_oc_plasma.sh 2>&1 | tee logs/oc_plasma_run.log

set -e
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
PY=/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin/python

echo "=== Step 4: 预处理 ===" && $PY src/scripts/run_analysis.py --plasma --step preprocess /mnt/f/meth/OC_plasma/

# 从状态文件读取 workdir（取最新的 temp/YYYYMMDD_* 目录）
WORKDIR=$(ls -dt temp/20[0-9]* | head -1)
echo "workdir: $WORKDIR"

echo "=== Step 5: 比对 ===" && $PY src/scripts/run_analysis.py --plasma --step align   --workdir $WORKDIR
echo "=== Step 6: 去重 ===" && $PY src/scripts/run_analysis.py --plasma --step dedup   --workdir $WORKDIR
echo "=== Step 7: On-target ===" && $PY src/scripts/run_analysis.py --plasma --step target  --workdir $WORKDIR
echo "=== Step 7.5: M-bias ===" && $PY src/scripts/run_analysis.py --plasma --step mbias   --workdir $WORKDIR

echo "=== 请检查 M-bias SVG 后按 Enter 继续，或 Ctrl+C 中止调整参数 ==="
read

echo "=== Step 8: 甲基化提取 ===" && $PY src/scripts/run_analysis.py --plasma --step extract --workdir $WORKDIR
echo "=== Step 9: 过滤 ===" && $PY src/scripts/run_analysis.py --plasma --step filter  --workdir $WORKDIR
echo "=== Step 10: 差异分析 ===" && $PY src/scripts/run_analysis.py --plasma --step diff    --workdir $WORKDIR
echo "=== Step 11: 可视化 ===" && $PY src/scripts/run_analysis.py --plasma --step visualize --workdir $WORKDIR
echo "=== Step 12: 报告 ===" && $PY src/scripts/run_analysis.py --plasma --step report  --workdir $WORKDIR

echo "=== 全流程完成！结果目录：result/$WORKDIR ==="
```

---

### 10.4 断点续跑说明

每一步完成后状态自动记录到 `$WORKDIR/pipeline_state.json`。若中途中断（如断电、OOM），直接重跑该步骤即可：

```bash
# 查看已完成步骤
python src/scripts/run_analysis.py --list-steps --workdir $WORKDIR

# 从中断步骤继续（以 align 为例）
python src/scripts/run_analysis.py --plasma --step align --workdir $WORKDIR
```

---

### 10.5 完成后清理建议

```bash
# 删除大型中间文件（align + dedup BAM），保留 targeted BAM 和结果
rm $WORKDIR/*.bwameth_pe.sorted.bam
rm $WORKDIR/*.bwameth_pe.sorted.bam.bai
rm $WORKDIR/*.bwameth_pe.dedup.bam
rm $WORKDIR/*.bwameth_pe.dedup.bam.bai
# 可选：删除 clean fastq（原始 fastq 仍在 /mnt/f/meth/OC_plasma/）
rm $WORKDIR/*.clean.R*.fq.gz

# 清理后保留约 6–8 GB（targeted BAM × 10 + 甲基化/报告文件）
```

---

## 十一、特征降维与诊断分类器（已实现）

> **状态**：已完成开发（2026-04-03）。当前 10 样本（5良性/5恶性）初步分析已发现 71 个 DMC，|Δβ| 均值 0.481。53 恶性 / 20 良性样本全队列 pipeline 正在运行中，完成后可直接进入 Step 14。

### 11.1 目标

在更大队列中验证现有 DMC，通过机器学习精简出核心 CpG 面板，并训练恶性概率分类器，支持未来临床液体活检诊断。

---

### 11.2 阶段一：扩大队列，重新验证 DMC（Step 13）

**背景**：当前 10 个样本（5/5）发现 71 个 DMC，统计效力有限，需在 `data/patient3.csv` 的 83 个患者中重新运行全流程，得到更可靠的差异位点集合。

**实施方案**：

- 运行完整流程（Steps 4–12），输入为 `patient3.csv` 中所有可用样本
- 调整差异分析阈值（`diff_p_value` / `diff_delta_beta`），对比 FDR 校正前后的 DMC 数量
- 输出：稳定的候选 DMC 列表（预期 20–100 个）

```bash
python src/scripts/run_analysis.py --plasma \
    --patient-csv data/patient3.csv \
    --workdir temp/patient3_cohort
```

---

### 11.3 阶段二：特征降维，精简核心 CpG 面板（Step 14）

**目标**：从验证后的 DMC 集合中剔除冗余，筛出最具区分度的核心位点。

**实现方案（LASSO + RF 并集，可通过 `config.toml [ml]` 调整）**：

| 方法 | 实现 | 说明 |
| --- | --- | --- |
| **L1 逻辑回归（LASSO）** | `LogisticRegressionCV(penalty="l1", solver="liblinear")` | 交叉验证选最优 C，非零系数特征入选 |
| **随机森林特征重要性** | `RandomForestClassifier(n_estimators=500)` | Gini 不纯度重要性 > 0 的特征入选 |
| **并集（默认）** | `feature_selection_mode = "union"` | LASSO 或 RF 选中任一即入选，可改 "lasso"/"rf" |

**覆盖率过滤**：DMC 位点在 < `min_feature_coverage`（默认 50%）样本中无覆盖时丢弃，避免大量 NaN 特征进入模型。

**输出（保存至 `result/$WORKDIR/`）**：

- `feature_importance.csv`：每个 DMC 的 LASSO 系数 + RF 重要性 + 入选标记
- `panel_cpg.csv`：最终入选的特征位点列表（`feature` 列，格式 `chr:pos`）
- `feature_matrix.npz`：样本×特征矩阵缓存（供 Step 15 复用，避免重读 bedGraph）

```bash
python src/scripts/run_analysis.py --plasma --step feature_select --workdir $WORKDIR
```

---

### 11.4 阶段三：构建诊断分类器（Step 15）

**目标**：以核心 CpG 面板甲基化水平为输入，LOO-CV 评估恶性诊断性能。

**实现方案**：

| 模型 | 实现 | 参数 |
| --- | --- | --- |
| **逻辑回归（LR）** | `LogisticRegression(class_weight="balanced")` | `lr_C`（默认 1.0） |
| **SVM（RBF 核）** | `SVC(kernel="rbf", probability=True, class_weight="balanced")` | `svm_C`（默认 1.0） |
| **随机森林** | `RandomForestClassifier(class_weight="balanced")` | `rf_n_estimators`（默认 200） |

**LOO-CV 实现细节**：

- 每轮在训练集上独立 `fit` StandardScaler，再 `transform` 测试样本，严格防止数据泄漏
- NaN 统一按全局列中位数插补（在 LOO 循环外，可视为已知先验分布，不影响泄漏判断）
- `class_weight="balanced"` 应对 53 恶性 / 20 良性的类别不均衡

**输出（保存至 `result/$WORKDIR/`）**：

- `roc_curve.png`：三分类器 ROC 曲线对比（含各自 AUC）
- `confusion_matrix.png`：各分类器混淆矩阵（2×2）
- `classifier_predictions.csv`：每个样本的预测概率和标签（每分类器一行）
- `classifier_summary.csv`：AUC / 准确率 / 灵敏度 / 特异度汇总

```bash
python src/scripts/run_analysis.py --plasma --step classify --workdir $WORKDIR
```

---

### 11.5 实施路线图

```text
当前状态（10样本，5/5，71个DMC）
  │
  ▼  Step 13：扩大队列（patient3.csv，53恶性+20良性）
重新验证 DMC（更大统计效力，FDR校正）
  │
  ▼  Step 14：LASSO + 随机森林特征降维   ← 已实现
核心 CpG 面板（panel_cpg.csv）
  │
  ▼  Step 15：LOO-CV 评估 LR / SVM / RF   ← 已实现
ROC 曲线、混淆矩阵、分类器汇总
  │
  ▼  （未来）独立验证队列评估
AUC / 敏感性 / 特异性报告
```

---

### 11.6 新增步骤对流程的影响

| 步骤 | 状态 | 涉及文件 | 依赖 |
| --- | --- | --- | --- |
| Step 13 | 复用现有流程 | `patient3.csv` | Steps 4–12 |
| Step 14 `feature_select` | **已实现** | `src/scripts/modules/feature_utils.py` | diff 结果（Step 10）+ bedGraph（Step 8） |
| Step 15 `classify` | **已实现** | `src/scripts/modules/classifier_utils.py` | feature_select 结果（Step 14）可选 |

**`config.toml [ml]` 配置节（已添加）**：

```toml
[ml]
# 特征选择模式: "union"（LASSO和RF并集）/ "lasso" / "rf"
feature_selection_mode = "union"
# DMC位点最低覆盖率（低于此比例的样本无覆盖则丢弃该特征）
min_feature_coverage = 0.5
# LASSO 交叉验证折数
lasso_cv_folds = 5
# 随机森林树数量（特征选择 + 分类均用此值）
rf_n_estimators = 500
# RF 最小特征重要性阈值（0 = 保留所有重要性 > 0 的特征）
rf_min_importance = 0.0
# 逻辑回归正则化参数 C
lr_C = 1.0
# SVM 正则化参数 C
svm_C = 1.0
```

**新增 Python 依赖**（已更新 `requirements.txt` / `pyproject.toml`）：

```text
scikit-learn>=1.3.0
```

安装：

```bash
pip install scikit-learn
# 或
conda install -n meth scikit-learn
```

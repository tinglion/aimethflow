#!/usr/bin/env python3
"""
run_analysis.py - 执行完整的甲基化分析流程

支持两种运行模式：
  【全流程模式】一次性运行所有步骤（默认）
  【分步模式】  指定单个步骤运行，通过状态文件在步骤间传递数据

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
基本用法:
  python src/scripts/run_analysis.py --plasma <文件/文件夹>

分步执行（血浆流程全部 8 步）:
  # 第1步：fastp + UMI提取 → 生成 workdir，打印路径
  python src/scripts/run_analysis.py --plasma --step preprocess <文件夹>

  # 后续各步骤：通过 --workdir 指定上一步生成的目录
  python src/scripts/run_analysis.py --plasma --step align     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step dedup     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step target    --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step mbias     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step extract   --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step filter    --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step diff      --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step visualize --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step report    --workdir temp/20240101_120000

  # 查看步骤定义和当前状态
  python src/scripts/run_analysis.py --list-steps
  python src/scripts/run_analysis.py --list-steps --workdir temp/20240101_120000

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
步骤说明:
  preprocess  4   fastp质控 + UMI提取（血浆：--umi --trim_front）
  align       5   Bismark双端比对（bismark --genome）
  dedup       6   UMI去重（umi_tools dedup --paired）
  target      7   On-target过滤（bedtools intersect -abam panel.bed）
  mbias       7.5 M-bias检查（MethylDackel mbias），确认tail2裁剪参数
  extract     8   甲基化信息提取（bismark_methylation_extractor）
  filter      9   甲基化数据过滤 + 候选区域过滤
  diff        10  差异甲基化分析
  visualize   11  可视化图表生成（密度图/散点图/热图）
  report      12  生成分析报告
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import os
import time
import argparse
import logging
import shutil

# 加载配置（早于其他模块，以便注入 PATH）
# 注意：load_config 在此处尚未 import，使用 toml 直接加载
def _inject_conda_path():
    """将 conda meth 环境的 bin 目录注入到 PATH，确保子进程能找到工具。"""
    import toml
    _cfg_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "config.toml")
    try:
        _cfg = toml.load(_cfg_path)
        _bin = _cfg.get("tools", {}).get("conda_env_bin", "")
        if _bin and os.path.isdir(_bin) and _bin not in os.environ.get("PATH", ""):
            os.environ["PATH"] = _bin + ":" + os.environ.get("PATH", "")
    except Exception:
        pass

_inject_conda_path()

# 导入各个模块的函数
from utils.logger import setup_logger
from utils.config import load_config, get_output_dirs, create_output_dirs
from scripts.modules import (
    find_files,
    preprocess_fastq,
    align_methylation,
    deduplicate_bam,
    extract_target_regions,
    check_mbias,
    extract_methylation_methyldackel,
    extract_methylation_info,
    merge_methylation_files,
    filter_methylation_data,
    filter_to_candidate_regions_wrapper,
    perform_diff_analysis,
    visualize_results,
    generate_report,
)
from methylation.sample_grouping import resolve_patient_files

# 加载配置文件
config = load_config()

# 生成日志文件路径
log_dir = config["directories"]["log_dir"]
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/analysis_{time.strftime('%Y%m%d_%H%M%S')}.log"
logger = setup_logger(log_file=log_file, console_level=logging.INFO, file_level=logging.DEBUG)


def check_presets():
    """
    检查分析前的预设信息

    Returns:
        bool: 预设信息检查是否通过
    """
    logger.info("开始检查预设信息")

    # 1. 检查探针区域文件
    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if probe_panel:
        if not os.path.exists(probe_panel):
            logger.error(f"探针区域文件不存在: {probe_panel}")
            logger.warning("将跳过靶向区域提取步骤")
        else:
            logger.info(f"探针区域文件: {probe_panel}")
    else:
        logger.warning("未配置探针区域文件，跳过靶向区域提取")

    # 2. 检查参考基因组
    reference_genome = config.get("extraction", {}).get("reference_genome")
    if not reference_genome or not os.path.exists(reference_genome):
        logger.error(f"参考基因组不存在: {reference_genome}")
        return False
    logger.info(f"参考基因组: {reference_genome}")

    # 3. 检查Bismark索引（bismark_genome_preparation 生成 Bisulfite_Genome/CT_conversion/BS_CT.*.bt2）
    bismark_index = os.path.dirname(reference_genome)
    bismark_index_file = os.path.join(
        bismark_index, "Bisulfite_Genome", "CT_conversion", "BS_CT.1.bt2"
    )
    if not os.path.exists(bismark_index_file):
        logger.warning(f"Bismark索引未找到: {bismark_index_file}")
        logger.warning("请先运行: bismark_genome_preparation --bowtie2 " + bismark_index)
    else:
        logger.info(f"Bismark索引: {bismark_index}/Bisulfite_Genome/")

    # 4. 检查必要的工具（血浆模式必须: fastp, bismark, deduplicate_bismark, samtools, bedtools, MethylDackel, umi_tools）
    required_tools = ["fastp", "bismark", "deduplicate_bismark", "samtools", "bedtools", "MethylDackel"]
    # umi_tools 可能在 meth env 外，单独检测
    umi_tools_bin = shutil.which("umi_tools") or config.get("tools", {}).get("umi_tools_fallback", "")
    if umi_tools_bin and os.path.isfile(umi_tools_bin):
        logger.info(f"umi_tools: {umi_tools_bin}")
        # 确保 umi_tools 所在目录也在 PATH 中
        umi_dir = os.path.dirname(umi_tools_bin)
        if umi_dir and umi_dir not in os.environ.get("PATH", ""):
            os.environ["PATH"] = umi_dir + ":" + os.environ.get("PATH", "")
    else:
        logger.warning("umi_tools 未找到 — 血浆模式去重将回退到 deduplicate_bismark（不推荐）")

    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool):
            logger.info(f"工具 {tool}: {shutil.which(tool)}")
        else:
            logger.warning(f"工具 {tool} 未安装或不在PATH中")
            missing_tools.append(tool)

    if missing_tools:
        logger.warning(f"以下工具未找到: {', '.join(missing_tools)}")

    logger.info("预设信息检查完成")
    return True


# ──────────────────────────────────────────────
# 分步执行：步骤定义与状态文件
# ──────────────────────────────────────────────

# 步骤有序列表：(step_key, display_num, 中文说明)
PIPELINE_STEPS = [
    ("preprocess", "4",   "fastp质控 + UMI提取"),
    ("align",      "5",   "Bismark双端比对"),
    ("dedup",      "6",   "UMI去重"),
    ("target",     "7",   "On-target过滤"),
    ("mbias",      "7.5", "M-bias检查"),
    ("extract",    "8",   "甲基化信息提取"),
    ("filter",     "9",   "甲基化数据过滤 + 候选区域"),
    ("diff",       "10",  "差异甲基化分析"),
    ("visualize",  "11",  "可视化图表生成"),
    ("report",     "12",  "生成分析报告"),
]
STEP_KEYS = [s[0] for s in PIPELINE_STEPS]
_STATE_FILE = "pipeline_state.json"


def _state_path(workdir):
    return os.path.join(workdir, _STATE_FILE)


def _load_state(workdir):
    """从 workdir 下的 pipeline_state.json 加载上一步的输出。"""
    import json
    path = _state_path(workdir)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_state(workdir, state):
    """将当前状态写入 workdir/pipeline_state.json（追加更新，不覆盖已有字段）。"""
    import json
    path = _state_path(workdir)
    existing = _load_state(workdir)
    existing.update(state)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(existing, f, ensure_ascii=False, indent=2)


def _print_step_list(workdir=None):
    """打印步骤列表及各步完成状态（若提供 workdir）。"""
    state = _load_state(workdir) if workdir else {}
    completed = state.get("completed_steps", [])
    print("\n┌─────────────────────────────────────────────────────────┐")
    print("│              血浆甲基化流程 — 步骤一览                  │")
    print("├──────────────┬──────┬─────────────────────────────────┤")
    print("│ --step 参数  │ 编号 │ 说明                            │")
    print("├──────────────┼──────┼─────────────────────────────────┤")
    for key, num, desc in PIPELINE_STEPS:
        mark = "✓" if key in completed else " "
        print(f"│ [{mark}] {key:<10}│  {num:<4}│ {desc:<31} │")
    print("└──────────────┴──────┴─────────────────────────────────┘")
    if workdir:
        print(f"\n工作目录: {workdir}")
        print(f"状态文件: {_state_path(workdir)}")
    print()


def main():
    """
    主函数
    """
    # ── 解析命令行参数 ──
    parser = argparse.ArgumentParser(
        description="执行完整的甲基化分析流程（支持全流程或分步执行）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "files", nargs="*",
        help="FastQ文件或文件夹路径列表（分步模式且已有workdir时可省略）",
    )

    # ── 分步执行选项 ──
    parser.add_argument(
        "--step", metavar="STEP",
        choices=STEP_KEYS,
        help=f"只执行指定步骤（可选: {', '.join(STEP_KEYS)}）",
    )
    parser.add_argument(
        "--workdir", metavar="DIR",
        help="指定已有工作目录以恢复分步执行（会读取该目录下的 pipeline_state.json）",
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="列出所有步骤说明及完成状态后退出",
    )

    # ── 模式开关 ──
    parser.add_argument(
        "--patient-csv", metavar="CSV",
        help="患者ID列表文件（每行一个ID），从Excel自动查找对应R1/R2文件，与 files 互斥",
    )
    parser.add_argument(
        "--data-dir", metavar="DIR",
        help="FastQ数据目录，与 --patient-csv 配合使用，默认读取 config.toml [data] plasma_dir",
    )
    parser.add_argument(
        "--temp-dir", metavar="DIR",
        help="临时文件存放目录（覆盖 config.toml [directories] temp_dir），例如 /mnt/g/temp",
    )
    parser.add_argument(
        "--plasma", action="store_true", default=True,
        help="启用血浆cfDNA模式（默认已开启）",
    )
    parser.add_argument(
        "--group-method", choices=["v1", "v2"], default="v2",
        help="分组方法: v1（基于索引）或 v2（基于Excel文件）",
    )

    # ── 全流程模式的 skip 开关（分步模式下这些由 --step 隐式控制） ──
    parser.add_argument("--skip-quality",           action="store_true", default=True)
    parser.add_argument("--skip-preprocessing",     action="store_true", default=False)
    parser.add_argument("--skip-alignment",         action="store_true", default=False)
    parser.add_argument("--skip-deduplication",     action="store_true", default=False)
    parser.add_argument("--skip-target-extraction", action="store_true", default=False)
    parser.add_argument("--skip-extraction",        action="store_true", default=False)
    parser.add_argument("--skip-filtering",         action="store_true", default=False)
    parser.add_argument("--skip-candidate-regions", action="store_true", default=False)
    parser.add_argument("--skip-diff",              action="store_true", default=False)
    parser.add_argument("--skip-visualization",     action="store_true", default=False)
    parser.add_argument("--skip-report",            action="store_true", default=False)
    parser.add_argument("--clean-temp",             action="store_true", default=False)
    parser.add_argument("--config",  type=str, help="使用指定的配置文件")
    parser.add_argument("--verbose", action="store_true", default=False)

    args = parser.parse_args()

    # ── --list-steps ──
    if args.list_steps:
        _print_step_list(args.workdir)
        return

    # ── 自定义配置文件 ──
    global config
    if args.config and os.path.exists(args.config):
        from utils.config import load_config as _load_cfg
        config = _load_cfg(args.config)

    # ── --temp-dir 覆盖临时目录 ──
    if args.temp_dir:
        config["directories"]["temp_dir"] = args.temp_dir
        os.makedirs(args.temp_dir, exist_ok=True)
        logger.info(f"临时目录已覆盖为: {args.temp_dir}")

    # ── 详细输出 ──
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        for h in logger.handlers:
            if getattr(h, "name", None) == "console":
                h.setLevel(logging.DEBUG)

    # ── 预设检查 ──
    if not check_presets():
        logger.error("预设信息检查失败，分析终止")
        return

    # ── 确定工作目录 ──
    if args.workdir:
        workdir = args.workdir
        # 目录不存在则创建（首次 preprocess 即可指定 --workdir）
        os.makedirs(workdir, exist_ok=True)
        # 复用 workdir 的时间戳，使 result_dir 与 workdir 对应同一次运行
        timestamp = os.path.basename(os.path.normpath(workdir))
        output_dirs = get_output_dirs(config, timestamp=timestamp)
        result_dir = output_dirs["result"]
        os.makedirs(result_dir, exist_ok=True)
    else:
        # 新建目录
        output_dirs = get_output_dirs(config)
        create_output_dirs(output_dirs)
        workdir = output_dirs["temp"]
        result_dir = output_dirs["result"]

    logger.info(f"工作目录: {workdir}")

    # ── 加载状态（分步恢复） ──
    state = _load_state(workdir)

    # ── 查找输入文件（优先 --patient-csv，其次 files 参数，最后从状态恢复） ──
    if getattr(args, "patient_csv", None):
        # 从患者ID CSV + Excel 自动解析文件列表
        data_dir = getattr(args, "data_dir", None)
        fastq_files = resolve_patient_files(args.patient_csv, data_dir=data_dir)
        bam_files_input = []
        methylation_files_input = []
        if not fastq_files:
            logger.error(f"--patient-csv 未解析到任何文件，请检查 {args.patient_csv} 及数据目录")
            return
        _save_state(workdir, {"fastq_files": fastq_files})
    elif args.files:
        found_files = find_files(args.files)
        fastq_files = found_files["fastq"]
        bam_files_input = found_files["bam"]
        methylation_files_input = found_files["methylation"]
        # 更新状态中的原始输入列表
        if fastq_files:
            state["fastq_files"] = fastq_files
        _save_state(workdir, {"fastq_files": fastq_files})
    else:
        fastq_files = state.get("fastq_files", [])
        bam_files_input = []
        methylation_files_input = []
        if not fastq_files:
            logger.error("未提供输入文件，且 workdir 中也无已记录的 fastq_files")
            return

    logger.info(f"输入 FastQ 文件 ({len(fastq_files)} 个): {fastq_files}")

    # ── 从状态恢复各步骤输出 ──
    processed_files       = state.get("processed_files", fastq_files)
    bam_files             = state.get("bam_files", bam_files_input)
    deduplicated_files    = state.get("deduplicated_files", bam_files)
    target_files          = state.get("target_files", deduplicated_files)
    sample1_methylation   = state.get("sample1_methylation")
    sample2_methylation   = state.get("sample2_methylation")
    filtered_sample1      = state.get("filtered_sample1")
    filtered_sample2      = state.get("filtered_sample2")
    candidate_sample1     = state.get("candidate_sample1")
    candidate_sample2     = state.get("candidate_sample2")
    benign_samples        = state.get("benign_samples", [])
    malignant_samples     = state.get("malignant_samples", [])

    # ── 判断哪些步骤需要执行 ──
    # 分步模式：只执行指定步骤；全流程模式：根据 skip_* 标志
    single_step = args.step  # None → 全流程

    def _should_run(step_key, skip_flag):
        if single_step:
            return step_key == single_step
        return not skip_flag

    # =========================================================
    # 步骤 4: preprocess
    # =========================================================
    if _should_run("preprocess", args.skip_preprocessing) and fastq_files:
        logger.info("━━ 步骤 4: fastp质控 + UMI提取 ━━")
        processed_files = preprocess_fastq(fastq_files, workdir, plasma=args.plasma, config=config)
        if processed_files:
            _save_state(workdir, {"processed_files": processed_files,
                                  "completed_steps": state.get("completed_steps", []) + ["preprocess"]})
            state["processed_files"] = processed_files
            logger.info(f"  → 预处理输出: {processed_files}")
        else:
            logger.error("预处理失败，processed_files 为空")
            if single_step:
                return
    elif not single_step:
        logger.info("4. 跳过数据预处理")

    # =========================================================
    # 步骤 5: align
    # =========================================================
    if _should_run("align", args.skip_alignment) and processed_files and not bam_files:
        logger.info("━━ 步骤 5: bwameth 双端比对 ━━")
        bam_files = align_methylation(processed_files, workdir, config, plasma=args.plasma)
        if bam_files:
            _save_state(workdir, {"bam_files": bam_files,
                                  "completed_steps": state.get("completed_steps", []) + ["align"]})
            state["bam_files"] = bam_files
            logger.info(f"  → 比对输出: {bam_files}")
        else:
            logger.error("比对失败，bam_files 为空")
            if single_step:
                return
    elif single_step and bam_files:
        # 单步模式：bam_files 已从 pipeline_state.json 恢复，说明 align 已完成，跳过
        logger.info(f"✓ align 已完成，跳过重复比对（共 {len(bam_files)} 个BAM已记录在状态文件中）")
    elif not single_step:
        logger.info("5. 跳过甲基化序列比对" if not bam_files else "5. 使用已有BAM文件")

    # =========================================================
    # 步骤 6: dedup
    # =========================================================
    if _should_run("dedup", args.skip_deduplication) and bam_files:
        logger.info("━━ 步骤 6: UMI去重 ━━")
        deduplicated_files = deduplicate_bam(bam_files, workdir, plasma=args.plasma)
        if deduplicated_files:
            _save_state(workdir, {"deduplicated_files": deduplicated_files,
                                  "completed_steps": state.get("completed_steps", []) + ["dedup"]})
            state["deduplicated_files"] = deduplicated_files
            logger.info(f"  → 去重输出: {deduplicated_files}")
        else:
            logger.error("去重失败，deduplicated_files 为空")
            if single_step:
                return
    elif not single_step:
        logger.info("6. 跳过PCR重复去除")

    # =========================================================
    # 步骤 7: target
    # =========================================================
    if _should_run("target", args.skip_target_extraction) and deduplicated_files:
        logger.info("━━ 步骤 7: On-target过滤 ━━")
        target_files = extract_target_regions(deduplicated_files, workdir, config)
        if target_files:
            _save_state(workdir, {"target_files": target_files,
                                  "completed_steps": state.get("completed_steps", []) + ["target"]})
            state["target_files"] = target_files
            logger.info(f"  → On-target输出: {target_files}")
        else:
            logger.warning("on-target过滤未生成输出，后续步骤使用去重后BAM")
    elif not single_step:
        logger.info("7. 跳过靶向区域提取")

    # =========================================================
    # 步骤 7.5: mbias
    # =========================================================
    if _should_run("mbias", False) and args.plasma:
        mbias_input = target_files or deduplicated_files
        if mbias_input:
            reference_genome = config.get("extraction", {}).get("reference_genome")
            if reference_genome and os.path.exists(reference_genome):
                logger.info("━━ 步骤 7.5: M-bias检查 ━━")
                check_mbias(mbias_input, workdir, reference_genome, config=config)
                _save_state(workdir, {"mbias_done": True,
                                      "completed_steps": state.get("completed_steps", []) + ["mbias"]})
                logger.info("  → M-bias SVG 已保存至 workdir，请检查 R2 3'端曲线")
            else:
                logger.warning("7.5. 跳过M-bias检查（参考基因组不存在）")
        else:
            logger.warning("7.5. 跳过M-bias检查（无可用BAM）")

    # =========================================================
    # 步骤 8: extract
    # =========================================================
    if _should_run("extract", args.skip_extraction) and not methylation_files_input:
        bam_input = target_files or deduplicated_files
        if bam_input:
            logger.info("━━ 步骤 8: 甲基化信息提取 ━━")
            aligner = config.get("plasma", {}).get("aligner", "bismark") if args.plasma else "bismark"

            if args.plasma and aligner == "bwameth":
                # ── bwameth 路径：MethylDackel extract → 按 Excel 分组 → 合并同组 ──
                logger.info(f"  [提取工具: MethylDackel（对应 aligner=bwameth）]")
                meth_files = extract_methylation_methyldackel(bam_input, workdir, config)
                if meth_files:
                    from methylation.sample_grouping import group_methylation_files_v2
                    benign_files, malignant_files = group_methylation_files_v2(meth_files, fastq_files)
                    logger.info(f"  → 分组：良性 {len(benign_files)} 个，恶性 {len(malignant_files)} 个")

                    s1_path = os.path.join(workdir, "sample1_methylation.bedGraph")
                    s2_path = os.path.join(workdir, "sample2_methylation.bedGraph")

                    def _sample_id(path):
                        """从 bedGraph 路径提取纯样本 ID，如 B10002_CpG.bedGraph → B10002。"""
                        return os.path.basename(path).split("_CpG")[0].split(".")[0].split("_")[0]

                    if benign_files:
                        merge_methylation_files(benign_files, s1_path)
                        sample1_methylation = s1_path
                        benign_samples = [_sample_id(f) for f in benign_files]
                    else:
                        sample1_methylation = meth_files[0]
                        benign_samples = [_sample_id(meth_files[0])]
                        logger.warning("  → 未找到良性样本，以第1个文件作为 sample1")

                    if malignant_files:
                        merge_methylation_files(malignant_files, s2_path)
                        sample2_methylation = s2_path
                        malignant_samples = [_sample_id(f) for f in malignant_files]
                    else:
                        sample2_methylation = meth_files[1] if len(meth_files) > 1 else meth_files[0]
                        malignant_samples = [_sample_id(sample2_methylation)]
                        logger.warning("  → 未找到恶性样本，以第2个文件作为 sample2")
            else:
                # ── bismark 路径：bismark_methylation_extractor（经 extract_methylation_info 封装）──
                logger.info(f"  [提取工具: bismark_methylation_extractor（对应 aligner=bismark）]")
                precomputed_bams = bam_input if args.plasma else None
                sample1_methylation, sample2_methylation = extract_methylation_info(
                    workdir, fastq_files,
                    group_method=args.group_method,
                    config=config,
                    bam_files=precomputed_bams,
                )

            _save_state(workdir, {
                "sample1_methylation": sample1_methylation,
                "sample2_methylation": sample2_methylation,
                "benign_samples": benign_samples,
                "malignant_samples": malignant_samples,
                "completed_steps": state.get("completed_steps", []) + ["extract"],
            })
            state["sample1_methylation"] = sample1_methylation
            state["sample2_methylation"] = sample2_methylation
            state["benign_samples"] = benign_samples
            state["malignant_samples"] = malignant_samples
            logger.info(f"  → 甲基化文件: {sample1_methylation}, {sample2_methylation}")
        else:
            logger.error("无可用BAM文件，无法提取甲基化信息")
            if single_step:
                return
    elif methylation_files_input and not single_step:
        logger.info("8. 使用提供的甲基化数据文件")
        sample1_methylation = methylation_files_input[0] if methylation_files_input else None
        sample2_methylation = methylation_files_input[1] if len(methylation_files_input) > 1 else sample1_methylation
    elif not single_step:
        logger.info("8. 跳过甲基化信息提取")

    # =========================================================
    # 步骤 9: filter
    # =========================================================
    if _should_run("filter", args.skip_filtering) and (sample1_methylation or sample2_methylation):
        logger.info("━━ 步骤 9: 甲基化数据过滤 + 候选区域 ━━")
        if sample1_methylation:
            filtered_sample1 = filter_methylation_data(sample1_methylation, workdir, config)
        if sample2_methylation:
            filtered_sample2 = filter_methylation_data(sample2_methylation, workdir, config)
        if not args.skip_candidate_regions:
            if filtered_sample1:
                candidate_sample1 = filter_to_candidate_regions_wrapper(filtered_sample1, workdir, config)
            if filtered_sample2:
                candidate_sample2 = filter_to_candidate_regions_wrapper(filtered_sample2, workdir, config)
        _save_state(workdir, {
            "filtered_sample1": filtered_sample1, "filtered_sample2": filtered_sample2,
            "candidate_sample1": candidate_sample1, "candidate_sample2": candidate_sample2,
            "completed_steps": state.get("completed_steps", []) + ["filter"],
        })
        state.update({
            "filtered_sample1": filtered_sample1, "filtered_sample2": filtered_sample2,
            "candidate_sample1": candidate_sample1, "candidate_sample2": candidate_sample2,
        })
        logger.info(f"  → 候选区域文件: {candidate_sample1}, {candidate_sample2}")
    elif not single_step:
        logger.info("9. 跳过甲基化数据过滤")

    # =========================================================
    # 步骤 10: diff
    # =========================================================
    if _should_run("diff", args.skip_diff) and candidate_sample1 and candidate_sample2:
        logger.info("━━ 步骤 10: 差异甲基化分析 ━━")
        perform_diff_analysis(candidate_sample1, candidate_sample2, result_dir)
        _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["diff"]})
    elif not single_step:
        logger.info("10. 跳过差异甲基化分析")

    # =========================================================
    # 步骤 11: visualize
    # =========================================================
    if _should_run("visualize", args.skip_visualization):
        vis_input = candidate_sample1 or sample1_methylation
        if vis_input and os.path.exists(vis_input):
            logger.info("━━ 步骤 11: 可视化图表生成 ━━")
            visualize_results(vis_input, result_dir)
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["visualize"]})
        else:
            logger.warning("11. 跳过可视化（无可用甲基化文件）")
    elif not single_step:
        logger.info("11. 跳过可视化图表生成")

    # =========================================================
    # 步骤 12: report
    # =========================================================
    if _should_run("report", args.skip_report):
        logger.info("━━ 步骤 12: 生成分析报告 ━━")
        generate_report(workdir, result_dir,
                        sample1_methylation, sample2_methylation,
                        candidate_sample1, candidate_sample2,
                        benign_samples=benign_samples,
                        malignant_samples=malignant_samples)
        _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["report"]})

    # ── 完成提示 ──
    if single_step:
        print(f"\n✓  步骤 [{single_step}] 完成")
        print(f"   工作目录（下一步传入 --workdir）: {workdir}")
        _print_step_list(workdir)
    else:
        print(f"\n✓  全流程分析完成！结果保存在: {result_dir}")

    # ── 清理 ──
    if args.clean_temp and not single_step:
        try:
            shutil.rmtree(workdir)
            logger.info(f"临时目录已清理: {workdir}")
        except Exception as e:
            logger.error(f"清理临时文件时出错: {e}")


if __name__ == "__main__":
    main()


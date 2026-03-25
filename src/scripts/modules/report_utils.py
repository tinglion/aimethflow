#!/usr/bin/env python3
"""
report_utils.py - 报告生成模块

该模块包含分析报告生成相关的函数
"""

import os
import json
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


def generate_report(
    temp_dir,
    result_dir,
    sample1_methylation,
    sample2_methylation,
    candidate_sample1,
    candidate_sample2,
):
    """
    生成分析报告

    Args:
        temp_dir: 临时目录
        result_dir: 结果目录
        sample1_methylation: 样本1甲基化文件路径
        sample2_methylation: 样本2甲基化文件路径
        candidate_sample1: 样本1候选区域文件路径
        candidate_sample2: 样本2候选区域文件路径
    """
    logger.info("开始生成分析报告")

    # 收集分析信息
    report_data = {
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "sample_info": {
            "sample1": sample1_methylation,
            "sample2": sample2_methylation,
            "candidate_sample1": candidate_sample1,
            "candidate_sample2": candidate_sample2,
        },
        "directories": {"temp": temp_dir, "result": result_dir},
    }

    # 生成JSON报告
    json_report = os.path.join(result_dir, "analysis_report.json")
    with open(json_report, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2, ensure_ascii=False)

    # 生成HTML报告
    html_report = os.path.join(result_dir, "analysis_report.html")
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>甲基化分析报告</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1, h2, h3 {{ color: #333; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .info-box {{ background-color: #f9f9f9; padding: 15px; margin: 10px 0; border-radius: 5px; }}
        </style>
    </head>
    <body>
        <h1>甲基化分析报告</h1>
        <div class="info-box">
            <h2>分析信息</h2>
            <p>分析日期: {report_data['analysis_date']}</p>
            <p>临时目录: {report_data['directories']['temp']}</p>
            <p>结果目录: {report_data['directories']['result']}</p>
        </div>
        <div class="info-box">
            <h2>样本信息</h2>
            <table>
                <tr>
                    <th>样本</th>
                    <th>文件路径</th>
                </tr>
                <tr>
                    <td>样本1甲基化数据</td>
                    <td>{report_data['sample_info']['sample1']}</td>
                </tr>
                <tr>
                    <td>样本2甲基化数据</td>
                    <td>{report_data['sample_info']['sample2']}</td>
                </tr>
                <tr>
                    <td>样本1候选区域数据</td>
                    <td>{report_data['sample_info']['candidate_sample1']}</td>
                </tr>
                <tr>
                    <td>样本2候选区域数据</td>
                    <td>{report_data['sample_info']['candidate_sample2']}</td>
                </tr>
            </table>
        </div>
        <div class="info-box">
            <h2>分析流程</h2>
            <ol>
                <li>数据预处理（质量控制和接头去除）</li>
                <li>甲基化序列比对</li>
                <li>PCR重复去除</li>
                <li>靶向区域提取</li>
                <li>提取甲基化信息</li>
                <li>过滤甲基化数据</li>
                <li>过滤到候选区域</li>
                <li>进行差异甲基化分析</li>
                <li>可视化结果</li>
            </ol>
        </div>
    </body>
    </html>
    """

    with open(html_report, "w", encoding="utf-8") as f:
        f.write(html_content)

    logger.info(f"分析报告生成完成: {html_report}")
    return html_report

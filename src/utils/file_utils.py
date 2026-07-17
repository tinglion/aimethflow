#!/usr/bin/env python3
"""
file_utils.py - 文件处理工具模块

包含处理FastQ文件的通用函数
"""

import gzip
from Bio import SeqIO


def read_fastq(file_path, start=0, count=None):
    """
    读取FastQ文件，支持指定开始位置和读取数量
    
    Args:
        file_path: FastQ文件路径
        start: 开始位置（从0开始）
        count: 读取数量，None表示读取所有
    
    Returns:
        list: 读取的序列记录列表
        int: 总序列数量
    """
    records = []
    total_records = 0
    
    with gzip.open(file_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            total_records += 1
            
            # 跳过直到开始位置
            if total_records <= start + 1:
                continue
            
            # 收集指定数量的序列
            records.append(record)
            if count is not None and len(records) >= count:
                break
    
    return records, total_records


def get_fastq_stats(file_path):
    """
    获取FastQ文件的统计信息
    
    Args:
        file_path: FastQ文件路径
    
    Returns:
        dict: 统计信息字典
    """
    total_records = 0
    total_bases = 0
    lengths = []
    quality_scores = []
    
    with gzip.open(file_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            total_records += 1
            length = len(record)
            total_bases += length
            lengths.append(length)
            
            # 计算平均质量分数
            if record.letter_annotations.get('phred_quality'):
                avg_quality = sum(record.letter_annotations['phred_quality']) / length
                quality_scores.append(avg_quality)
    
    stats = {
        'total_records': total_records,
        'total_bases': total_bases,
        'lengths': lengths,
        'quality_scores': quality_scores
    }
    
    if lengths:
        stats['min_length'] = min(lengths)
        stats['max_length'] = max(lengths)
        stats['avg_length'] = total_bases / total_records
    
    if quality_scores:
        stats['avg_quality'] = sum(quality_scores) / len(quality_scores)
    
    return stats

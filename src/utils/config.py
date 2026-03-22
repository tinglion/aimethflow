#!/usr/bin/env python3
"""
config.py - 配置文件加载模块

用于加载和解析TOML格式的配置文件
"""

import os
import toml


def load_config(config_file=None):
    """
    加载配置文件
    
    Args:
        config_file: 配置文件路径，如果为None则使用默认路径
    
    Returns:
        dict: 配置字典
    """
    if config_file is None:
        # 查找默认配置文件
        current_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        config_file = os.path.join(current_dir, 'config.toml')
    
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"配置文件不存在: {config_file}")
    
    with open(config_file, 'r', encoding='utf-8') as f:
        config = toml.load(f)
    
    return config


def get_data_dir(config):
    """
    获取数据目录
    
    Args:
        config: 配置字典
    
    Returns:
        dict: 数据目录字典
    """
    return {
        'root': config['data']['root_dir'],
        'tissue': config['data']['tissue_dir'],
        'plasma': config['data']['plasma_dir']
    }


def get_output_dirs(config, timestamp=None):
    """
    获取输出目录
    
    Args:
        config: 配置字典
        timestamp: 时间戳，如果为None则生成新的
    
    Returns:
        dict: 输出目录字典
    """
    import time
    
    if timestamp is None:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
    
    return {
        'temp': os.path.join(config['directories']['temp_dir'], timestamp),
        'result': os.path.join(config['directories']['result_dir'], timestamp),
        'log': os.path.join(config['directories']['log_dir'], timestamp)
    }


def create_output_dirs(output_dirs):
    """
    创建输出目录
    
    Args:
        output_dirs: 输出目录字典
    """
    for dir_path in output_dirs.values():
        os.makedirs(dir_path, exist_ok=True)

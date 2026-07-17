import logging
import os


# 创建一个函数来获取当前文件的相对路径
def get_relative_path(file_path):
    """
    获取文件的相对路径
    
    Args:
        file_path: 文件路径
    
    Returns:
        str: 相对于当前工作目录的路径
    """
    return os.path.relpath(file_path, os.getcwd())


# 由于Formatter没有内置的relative_path选项，我们需要通过一个小技巧来添加它
# 我们可以使用Filter来实现这一点
class RelativePathFilter(logging.Filter):
    """
    相对路径过滤器，用于在日志中添加相对路径信息
    """
    def filter(self, record):
        record.relative_path = get_relative_path(record.pathname)
        return True


# 全局日志记录器
logger = None


def setup_logger(log_file=None, console_level=logging.INFO, file_level=logging.DEBUG):
    """
    设置日志记录器
    
    Args:
        log_file: 日志文件路径，如果为None则不记录到文件
        console_level: 控制台日志级别
        file_level: 文件日志级别
    
    Returns:
        logging.Logger: 配置好的日志记录器
    """
    global logger
    
    # 如果已经设置过，直接返回
    if logger is not None:
        return logger
    
    # 创建日志记录器
    logger = logging.getLogger(__name__)
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(logging.DEBUG)
    
    # 添加过滤器到日志记录器
    logger.addFilter(RelativePathFilter())
    
    # 创建流处理器（console handler）
    ch = logging.StreamHandler()
    ch.setLevel(console_level)
    ch.setFormatter(
        logging.Formatter(
            "%(asctime)s - %(levelname)s - %(relative_path)s %(funcName)s:%(lineno)d - %(message)s"
        )
    )
    logger.addHandler(ch)
    
    # 如果指定了日志文件，添加文件处理器
    if log_file:
        # 确保日志目录存在
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
        
        fh = logging.FileHandler(log_file, encoding='utf-8')
        fh.setLevel(file_level)
        fh.setFormatter(
            logging.Formatter(
                "%(asctime)s - %(levelname)s - %(relative_path)s %(funcName)s:%(lineno)d - %(message)s"
            )
        )
        logger.addHandler(fh)
    
    return logger


# 默认设置一个基本的logger，只输出到控制台
setup_logger()

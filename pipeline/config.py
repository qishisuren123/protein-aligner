"""
配置加载工具
"""
import yaml
import os

def load_config(config_path=None):
    """加载配置文件"""
    if config_path is None:
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "configs", "default.yaml"
        )
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

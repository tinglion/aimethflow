import os
import sys
from utils.config import load_config
from scripts.modules.classifier_utils import _build_sample_date_map

config = load_config()
s = _build_sample_date_map(config)
print(s)
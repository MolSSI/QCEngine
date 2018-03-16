"""
A set of scripts to setup testing
"""

import os
os.environ["DQM_CONFIG_PATH"] = os.path.dirname(os.path.abspath(__file__))
os.environ["TMPDIR"] = "something_scratch"

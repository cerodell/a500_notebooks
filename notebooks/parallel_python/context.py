# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: all,-language_info,-toc,-latex_envs
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
# ---

# %%
"""
http://docs.python-guide.org/en/latest/writing/structure

import context

will add the directory that contains this file to the front of sys.path
This allows notebooks to work with packages stored in subdirectories of the
notebook folder

"""
import sys
import site
from pathlib import Path

# %%
this_dir = Path(__file__).resolve().parent
data_dir = this_dir.parent / 'data'
sys.path.insert(0, str(this_dir))
sep = "*" * 30
print(f"{sep}\ncontext imported. Front of path:\n{sys.path[0]}\n{sys.path[1]}\n{sep}\n")

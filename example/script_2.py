# ============================================================
# example/script_2.py  (Python)
# Purpose: minimal example of a runnable pipeline file that depends on
#          script_1 (see batch_submit.sh — script_2 is submitted with
#          --depend=afterok on script_1's job).
# Inputs:  none
# Outputs: none (prints to stdout/err captured by SLURM)
#
# batch_controller.sh runs this as:  python -u example/script_2.py \
#   <SLURM_JOB_ID> <SLURM_ARRAY_TASK_ID> <extra_args>
# so argv[1:] may carry SLURM identifiers — ignore them if unused.
#
# Standard preamble pulls project paths from project_strings.py:
#   work_dir = os.environ.get('work_dir')
#   sys.path.insert(0, work_dir + '/big_cy')
#   from project_strings import *
# (omitted here so the example runs even before project_strings is filled in)
# ============================================================

print("Hello world from script_2.py")

# A real script would, e.g.:
#   import pandas as pd
#   df = pd.read_stata(output / "some_result.dta")
#   ... do work ...
#   df.to_stata(output / "final.dta")

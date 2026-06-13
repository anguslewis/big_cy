#!/bin/bash
############################################################################
# batch_controller.sh — runs one pipeline file on Sherlock. Invoked by
# batch_submit.sh via sbatch:
#     batch_controller.sh <folder> <file> [<extra_arg>]
# Detects the file's extension and dispatches to the right runtime
# (.do -> Stata, .py -> Python, .sh -> bash, .sif -> Apptainer), then
# checks the log for errors so a failed step fails the SLURM job.
############################################################################

# Repo code root and log root. These two lines are the only repo-specific
# edits needed if you reuse this controller elsewhere.
base_folder="${HOME}/big_cy"
logs="${big_cy}/logs"

# Per-user parameters (paths, env names).
source ~/shell_profile.sh

# SLURM context
echo "SLURM_JOB_ID="$SLURM_JOB_ID
echo "SLURM_ARRAY_JOB_ID="$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID
echo "Build folder    = ${base_folder}"
echo "Build subfolder = ${1}"
echo "Build code file = ${2}"
echo "Additional args = ${3}"

echo "Beginning step: ${1} ${2} ${SLURM_ARRAY_TASK_ID}"

############################################################################
# Locate the code file and capture its extension
############################################################################
if [ -f "${base_folder}/${1}/${2}."* ]; then
    name_w_extension=`basename $(ls ${base_folder}/${1}/${2}.*)`
    extension="${name_w_extension##*.}"
elif [ -f "${big_cy}/${1}/${2}."* ]; then
    name_w_extension=`basename $(ls ${big_cy}/${1}/${2}.*)`
    extension="${name_w_extension##*.}"
else
    echo "No code file found at ${base_folder}/${1}/${2}.* or ${big_cy}/${1}/${2}.*"
    echo "Exiting..."
    exit 1
fi
echo "Extension       = ${extension}"

############################################################################
# Dispatch by extension
############################################################################

# Stata
if [ "$extension" == "do" ]; then
    source ~/shell_modules.sh stata
    umask 007
    stata-mp -b "${base_folder}/batch_controller.do" ${1} ${2} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}
    rm -f batch_controller.log
    if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
        logfile="${logs}/${1}/${2}_${SLURM_JOB_ID}.log"
    else
        logfile="${logs}/${1}/${2}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
    fi
    echo "Checking Stata log for errors at $logfile ..."
    if ! test -f "$logfile"; then
        echo "No Stata log file found at $logfile"; exit 1
    fi
    if egrep --before-context=2 --max-count=1 "^r\([0-9]+\);$" "$logfile"; then
        echo "Stata error found."; exit 1
    fi
    echo "No errors found in Stata log."
fi

# Python
if [ "$extension" == "py" ]; then
    eval "$(conda shell.bash hook)"
    # Activate this repo's conda environment (klrep model-solve env).
    conda activate big_cy_klrep_env
    umask 007
    python -u "${base_folder}/${1}/${2}.py" ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID} ${3}
    if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
        logfile="${logs}/${1}/${2}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
    else
        logfile="${logs}/${1}/${2}_${SLURM_JOB_ID}.err"
    fi
    echo "Checking Python error log at $logfile ..."
    if ! test -f "$logfile"; then
        echo "No error log file found at $logfile"; exit 1
    fi
    for error_str in Traceback SyntaxError; do
        if grep -q "${error_str}" "$logfile"; then
            echo "Python error in log: $logfile"
            grep -B 2 -m 1 "${error_str}" "$logfile"
            exit 1
        fi
    done
    echo "No Python errors found in error log."
fi

# Bash
if [ "$extension" == "sh" ]; then
    source "${base_folder}/${1}/${2}.sh"
fi

# Apptainer / Singularity image
if [ "$extension" == "sif" ]; then
    unset MODULEPATH MODULESHOME
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
    apptainer run --nv "${big_cy}/${1}/${2}.sif" ${3}
fi

echo "Finished Step: ${1} ${2}"
exit

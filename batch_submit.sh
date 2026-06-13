#!/bin/bash
############################################################################
# batch_submit.sh — submit big_cy jobs to Sherlock with dependencies.
#
# HOW THIS FILE WORKS
#   * One block per runnable code file. Each block is guarded by:
#         if [[ ${to_run} == *" <file> "* ]]; then ... fi
#     so a file is submitted only when its name appears in `to_run`.
#   * `to_run` is a SPACE-PADDED list naming which files to submit THIS
#     run. It accepts ANY SUBSET of the files defined below, in any
#     combination — only the named files fire, the rest are skipped.
#     (The leading/trailing spaces matter: the guard matches " name ".)
#   * DEPENDENCIES use SLURM `--depend=afterok:<job_id>`. A block captures
#     its `job_id`; a downstream block passes `--depend=afterok:<id>` so it
#     starts only after the upstream job succeeds. If the upstream file is
#     NOT in `to_run` this run, its depend-var stays empty and the
#     downstream job starts immediately. So you can submit any subset and
#     the dependency wiring still does the right thing.
#
# HOW TO SUBMIT (do not normally edit `to_run` in this file)
#   sherlock-agent-com big_cy submit "script_1 script_2"   # both, chained
#   sherlock-agent-com big_cy submit "script_1"            # just the first
#   The wrapper sets `to_run` for that one invocation via the environment.
############################################################################

# Project path + SLURM globals — sets $big_cy, $big_cy_code, $mailuser,
# $mailtype, $system_part, etc. (self-contained; no ~/shell_profile.sh).
source "$HOME/big_cy/project_variables.sh"

# Sanity check: the project-root variable must be set by project_variables.sh.
if [[ -z "$big_cy" ]]; then
   echo "Empty big_cy variable: Exiting..."
   exit 1
fi
echo "Submitting big_cy"

# `to_run`: which files to submit this run. NOTE: this `sherlock-agent-com`
# version runs `bash ./batch_submit.sh` with NO job argument, so `to_run` is the
# committed default below (space-padded; the guards match " name "). Edit it to
# pick what a `submit` fires.
#   Env built OK (job 29446396) -> drop setup_klrep_env_sherlock; just solve.
#   PHASE 1 (now): to_run=" run_klrep " with specs=1 (baseline).
#   PHASE 2: set specs=1-9 (above) and keep to_run=" run_klrep ".
#   (To rebuild the env, add ' setup_klrep_env_sherlock ' back.)
to_run="${to_run:- run_klrep }"

############################################################################
# Shared SLURM settings
############################################################################
shared_settings="--partition=${system_part} --nodes=1 --requeue --mail-user=${mailuser} --mail-type=${mailtype}"
shared_settings_gpu="--partition=maggiori --nodes=1 --requeue --mail-user=${mailuser} --mail-type=${mailtype}"

############################################################################
# klrep specifications to run as a GPU array job. Each spec <-> a SLURM array
# index <-> param_file_<index> (index 1 = benchmark). The run_klrep block
# below derives `--array=${specs}` from this. Override per submit via the env
# (e.g.  specs=1  for Phase-1 baseline;  specs=1-9  for the full set).
#   PHASE 1 (baseline): specs="1"
#   PHASE 2 (all 9):    specs="1-9"
specs="${specs:-1}"

############################################################################
# EXAMPLE pipeline (delete once you add real jobs).
#
# script_1.do  ->  script_2.py
# script_2 depends on script_1: it only starts after script_1 succeeds,
# but only if script_1 was actually submitted this run.
############################################################################

folder="example"
logs="${big_cy}/logs/${folder}"
mkdir -p "${logs}"

# --- Step 1: script_1.do (no upstream dependency) ------------------------
file="script_1"
depend_script_2=""                       # reset; set only if step 1 fires
if [[ ${to_run} == *" $file "* ]]; then
    job_id=`sbatch ${shared_settings} --time=0-0:10:00 --ntasks=1 --mem=2G \
                --job-name=$file --output="${logs}/${file}_%A.out" --error="${logs}/${file}_%A.err" \
                "${big_cy_code}/batch_controller.sh" $folder $file | awk '{print $NF}'`
    echo "Submitted $folder/$file: ${job_id}"
    sleep 1
    depend_script_2="--depend=afterok:${job_id}"   # script_2 waits on this
fi

# --- Step 2: script_2.py (depends on script_1 if it ran this run) --------
file="script_2"
if [[ ${to_run} == *" $file "* ]]; then
    job_id=`sbatch ${shared_settings} ${depend_script_2} --time=0-0:10:00 --ntasks=1 --mem=2G \
                --job-name=$file --output="${logs}/${file}_%A.out" --error="${logs}/${file}_%A.err" \
                "${big_cy_code}/batch_controller.sh" $folder $file | awk '{print $NF}'`
    echo "Submitted $folder/$file: ${job_id}"
    sleep 1
fi

############################################################################
# ADD REAL JOBS BELOW.
# Pattern for each file:
#   1. set  folder="<subdir>"  and  file="<basename without extension>"
#   2. ensure  logs="${big_cy}/logs/${folder}"; mkdir -p "${logs}"
#   3. guard with  if [[ ${to_run} == *" $file "* ]]; then ... fi
#   4. inside: sbatch <settings> [<--depend=afterok:...>] \
#        --job-name=$file --output/--error in ${logs} \
#        "${big_cy_code}/batch_controller.sh" $folder $file
#   5. capture job_id and expose a depend var for any downstream file.
# Use ${shared_settings_gpu} and add --gpus=1 for GPU jobs.
############################################################################

############################################################################
# klrep: Kekre-Lenel model solve (GPU array-per-specification).
#
#   setup_klrep_env_sherlock.sh  ->  run_klrep.py
#
# setup_klrep_env_sherlock builds the conda env big_cy_klrep_env on the cluster
# (one-time; include in to_run only when (re)building the env). run_klrep is a
# GPU ARRAY job: --array=${specs} maps each array index N -> param_file_N,
# solves that calibration on the full grid on GPU, validates vs the steady
# state, and writes outputs under $big_cy/data/output/klreplication/.
############################################################################

folder="klreplication"
logs="${big_cy}/logs/${folder}"
mkdir -p "${logs}"

# --- one-time conda env build (maggiori GPU node) ------------------------
# Run on the maggiori GPU nodes (sh04-06n05 H100/80G etc., ~1000G RAM) so the big
# conda pytorch+cuda solve has ample memory (16G OOM'd at job 29440076), AND so
# the build's `torch.cuda.is_available()` verification actually exercises the GPU
# (confirms the CUDA build + float64 on the H100 at build time).
file="setup_klrep_env_sherlock"
depend_run_klrep=""                       # reset; set only if env build fires
if [[ ${to_run} == *" $file "* ]]; then
    job_id=`sbatch ${shared_settings_gpu} --gpus=1 --time=0-2:00:00 --ntasks=1 --cpus-per-task=8 --mem=96G \
                --job-name=$file --output="${logs}/${file}_%A.out" --error="${logs}/${file}_%A.err" \
                "${big_cy_code}/batch_controller.sh" $folder $file | awk '{print $NF}'`
    echo "Submitted $folder/$file: ${job_id}"
    sleep 1
    depend_run_klrep="--depend=afterok:${job_id}"   # run_klrep waits on the env
fi

# --- GPU array solve: one array index per specification ------------------
# maggiori 80G GPU (--gpus=1); the node has ~1000G RAM so mem is generous.
file="run_klrep"
if [[ ${to_run} == *" $file "* ]]; then
    # ~7-10s/iter on A100; conv to 1e-8 can take ~1-1.5k dampened iters. 7h wall
    # (maggiori allows 7 days) so the per-spec solve never TIMEOUTs mid-convergence.
    job_id=`sbatch ${shared_settings_gpu} ${depend_run_klrep} --gpus=1 --array=${specs} \
                --time=0-07:00:00 --ntasks=1 --cpus-per-task=8 --mem=96G \
                --job-name=$file --output="${logs}/${file}_%A_%a.out" --error="${logs}/${file}_%A_%a.err" \
                "${big_cy_code}/batch_controller.sh" $folder $file | awk '{print $NF}'`
    echo "Submitted $folder/$file (array=${specs}): ${job_id}"
    sleep 1
fi

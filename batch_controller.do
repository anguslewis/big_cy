* ============================================================
* batch_controller.do — Stata entry point invoked by batch_controller.sh
* for every .do job. Sets globals, opens a per-job log, runs the file.
*
* Args (positional, from batch_controller.sh):
*   `1' = folder (subdir under the repo)
*   `2' = file   (basename, no extension)
*   `3' = SLURM_JOB_ID
*   `4' = SLURM_ARRAY_TASK_ID  (empty for non-array jobs)
* ============================================================

about
clear
set more off
global sherlock = 1
qui do ~/big_cy/project_globals.do

* Print arguments
di "Folder:   `1'"
di "File:     `2'"
di "Job ID:   `3'"
di "Array ID: `4'"

* Open log file (per-job, under $logs/<folder>/)
cap log close
cap mkdir "$logs/`1'"
if "`4'"=="" log using "$logs/`1'/`2'_`3'.log" , replace
if "`4'"!="" log using "$logs/`1'/`2'_`3'_`4'.log" , replace

* Run the requested file (array id passed through as its argument)
do "~/big_cy/`1'/`2'.do" `4'

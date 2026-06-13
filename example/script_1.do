* ============================================================
* example/script_1.do  (Stata)
* Purpose: minimal example of a runnable pipeline file.
* Inputs:  none
* Outputs: none (prints to log)
*
* This file is launched by batch_controller.sh (extension .do ->
* stata-mp -b batch_controller.do, which sets globals and opens a log).
* By the time this runs, $output / $temp / $raw / $logs are defined.
* ============================================================

display "Hello world from script_1.do"

* A real script would, e.g.:
*   use "$temp/some_input.dta", clear
*   ... do work ...
*   save "$output/some_result.dta", replace

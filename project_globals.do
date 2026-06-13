* project_globals.do — canonical path globals for big_cy (Stata side).
* Imported at the top of every .do file:  qui do $work_dir/big_cy/project_globals.do
*
* Remote (/oak) root uses Lab subdir "lewis"; local uses Dropbox with no subdir.

set more off , perm
set varabbrev off

****************************************************************************
* PATH GLOBALS
****************************************************************************

* $sherlock is set to 1 in batch_controller.do on the cluster; default local.
if "$sherlock" == "" global sherlock = 0

if $sherlock==1 {
    global gcap_lab  = "/oak/stanford/groups/maggiori/Lab"
    global gcap_data = "/oak/stanford/groups/maggiori/GCAP/data"
    global big_cy      = "$gcap_lab/lewis/big_cy"
    global big_cy_code = "~/big_cy"
}
else {
    global gcap_lab  = "/Volumes/maggiori/Lab"
    global gcap_data = "/Volumes/maggiori/GCAP/data"
    global big_cy      = "$droppath/big_cy"
    global big_cy_code = "$work_dir/big_cy"
}

* Standard project paths (do not change the names).
global raw    = "$big_cy/data/raw"
global temp   = "$big_cy/data/temp"
global output = "$big_cy/data/output"
global figtab = "$big_cy/figtab"
global logs   = "$big_cy/logs"

****************************************************************************
* PROJECT-SPECIFIC GLOBALS / ADO PATH / SCHEMES go below
****************************************************************************
* e.g. data-source locations, draft vintage, adopath additions, scheme imports.

#!/bin/bash
############################################################################
# project_variables.sh — big_cy project path + SLURM globals for Sherlock.
# Sourced by batch_submit.sh / batch_controller.sh (after ~/shell_profile.sh).
# Mirrors rier/project_variables.sh. Data root on oak matches the
# sherlock-agent-com remote root (/oak/.../Lab/lewis) + repo name.
############################################################################

############################################################################
# PATH GLOBALS
############################################################################

export big_cy="/oak/stanford/groups/maggiori/Lab/lewis/big_cy"
export raw="$big_cy/data/raw"
export temp="$big_cy/data/temp"
export output="$big_cy/data/output"
export figtab="$big_cy/figtab"

export big_cy_code="$HOME/big_cy"

export gcap_data="/oak/stanford/groups/maggiori/GCAP/data"

export system_part="maggiori,normal,owners,gsb"

export whoami="angusl"
export mailuser="angusl@stanford.edu"
export mailtype="ALL"

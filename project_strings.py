# project_strings.py — canonical path constants for big_cy (Python side).
#
# ONE-TIME SETUP (so `from project_strings import *` works anywhere):
#   SITEDIR=$(python -m site --user-site)
#   mkdir -p "$SITEDIR"
#   echo "$HOME/big_cy" > "$SITEDIR/big_cy.pth"   # (on Sherlock)
#
# Per-user path branches: Sherlock /oak uses Lab subdir "lewis"; local Dropbox
# has no subdir.

import getpass
import os
from pathlib import Path
import sys

username = getpass.getuser()
print("Username: ", username)

# --- Per-user roots -------------------------------------------------------
# Angus, Sherlock
if username == "angusl":
    big_cy = Path("/oak/stanford/groups/maggiori/Lab/lewis/big_cy")
    big_cy_code = Path("~/big_cy").expanduser()

# Angus, local (Dropbox). To use SMB-mounted Sherlock data instead, set
#   import builtins; builtins.use_sherlock_data = True   (before import)
import builtins
use_sherlock_data = getattr(builtins, "use_sherlock_data", False)

if username == "angus" and use_sherlock_data:
    big_cy = Path("/Volumes/maggiori/Lab/lewis/big_cy")
    big_cy_code = Path("~/code/big_cy").expanduser()
if username == "angus" and not use_sherlock_data:
    big_cy = Path("/Users/angus/CBS Dropbox/Angus Lewis/big_cy")
    big_cy_code = Path("~/code/big_cy").expanduser()

# Refuse to silently NameError downstream if the user isn't recognized.
if "big_cy" not in dir():
    raise RuntimeError(
        f"Unknown user {username!r}; configure paths in project_strings.py."
    )

# --- Standard project paths (do not change the names) ---------------------
raw = big_cy / "data/raw"
temp = big_cy / "data/temp"
output = big_cy / "data/output"
figtab = big_cy / "figtab"
logs = big_cy / "logs"

# Scratch-backed temp (fast I/O on Sherlock; same as temp locally).
if username == "angusl":
    temp_scr = Path(os.environ.get("SCRATCH", "/scratch/users/angusl")) / "big_cy/data/temp"
else:
    temp_scr = temp

# --- Project-specific data sources go below -------------------------------
# e.g. gcap_data, shared outputs, vintage strings, year lists, etc.

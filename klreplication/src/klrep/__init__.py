"""klrep — tensorized-Python replication of the Kekre-Lenel (2024 AER)
quantitative model (big_cy task Q0).

See ../../PLAN.md and ../../reference/ for the design and the Fortran/MATLAB maps.
Pure-numerical leaf modules do not import project_strings (they take tensors, not
paths); path-touching entry/IO scripts (run/, post/, params loader) do.
"""

# Checked Example 2 reference outputs

These outputs come from the corrected full-size runs validated on 2026-10-08. They are optional references; the scripts solve the problems afresh and do not load these files.

- `line_results.mat`, `shared_results.mat`: saved input data and both capacity cases.
- `line_diagnostics.txt`, `shared_diagnostics.txt`: parameters, initialization, convergence, conservation, and costs.
- `VALIDATION.txt`: the original validation record for the corrected package.
- `figures/line_graph.*`, `figures/simple_graph.*`: the two main comparisons.
- `figures/*_capacity_profiles.*`: accompanying capacity-profile plots.

New runs write to `../results/`, preserving these reference files. The validation record describes the earlier code-correction stage; the manuscript and repository documentation were edited subsequently. Reported timings are reference timings and are not controlled replacements for the Example 3 performance measurements.

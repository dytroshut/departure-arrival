# Repository organization and migration

Prepared from `dytroshut/departure-arrival` commit `6a35b0b7c0e9d25d621b3da26bafdabf6a64b74c` on 2026-10-08.

## File migration

The old root-level copies should be removed when the corresponding new folders are installed. Their contents are retained at the locations below.

| Old root file | New location | Action |
| --- | --- | --- |
| `coupled_DA_new_3d.m` | `examples/example1_single_node/coupled_DA_new_3d.m` | Move unchanged |
| `line_graph_revision.m` | `examples/example2_nodal_capacities/line_graph_revision.m` | Replace with corrected entry script; original retained in `legacy/example2/` |
| `simple_sinkhorn.m` | `examples/example2_nodal_capacities/simple_sinkhorn.m` | Replace with corrected entry script; original retained in `legacy/example2/` |
| `int_DAline_2d.m` | `legacy/example2/int_DAline_2d.m` | Archive different CVX experiment |
| `path_sinkhorn_2dplot.m` | `legacy/example2/path_sinkhorn_2dplot.m` | Archive earlier Sinkhorn experiment |
| `complex_graph_final_v10.m` | `examples/example3_large_network/complex_graph_final_v10.m` | Move unchanged |
| `da_table_final.m` | `examples/example3_large_network/da_table_final.m` | Move unchanged |
| `complex_graph.m` | `legacy/example3/complex_graph.m` | Archive earlier visualization |
| `comparision_table.m` | `legacy/example3/comparision_table.m` | Archive earlier comparison |
| `convergence.eps`, `convergence.svg`, `snapshot.eps`, `snapshot.svg`, `graph_topology.svg` | `examples/example3_large_network/reference/figures/` | Move unchanged |
| `README.md` | `README.md` | Replace the incomplete introduction with the repository guide |

The corrected Example 2 adds seven helpers under `examples/example2_nodal_capacities/helpers/`, a test entry point under `tests/`, checked reference outputs, and documentation. `.gitignore` excludes newly generated runs and local editor files. `.gitattributes` preserves exported figures and MATLAB data and identifies generated figure files for review.

## Applying the organized package

1. Add the `examples`, `tests`, `docs`, and `legacy` folders together with the new `README.md`, `.gitignore`, and `.gitattributes`.
2. Remove the old root files listed in the migration table after their new locations are present. Keeping both root and folder copies makes it easy to run the wrong version.
3. Run the small Example 2 verification command from the main README.
4. Commit the moves, corrected Example 2 code, and documentation together.

The ZIP is a complete organized working tree without Git metadata. It contains files to install into the existing repository; uploading the ZIP alone does not create the folder structure. The companion patch, when supplied, applies the same changes to the source commit without rewriting Git history.

## Scope and validation

This reorganization preserves all original numerical scripts: existing Example 1 and Example 3 entry scripts are relocated unchanged, and superseded alternatives are archived unchanged. For Example 2, the corrected numerical helpers are retained; folder-aware entry points and the test path are the only integration changes, together with the default output directory of `example2_compare`.

The supplied full-size Example 2 results were validated before reorganization. The small regression suite is rerun after helper relocation. See `organization_validation.txt` for the checks completed on this organized tree.

The existing Example 1 repository coverage is coupled-only. Its README records a pre-existing tensor-ordering concern; the independent DA figure scripts still need a separate final selection. This task does not silently substitute a different local Example 1 experiment.

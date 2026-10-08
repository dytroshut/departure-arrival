# Example 3 · Large network

The current repository scripts are retained with their numerical implementations unchanged:

| Script | Purpose | Requirement |
| --- | --- | --- |
| `complex_graph_final_v10.m` | Network schedules, topology, and convergence figures | MATLAB |
| `da_table_final.m` | Comparison and scalability experiments | MATLAB with Optimization Toolbox (`linprog`) |

From the repository root:

```matlab
run(fullfile('examples','example3_large_network','complex_graph_final_v10.m'));
```

For the benchmark experiment, inspect the switches and instance sizes near the top of the table script before running:

```matlab
run(fullfile('examples','example3_large_network','da_table_final.m'));
```

The benchmark may be lengthy. It writes `da_table_final_results.mat` and `da_table_final_results.csv` to its working directory; those generated files are ignored by Git. The visualization opens MATLAB figure windows.

The existing exported figures have moved to `reference/figures/`. They were carried over from the repository without regeneration. This organization task does not rerun the Example 3 benchmarks or certify new runtime measurements.

The earlier `complex_graph.m` and `comparision_table.m` are in [`legacy/example3/`](../../legacy/example3/). The old comparison uses different formulations/settings and is not the current table entry point.

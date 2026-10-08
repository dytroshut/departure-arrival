# Example 2 · Time-varying and shared nodal capacities

Two scripts generate the corrected comparisons:

| Script | Comparison |
| --- | --- |
| `line_graph_revision.m` | One line graph, six intermediate nodes, seven edges |
| `simple_sinkhorn.m` | Three selected paths with shared intermediate-node capacities |

Both use the functions in `helpers/`. Run the scripts from MATLAB or use the commands in the [main README](../../README.md#quick-start). Data, diagnostics, and figures are written to `results/`, which is ignored by Git. The supplied [reference outputs](reference/README.md) are kept separately.

## Model and parameters

- Time horizon `[0,1]`, `200` temporal grid points, and final `epsilon = 0.065`.
- Unit edge weights and cost `1/(t-s)` on the strict causal support `t > s`.
- Noncausal transitions are excluded. There is no row normalization of the Gibbs kernel.
- Original Gaussian-mixture boundary shapes are retained. Only cells that cannot belong to a complete causal path are removed, followed by renormalization. Deleted mass is recorded in the diagnostics.
- Marginals and capacities are cell masses, not continuous-time density values.
- The mean cell capacity is **`0.0105` for both line cases** and **`0.009` for both shared-network cases**.

For each intermediate node in display order, the varying profile is proportional to

```text
1 + 0.65*sin(2*pi*t/0.60 + (k-1)*pi/2),  k = 1,...,6,
```

and is rescaled to the stated mean. The original line mean `0.0100` is infeasible for the varying profiles under strict causality. Independent feasibility checks gave a minimum multiplier of approximately `1.01648061003`; the corrected mean uses a multiplier of `1.05`.

The shared paths are:

```text
p1: v0 -> v1 -> v3 -> v4 -> v5 -> vT
p2: v0 -> v2 -> v3 -> v4 -> v6 -> vT
p3: v0 -> v1 -> v3 -> v4 -> v6 -> vT
```

A node's capacity applies to the sum of all path marginals at that node.

## Solver and initialization

The solver uses exact cyclic Gauss–Seidel updates in log coordinates, with updated marginals after each block. Capacity multipliers can both activate and release bounds. No message floors, boundary-only polishing, or display-time offsets are used.

The default higher-epsilon initialization sequence is `(1.04, 0.26, 0.065)`. The more restrictive time-varying line case uses `(16.64, 4.16, 1.04, 0.26, 0.065)`, with gauge-aligned extrapolation between initialization phases. Before its final solve, a safeguarded dual Newton initializer refines the starting potentials. This initializer uses MATLAB only.

All reported plans come from the final exact cyclic solve at **`epsilon = 0.065`**. Earlier phases provide starting scalings. Their settings and iteration counts are recorded in the diagnostics.

For a fixed-epsilon solve from all-ones scalings, after adding `helpers/` to the MATLAB path:

```matlab
S = example2_data('line');
R = example2_solve(S, S.uniform, struct( ...
    'epsilonSchedule', S.epsilon, ...
    'useNewtonInitialization', false, ...
    'maxIter', 1000000));
```

This can take substantially more cycles. Final stopping checks departure/arrival L1 errors, total capacity excess, and a scaled complementarity residual at tolerance `1e-8`. Adjacent edge plans are reconstructed to check conservation, marginal consistency, and zero noncausal mass. A failed solve raises an error before exporting figures.

## Files

| Helper | Purpose |
| --- | --- |
| `example2_data.m` | Paths, boundaries, costs, capacities |
| `example2_compare.m` | Run the two capacity cases |
| `example2_solve.m` | Initialization and final solve |
| `example2_sinkhorn.m` | Log-domain cyclic solver and consistency checks |
| `example2_newton_initial.m` | Safeguarded dual initialization |
| `example2_save.m` | Checked data and diagnostic export |
| `example2_plot.m` | Figure generation |

The verification entry point is [`tests/verify_example2.m`](../../tests/verify_example2.m). The numerical algorithms are the same as the corrected package validated on 2026-10-08; only helper lookup and output locations changed during folder organization.

## Paper figures

Use `results/figures/line_graph.eps` and `results/figures/simple_graph.eps` after running the scripts, or the supplied checked versions in `reference/figures/`. The new EPS files are tightly cropped. Remove the old LaTeX `trim`/`clip` options when replacing the manuscript figures.

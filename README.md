# Departure–arrival constrained optimal transport

**Transport schedules on networks with flexible departure and arrival times, and capacity-limited intermediate nodes.**

MATLAB code accompanying *Temporally Flexible Transport Scheduling on Networks with Departure–Arrival Constraints and Nodal Capacity Limits* by Anqi Dong, Karl H. Johansson, and Johan Karlsson.

[Examples](#examples) · [Quick start](#quick-start) · [Requirements](#requirements) · [Validation](#validation) · [Repository guide](docs/repository_guide.md)

<p align="center">
  <img src="examples/example2_nodal_capacities/reference/figures/line_graph.png" width="48%" alt="Example 2: line graph with uniform and time-varying nodal capacities">
  <img src="examples/example2_nodal_capacities/reference/figures/simple_graph.png" width="48%" alt="Example 2: selected paths sharing nodal capacities">
</p>

<p align="center"><em>Example 2: capacity profiles shape when mass crosses each node; shared capacities aggregate flow across paths.</em></p>

## Examples

| Example | Contents | Start here |
| --- | --- | --- |
| **1 · Single intermediate node** | Existing coupled departure–arrival illustration | [Example 1](examples/example1_single_node/README.md) |
| **2 · Nodal capacities** | Corrected line-graph and shared-network comparisons with uniform and time-varying capacities | [Example 2](examples/example2_nodal_capacities/README.md) |
| **3 · Large network** | Network visualization, convergence plots, and comparison/scalability experiments | [Example 3](examples/example3_large_network/README.md) |

Example 1 currently contains the coupled illustration only. The complete independent/coupled figure set is not yet packaged here. The [legacy folder](legacy/README.md) preserves earlier experiments separately from the current entry points.

## Quick start

Download or clone the repository, then set MATLAB's Current Folder to the repository root.

Run the small Example 2 verification suite:

```matlab
addpath('tests');
verify_example2;
```

Generate the line-graph comparison:

```matlab
run(fullfile('examples','example2_nodal_capacities','line_graph_revision.m'));
```

Generate the shared-network comparison:

```matlab
run(fullfile('examples','example2_nodal_capacities','simple_sinkhorn.m'));
```

Each script compares uniform and time-varying capacities. New data, diagnostics, and figures go into `examples/example2_nodal_capacities/results/`. Checked outputs supplied with the repository are in the separate [reference folder](examples/example2_nodal_capacities/reference/README.md).

The full Example 2 solves use unit reciprocal edge costs and can take substantially longer than the small verification suite, particularly for the time-varying line case. Initialization settings are documented in the [example notes](examples/example2_nodal_capacities/README.md#solver-and-initialization).

## Requirements

| Component | Requirements |
| --- | --- |
| Example 2 and its tests | MATLAB; checked with R2025b Update 1. No CVX or Optimization Toolbox required. |
| Existing Example 1 illustration | MATLAB and an installed, initialized CVX environment. |
| Example 3 visualization | MATLAB. |
| Example 3 comparison tables | MATLAB with Optimization Toolbox (`linprog`). |

## Example 2 settings

| Setting | Line graph | Shared network |
| --- | --- | --- |
| Time horizon | `[0,1]` | `[0,1]` |
| Temporal grid points | `200` | `200` |
| Final regularization parameter | `0.065` | `0.065` |
| Mean cell capacity, both cases | **`0.0105`** | **`0.009`** |
| Edge cost | `1/(t-s)` for `t > s` | `1/(t-s)` for `t > s` |
| Intermediate nodes | Six, along one path | Six, shared by three selected paths |

All capacities above are **cell masses**. Noncausal transitions are excluded. The corrected line capacity replaces the earlier `0.0100` setting; its time-varying case was infeasible under strict causality. See the [Example 2 notes](examples/example2_nodal_capacities/README.md) for the sinusoidal profiles, support handling, and initialization.

## Validation

The corrected Example 2 package includes six small regression checks and checked outputs for all four full-size comparisons. Boundary errors and total capacity excess are at most `1e-8` in the supplied runs; edge-marginal consistency is checked separately. See the [validation record](examples/example2_nodal_capacities/reference/VALIDATION.txt) and per-case diagnostics.

The existing Example 1 and Example 3 numerical scripts are preserved during this reorganization. Example 2 verification does not validate those other experiments. Their own READMEs describe the available scripts and requirements.

## Repository layout

```text
examples/
  example1_single_node/       Existing coupled DA illustration
  example2_nodal_capacities/  Two entry scripts, helpers, checked references
  example3_large_network/    Visualization and benchmark scripts
tests/                       Small Example 2 verification suite
docs/                        File migration and repository guide
legacy/                      Earlier scripts retained for provenance
```

To find an old filename or understand which files replaced it, use the [migration table](docs/repository_guide.md#file-migration).

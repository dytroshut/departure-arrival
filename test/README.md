# Verification

From the repository root in MATLAB:

```matlab
addpath('tests');
verify_example2;
```

The test function locates the Example 2 helper folder relative to its own file. It requires MATLAB only and checks:

1. A capped single-node solution against the clipped Gibbs law.
2. Newton initialization followed by cyclic scaling against the same solution.
3. Shared-node marginals and costs against explicitly enumerated tensors.
4. Release of previously active capacity multipliers.
5. Active shared capacities against independent full-tensor scaling.
6. Rejection of the original infeasible varying line capacity.

These small checks do not run the full paper experiments. The [reference validation record](../examples/example2_nodal_capacities/reference/VALIDATION.txt) documents the full-size Example 2 runs.

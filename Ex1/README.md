# Example 1 · Single intermediate node

`coupled_DA_new_3d.m` is the existing coupled departure–arrival illustration from the repository. It solves a fixed joint endpoint problem using CVX and opens figures in MATLAB.

From the repository root, with CVX installed and initialized:

```matlab
run(fullfile('examples','example1_single_node','coupled_DA_new_3d.m'));
```

This file has been relocated without changing its numerical implementation. It has not been revalidated as part of the Example 2 correction. The independent DA script and the complete final Example 1 figure set are not included in this folder yet.

Before treating this existing script as the final paper reproducer, review its cost-tensor axis ordering: the script permutes the first two axes of `C`, while `Pi` is declared in departure–arrival order. This pre-existing issue has been recorded rather than silently changing a different experiment during repository organization.

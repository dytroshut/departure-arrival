%% baseline_comparison.m
%
%  Compares 5 methods on the transport scheduling problem.
%  Uses the SAME 100-node series-parallel DAG (10 stages × 10 nodes).
%  Varying: both P (paths) and T (time grid) to reveal curse of dimensionality.
%
%  Methods:
%    1. Direct LP          (CVX + MOSEK; marginal constraints, no capacity)
%    2. Network simplex    (MATLAB linprog dual-simplex; same LP)
%    3. Column generation  (Dantzig-Wolfe; starts with 3 seed paths)
%    4. Standard Sinkhorn  (entropic OT, aggregate kernel, Dykstra capacity)
%    5. Graphical Sinkhorn (entropic OT, path-wise, WITH nodal capacity)
%
%  NOTE: Methods 1-3 solve the UNCAPACITATED problem.  Methods 4-5 enforce
%  nodal capacity.
%
%  KEY SCALABILITY POINT:
%    LP  vars  ∝  P × T²   →  grows super-linearly; OOM for large instances.
%    CG  pricing ∝ P × T²  per iteration  →  slow for large P.
%    Graphical Sinkhorn ∝  L × P × T  per iteration  →  scales linearly.
%
%  Memory strategy:
%    - All paths share the same kernel K^L (homogeneous network) → stored ONCE
%      using MATLAB copy-on-write; no P × T² kernel copies.
%    - LP is built lazily per-method only when mem_est ≤ MEM_LIM.
%    - skip_mat tracks memory-skipped entries (shown as "---" in table)
%      vs genuine time-limit hits (shown as ">T s").

clear; clc;
rng(42);

%% ═══════════════════════════════════════════════════════════════════════════
%% Network topology
%% ═══════════════════════════════════════════════════════════════════════════
nStages = 10;  nWidth  = 10;
nInner  = nStages * nWidth;
nNodes  = nInner + 2;
src     = nInner + 1;
snk     = nInner + 2;
node_of = @(s,pp) (s-1)*nWidth + pp;

adj = cell(nNodes,1);
adj{src} = arrayfun(@(pp) node_of(1,pp), 1:nWidth);
for s = 1:nStages-1
    for pp = 1:nWidth
        adj{node_of(s,pp)} = node_of(s+1, randperm(nWidth,3));
    end
end
for pp = 1:nWidth,  adj{node_of(nStages,pp)} = snk;  end
adj{snk} = [];

%% ═══════════════════════════════════════════════════════════════════════════
%% Path pool  (up to 2000 unique random walks, generated once)
%% ═══════════════════════════════════════════════════════════════════════════
max_pool = 2000;
all_paths = {};
seen      = containers.Map('KeyType','char','ValueType','logical');
ntry      = 0;
while numel(all_paths) < max_pool && ntry < 2000000
    ntry = ntry+1;
    kk = src;  pth = kk;  ok = true;
    while kk ~= snk
        nb = adj{kk};
        if isempty(nb),  ok = false;  break;  end
        kk = nb(randi(numel(nb)));
        pth(end+1) = kk;  %#ok<AGROW>
    end
    if ok
        key = sprintf('%d,', pth);
        if ~isKey(seen, key)
            seen(key) = true;  all_paths{end+1} = pth;  %#ok<AGROW>
        end
    end
end
pool_size = numel(all_paths);
fprintf('Path pool: %d unique paths generated.\n', pool_size);

%% ═══════════════════════════════════════════════════════════════════════════
%% Fixed problem parameters
%% ═══════════════════════════════════════════════════════════════════════════
eps_ent  = 0.04;
dt_min   = 0.035;
cap_sc   = 0.005;   % nodal capacity (Methods 4 & 5)
TIME_LIM = 300;     % seconds – wall-clock limit per solve
%  Memory limit for LP Aeq construction (bytes).
%  Estimate: 2*nVars non-zeros × ~20 bytes each = nVars*40.
%  Set to 4 GB — lets (P=1000,T=500) try; skips (P=2000,T=600).
MEM_LIM  = 4e9;

sinkhorn_opts = struct('eps0',1e-8,'damp',0.7, ...
    'clip_min',1e-8,'clip_max',1e8,'maxIter',300,'tol',1e-3);

%% ═══════════════════════════════════════════════════════════════════════════
%% Problem instances
%% ═══════════════════════════════════════════════════════════════════════════
problem_instances = [
    100,  200;   % baseline
    200,  400;   % T scaling
    500,  400;   % P scaling
   1000,  500;   % large P: LP ~47M vars → very slow / timeout
   2000,  600;   % curse of dimensionality: LP ~137M vars → OOM
];
nInst    = size(problem_instances, 1);
nMethods = 5;

% Result matrices  (rows = methods, cols = instances)
rt_mat   = nan(nMethods, nInst);
mem_mat  = nan(nMethods, nInst);
gap_mat  = nan(nMethods, nInst);
obj_mat  = nan(nMethods, nInst);
iter_mat = nan(nMethods, nInst);
nvar_mat = nan(1, nInst);        % LP vars shared across methods 1-3
skip_mat = false(nMethods, nInst); % true = memory-skipped (show "---")

%% ═══════════════════════════════════════════════════════════════════════════
%% Main benchmark loop
%% ═══════════════════════════════════════════════════════════════════════════
for si = 1:nInst
    nPaths = min(problem_instances(si,1), pool_size);
    nt     = problem_instances(si,2);
    fprintf('\n\n══════ P=%d, T=%d ══════\n', nPaths, nt);

    paths   = all_paths(1:nPaths);
    middles = unique([paths{:}]);
    middles = middles(middles ~= src & middles ~= snk);

    % ── Time grid & marginals ──────────────────────────────────────────────
    t   = linspace(0,1,nt)';
    dt  = t(2)-t(1);
    G   = @(m,s) exp(-0.5*((t-m)/s).^2);
    mu0 = 0.9*G(0.10,0.05) + 0.6*G(0.22,0.05);  mu0 = mu0/sum(mu0);
    muT = 0.8*G(0.76,0.05) + 0.7*G(0.88,0.05);  muT = muT/sum(muT);

    % ── Edge kernel ───────────────────────────────────────────────────────
    [U,V]  = ndgrid(t,t);
    Kedge  = zeros(nt);
    kmask  = (V - U) >= dt_min;
    Kedge(kmask) = exp(-(V(kmask)-U(kmask)) / eps_ent);
    Cost   = max(0, V - U);

    % ── Path kernel: K_common = K_edge^L  (same for all paths) ────────────
    path_lengths = cellfun(@(p) numel(p)-1, paths);
    Lmax = max(path_lengths);
    Kpow = cell(Lmax+1,1);
    Kpow{1} = eye(nt);
    for ll = 1:Lmax
        Kpow{ll+1} = Kpow{ll} * Kedge;
    end
    % All paths have the same length in this homogeneous DAG.
    % Store ONE shared kernel; Kpath cell uses MATLAB copy-on-write → O(T²) memory.
    K_common = Kpow{Lmax+1};
    Kpath = repmat({K_common}, nPaths, 1);   % COW: no copies until write

    % ── Truncate marginals to kernel support ──────────────────────────────
    src_reachable  = sum(K_common, 2) > 0;
    snk_reachable  = sum(K_common, 1)' > 0;
    mu0 = mu0 .* src_reachable;  if sum(mu0)>0, mu0=mu0/sum(mu0); else, mu0=ones(nt,1)/nt; end
    muT = muT .* snk_reachable;  if sum(muT)>0, muT=muT/sum(muT); else, muT=ones(nt,1)/nt; end

    % ── LP variable count (lightweight — just count support) ──────────────
    supp_per_path = nnz(K_common > 1e-300);
    nVars         = nPaths * supp_per_path;
    mem_est_lp    = nVars * 40;   % bytes: 2*nVars nnz × ~20 bytes each
    nvar_mat(si)  = nVars;
    fprintf('  LP: %d variables (%.1fM), %.2f GB est  |  %d constraints\n', ...
        nVars, nVars/1e6, mem_est_lp/1e9, 2*nt);

    lp_obj_ref = NaN;
    lp_built   = false;
    Aeq = []; beq = []; c_lp = []; lb_lp = [];

    %% ── Method 1: Direct LP  (CVX + MOSEK) ──────────────────────────────
    fprintf('  [1] Direct LP (CVX+MOSEK) ...');
    if mem_est_lp > MEM_LIM
        skip_mat(1,si) = true;
        fprintf(' SKIP – mem %.1f GB > %.0f GB limit\n', mem_est_lp/1e9, MEM_LIM/1e9);
    else
        try
            if ~lp_built
                [Aeq,beq,c_lp,lb_lp,~,~] = build_lp_system(paths,Kpath,Cost,mu0,muT,nt);
                lp_built = true;
            end
            mb_before = get_mem_mb();
            tm = tic;
            cvx_begin quiet
                cvx_solver mosek
                variable x_lp(nVars)
                minimize( c_lp' * x_lp )
                subject to
                    Aeq * x_lp == beq
                    x_lp >= lb_lp
            cvx_end
            rt = toc(tm);
            mb_after = get_mem_mb();
            if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
                lp_obj_ref    = cvx_optval;
                rt_mat(1,si)  = rt;
                obj_mat(1,si) = lp_obj_ref;
                gap_mat(1,si) = 0;
                mem_mat(1,si) = max(mb_after-mb_before, 0);
                fprintf(' %.2f s  (obj=%.4f)\n', rt, lp_obj_ref);
            else
                fprintf(' CVX status: %s\n', cvx_status);
            end
        catch ME
            fprintf(' ERROR: %s\n', ME.message);
        end
    end

    %% ── Method 2: Network simplex  (linprog dual-simplex) ────────────────
    fprintf('  [2] Network simplex (linprog) ...');
    if mem_est_lp > MEM_LIM
        skip_mat(2,si) = true;
        fprintf(' SKIP – mem %.1f GB\n', mem_est_lp/1e9);
    else
        try
            if ~lp_built
                [Aeq,beq,c_lp,lb_lp,~,~] = build_lp_system(paths,Kpath,Cost,mu0,muT,nt);
                lp_built = true;
            end
            lp_opts = optimoptions('linprog','Algorithm','dual-simplex', ...
                'Display','off','MaxIterations',5e6,'MaxTime',TIME_LIM);
            mb_before = get_mem_mb();
            tm = tic;
            [~, obj_ns, flag] = linprog(c_lp,[],[],Aeq,beq,lb_lp,[],lp_opts);
            rt = toc(tm);
            mb_after = get_mem_mb();
            if flag == 1
                rt_mat(2,si)  = rt;
                obj_mat(2,si) = obj_ns;
                mem_mat(2,si) = max(mb_after-mb_before, 0);
                ref = lp_obj_ref; if isnan(ref), ref = obj_ns; end
                gap_mat(2,si) = abs(obj_ns-ref)/(abs(ref)+1e-16);
                fprintf(' %.2f s  (gap=%.2e)\n', rt, gap_mat(2,si));
            elseif flag == -9
                fprintf(' time limit (%.0f s)\n', TIME_LIM);
            else
                fprintf(' linprog flag %d\n', flag);
            end
        catch ME
            fprintf(' ERROR: %s\n', ME.message);
        end
    end

    %% ── Method 3: Column generation  (3-path seed, full pricing) ─────────
    %   Pricing scans ALL P inactive paths → O(P×T²) per CG iter.
    %   CG does NOT need the full LP (only builds small restricted LP).
    fprintf('  [3] Column generation (3-path seed) ...');
    try
        mb_before = get_mem_mb();
        tm = tic;
        [obj_cg, nCols, nCGiter, flag_cg] = col_gen_solve( ...
            paths, Kpath, Cost, mu0, muT, nt, TIME_LIM);
        rt = toc(tm);
        mb_after = get_mem_mb();
        if ~isnan(obj_cg)
            rt_mat(3,si)   = rt;
            obj_mat(3,si)  = obj_cg;
            iter_mat(3,si) = nCGiter;
            mem_mat(3,si)  = max(mb_after-mb_before, 0);
            ref = lp_obj_ref; if isnan(ref), ref = obj_cg; end
            gap_mat(3,si) = abs(obj_cg-ref)/(abs(ref)+1e-16);
            fprintf(' %.2f s  (%d cols, %d CG iters, gap=%.2e)\n', ...
                rt, nCols, nCGiter, gap_mat(3,si));
        elseif flag_cg == -1
            fprintf(' time limit\n');
        else
            fprintf(' did not converge\n');
        end
    catch ME
        fprintf(' ERROR: %s\n', ME.message);
    end

    %% ── Method 4: Standard Sinkhorn  (aggregate kernel + Dykstra) ────────
    %   K_agg = K_common (all paths same kernel).
    %   Dykstra big-projection loop: O(T³) per iter → slow for large T.
    fprintf('  [4] Standard Sinkhorn + Dykstra (cap) ...');
    try
        K_agg = K_common;       % all Kpath{pp} identical → K_agg = K_common
        r_agg = numel(middles) * cap_sc;

        mb_before = get_mem_mb();
        tm = tic;
        [~, obj_ss, nit_ss] = std_sinkhorn_cap( ...
            K_agg, Cost, mu0, muT, r_agg, nt, sinkhorn_opts, TIME_LIM);
        rt = toc(tm);
        mb_after = get_mem_mb();

        rt_mat(4,si)   = rt;
        obj_mat(4,si)  = obj_ss;
        iter_mat(4,si) = nit_ss;      % negative = time limit hit
        mem_mat(4,si)  = max(mb_after-mb_before, 0);
        ref = lp_obj_ref; if isnan(ref), ref = obj_ss; end
        gap_mat(4,si)  = abs(obj_ss-ref)/(abs(ref)+1e-16);
        if nit_ss < 0
            fprintf(' %.1f s  TIME LIMIT (%d iter)\n', rt, -nit_ss);
        else
            fprintf(' %.2f s  (%d iter, gap=%.2e)\n', rt, nit_ss, gap_mat(4,si));
        end
    catch ME
        fprintf(' ERROR: %s\n', ME.message);
    end

    %% ── Method 5: Graphical Sinkhorn  (path-wise, with capacity) ─────────
    fprintf('  [5] Graphical Sinkhorn (cap=%.3f) ...', cap_sc);
    try
        r_cap = containers.Map('KeyType','double','ValueType','any');
        for vtx = middles,  r_cap(vtx) = cap_sc * ones(nt,1);  end

        mb_before = get_mem_mb();
        tm = tic;
        [~, m0gs, mTgs, nit_gs] = runPathSinkhorn( ...
            Kedge, mu0, muT, paths, middles, nt, sinkhorn_opts, r_cap);
        rt = toc(tm);
        mb_after = get_mem_mb();

        obj_gs = (t'*mTgs - t'*m0gs) * dt;

        rt_mat(5,si)   = rt;
        obj_mat(5,si)  = obj_gs;
        iter_mat(5,si) = nit_gs;
        mem_mat(5,si)  = max(mb_after-mb_before, 0);
        ref = lp_obj_ref; if isnan(ref), ref = obj_gs; end
        gap_mat(5,si)  = abs(obj_gs-ref)/(abs(ref)+1e-16);
        fprintf(' %.2f s  (%d iter, gap=%.2e)\n', rt, nit_gs, gap_mat(5,si));
    catch ME
        fprintf(' ERROR: %s\n', ME.message);
    end
end

%% ═══════════════════════════════════════════════════════════════════════════
%% Console summary
%% ═══════════════════════════════════════════════════════════════════════════
method_labels = {'Direct LP (CVX+MOSEK)', 'Network simplex (linprog)', ...
                 'Column generation', 'Std Sinkhorn+Dykstra', 'Graphical Sinkhorn'};

W = 15;
sep = repmat('─', 1, 6 + nInst*W);
fprintf('\n\n%s\n', sep);
fprintf('%-26s', 'Method');
for si = 1:nInst
    fprintf('  P=%-4d T=%-4d', problem_instances(si,1), problem_instances(si,2));
end
fprintf('\n%s\n', sep);
for m = 1:nMethods
    fprintf('%-26s', method_labels{m});
    for si = 1:nInst
        rt = rt_mat(m,si);
        if isnan(rt) && skip_mat(m,si)
            fprintf('%*s', W, '---');
        elseif isnan(rt)
            fprintf('%*s', W, '>300s');
        else
            fprintf('%*.2fs', W-1, rt);
        end
    end
    fprintf('\n');
end
fprintf('%s\n', sep);
fprintf('%-26s', '#LP vars');
for si = 1:nInst
    nv = nvar_mat(si);
    if isnan(nv),         fprintf('%*s', W, '--');
    elseif nv >= 1e6,     fprintf('%*s', W, sprintf('%.1fM', nv/1e6));
    else,                 fprintf('%*s', W, sprintf('%.0fk', nv/1e3));
    end
end
fprintf('\n%s\n', sep);

%% ═══════════════════════════════════════════════════════════════════════════
%% LaTeX table
%% ═══════════════════════════════════════════════════════════════════════════
tex_labels = {'Direct LP (CVX/MOSEK)', 'Network simplex', 'Col.\ generation$^{\dagger}$', ...
              'Std.\ Sinkhorn + Dykstra', '\textbf{Graphical Sinkhorn (ours)}'};

fprintf('\n\n%%  ── LaTeX table ──\n');
fprintf('\\begin{tabular}{l%s}\n', repmat('r',1,nInst));
fprintf('\\toprule\n');
fprintf('& \\multicolumn{%d}{c}{$(P,\\,T)$}\\\\\n', nInst);
fprintf('\\cmidrule(lr){2-%d}\n', nInst+1);
fprintf('Method');
for si = 1:nInst
    fprintf(' & $(%d,\\;%d)$', problem_instances(si,1), problem_instances(si,2));
end
fprintf(' \\\\\n');
fprintf('\\#var');
for si = 1:nInst
    nv = nvar_mat(si);
    if isnan(nv),       fprintf(' & --');
    elseif nv >= 1e6,   fprintf(' & $%.1f$M', nv/1e6);
    else,               fprintf(' & $%.0f$k', nv/1e3);
    end
end
fprintf(' \\\\\n\\midrule\n');

for m = 1:nMethods
    fprintf('%s', tex_labels{m});
    for si = 1:nInst
        rt  = rt_mat(m,si);
        nit = iter_mat(m,si);
        if isnan(rt) && skip_mat(m,si)
            fprintf(' & ---');                            % memory limit
        elseif isnan(rt)
            fprintf(' & $>\\!%d$\\,s', TIME_LIM);        % time limit
        elseif ~isnan(nit) && nit < 0                    % Sinkhorn TL
            fprintf(' & $>\\!%d$\\,s', TIME_LIM);
        elseif m == 3 && ~isnan(nit)                     % CG: show iters
            fprintf(' & %.2f\\,s (%d\\,it)', rt, nit);
        elseif m >= 4 && ~isnan(nit)                     % Sinkhorn
            fprintf(' & %.2f\\,s (%d\\,it)', rt, nit);
        else
            fprintf(' & %.2f\\,s', rt);
        end
    end
    fprintf(' \\\\\n');
    if m == 3,  fprintf('\\midrule\n');  end  % separator before Sinkhorn block
end
fprintf('\\bottomrule\n\\end{tabular}\n');
fprintf(['%%\n%% LP vars \\propto P T^2;  GS \\propto L P T per iter.\n' ...
         '%% --- = memory limit (%.0f GB);  $>$%ds = time limit.\n' ...
         '%% ^{\\dagger} CG converges in 1 pricing iter (homogeneous network);\n' ...
         '%%   pricing oracle still costs O(P T^2).\n'], MEM_LIM/1e9, TIME_LIM);

%% ═══════════════════════════════════════════════════════════════════════════
%% Local functions
%% ═══════════════════════════════════════════════════════════════════════════

function [obj, nCols, nCGiter, flag] = col_gen_solve( ...
        paths, Kpath, Cost, mu0, muT, nt, time_lim)
%COL_GEN_SOLVE  Dantzig-Wolfe column generation (3-path seed).
%   Pricing scans ALL inactive paths → O(P×T²) per CG iteration.

    nP      = numel(paths);
    obj     = NaN;
    flag    = 0;
    nCGiter = 0;
    t0      = tic;

    n_seed  = min(3, nP);
    active  = 1:n_seed;
    obj_r   = NaN;

    lp_opts = optimoptions('linprog','Algorithm','dual-simplex', ...
        'Display','off','MaxIterations',3e6,'MaxTime',time_lim*0.8);

    for iter_cg = 1:nP
        nCGiter = iter_cg;

        % Restricted LP
        [Aeq_r, beq_r, c_r, lb_r, ~, ~] = build_lp_system( ...
            paths(active), Kpath(active), Cost, mu0, muT, nt);
        [~, obj_r_raw, flag_r, ~, lam] = linprog( ...
            c_r,[],[],Aeq_r,beq_r,lb_r,[],lp_opts);
        if isscalar(obj_r_raw) && ~isnan(obj_r_raw),  obj_r = obj_r_raw;  end

        if flag_r ~= 1
            if flag_r == -9,  flag = -1;  end
            break
        end

        % Pricing: scan all inactive paths  ← O(P × T²)
        lambda0 = lam.eqlin(1:nt);
        lambdaT = lam.eqlin(nt+1:2*nt);

        inactive = setdiff(1:nP, active);
        best_rc  = -1e-8;
        best_p   = -1;
        for pp = inactive
            [rs, cs] = find(Kpath{pp} > 1e-300);
            if isempty(rs),  continue;  end
            c_p = Cost(sub2ind([nt,nt], rs, cs));
            rc  = min(c_p - lambda0(rs) - lambdaT(cs));
            if rc < best_rc,  best_rc = rc;  best_p = pp;  end
        end

        if best_p < 0
            obj   = obj_r;
            nCols = numel(active);
            flag  = 1;
            return
        end
        active = [active, best_p];  %#ok<AGROW>
        if toc(t0) > time_lim,  flag = -1;  break;  end
    end

    if flag == 0 && ~isnan(obj_r)
        obj   = obj_r;
        nCols = numel(active);
        flag  = 1;
    end
    if isnan(obj),  nCols = numel(active);  end
end


function [Aeq, beq, c_obj, lb, nVars, var_info] = build_lp_system( ...
        paths, Kpath, Cost, mu0, muT, nt)
%BUILD_LP_SYSTEM  Sparse LP for path-flow OT (marginal constraints only).

    nP     = numel(paths);
    THRESH = 1e-300;

    supp  = cell(nP,1);
    nVpP  = zeros(nP,1);
    for pp = 1:nP
        [rs, cs] = find(Kpath{pp} > THRESH);
        supp{pp} = [rs, cs];
        nVpP(pp) = numel(rs);
    end
    nVars    = sum(nVpP);
    var_info = struct('nVarsPerPath', nVpP, 'supp', {supp});

    nEq    = 2*nt;
    beq    = [mu0; muT];
    nnz_est = 2 * nVars;
    I = zeros(nnz_est,1);  J = zeros(nnz_est,1);  V = zeros(nnz_est,1);
    c_obj = zeros(nVars,1);
    lb    = zeros(nVars,1);

    offset = 0;  ptr = 0;
    for pp = 1:nP
        S  = supp{pp};
        n  = size(S,1);
        rs = S(:,1);  cs = S(:,2);
        jj = offset + (1:n)';

        I(ptr+1:ptr+n) = rs;
        J(ptr+1:ptr+n) = jj;
        V(ptr+1:ptr+n) = 1;
        ptr = ptr + n;

        I(ptr+1:ptr+n) = nt + cs;
        J(ptr+1:ptr+n) = jj;
        V(ptr+1:ptr+n) = 1;
        ptr = ptr + n;

        c_obj(jj) = Cost(sub2ind([nt,nt], rs, cs));
        offset = offset + n;
    end
    Aeq = sparse(I(1:ptr), J(1:ptr), V(1:ptr), nEq, nVars);
end


function [Pi, obj, nIter] = std_sinkhorn_cap( ...
        K_agg, Cost, mu0, muT, r_agg, nt, opts, time_lim)
%STD_SINKHORN_CAP  Sinkhorn + sequential Dykstra capacity projections.
%   Big-projection τ-loop: O(nt³) per Sinkhorn iter — bottleneck for large T.
%   nIter < 0 signals time limit hit (|nIter| = iterations completed).

    if nargin < 8 || isempty(time_lim),  time_lim = Inf;  end
    eps0    = opts.eps0;
    maxIter = opts.maxIter;
    u = ones(nt,1);  v = ones(nt,1);
    nIter = maxIter;
    t0 = tic;

    for it = 1:maxIter
        v = muT ./ max(K_agg' * u, eps0);
        u = mu0 ./ max(K_agg  * v, eps0);

        Pi = u .* K_agg .* v';
        any_proj = false;
        for tau = 1:nt
            sub   = Pi(1:tau, tau:nt);
            m_tau = sum(sub(:));
            if m_tau > r_agg + 1e-12
                Pi(1:tau, tau:nt) = sub * (r_agg / m_tau);
                any_proj = true;
            end
        end

        if any_proj
            row_s = sum(Pi, 2) + eps0;
            col_s = sum(Pi, 1)' + eps0;
            u = u .* (mu0 ./ row_s);
            v = v .* (muT ./ col_s);
            Pi = u .* K_agg .* v';
        end

        res = norm(sum(Pi,2)-mu0,1) + norm(sum(Pi,1)'-muT,1);
        if res < 1e-6,  nIter = it;  break;  end

        if toc(t0) > time_lim,  nIter = -it;  break;  end
    end

    Pi  = u .* K_agg .* v';
    obj = sum(Cost(:) .* Pi(:));
end


function [ribbons, m0, mT, nIter] = runPathSinkhorn( ...
        Kedge, p, q, paths, middles, nt, opts, r_cap)
%RUNPATHSINKHORN  Path-wise Sinkhorn — batched BLAS implementation.

    eps0 = opts.eps0;  damp = opts.damp;
    cmin = opts.clip_min;  cmax = opts.clip_max;
    maxIter = opts.maxIter;
    if isfield(opts,'tol'), tol = opts.tol; else, tol = 1e-3; end
    KT = Kedge';  nP = numel(paths);

    L = numel(paths{1}) - 2;
    stage_node = zeros(L, nP, 'double');
    for pp = 1:nP,  stage_node(:,pp) = paths{pp}(2:end-1)';  end
    middles_v = middles(:)';

    max_node  = max(stage_node(:));
    w_mat     = ones(nt, max_node);
    r_cap_mat = zeros(nt, max_node);
    for vtx = middles_v,  r_cap_mat(:,vtx) = r_cap(vtx);  end
    mid_idx   = middles_v;

    uniq_at = cell(L,1);
    for k = 1:L,  uniq_at{k} = unique(stage_node(k,:));  end

    u = ones(nt,1);  v = ones(nt,1);
    nIter = maxIter;
    B = cell(L,1);  F = cell(L,1);

    for it = 1:maxIter
        W = cell(L,1);
        for k = 1:L,  W{k} = w_mat(:, stage_node(k,:));  end

        Kv = Kedge * v;
        B{L} = max(Kv(:, ones(1,nP)), eps0);
        for k = L-1:-1:1
            B{k} = max(Kedge * (W{k+1} .* B{k+1}), eps0);
        end

        KTu = KT * u;
        F{1} = max(KTu(:, ones(1,nP)), eps0);
        for k = 2:L
            F{k} = max(KT * (W{k-1} .* F{k-1}), eps0);
        end

        A0 = Kedge * (W{1} .* B{1});
        AT = KT    * (W{L} .* F{L});

        m0 = u .* max(sum(A0,2), eps0);
        mT = v .* max(sum(AT,2), eps0);

        mn_raw = zeros(nt, max_node);
        for k = 1:L
            FkBk = F{k} .* B{k};
            for vn = uniq_at{k}
                pmask = stage_node(k,:) == vn;
                mn_raw(:,vn) = mn_raw(:,vn) + sum(FkBk(:,pmask), 2);
            end
        end

        mv_mat  = w_mat(:,mid_idx) .* max(mn_raw(:,mid_idx), eps0);
        r_sub   = r_cap_mat(:,mid_idx);
        cap_res = sum(sum(max(mv_mat - r_sub, 0)));

        u = min(max(u .* (p  ./ max(m0,    eps0)).^damp, cmin), cmax);
        v = min(max(v .* (q  ./ max(mT,    eps0)).^damp, cmin), cmax);
        w_mat(:,mid_idx) = min(max( ...
            w_mat(:,mid_idx) .* min(1, r_sub ./ max(mv_mat,eps0)).^damp, ...
            cmin), cmax);

        bnd_res = sum(abs(m0-p)) + sum(abs(mT-q));
        if bnd_res < tol && cap_res < tol * numel(mid_idx) * nt
            nIter = it;  break
        end
    end

    % Light polish: 6 Sinkhorn steps with w frozen
    for kp = 1:6  %#ok<FORPF>
        for k = 1:L,  W{k} = w_mat(:, stage_node(k,:));  end
        Kv = Kedge*v;
        B{L} = max(Kv(:,ones(1,nP)), eps0);
        for k = L-1:-1:1,  B{k} = max(Kedge*(W{k+1}.*B{k+1}), eps0);  end
        A0 = Kedge*(W{1}.*B{1});
        u = u .* (p ./ max(u .* max(sum(A0,2),eps0), eps0));
        KTu = KT*u;
        F{1} = max(KTu(:,ones(1,nP)), eps0);
        for k = 2:L,  F{k} = max(KT*(W{k-1}.*F{k-1}), eps0);  end
        AT = KT*(W{L}.*F{L});
        v = v .* (q ./ max(v .* max(sum(AT,2),eps0), eps0));
    end
    m0 = u .* max(sum(A0,2), eps0);
    mT = v .* max(sum(AT,2), eps0);

    ribbons = containers.Map('KeyType','double','ValueType','any');
    for vtx = mid_idx
        rv = zeros(nt,1);
        for k = 1:L
            pmask = stage_node(k,:) == vtx;
            if any(pmask)
                rv = rv + w_mat(:,vtx) .* sum(F{k}(:,pmask).*B{k}(:,pmask), 2);
            end
        end
        ribbons(vtx) = rv;
    end
end


function mb = get_mem_mb()
%GET_MEM_MB  Cross-platform resident memory of current process (MB).
    try
        [~,sys] = memory();
        mb = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e6;
    catch
        try
            pid = feature('getpid');
            [~, out] = system(sprintf('ps -o rss= -p %d', pid));
            mb = str2double(strtrim(out)) / 1024;
        catch
            mb = NaN;
        end
    end
end





% %% baseline_comparison.m
% %
% %  Compares 5 methods on the transport scheduling problem.
% %  Uses the SAME 100-node series-parallel DAG (10 stages × 10 nodes) as complex_graph.m.
% %  Varying: both P (paths) and T (time grid) to reveal curse of dimensionality.
% %
% %  Methods:
% %    1. Direct LP          (CVX + MOSEK; marginal constraints, no capacity)
% %    2. Network simplex    (MATLAB linprog dual-simplex; same LP)
% %    3. Column generation  (Dantzig-Wolfe; starts with 3 seed paths)
% %    4. Standard Sinkhorn  (entropic OT, aggregate kernel, Dykstra capacity)
% %    5. Graphical Sinkhorn (entropic OT, path-wise, WITH nodal capacity)
% %
% %  NOTE: Methods 1-3 solve the UNCAPACITATED problem.  Methods 4-5 enforce
% %  nodal capacity constraints.
% %
% %  KEY SCALABILITY POINT:
% %    LP  vars  ∝  P × T²   →  grows super-linearly; skipped above VAR_LIM.
% %    CG  pricing ∝ P × T²  per iteration  →  slow for large P even if few
% %                           columns are added (pricing oracle dominates).
% %    Graphical Sinkhorn ∝  L × P × T  per iteration  →  scales linearly.
% %
% %  Outputs:
% %    - Console table: runtime (s), peak memory (MB), LP variable count,
% %                     iterations (Sinkhorn only), obj gap vs Direct LP
% %    - LaTeX table fragment printed to console
% %
% 
% clear; clc;
% rng(42);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Network topology  (identical to complex_graph.m)
% %% ═══════════════════════════════════════════════════════════════════════════
% nStages = 10;  nWidth  = 10;
% nInner  = nStages * nWidth;
% nNodes  = nInner + 2;
% src     = nInner + 1;
% snk     = nInner + 2;
% node_of = @(s,pp) (s-1)*nWidth + pp;
% 
% adj = cell(nNodes,1);
% adj{src} = arrayfun(@(pp) node_of(1,pp), 1:nWidth);
% for s = 1:nStages-1
%     for pp = 1:nWidth
%         adj{node_of(s,pp)} = node_of(s+1, randperm(nWidth,3));
%     end
% end
% for pp = 1:nWidth,  adj{node_of(nStages,pp)} = snk;  end
% adj{snk} = [];
% 
% % Generate up to 2000 unique random paths (pool used by all experiments)
% max_pool = 2000;
% all_paths = {};
% seen      = containers.Map('KeyType','char','ValueType','logical');
% ntry      = 0;
% while numel(all_paths) < max_pool && ntry < 2000000
%     ntry = ntry+1;
%     kk = src;  pth = kk;  ok = true;
%     while kk ~= snk
%         nb = adj{kk};
%         if isempty(nb),  ok = false;  break;  end
%         kk = nb(randi(numel(nb)));
%         pth(end+1) = kk;  %#ok<AGROW>
%     end
%     if ok
%         key = sprintf('%d,', pth);
%         if ~isKey(seen, key)
%             seen(key) = true;  all_paths{end+1} = pth;  %#ok<AGROW>
%         end
%     end
% end
% pool_size = numel(all_paths);
% fprintf('Path pool: %d unique paths generated.\n', pool_size);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Fixed problem parameters
% %% ═══════════════════════════════════════════════════════════════════════════
% eps_ent  = 0.04;
% dt_min   = 0.035;
% cap_sc   = 0.005;   % nodal capacity (Methods 3 & 4)
% TIME_LIM = 300;     % seconds – wall-clock limit per solve
% MEM_LIM  = 8e9;     % bytes  – skip LP if estimated memory exceeds this
% VAR_LIM  = 20e6;    % skip LP if variable count exceeds this (CVX overhead)
% 
% sinkhorn_opts = struct('eps0',1e-8,'damp',0.7, ...
%     'clip_min',1e-8,'clip_max',1e8,'maxIter',300,'tol',1e-3);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Problem instances  — vary BOTH nPaths (P) and nt (T) to reveal the
% %% curse of dimensionality:  LP variables ∝ P × T²,  Sinkhorn ∝ L × P × T
% %% ═══════════════════════════════════════════════════════════════════════════
% %   Each row: [nPaths, nt]   (P = paths,  T = time grid)
% problem_instances = [
%     100,  200;   % baseline
%     200,  400;   % T scaling: LP vars ~100× larger than baseline
%     500,  400;   % P scaling
%    1000,  500;   % large P: LP O(P·T²) = 80M vars → very slow
%    2000,  600;   % curse of dimensionality: LP ~500M vars → timeout/OOM
% ];
% nInst    = size(problem_instances, 1);
% nMethods = 5;
% 
% % Result matrices  (rows = methods, cols = instances)
% rt_mat   = nan(nMethods, nInst);
% mem_mat  = nan(nMethods, nInst);
% gap_mat  = nan(nMethods, nInst);
% obj_mat  = nan(nMethods, nInst);
% iter_mat = nan(nMethods, nInst);
% nvar_mat = nan(nMethods, nInst);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Main benchmark loop
% %% ═══════════════════════════════════════════════════════════════════════════
% for si = 1:nInst
%     nPaths = min(problem_instances(si,1), pool_size);
%     nt     = problem_instances(si,2);
%     fprintf('\n\n══════ P=%d, T=%d ══════\n', nPaths, nt);
% 
%     % Subset paths and recompute middles for this instance
%     paths   = all_paths(1:nPaths);
%     middles = unique([paths{:}]);
%     middles = middles(middles ~= src & middles ~= snk);
% 
%     % ── Time grid & marginals ──────────────────────────────────────────────
%     t   = linspace(0,1,nt)';
%     dt  = t(2)-t(1);
%     G   = @(m,s) exp(-0.5*((t-m)/s).^2);
%     % Peaks well-separated so min-travel-time (11*dt_min=0.385) leaves
%     % feasible mass at every support point — avoids LP infeasibility.
%     mu0 = 0.9*G(0.10,0.05) + 0.6*G(0.22,0.05);  mu0 = mu0/sum(mu0);
%     muT = 0.8*G(0.76,0.05) + 0.7*G(0.88,0.05);  muT = muT/sum(muT);
% 
%     % ── Kernel ────────────────────────────────────────────────────────────
%     [U,V]  = ndgrid(t,t);
%     Kedge  = zeros(nt);
%     kmask  = (V - U) >= dt_min;
%     Kedge(kmask) = exp(-(V(kmask)-U(kmask)) / eps_ent);
%     Cost   = max(0, V - U);              % C(s,t) = t-s  (normalised ∈[0,1])
% 
%     % ── Path kernels  K_path{p} = K_edge^{L_p} ────────────────────────────
%     path_lengths = cellfun(@(p) numel(p)-1, paths);   % recompute per instance
%     Lmax = max(path_lengths);
%     Kpow = cell(Lmax+1,1);
%     Kpow{1} = eye(nt);
%     for ll = 1:Lmax
%         Kpow{ll+1} = Kpow{ll} * Kedge;
%     end
%     Kpath = cell(nPaths,1);
%     for pp = 1:nPaths
%         Kpath{pp} = Kpow{path_lengths(pp)+1};
%     end
% 
%     % ── Truncate mu0/muT to kernel support (prevents LP infeasibility) ────
%     % Gaussian tails extend past the feasible time window (min travel 0.385).
%     % Any s with zero reachable t (or t with zero reachable s) must be zeroed.
%     K_agg_support = Kpath{1};  % all same length
%     src_reachable  = sum(K_agg_support, 2) > 0;   % departure indices with ≥1 feasible arrival
%     snk_reachable  = sum(K_agg_support, 1)' > 0;  % arrival  indices with ≥1 feasible departure
%     mu0 = mu0 .* src_reachable;  if sum(mu0)>0, mu0 = mu0/sum(mu0); else, mu0 = ones(nt,1)/nt; end
%     muT = muT .* snk_reachable;  if sum(muT)>0, muT = muT/sum(muT); else, muT = ones(nt,1)/nt; end
% 
%     % ── LP structure (shared by methods 1-2) ──────────────────────────────
%     [Aeq, beq, c_lp, lb_lp, nVars, var_info] = build_lp_system( ...
%         paths, Kpath, Cost, mu0, muT, nt);
%     fprintf('  LP: %d variables, %d constraints\n', nVars, size(Aeq,1));
%     nvar_mat([1 2 3], si) = nVars;
% 
%     lp_obj_ref = NaN;   % reference value for gap computation
% 
%     %% ── Method 1: Direct LP  (CVX + MOSEK) ──────────────────────────────
%     fprintf('  [1] Direct LP (CVX+MOSEK) ...');
%     if nVars * 8 > MEM_LIM || nVars > VAR_LIM
%         fprintf(' skipped – %d vars (limit %.0fM)\n', nVars, VAR_LIM/1e6);
%     else
%         try
%             mb_before = get_mem_mb();
%             tm = tic;
%             cvx_begin quiet
%                 cvx_solver mosek
%                 variable x_lp(nVars)
%                 minimize( c_lp' * x_lp )
%                 subject to
%                     Aeq * x_lp == beq
%                     x_lp >= lb_lp
%             cvx_end
%             rt = toc(tm);
%             mb_after = get_mem_mb();
% 
%             if strcmp(cvx_status,'Solved') || contains(cvx_status,'Inaccurate')
%                 lp_obj_ref  = cvx_optval;
%                 rt_mat(1,si)  = rt;
%                 obj_mat(1,si) = lp_obj_ref;
%                 gap_mat(1,si) = 0;
%                 mem_mat(1,si) = max(mb_after - mb_before, 0);
%                 fprintf(' %.2f s  (obj=%.4f)\n', rt, lp_obj_ref);
%             else
%                 fprintf(' CVX status: %s\n', cvx_status);
%             end
%         catch ME
%             fprintf(' ERROR: %s\n', ME.message);
%         end
%     end
% 
%     %% ── Method 2: Network simplex  (linprog dual-simplex) ────────────────
%     fprintf('  [2] Network simplex (linprog) ...');
%     if nVars * 8 > MEM_LIM || nVars > VAR_LIM
%         fprintf(' skipped – %d vars\n', nVars);
%     else
%         try
%             lp_opts = optimoptions('linprog', ...
%                 'Algorithm','dual-simplex','Display','off', ...
%                 'MaxIterations',5e6,'MaxTime',TIME_LIM);
%             mb_before = get_mem_mb();
%             tm = tic;
%             [~, obj_ns, flag] = linprog(c_lp,[],[],Aeq,beq,lb_lp,[],lp_opts);
%             rt = toc(tm);
%             mb_after = get_mem_mb();
%             if flag == 1
%                 rt_mat(2,si)  = rt;
%                 obj_mat(2,si) = obj_ns;
%                 mem_mat(2,si) = max(mb_after - mb_before, 0);
%                 ref = lp_obj_ref; if isnan(ref), ref = obj_ns; end
%                 gap_mat(2,si) = abs(obj_ns - ref) / (abs(ref)+1e-16);
%                 fprintf(' %.2f s  (gap=%.2e)\n', rt, gap_mat(2,si));
%             elseif flag == -9
%                 fprintf(' time limit reached\n');
%             else
%                 fprintf(' linprog flag %d\n', flag);
%             end
%         catch ME
%             fprintf(' ERROR: %s\n', ME.message);
%         end
%     end
% 
%     %% ── Method 3: Column generation  (Dantzig-Wolfe, 3-path seed) ───────────
%     %   Starts with 3 seed paths (not 15) so the pricing oracle runs at least
%     %   once over all P inactive paths.  Pricing cost = O(P × T²) per CG iter
%     %   — this dominates at large P and shows CG does not scale in P.
%     fprintf('  [3] Column generation (3-path seed) ...');
%     if nVars * 8 > MEM_LIM || nVars > VAR_LIM
%         fprintf(' skipped\n');
%     else
%         try
%             mb_before = get_mem_mb();
%             tm = tic;
%             [obj_cg, nCols, nCGiter, flag_cg] = col_gen_solve( ...
%                 paths, Kpath, Cost, mu0, muT, nt, TIME_LIM);
%             rt = toc(tm);
%             mb_after = get_mem_mb();
%             if ~isnan(obj_cg)
%                 rt_mat(3,si)  = rt;
%                 obj_mat(3,si) = obj_cg;
%                 iter_mat(3,si) = nCGiter;
%                 mem_mat(3,si) = max(mb_after - mb_before, 0);
%                 ref = lp_obj_ref; if isnan(ref), ref = obj_cg; end
%                 gap_mat(3,si) = abs(obj_cg - ref) / (abs(ref)+1e-16);
%                 fprintf(' %.2f s  (%d cols, %d CG iters, gap=%.2e)\n', ...
%                     rt, nCols, nCGiter, gap_mat(3,si));
%             elseif flag_cg == -1
%                 fprintf(' time limit\n');
%             else
%                 fprintf(' did not converge\n');
%             end
%         catch ME
%             fprintf(' ERROR: %s\n', ME.message);
%         end
%     end
% 
%     %% ── Method 4: Standard Sinkhorn  (aggregate kernel + Dykstra capacity) ─
%     %   Solves with CAPACITY via a "big projection" loop over all nt time
%     %   steps at each Sinkhorn iteration → O(nt^3) per iter, slow for large nt.
%     fprintf('  [4] Standard Sinkhorn + Dykstra (cap) ...');
%     % Note: K_agg is nt×nt (fine), but the Dykstra τ-loop is O(nt³) per iter.
%     % For large nt this is the dominant bottleneck demonstrating scalability.
%     try
%         % Aggregate kernel: uniform mixture of path kernels
%         K_agg = zeros(nt);
%         for pp = 1:nPaths,  K_agg = K_agg + Kpath{pp};  end
%         K_agg = K_agg / nPaths;
% 
%         % Aggregate transit-capacity = sum of all nodal capacities
%         r_agg = numel(middles) * cap_sc;
% 
%         mb_before = get_mem_mb();
%         tm = tic;
%         [~, obj_ss, nit_ss] = std_sinkhorn_cap( ...
%             K_agg, Cost, mu0, muT, r_agg, nt, sinkhorn_opts, TIME_LIM);
%         rt = toc(tm);
%         mb_after = get_mem_mb();
% 
%         rt_mat(4,si)   = rt;
%         obj_mat(4,si)  = obj_ss;
%         iter_mat(4,si) = nit_ss;   % negative = time limit hit
%         mem_mat(4,si)  = max(mb_after - mb_before, 0);
%         ref = lp_obj_ref; if isnan(ref), ref = obj_ss; end
%         gap_mat(4,si)  = abs(obj_ss - ref) / (abs(ref)+1e-16);
%         if nit_ss < 0
%             fprintf(' %.1f s  TIME LIMIT (%d iter completed)\n', rt, -nit_ss);
%         else
%             fprintf(' %.2f s  (%d iter, gap=%.2e)\n', rt, nit_ss, gap_mat(4,si));
%         end
%     catch ME
%         fprintf(' ERROR: %s\n', ME.message);
%     end
% 
%     %% ── Method 5: Graphical Sinkhorn  (path-wise, WITH capacity) ──────────
%     fprintf('  [5] Graphical Sinkhorn (cap=%.3f) ...', cap_sc);
%     try
%         r_cap = containers.Map('KeyType','double','ValueType','any');
%         for vtx = middles,  r_cap(vtx) = cap_sc * ones(nt,1);  end
% 
%         mb_before = get_mem_mb();
%         tm = tic;
%         [~, m0gs, mTgs, nit_gs] = runPathSinkhorn( ...
%             Kedge, mu0, muT, paths, middles, nt, sinkhorn_opts, r_cap);
%         rt = toc(tm);
%         mb_after = get_mem_mb();
% 
%         % Objective: E[T-S] = int t*mT dt - int s*m0 ds
%         obj_gs = (t'*mTgs - t'*m0gs) * dt;
% 
%         rt_mat(5,si)   = rt;
%         obj_mat(5,si)  = obj_gs;
%         iter_mat(5,si) = nit_gs;
%         mem_mat(5,si)  = max(mb_after - mb_before, 0);
%         ref = lp_obj_ref; if isnan(ref), ref = obj_gs; end
%         gap_mat(5,si)  = abs(obj_gs - ref) / (abs(ref)+1e-16);
%         fprintf(' %.2f s  (%d iter, gap=%.2e)\n', rt, nit_gs, gap_mat(5,si));
%     catch ME
%         fprintf(' ERROR: %s\n', ME.message);
%     end
% end
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Console summary
% %% ═══════════════════════════════════════════════════════════════════════════
% method_labels = {'Direct LP (CVX+MOSEK)', 'Network simplex (linprog)', ...
%                  'Column generation', 'Std Sinkhorn+Dykstra', 'Graphical Sinkhorn'};
% 
% W = 14;
% sep = repmat('─',1,6 + nInst*W);
% fprintf('\n\n%s\n', sep);
% fprintf('%-24s', 'Method');
% for si = 1:nInst
%     fprintf('  P=%-3d T=%-3d', problem_instances(si,1), problem_instances(si,2));
% end
% fprintf('\n%s\n', sep);
% for m = 1:nMethods
%     fprintf('%-24s', method_labels{m});
%     for si = 1:nInst
%         rt = rt_mat(m,si);
%         if isnan(rt), fprintf('%*s',W,'--'); else, fprintf('%*.2fs',W-1,rt); end
%     end
%     fprintf('\n');
% end
% fprintf('%s\n', sep);
% fprintf('%-24s', '#LP vars');
% for si = 1:nInst
%     nv = nvar_mat(1,si);
%     if isnan(nv)
%         fprintf('%*s',W,'--');
%     elseif nv >= 1e6
%         fprintf('%*s',W,sprintf('%.1fM',nv/1e6));
%     else
%         fprintf('%*s',W,sprintf('%.0fk',nv/1e3));
%     end
% end
% fprintf('\n%s\n', sep);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% LaTeX table — combined (P, T) columns, showing curse of dimensionality
% %%
% %%   Column header:  $P=30,\,T=50$  etc.
% %%   Each cell:      runtime (s) with gap% in parens where applicable
% %%   Extra row:      \#LP vars  shows O(P \times T^2) growth for LP methods
% %%   Caption note:   Sinkhorn methods include nodal capacity (r=0.005).
% %%                   LP methods solve the uncapacitated relaxation.
% %% ═══════════════════════════════════════════════════════════════════════════
% fprintf('\n\n%%  ── LaTeX table (combined P×T scaling) ──\n');
% 
% col_hdr = cell(1,nInst);
% for si = 1:nInst
%     col_hdr{si} = sprintf('$P{=}%d,\\,T{=}%d$', ...
%         problem_instances(si,1), problem_instances(si,2));
% end
% 
% tex_labels = {'Direct LP (CVX/MOSEK)', 'Network simplex', 'Column generation', ...
%               'Std.\ Sinkhorn + Dykstra', '\textbf{Graphical Sinkhorn (ours)}'};
% 
% fprintf('\\begin{tabular}{l%s}\n', repmat('c',1,nInst));
% fprintf('\\toprule\n');
% fprintf('Method');
% for si = 1:nInst,  fprintf(' & %s', col_hdr{si});  end
% fprintf(' \\\\\n');
% % Sub-header: LP variable count
% fprintf('\\# LP vars');
% for si = 1:nInst
%     nv = nvar_mat(1,si);
%     if isnan(nv), fprintf(' & --');
%     else
%         if nv >= 1e6, fprintf(' & $%.1f$M', nv/1e6);
%         else,         fprintf(' & $%.0f$k', nv/1e3);
%         end
%     end
% end
% fprintf(' \\\\\n\\midrule\n');
% 
% for m = 1:nMethods
%     fprintf('%s', tex_labels{m});
%     for si = 1:nInst
%         rt  = rt_mat(m,si);
%         nit = iter_mat(m,si);
%         if isnan(rt)
%             fprintf(' & $>$%d\\,s', TIME_LIM);    % skipped / OOM
%         elseif ~isnan(nit) && nit < 0             % Sinkhorn time limit
%             fprintf(' & $>$%d\\,s', TIME_LIM);
%         elseif rt >= TIME_LIM * 0.95              % LP near time limit
%             fprintf(' & $>$%d\\,s', TIME_LIM);
%         elseif m == 3 && ~isnan(nit)              % CG: show iters (= #CG iters)
%             fprintf(' & %.2f\\,s (%d\\,it)', rt, nit);
%         elseif m >= 4 && ~isnan(nit)              % Sinkhorn: show iters
%             fprintf(' & %.2f\\,s (%d\\,it)', rt, nit);
%         else
%             fprintf(' & %.2f\\,s', rt);
%         end
%     end
%     fprintf(' \\\\\n');
% end
% fprintf('\\bottomrule\n\\end{tabular}\n');
% fprintf(['%%\n%% Columns: P = number of paths,  T = time-grid size.\n' ...
%          '%% LP vars \\propto P \\times T^2;   ' ...
%          'Graphical Sinkhorn \\propto L \\times P \\times T  per iteration.\n' ...
%          '%% Entries marked $>$%.0fs hit the wall-clock time limit.\n'], TIME_LIM);
% 
% %% ═══════════════════════════════════════════════════════════════════════════
% %% Local functions
% %% ═══════════════════════════════════════════════════════════════════════════
% 
% function [obj, nCols, nCGiter, flag] = col_gen_solve( ...
%         paths, Kpath, Cost, mu0, muT, nt, time_lim)
% %COL_GEN_SOLVE  Dantzig-Wolfe column generation for the path-flow LP.
% %
% %   Starts with 3 seed paths (minimum to get a well-conditioned restricted LP).
% %   Each CG iteration:
% %     1. Solve restricted LP over active columns  →  get dual variables λ0, λT
% %     2. PRICING: for every inactive path, compute min reduced cost           ← O(P×T²)
% %     3. If min reduced cost < -ε, add that path and repeat
% %   The pricing step dominates for large P: it evaluates O(P) paths, each
% %   requiring O(T²) work (full kernel support scan) → O(P×T²) per CG iteration.
% %   This matches the LP variable-count scaling and explains why CG does not
% %   scale better than direct LP as P grows.
% %
% %   Returns:
% %     obj      = optimal objective (NaN if failed / timed out)
% %     nCols    = number of active columns at termination
% %     nCGiter  = number of CG iterations performed
% %     flag     = 1 (optimal), -1 (time limit), 0 (did not converge)
% 
%     nP       = numel(paths);
%     obj      = NaN;
%     flag     = 0;
%     nCGiter  = 0;
%     t0       = tic;
% 
%     % Seed: 3 paths (small enough to force real pricing; enough for feasibility)
%     n_seed = min(3, nP);
%     active = 1:n_seed;
%     obj_r  = NaN;
% 
%     lp_opts = optimoptions('linprog','Algorithm','dual-simplex', ...
%         'Display','off','MaxIterations',3e6,'MaxTime',time_lim*0.8);
% 
%     for iter_cg = 1:nP
% 
%         nCGiter = iter_cg;
% 
%         % ── 1. Solve restricted LP ────────────────────────────────────────
%         [Aeq_r, beq_r, c_r, lb_r, ~, ~] = build_lp_system( ...
%             paths(active), Kpath(active), Cost, mu0, muT, nt);
% 
%         [~, obj_r_raw, flag_r, ~, lam] = linprog( ...
%             c_r,[],[],Aeq_r,beq_r,lb_r,[],lp_opts);
%         if isscalar(obj_r_raw) && ~isnan(obj_r_raw)
%             obj_r = obj_r_raw;
%         end
% 
%         if flag_r ~= 1
%             if flag_r == -9,  flag = -1;  end
%             break
%         end
% 
%         % ── 2. Pricing: scan ALL inactive paths  ← O(P × T²) ─────────────
%         lambda0 = lam.eqlin(1:nt);
%         lambdaT = lam.eqlin(nt+1:2*nt);
% 
%         inactive = setdiff(1:nP, active);
%         best_rc  = -1e-8;
%         best_p   = -1;
% 
%         for pp = inactive
%             [rs, cs] = find(Kpath{pp} > 1e-300);
%             if isempty(rs),  continue;  end
%             c_p = Cost(sub2ind([nt,nt], rs, cs));
%             rc  = min(c_p - lambda0(rs) - lambdaT(cs));
%             if rc < best_rc
%                 best_rc = rc;  best_p = pp;
%             end
%         end
% 
%         % ── 3. Termination or column addition ────────────────────────────
%         if best_p < 0
%             obj   = obj_r;
%             nCols = numel(active);
%             flag  = 1;
%             return
%         end
% 
%         active = [active, best_p];  %#ok<AGROW>
% 
%         if toc(t0) > time_lim
%             flag = -1;  break
%         end
%     end
% 
%     % Fell through (all paths added) — use last LP value
%     if flag == 0 && ~isnan(obj_r)
%         obj   = obj_r;
%         nCols = numel(active);
%         flag  = 1;
%     end
%     if isnan(obj), nCols = numel(active); end
% end
% 
% function [Aeq, beq, c_obj, lb, nVars, var_info] = build_lp_system( ...
%         paths, Kpath, Cost, mu0, muT, nt)
% %BUILD_LP_SYSTEM  Construct sparse equality-constrained LP for path-flow OT.
% %   Variables: π_p(s,t) for each path p and feasible (s,t) pair.
% %   Constraints: source marginal  Σ_p Σ_t π_p(s,t) = mu0(s)  ∀s
% %                sink   marginal  Σ_p Σ_s π_p(s,t) = muT(t)  ∀t
% %   No capacity constraints (LP baselines solve the uncapacitated problem).
% 
%     nP = numel(paths);
%     THRESH = 1e-300;
% 
%     % ── Count feasible variables and collect their (path, s, t) indices ──
%     supp  = cell(nP,1);
%     nVpP  = zeros(nP,1);
%     for pp = 1:nP
%         [rs, cs] = find(Kpath{pp} > THRESH);
%         supp{pp} = [rs, cs];
%         nVpP(pp) = numel(rs);
%     end
%     nVars    = sum(nVpP);
%     var_info = struct('nVarsPerPath', nVpP, 'supp', {supp});
% 
%     % ── Sparse Aeq assembly ───────────────────────────────────────────────
%     %   Rows 1..nt   : source marginal (sum over t and paths)
%     %   Rows nt+1..2nt: sink   marginal (sum over s and paths)
%     nEq    = 2*nt;
%     beq    = [mu0; muT];
% 
%     nnz_est = 2 * nVars;
%     I = zeros(nnz_est,1);  J = zeros(nnz_est,1);  V = zeros(nnz_est,1);
%     c_obj = zeros(nVars,1);
%     lb    = zeros(nVars,1);
% 
%     offset = 0;  ptr = 0;
%     for pp = 1:nP
%         S  = supp{pp};
%         n  = size(S,1);
%         rs = S(:,1);   % departure time indices
%         cs = S(:,2);   % arrival   time indices
%         jj = offset + (1:n)';
% 
%         % Source marginal rows
%         I(ptr+1:ptr+n) = rs;
%         J(ptr+1:ptr+n) = jj;
%         V(ptr+1:ptr+n) = 1;
%         ptr = ptr + n;
% 
%         % Sink marginal rows
%         I(ptr+1:ptr+n) = nt + cs;
%         J(ptr+1:ptr+n) = jj;
%         V(ptr+1:ptr+n) = 1;
%         ptr = ptr + n;
% 
%         % Objective: c(s,t) = t - s  (normalised)
%         c_obj(jj) = Cost(sub2ind([nt,nt], rs, cs));
% 
%         offset = offset + n;
%     end
%     Aeq = sparse(I(1:ptr), J(1:ptr), V(1:ptr), nEq, nVars);
% end
% 
% 
% function [Pi, obj, nIter] = std_sinkhorn_cap(K_agg, Cost, mu0, muT, r_agg, nt, opts, time_lim)
% %STD_SINKHORN_CAP  Sinkhorn + sequential Dykstra capacity projections.
% %
% %   At each iteration:
% %     1. Standard Sinkhorn-Knopp update (marginals)
% %     2. "Big projection": loop over τ = 1..nt and enforce that the total
% %        mass in transit at time τ does not exceed r_agg.  Each projection
% %        step re-slices the nt×nt plan matrix → O(nt^3) per iteration.
% %        This is the bottleneck that makes Standard Sinkhorn slow for
% %        large nt compared with the path-structured Graphical Sinkhorn.
% %
% %   r_agg    = total aggregate transit capacity (= n_nodes × r_cap_per_node).
% %   time_lim = (optional) wall-clock time limit in seconds; nIter = -it on exit.
% 
%     if nargin < 8 || isempty(time_lim), time_lim = Inf; end
%     eps0    = opts.eps0;
%     maxIter = opts.maxIter;
%     u = ones(nt,1);
%     v = ones(nt,1);
%     nIter = maxIter;
%     t0 = tic;
% 
%     for it = 1:maxIter
%         % ── 1. Sinkhorn-Knopp marginal update ────────────────────────────
%         v = muT ./ max(K_agg' * u, eps0);
%         u = mu0 ./ max(K_agg  * v, eps0);
% 
%         % ── 2. "Big projection": enforce transit capacity at each τ ──────
%         %   Cost: the τ-loop runs nt times; each body computes a triangular
%         %   submatrix sum and a conditional rescale → O(nt^2) per τ.
%         Pi = u .* K_agg .* v';
%         any_proj = false;
%         for tau = 1:nt
%             % mass in transit = entries where departure ≤ τ ≤ arrival
%             sub = Pi(1:tau, tau:nt);
%             m_tau = sum(sub(:));
%             if m_tau > r_agg + 1e-12
%                 Pi(1:tau, tau:nt) = sub * (r_agg / m_tau);
%                 any_proj = true;
%             end
%         end
% 
%         % ── 3. If projection fired, re-extract u,v from projected plan ───
%         if any_proj
%             row_s = sum(Pi, 2) + eps0;
%             col_s = sum(Pi, 1)' + eps0;
%             u = u .* (mu0 ./ row_s);
%             v = v .* (muT ./ col_s);
%             Pi = u .* K_agg .* v';
%         end
% 
%         % ── 4. Convergence on marginal residuals ─────────────────────────
%         res = norm(sum(Pi,2) - mu0, 1) + norm(sum(Pi,1)' - muT, 1);
%         if res < 1e-6
%             nIter = it;  break
%         end
%         % ── 5. Wall-clock time limit ──────────────────────────────────────
%         if toc(t0) > time_lim
%             nIter = -it;   % negative signals time limit hit
%             break
%         end
%     end
% 
%     Pi  = u .* K_agg .* v';
%     obj = sum(Cost(:) .* Pi(:));
% end
% 
% 
% function [ribbons, m0, mT, nIter] = runPathSinkhorn( ...
%         Kedge, p, q, paths, middles, nt, opts, r_cap)
% %RUNPATHSINKHORN  Path-wise Sinkhorn — BATCHED matrix-matrix implementation.
% %
% %   Replaces the per-path for-loop with batched nt×(nt×nP) products so
% %   MATLAB's BLAS handles the heavy lifting.  For paths of equal length L:
% %     B{k} = Kedge * (W{k+1} .* B{k+1})    nt×nt  *  nt×nP  each stage
% %     F{k} = KT    * (W{k-1} .* F{k-1})    same
% %   Node marginals are accumulated via per-stage unique-node scatter.
% %   Complexity per iteration: O(L × nt² × nP) — same as before but ×100
% %   faster because BLAS vs. interpreted MATLAB loops.
% 
%     eps0 = opts.eps0;  damp = opts.damp;
%     cmin = opts.clip_min;  cmax = opts.clip_max;
%     maxIter = opts.maxIter;
%     if isfield(opts,'tol'), tol = opts.tol; else, tol = 1e-3; end
%     KT = Kedge';  nP = numel(paths);
% 
%     % ── Pre-compute stage_node: L × nP ───────────────────────────────────
%     L = numel(paths{1}) - 2;   % interior nodes per path (all same here)
%     stage_node = zeros(L, nP, 'double');
%     for pp = 1:nP
%         stage_node(:,pp) = paths{pp}(2:end-1)';
%     end
%     middles_v = middles(:)';   % row vector for iteration
% 
%     % ── Dense weight + capacity matrices (avoid Map overhead) ─────────────
%     max_node = max(stage_node(:));
%     w_mat     = ones(nt, max_node);           % current capacity multipliers
%     r_cap_mat = zeros(nt, max_node);          % per-node capacity limits
%     for vtx = middles_v
%         r_cap_mat(:,vtx) = r_cap(vtx);
%     end
%     mid_idx = middles_v;   % column indices into w_mat for middles
% 
%     % ── Per-stage unique-node list (precomputed for scatter) ──────────────
%     uniq_at = cell(L,1);
%     for k = 1:L
%         uniq_at{k} = unique(stage_node(k,:));
%     end
% 
%     u = ones(nt,1);  v = ones(nt,1);
%     nIter = maxIter;
%     B = cell(L,1);  F = cell(L,1);
% 
%     for it = 1:maxIter
% 
%         % ── Build W{k}: weight at each stage, nt × nP  (fast column index)
%         W = cell(L,1);
%         for k = 1:L
%             W{k} = w_mat(:, stage_node(k,:));
%         end
% 
%         % ── Batched backward pass ─────────────────────────────────────────
%         Kv = Kedge * v;
%         B{L} = max(Kv(:, ones(1,nP)), eps0);   % same Kv for all paths
%         for k = L-1:-1:1
%             B{k} = max(Kedge * (W{k+1} .* B{k+1}), eps0);  % nt×nt * nt×nP
%         end
% 
%         % ── Batched forward pass ──────────────────────────────────────────
%         KTu = KT * u;
%         F{1} = max(KTu(:, ones(1,nP)), eps0);
%         for k = 2:L
%             F{k} = max(KT * (W{k-1} .* F{k-1}), eps0);
%         end
% 
%         % ── A0, AT (one matmul each) ──────────────────────────────────────
%         A0 = Kedge * (W{1} .* B{1});   % nt×nP
%         AT = KT    * (W{L} .* F{L});
% 
%         % ── Marginals ─────────────────────────────────────────────────────
%         m0 = u .* max(sum(A0,2), eps0);
%         mT = v .* max(sum(AT,2), eps0);
% 
%         % ── Scatter node marginals: mn_raw(:,vtx) = Σ_{k,p} F{k}[:,p].*B{k}[:,p]
%         mn_raw = zeros(nt, max_node);
%         for k = 1:L
%             FkBk = F{k} .* B{k};   % nt × nP
%             for vn = uniq_at{k}
%                 pmask = stage_node(k,:) == vn;
%                 mn_raw(:,vn) = mn_raw(:,vn) + sum(FkBk(:,pmask), 2);
%             end
%         end
% 
%         % ── Capacity residual & updates (fully vectorised over middles) ───
%         mv_mat  = w_mat(:,mid_idx) .* max(mn_raw(:,mid_idx), eps0);
%         r_sub   = r_cap_mat(:,mid_idx);
%         cap_res = sum(sum(max(mv_mat - r_sub, 0)));
% 
%         u = min(max(u .* (p  ./ max(m0,    eps0)).^damp, cmin), cmax);
%         v = min(max(v .* (q  ./ max(mT,    eps0)).^damp, cmin), cmax);
%         w_mat(:,mid_idx) = min(max( ...
%             w_mat(:,mid_idx) .* min(1, r_sub ./ max(mv_mat,eps0)).^damp, ...
%             cmin), cmax);
% 
%         bnd_res = sum(abs(m0-p)) + sum(abs(mT-q));
%         if bnd_res < tol && cap_res < tol * numel(mid_idx) * nt
%             nIter = it;  break
%         end
%     end
% 
%     % ── Light polish: 6 Sinkhorn steps with w frozen ─────────────────────
%     for kp = 1:6  %#ok<FORPF>
%         for k = 1:L,  W{k} = w_mat(:, stage_node(k,:));  end
%         Kv = Kedge*v;
%         B{L} = max(Kv(:,ones(1,nP)), eps0);
%         for k = L-1:-1:1,  B{k} = max(Kedge*(W{k+1}.*B{k+1}), eps0);  end
%         A0 = Kedge*(W{1}.*B{1});
%         u = u .* (p ./ max(u .* max(sum(A0,2),eps0), eps0));
%         KTu = KT*u;
%         F{1} = max(KTu(:,ones(1,nP)), eps0);
%         for k = 2:L,  F{k} = max(KT*(W{k-1}.*F{k-1}), eps0);  end
%         AT = KT*(W{L}.*F{L});
%         v = v .* (q ./ max(v .* max(sum(AT,2),eps0), eps0));
%     end
%     m0 = u .* max(sum(A0,2), eps0);
%     mT = v .* max(sum(AT,2), eps0);
% 
%     % ── Ribbon marginals (per-node accumulated flow) ──────────────────────
%     ribbons = containers.Map('KeyType','double','ValueType','any');
%     for vtx = mid_idx
%         rv = zeros(nt,1);
%         for k = 1:L
%             pmask = stage_node(k,:) == vtx;
%             if any(pmask)
%                 rv = rv + w_mat(:,vtx) .* sum(F{k}(:,pmask).*B{k}(:,pmask), 2);
%             end
%         end
%         ribbons(vtx) = rv;
%     end
% end
% 
% 
% function mb = get_mem_mb()
% %GET_MEM_MB  Cross-platform resident memory of current process (MB).
% %   Uses MATLAB's memory() on Windows; falls back to 'ps' on macOS/Linux.
%     try
%         [~,sys] = memory();
%         mb = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e6;
%     catch
%         try
%             pid = feature('getpid');
%             [~, out] = system(sprintf('ps -o rss= -p %d', pid));
%             mb = str2double(strtrim(out)) / 1024;   % KB → MB
%         catch
%             mb = NaN;
%         end
%     end
% end




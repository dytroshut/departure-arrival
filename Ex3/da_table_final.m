%% da_table_final.m — final comparison and scalability tables
%
%  Companion to "Temporally Flexible Transport Scheduling on Networks with
%  Departure--Arrival Constraints and Nodal Capacity Limits".
%
%  ══════════════════════════════════════════════════════════════════════════
%  METHOD LIST
%  ══════════════════════════════════════════════════════════════════════════
%   R1  LP reference           path-restricted LP, industrial solver
%                              (linprog barrier / CVX+MOSEK).  NOT part of the
%                              textbook comparison -- it exists to give exact
%                              ground truth so the Sinkhorn bias is measurable.
%
%   M1  Time-expanded MILP     NOT SOLVED, and it does not need to be: the
%                              node-split time-expanded flow matrix is a
%                              node-arc incidence matrix with arc bounds,
%                              hence totally unimodular, so with integer
%                              supplies and capacities the LP relaxation is
%                              already integral (Ahuja-Magnanti-Orlin, Thm
%                              11.x / TU).  This row RUNS THAT CHECK and
%                              reports the largest deviation of any arc flow
%                              from an integer.  Reviewer answer = a
%                              certificate, not a benchmark.
%
%   M2  Network simplex        TEXTBOOK primal network simplex (AMO Ch. 11) on
%                              the node-split time-expanded graph.  Hand
%                              written: artificial-root basis, potentials from
%                              the tree, Dantzig entering rule with full arc
%                              scan, cycle by common ancestor, ratio test.
%                              Solves the AGGREGATED problem (all DAG paths),
%                              so it is a LOWER BOUND on Problem 3.
%
%   M3  Column generation      TEXTBOOK two-phase Dantzig-Wolfe on the
%                              path-restricted problem.  Columns are
%                              (path, causal time tuple); pricing is an exact
%                              O(L n_t^2) min-plus DP per path.  Hand written
%                              apart from the restricted master LP, which is
%                              what Dantzig-Wolfe is defined to delegate.
%
%   M4  Standard MM Sinkhorn   needs the full n_t^(L+1) path tensor -> not
%                              storable at any size in this table.
%
%   M5  Graphical Sinkhorn     OURS.  Algorithm 1: per-edge kernels, cyclic
%                              block updates, shared nodal multipliers,
%                              forward-backward message passing on each path's
%                              CHAIN-structured cost.  Structure exploitation
%                              is allowed here and only here -- that is the
%                              contribution.
%                              NB the interior-node blocks are updated ONE
%                              STAGE AT A TIME.  All-nodes-at-once is a Jacobi
%                              step across stages and oscillates as soon as
%                              the capacity binds; stage-wise converges at
%                              every capacity tested and in ~10x fewer cycles.
%
%  IMPLEMENTATION POLICY (state this in the paper):
%    baselines are textbook pseudocode with no tuning; the proposed method
%    exploits the chain structure of the path cost.  All in MATLAB, same
%    machine, same stopping rule, same wall-clock limit.  Commercial
%    network-simplex implementations are considerably faster than a textbook
%    one; the comparison isolates algorithmic structure, not engineering.
%    For that reason OPERATION COUNTS (pivots / rounds / cycles) are reported
%    beside wall-clock -- they are language independent.
%
%  ══════════════════════════════════════════════════════════════════════════
%  TWO FORMULATIONS, AND WHY BOTH APPEAR
%  ══════════════════════════════════════════════════════════════════════════
%  path-restricted (Problem 3): R1, M3, M5, M6.  Mass is assigned to P
%      prescribed paths and cannot switch route at a shared node.
%  aggregated flow:             M1, M2.  Any path in the DAG is allowed, so it
%      is a RELAXATION -- a valid lower bound, not the same problem.  Its size
%      does not depend on P at all.
%  The gap between them is the PRICE OF THE PATH RESTRICTION and is reported
%  explicitly.  Measured on a small instance it decreases with P and saturates
%  near 12%, i.e. the flow model is not a substitute for Problem 3.  That is
%  the answer to "why not just use min-cost flow", and it is a modelling
%  argument, not a computational one.
%
%  ══════════════════════════════════════════════════════════════════════════
%  VALIDATION ALREADY DONE (independent NumPy/HiGHS reimplementation)
%  ══════════════════════════════════════════════════════════════════════════
%   * network simplex vs HiGHS, 12 random capacitated min-cost flows:
%     exact agreement (rel. diff 0.00e+00), conservation and bounds checked
%   * network simplex on the node-split time-expanded graph:
%     1.147968 vs HiGHS 1.147968, rel. diff 1.9e-16, 272 pivots / 240 nodes
%   * two-phase Dantzig-Wolfe reproduces the exact LP optimum to 6.5e-16
%   * graphical Sinkhorn -> exact LP optimum from above as eps decreases:
%     +6.97% at eps=0.08, +4.07% at 0.04, +2.43% at 0.02
%   * MILP integrality: 218890 arcs, integer data, ZERO fractional flows
%   * a big-M phase 1 for DW destroys the pricing signal (duals saturate at
%     +/-M and the scheme cycles) -- hence the two-phase version
%
%  ══════════════════════════════════════════════════════════════════════════

clear; clc;
rng(42);

%% ═══════════════════════════════════════════════════════════════════════════
%% Topology (10 stages x 10 width, 3 random successors)
%% ═══════════════════════════════════════════════════════════════════════════
nStages = 10;  nWidth = 10;
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
L = nStages + 1;

%% ── Edge list and heterogeneous edge weights ─────────────────────────────
eu = zeros(1,0);  ev = zeros(1,0);
for u = 1:nNodes
    for v = adj{u},  eu(end+1) = u;  ev(end+1) = v;  end  %#ok<AGROW>
end
nE  = numel(eu);
eid = sparse(eu, ev, 1:nE, nNodes, nNodes);

nWclass = 8;                       % weights take one of nWclass values
w_ref   = 1/(L*L);
w_vals  = w_ref * linspace(0.6, 1.4, nWclass);
e_class = randi(nWclass, 1, nE);

%% ═══════════════════════════════════════════════════════════════════════════
%% Path pool  (fixed, generated BEFORE any solve; every method gets the same
%% prefix all_paths(1:P) -- no method ever sees another's solution)
%% ═══════════════════════════════════════════════════════════════════════════
max_pool  = 2000;
all_paths = {};
seen = containers.Map('KeyType','char','ValueType','logical');
ntry = 0;
while numel(all_paths) < max_pool && ntry < 2000000
    ntry = ntry + 1;
    kk = src;  pth = kk;
    while kk ~= snk
        nb = adj{kk};  kk = nb(randi(numel(nb)));  pth(end+1) = kk;  %#ok<AGROW>
    end
    key = sprintf('%d,', pth);
    if ~isKey(seen, key),  seen(key) = true;  all_paths{end+1} = pth;  end  %#ok<AGROW>
end
pool_size = numel(all_paths);
fprintf('Path pool: %d unique paths, L = %d edges each.\n', pool_size, L);

%% ═══════════════════════════════════════════════════════════════════════════
%% Parameters
%% ═══════════════════════════════════════════════════════════════════════════
eps_ent  = 0.01;%0.04;      % entropic regularisation
dt_min   = 0.035;     % Delta, minimum edge traversal time
TIME_LIM = 300;       % wall-clock limit, applied to EVERY method incl. ours
TOL      = 1e-6;      % shared stopping rule max{E_0,E_T,V} <= TOL
MAXIT    = 20000;     % iteration cap for the Sinkhorn methods
DAMP     = 1.0;       % Algorithm 1 has no damping
TOL_LAD  = [1e-2 1e-3 1e-4 1e-5 1e-6];  % time-to-accuracy ladder
%  ACCURACY/TIME FRONTIER (off by default).  Re-runs OUR method at several
%  eps and reports (time, error vs ground truth), so the claim can be stated
%  as "X times faster AT A GIVEN ACCURACY" rather than as a bare speed
%  number.  That is the direct answer to "your method is N% off".  Turn it on
%  when writing the response letter.
%  OURS_ONLY: skip every baseline and time only the proposed method.  Use it
%  on instances no baseline can reach, so the scalability claim rests on
%  MEASURED times rather than on projections in a comparison table.  A table
%  whose baseline cells are estimates invites the objection that they are not
%  measurements; a sentence saying "our method solves (2000,600) in X s, and
%  the largest instance any baseline completes is (200,50)" does not.
RUN_EPS_SWEEP = false;
%  The exact reference: at (10,50) the barrier LP took 758 s and column
%  generation reached the SAME optimum (1.256156, agreeing to 6 digits) in
%  10.6 s.  Both are exact, so CG is the cheaper ground truth.  With this flag
%  the industrial LP is run only where CG fails to converge, which saves most
%  of the wall clock without losing the reference.
LP_ONLY_IF_CG_FAILS = true;
%  Column generation carries a RIGOROUS Dantzig-Wolfe bound, so even when it
%  is cut off it brackets the optimum.  At (20,100) it returned obj 1.245318
%  with a 0.22% gap, i.e. the optimum is in [1.2426, 1.2453].  The quantity we
%  need it for is a ~10% entropic bias, so a 0.22% uncertainty in the
%  reference is irrelevant -- and chasing it with a barrier LP cost over an
%  hour on a 1.02M-variable problem.  Accept the bracket instead.
CG_GAP_OK  = 0.01;    % CG counts as the reference if its DW gap is <= this
LP_VAR_LIM = 4e5;     % and never attempt the industrial LP above this size:
                      % MATLAB honours MaxTime only for simplex, so a barrier
                      % solve cannot be time-capped and will run unbounded.
EPS_SWEEP = [0.08 0.04 0.02 0.01];

%  Capacity is a RATE.  A fixed per-cell value is not equally binding across
%  instances (it came out 8-11x slack, so the constraint was idle and our
%  method converged in a constant ~15 cycles at every size).  Calibrate it
%  per instance to a fixed fraction of the unconstrained peak nodal rate.
CAP_MODE = 'relative';   % 'relative' | 'rate'
%  KAPPA controls how hard the capacity binds.  Measured sweep at exactly the
%  instance that failed (P=10, n_t=50, eps=0.04, tol=1e-4):
%      kappa  cap        cycles  converged  active cells
%      1.00   0.086256       14  yes          1/3050   <- constraint IDLE
%      0.85   0.073318       15  yes         19/3050
%      0.70   0.060379      281  yes         38/3050
%      0.50   0.043128    20000  NO           136/3050  <- with the OLD all-at-once update
%  With the stage-wise update every capacity converges:
%      0.70 -> 29 cycles (38 active), 0.50 -> 38 (52), 0.35 -> 56 (83),
%      0.25 -> 60 cycles (114 active).
%  KAPPA=0.70 turned out to be nearly decorative in the real run: at (50,100)
%  it left 0/9500 cells active with peak/cap 0.906, and the capacity changed
%  the optimum by only +0.14% at (100,200).  0.35 binds properly and is
%  measured convergent.
%  0.70 binds genuinely (38 active cells) and converges in 281 cycles, i.e.
%  real work that scales with problem size.  0.50 is past the edge of what
%  this scheme handles at eps=0.04 and produced the E0=1.7 / NaN failure.
KAPPA    = 0.35;         % cap = KAPPA * unconstrained peak   ('relative')
RBAR     = 2.0;          % cap = RBAR * dt                    ('rate')

MEM_LIM  = 4e9;       % skip path LP construction above this estimated size
CVX_LIM  = 3e7;
FEAS_LP  = 3e5;       % exact feasibility LP only below this variable count.
                      % Above it the check came back 'inconclusive (flag 0)'
                      % after burning its whole budget -- the per-stage bound
                      % plus the V-plateau test carry it from there.
AGG_FLOP = 2e10;      % skip aggregate Sinkhorn above this setup cost
CG_ADD   = 50;        % columns added per CG round
CG_LAD   = [0.10 0.01 0.001];   % report WHEN CG reaches 10% / 1% / 0.1% gap
NS_ARCS  = 6e6;       % skip network simplex above this arc count
NS_MAXPIV= 500000;
%  Network simplex gets a LARGER budget than the other methods on purpose.
%  A textbook implementation on ~350k arcs needs tens of thousands of pivots
%  and a bare ">TIME_LIM" tells the reader nothing.  Being generous to the
%  baseline is the right direction to err; the true time is reported.
%  Calibrated from the measured run: 316.88 s / 28222 pivots at 346040 arcs
%  = 11.23 ms per pivot, and pivots ~ 2.8 x nodes.  So
%      est_seconds ~ 2.8 * nodes * 3.245e-8 * arcs
%  which reproduces the measured 5.3 min at n_t=50 exactly.  Projections:
%      n_t=100  ->  42 min      n_t=200  ->  5.6 h
%      n_t=400  ->  44 h        n_t=600  ->  149 h
%  Rather than burn NS_TIME producing nothing, PROJECT first and skip with the
%  projection printed.  NS_TIME=3000 lets n_t=100 finish, which gives a second
%  point and hence a scaling claim; drop it to 600 for one point and a much
%  shorter run.
NS_TIME  = 3000;   % largest planned NS solve is ~34 min at n_t=100
NS_PIV_PER_NODE = 2.8;          % measured
NS_SEC_PER_PIV_ARC = 3.245e-8;  % measured
MILP_SCALE = 1000;    % integer scaling for the integrality certificate
%  The certificate establishes a THEOREM (total unimodularity => integral LP
%  relaxation), so one instance settles it.  It already returned
%  max|flow-int| = 0.00e+00 over 346040 arcs.  Running it at every n_t costs
%  ~15 min and adds nothing, so do it once.
MILP_ONCE  = true;

%% ═══════════════════════════════════════════════════════════════════════════
%% Instances: five small (exact LP solvable => real comparison + ground truth)
%% then two large for the scalability claim.
%% ═══════════════════════════════════════════════════════════════════════════
%  ── TABLE A: EVERY CELL A MEASURED NUMBER ────────────────────────────────
%  What limits a COMPLETE table is the slowest baseline, and the two baselines
%  are limited by different things:
%
%    network simplex  costs ~ nodes x arcs, and BOTH are independent of P
%                     (nodes = 2*nI*n_t, arcs = nE*cells + nI*n_t).  So it is
%                     identical at P=5 and P=2000 -- it caps n_t, never P.
%                        n_t= 40 -> 2.2 min      n_t= 50 -> 4.4 min
%                        n_t=100 -> 34 min       n_t=200 -> 5.6 h
%    column generation costs ~ P x cells x rounds, i.e. it caps P.
%                        measured n_t=50: P=10 -> 8.9 s, P=20 -> 15.7 s
%                        so P=200 -> ~177 s;  at n_t=100, P=50 -> 249 s
%
%  Hence a staircase rather than a rectangle: larger P where n_t is small.
%  P spans 10 -> 200 (20x) and n_t spans 40 -> 100 (2.5x), ~55 min total, and
%  the largest entry is a 2.6M-variable path LP over a 346k-arc flow model --
%  not a toy.
problem_instances = [ 10,  40
                      50,  40
                     200,  40
                      10,  50
                      50,  50
                     200,  50
                      10, 100
                      50, 100 ];

%  ── TABLE B, IN THE SAME RUN ─────────────────────────────────────────────
%  Appended below the main list and run with the baselines skipped, so ONE
%  clean run produces both tables and nothing has to be stitched by hand.
extra_ours_only = [100, 200
                   500, 400
                  2000, 600];

%  ── (was: a separate second run) ─────────────────────────────────────────
%  These carry the scalability claim, which Table A structurally cannot: a
%  head-to-head table can only extend to where the slowest baseline finishes.
%  Swap this in for a second run; the baselines report calibrated projections
%  and our method reports measured time.
%  (superseded by extra_ours_only above)

%  NB both tables are needed and they answer different questions.  Table A:
%  "is it correct, and how much faster where a fair comparison exists".
%  Table B: "does it run where nothing else does".  Neither alone is an
%  argument; a baseline being unable to run IS the scalability result, it just
%  belongs in a different table from the one that compares numbers.
%  concatenate: the appended rows are timed for OUR method only
ours_only_row = [false(size(problem_instances,1),1); true(size(extra_ours_only,1),1)];
problem_instances = [problem_instances; extra_ours_only];
nInst = size(problem_instances,1);

method_labels = { ...
    'R1 LP reference (industrial)', ...
    'M1 Time-expanded MILP', ...
    'M2 Network simplex (textbook)', ...
    'M3 Column generation (textbook)', ...
    'M4 Standard MM Sinkhorn', ...
    'M5 Graphical Sinkhorn (ours)'};
nM = numel(method_labels);
IDX_LP = 1; IDX_MILP = 2; IDX_NS = 3; IDX_CG = 4; IDX_MM = 5; IDX_GS = 6;

rt   = nan(nM, nInst);      % wall clock
obj  = nan(nM, nInst);      % objective
ops  = nan(nM, nInst);      % pivots / rounds / cycles  (language independent)
gap  = nan(nM, nInst);      % vs the exact path-restricted LP
stat = repmat({''}, nM, nInst);
nvar = nan(1, nInst);       % path LP variables
narc = nan(1, nInst);       % flow arcs
flowlb = nan(1, nInst);     % aggregated flow optimum (lower bound)
pprice = nan(1, nInst);     % price of the path restriction
act    = nan(1, nInst);     % fraction of node-time cells at the bound
cggap  = nan(1, nInst);     % rigorous DW gap if cut off
milpdev= nan(1, nInst);     % max |flow - integer| in the certificate
lad    = nan(nInst, numel(TOL_LAD));      % time-to-accuracy, ours
epsT   = nan(nInst, numel(EPS_SWEEP));
epsB   = nan(nInst, numel(EPS_SWEEP));

%  ── caches ───────────────────────────────────────────────────────────────
%  The FLOW model (network simplex, MILP certificate) has a size that does not
%  depend on P at all, and with the capacity calibrated per n_t the flow
%  problem at a given n_t is literally the same problem.  Your run confirmed
%  it: 316.88 s / 28222 pivots at P=10 and 312.51 s / 28668 pivots at P=20.
%  So solve it ONCE per distinct n_t and reuse -- and calibrate the capacity
%  once per n_t too, which additionally makes a P sweep at fixed n_t a sweep
%  over ONE physical capacity constraint instead of a different one per column.
cap_cache  = containers.Map('KeyType','double','ValueType','double');
ns_cache   = containers.Map('KeyType','double','ValueType','any');
milp_cache = containers.Map('KeyType','double','ValueType','any');
unc_obj    = nan(1, nInst);   % unconstrained objective, per instance
capcost    = nan(1, nInst);   % change in transport-cost component
cg_lad     = nan(nInst, 3);   % CG time to 10% / 1% / 0.1% DW gap
nsproj     = nan(1, nInst);   % calibrated network-simplex projection (s)

%% ═══════════════════════════════════════════════════════════════════════════
%% Main loop
%% ═══════════════════════════════════════════════════════════════════════════
for si = 1:nInst
    P  = min(problem_instances(si,1), pool_size);
    nt = problem_instances(si,2);
    if si == 1
        fprintf('\n\n########## TABLE A: all methods, every cell measured ##########\n');
    elseif ours_only_row(si) && ~ours_only_row(max(si-1,1))
        fprintf('\n\n########## TABLE B: scalability, our method only ##########\n');
        fprintf('##########  baselines are out of reach at these sizes  ##########\n');
    end
    fprintf('\n\n═════════════ P = %d,  n_t = %d ═════════════\n', P, nt);

    S  = build_instance(all_paths(1:P), nt, L, dt_min, eps_ent, eid, e_class, w_vals);
    dt = S.t(2) - S.t(1);
    nV = P*L*S.nnzE;  nvar(si) = nV;
    mem_est = nV * 120;
    nArc = nE*S.nnzE + nInner*nt;  narc(si) = nArc;

    fprintf('  path LP  : %.2fM var, %d eq rows, %d cap rows (~%.2f GB)\n', ...
        nV/1e6, 2*nt + P*(L-1)*nt, S.nN*nt, mem_est/1e9);
    fprintf('  flow model: %d arcs, %d nodes  (independent of P)\n', ...
        nArc, nt + 2*nInner*nt + nt);

    % ── unconstrained solution for the current (P,n_t) instance ─────────
    % The capacity may be shared across a P-sweep at fixed n_t, but the
    % unconstrained objective depends on P and must be recomputed.
    R0 = graph_sinkhorn(S, Inf, TOL, MAXIT, TIME_LIM, DAMP, TOL_LAD);
    unc_obj(si) = R0.obj;

    % ── capacity calibration ─────────────────────────────────────────────
    if strcmp(CAP_MODE,'rate')
        cap = RBAR*dt;
        fprintf('  capacity  : fixed rate %.3f -> %.6f per cell\n', RBAR, cap);
    elseif isKey(cap_cache, nt)
        cap = cap_cache(nt);
        fprintf('  capacity  : %.6f per cell (rate %.3f), calibrated at n_t=%d\n', ...
                cap, cap/dt, nt);
    else
        cap = KAPPA * R0.peak;
        cap_cache(nt) = cap;
        fprintf('  capacity  : %.2f x unconstrained peak %.6f -> %.6f per cell (rate %.3f)\n', ...
                KAPPA, R0.peak, cap, cap/dt);
    end
    fprintf('              unconstrained objective %.6f\n', R0.obj);

    % ── feasibility pre-check ────────────────────────────────────────────
    [fv, fd] = feas_check(S, cap, nV <= FEAS_LP, min(60, TIME_LIM));
    fprintf('  feasibility: %s\n', fd);
    if strcmp(fv,'infeasible')
        for m = 1:nM,  stat{m,si} = 'infeas';  end
        fprintf('  >> INFEASIBLE - skipping all methods.  Raise KAPPA.\n');
        continue
    end

    lp_ref = NaN;  ref_src = 'none';  built = false;  LP = struct();
    cg_ok  = false;

    %% ── M3  textbook two-phase Dantzig-Wolfe ────────────────────────────
    fprintf('  [M3] Column generation (textbook) .. ');
    if ours_only_row(si)
        stat{IDX_CG,si} = 'skip';  fprintf('SKIP (scalability row: our method only)\n');
    else
    try
        tm = tic;
        [oCG, lbCG, nCol, nRnd, fCG, tlCG] = colgen(S, cap, TIME_LIM, CG_ADD, CG_LAD);
        cg_lad(si,:) = tlCG;
        e = toc(tm);
        ops(IDX_CG,si) = nRnd;
        if fCG == 1
            rt(IDX_CG,si) = e;  obj(IDX_CG,si) = oCG;  stat{IDX_CG,si} = 'ok';
            if isnan(lp_ref)
                lp_ref = oCG;  ref_src = 'exact optimum (CG)';  cg_ok = true;
            end
            gap(IDX_CG,si) = abs(oCG-lp_ref)/(abs(lp_ref)+1e-16);
            fprintf('%8.2f s   obj %.6f  (%d rounds, %d cols)\n', e, oCG, nRnd, nCol);
        elseif fCG == -2
            stat{IDX_CG,si} = 'infeas';
            fprintf('PHASE-1 STALL (%d rounds) - looks capacity-infeasible\n', nRnd);
        else
            stat{IDX_CG,si} = 'time';  rt(IDX_CG,si) = e;  obj(IDX_CG,si) = oCG;
            if isfinite(lbCG) && lbCG > 0
                cggap(si) = (oCG-lbCG)/abs(lbCG);
                fprintf('TIME LIMIT > %d s  (%d rounds, %d cols, obj %.6f, DW gap %.2f%%)\n', ...
                        TIME_LIM, nRnd, nCol, oCG, 100*cggap(si));
                %  a tight DW bracket is a perfectly good reference
                if cggap(si) <= CG_GAP_OK && isnan(lp_ref)
                    lp_ref = oCG;  cg_ok = true;
                    ref_src = sprintf('CG bound, optimum within %.2f%%', 100*cggap(si));
                    fprintf('        -> using it as the reference: optimum in [%.6f, %.6f]\n', ...
                            lbCG, oCG);
                end
            else
                fprintf('TIME LIMIT > %d s  (%d rounds, %d cols)\n', TIME_LIM, nRnd, nCol);
            end
        end
    catch ME
        stat{IDX_CG,si} = 'err';  fprintf('ERROR %s\n', ME.message);
    end
    end

    %% ── R1  exact path-restricted LP (industrial reference) ─────────────
    fprintf('  [R1] LP reference (barrier) ........ ');
    if ours_only_row(si)
        stat{IDX_LP,si} = 'skip';  fprintf('SKIP (scalability row: our method only)\n');
    elseif LP_ONLY_IF_CG_FAILS && cg_ok
        stat{IDX_LP,si} = 'skip';
        fprintf('SKIP  reference already available (%s): %.6f\n', ref_src, lp_ref);
    elseif nV > LP_VAR_LIM
        stat{IDX_LP,si} = 'mem';
        fprintf('SKIP  %.2fM var > %.2fM: a barrier solve cannot be time-capped\n', ...
                nV/1e6, LP_VAR_LIM/1e6);
    elseif mem_est > MEM_LIM
        stat{IDX_LP,si} = 'mem';  fprintf('SKIP  ~%.1f GB > %.0f GB\n', mem_est/1e9, MEM_LIM/1e9);
    else
        try
            LP = build_path_lp(S, cap);  built = true;
            o = optimoptions('linprog','Algorithm','interior-point', ...
                             'Display','off','MaxIterations',5000, ...
                             'OptimalityTolerance',1e-9);
            tm = tic;
            [~, oR, fR] = linprog(LP.c, LP.Aub, LP.bub, LP.Aeq, LP.beq, ...
                                  zeros(LP.nV,1), [], o);
            if fR ~= 1        % barrier hit its iteration cap -> try simplex
                o2 = optimoptions('linprog','Algorithm','dual-simplex', ...
                                  'Display','off','MaxTime',TIME_LIM, ...
                                  'MaxIterations',5e6);
                [~, oR, fR] = linprog(LP.c, LP.Aub, LP.bub, LP.Aeq, LP.beq, ...
                                      zeros(LP.nV,1), [], o2);
            end
            e = toc(tm);
            if fR == 1
                lp_ref = oR;  ref_src = 'exact LP';
                rt(IDX_LP,si) = e;  obj(IDX_LP,si) = oR;
                gap(IDX_LP,si) = 0;  stat{IDX_LP,si} = 'ok';
                fprintf('%8.2f s   obj %.6f\n', e, oR);
            else
                stat{IDX_LP,si} = 'err';  fprintf('linprog flag %d\n', fR);
            end
        catch ME
            stat{IDX_LP,si} = 'err';  fprintf('ERROR %s\n', ME.message);
        end
    end

    %% ── M2  network simplex on the node-split flow graph ────────────────
    %  Built before M1 because M1 reuses the same graph.
    F = [];
    fprintf('  [M2] Network simplex (textbook) .... ');
    if ours_only_row(si)
        stat{IDX_NS,si} = 'skip';  fprintf('SKIP (scalability row: our method only)\n');
    elseif isKey(ns_cache, nt)
        Z = ns_cache(nt);
        rt(IDX_NS,si) = Z.rt;  obj(IDX_NS,si) = Z.obj;  ops(IDX_NS,si) = Z.piv;
        stat{IDX_NS,si} = Z.stat;  flowlb(si) = Z.obj;
        fprintf('reuse (flow model is identical at n_t=%d): %.2f s, obj %.6f, %d pivots\n', ...
                nt, Z.rt, Z.obj, Z.piv);
        if ~isnan(lp_ref) && ~isnan(Z.obj)
            pprice(si) = (lp_ref - Z.obj)/abs(Z.obj);
            fprintf('        price of the path restriction: %+.1f%%\n', 100*pprice(si));
        end
    elseif nArc > NS_ARCS
        stat{IDX_NS,si} = 'mem';
        fprintf('SKIP  %d arcs > %.0e limit\n', nArc, NS_ARCS);
    elseif NS_PIV_PER_NODE*(nt+2*nInner*nt+nt)*NS_SEC_PER_PIV_ARC*nArc > NS_TIME
        est = NS_PIV_PER_NODE*(nt+2*nInner*nt+nt)*NS_SEC_PER_PIV_ARC*nArc;
        stat{IDX_NS,si} = 'proj';  nsproj(si) = est;
        fprintf('SKIP  projected %.1f h (%.0f pivots x %d arcs) > NS_TIME %.0f s\n', ...
                est/3600, NS_PIV_PER_NODE*(nt+2*nInner*nt+nt), nArc, NS_TIME);
        fprintf('        projection is calibrated on the measured n_t=50 solve\n');
    else
        try
            F = build_flow(S, adj, src, snk, nInner, nE, cap, eid, e_class, w_vals);
            tm = tic;
            [oNS, xNS, piv, stNS] = netsimplex(F.n, F.tail, F.head, F.cost, F.cap, ...
                                               F.b, NS_MAXPIV, 1e-10, NS_TIME);
            e = toc(tm);
            ops(IDX_NS,si) = piv;
            if strcmp(stNS,'optimal')
                rt(IDX_NS,si) = e;  obj(IDX_NS,si) = oNS;  stat{IDX_NS,si} = 'ok';
                flowlb(si) = oNS;
                ns_cache(nt) = struct('obj',oNS,'rt',e,'piv',piv,'stat','ok');
                fprintf('%8.2f s   obj %.6f  (%d pivots, %.2f/node)\n', ...
                        e, oNS, piv, piv/F.n);
                if ~isnan(lp_ref)
                    pprice(si) = (lp_ref - oNS)/abs(oNS);
                    fprintf('        price of the path restriction: %+.1f%% ', 100*pprice(si));
                    fprintf('(flow %.6f <= path %.6f)\n', oNS, lp_ref);
                end
            elseif strcmp(stNS,'timelimit')
                stat{IDX_NS,si} = 'time';
                ns_cache(nt) = struct('obj',NaN,'rt',e,'piv',piv,'stat','time');
                fprintf('TIME LIMIT > %d s  (%d pivots, %.0f pivots/s)\n', ...
                        NS_TIME, piv, piv/max(e,eps));
            else
                stat{IDX_NS,si} = 'err';
                fprintf('status %s after %d pivots\n', stNS, piv);
            end
        catch ME
            stat{IDX_NS,si} = 'err';  fprintf('ERROR %s\n', ME.message);
        end
    end

    %% ── M1  MILP integrality certificate ────────────────────────────────
    fprintf('  [M1] Time-expanded MILP ............ ');
    if ours_only_row(si)
        stat{IDX_MILP,si} = 'skip';  fprintf('SKIP (scalability row: our method only)\n');
    elseif isKey(milp_cache, nt)
        Z = milp_cache(nt);
        milpdev(si) = Z.dev;  rt(IDX_MILP,si) = Z.rt;  stat{IDX_MILP,si} = Z.stat;
        fprintf('reuse (n_t=%d): max |flow-int| = %.2e\n', nt, Z.dev);
    elseif MILP_ONCE && milp_cache.Count > 0
        stat{IDX_MILP,si} = 'skip';
        fprintf('SKIP  certificate already established (one instance suffices)\n');
    elseif isempty(F)
        stat{IDX_MILP,si} = 'n/a';
        fprintf('n/a  (flow graph not built at this size)\n');
    else
        try
            % integer supplies and integer arc capacity
            sup = floor(S.mu0*MILP_SCALE);  dem = floor(S.muT*MILP_SCALE);
            [~,imx] = max(sup);  sup(imx) = sup(imx) + (sum(dem)-sum(sup));
            bi = zeros(F.n,1);
            bi(1:nt) = sup;
            bi(F.off_snk+1:F.off_snk+nt) = -dem;
            capi = F.cap;  capi(F.split) = floor(cap*MILP_SCALE);
            Ai = sparse([F.tail; F.head], [1:F.m, 1:F.m]', ...
                        [ones(F.m,1); -ones(F.m,1)], F.n, F.m);
            %  The certificate MUST come from a VERTEX (basic) solution.
            %  Total unimodularity says an integral OPTIMAL VERTEX exists; it
            %  says nothing about interior points, and MATLAB's interior-point
            %  linprog does no crossover, so it returns a point in the relative
            %  interior of the optimal face -- typically fractional even when
            %  the theorem holds.  (Observed: max |flow-int| = 4.97e-01 with
            %  barrier vs exactly 0 with a vertex solver.)  Hence dual simplex.
            o = optimoptions('linprog','Algorithm','dual-simplex', ...
                             'Display','off','MaxTime',TIME_LIM, ...
                             'MaxIterations',5e6);
            tm = tic;
            [xi, oi, fi] = linprog(F.cost, [], [], Ai, bi, zeros(F.m,1), capi, o);
            e = toc(tm);
            if fi == 1
                dev = max(abs(xi - round(xi)));
                milpdev(si) = dev;
                stat{IDX_MILP,si} = 'ok';  rt(IDX_MILP,si) = e;
                milp_cache(nt) = struct('dev',dev,'rt',e,'stat','ok');
                fprintf('LP relaxation optimum %.4f  ', oi);
                fprintf('max |flow-int| = %.2e over %d arcs\n', dev, F.m);
                if dev < 1e-6
                    fprintf('        => integral vertex: matrix is totally unimodular,\n');
                    fprintf('           a MILP formulation adds nothing but branching.\n');
                else
                    fprintf('        => NOT integral - check the solver returned a VERTEX\n');
                    fprintf('           (barrier without crossover will not).\n');
                end
            else
                stat{IDX_MILP,si} = 'err';  fprintf('certificate LP flag %d\n', fi);
            end
        catch ME
            stat{IDX_MILP,si} = 'err';  fprintf('ERROR %s\n', ME.message);
        end
    end
    F = [];   %#ok<NASGU>  free the flow graph

    %% ── M4  standard multi-marginal Sinkhorn ────────────────────────────
    fprintf('  [M4] Standard MM Sinkhorn .......... ');
    stat{IDX_MM,si} = 'n/a';
    fprintf('n/a  full tensor n_t^%d = 10^%.0f entries (10^%.0f bytes)\n', ...
            L+1, (L+1)*log10(nt), (L+1)*log10(nt)+log10(8));

    %% ── M5  graphical Sinkhorn (ours) ───────────────────────────────────
    fprintf('  [M5] Graphical Sinkhorn (ours) ..... ');
    try
        tm = tic;
        R = graph_sinkhorn(S, cap, TOL, MAXIT, TIME_LIM, DAMP, TOL_LAD);
        e = toc(tm);
        rt(IDX_GS,si) = e;  obj(IDX_GS,si) = R.obj;  ops(IDX_GS,si) = R.it;
        lad(si,:) = R.tl;  act(si) = R.nAct / R.nCell;
        if R.flag == 1
            stat{IDX_GS,si} = 'ok';
            fprintf('%8.2f s   obj %.6f  (%d cyc, E0 %.1e ET %.1e V %.1e)\n', ...
                    e, R.obj, R.it, R.E0, R.ET, R.V);
            fprintf('        capacity active at %d/%d cells (%.1f%%), peak/cap %.3f\n', ...
                    R.nAct, R.nCell, 100*act(si), R.util);
            if act(si) < 0.005
                fprintf('        >> constraint essentially IDLE: lower KAPPA, otherwise\n');
                fprintf('        >> this is a chain-Sinkhorn benchmark, not the paper''s problem\n');
            end
        elseif R.flag == -3
            stat{IDX_GS,si} = 'brkdn';
            fprintf('NUMERICAL BREAKDOWN (%d cyc, %.2f s) E0 %.1e ET %.1e V %.1e\n', ...
                    R.it, e, R.E0, R.ET, R.V);
            fprintf('        >> multipliers left the representable range: the iteration is\n');
            fprintf('        >> DIVERGING at this capacity, not underflowing.  Raise KAPPA\n');
            fprintf('        >> (bind less hard), or move to log-domain updates.\n');
        elseif R.flag == -1
            stat{IDX_GS,si} = 'time';
            fprintf('TIME LIMIT > %d s (%d cyc)\n', TIME_LIM, R.it);
        else
            stat{IDX_GS,si} = 'iter';
            fprintf('ITERATION CAP (%d cyc, %.2f s) E0 %.1e ET %.1e V %.1e\n', ...
                    R.it, e, R.E0, R.ET, R.V);
            if R.Vstall
                fprintf('        >> V plateaued with small boundary residuals:\n');
                fprintf('        >> capacity-INFEASIBLE, not slow convergence\n');
            end
        end
        if ~isnan(unc_obj(si)) && isfinite(R.obj)
            capcost(si) = (R.obj - unc_obj(si))/abs(unc_obj(si));
            fprintf('        transport-cost change: constrained %.6f vs unconstrained %.6f = %+.2f%%\n', ...
                    R.obj, unc_obj(si), 100*capcost(si));
        end
        if ~isnan(lp_ref)
            gap(IDX_GS,si) = (R.obj-lp_ref)/abs(lp_ref);
            if ~isnan(cggap(si)) && cggap(si) > 1e-9
                fprintf('        entropic bias %+.3f %% vs %s\n', ...
                        100*gap(IDX_GS,si), ref_src);
                fprintf('        (reference itself known only to %.2f%%, so the bias is\n', ...
                        100*cggap(si));
                fprintf('         %+.3f %% +/- %.2f %% -- immaterial at this scale)\n', ...
                        100*gap(IDX_GS,si), 100*cggap(si));
            else
                fprintf('        entropic bias vs %s: %+.3f %%\n', ref_src, 100*gap(IDX_GS,si));
            end
        end
    catch ME
        stat{IDX_GS,si} = 'err';  fprintf('ERROR %s\n', ME.message);
    end

    %% ── eps frontier where ground truth exists ──────────────────────────
    if RUN_EPS_SWEEP && ~isnan(lp_ref) && numel(EPS_SWEEP) > 1
        fprintf('  [M5b] eps frontier (reference = %.6f, %s):\n', lp_ref, ref_src);
        for q = 1:numel(EPS_SWEEP)
            try
                Sq = build_instance(all_paths(1:P), nt, L, dt_min, EPS_SWEEP(q), ...
                                    eid, e_class, w_vals);
                tm = tic;
                Rq = graph_sinkhorn(Sq, cap, TOL, MAXIT, TIME_LIM, DAMP, TOL_LAD);
                epsT(si,q) = toc(tm);
                epsB(si,q) = (Rq.obj-lp_ref)/abs(lp_ref);
                fprintf('        eps %-6.3f %7.2f s (%5d cyc)  obj %.6f  bias %+7.3f %%  %s\n', ...
                    EPS_SWEEP(q), epsT(si,q), Rq.it, Rq.obj, 100*epsB(si,q), ...
                    tern(Rq.flag==1,'converged','NOT conv'));
            catch ME
                fprintf('        eps %-6.3f ERROR %s\n', EPS_SWEEP(q), ME.message);
            end
        end
    end
end

%% ═══════════════════════════════════════════════════════════════════════════
%% Report
%% ═══════════════════════════════════════════════════════════════════════════
report(problem_instances, method_labels, rt, obj, ops, gap, stat, ...
       nvar, narc, flowlb, pprice, act, cggap, milpdev, lad, ...
       epsT, epsB, TOL_LAD, EPS_SWEEP, TIME_LIM, TOL, MEM_LIM, capcost, cg_lad, nsproj, CG_LAD);

save('da_table_final_results.mat', 'problem_instances', 'method_labels', ...
     'rt', 'obj', 'ops', 'gap', 'stat', 'nvar', 'narc', 'flowlb', 'pprice', ...
     'act', 'cggap', 'milpdev', 'lad', 'epsT', 'epsB', 'unc_obj', 'capcost', ...
     'cg_lad', 'nsproj');
fprintf('\nSaved raw results to da_table_final_results.mat\n');

%  Flat CSV: one row per (method, instance).  Build the paper table from this
%  rather than by scraping the console.
fid = fopen('da_table_final_results.csv','w');
fprintf(fid, 'table,P,n_t,method,status,seconds,ops,objective,gap_vs_ref,%s\n', ...
        'flow_lower_bound,path_restriction_price,capacity_cost,cap_active_frac');
for si = 1:nInst
    for m = 1:nM
        fprintf(fid, '%s,%d,%d,%s,%s,', tern(ours_only_row(si),'B','A'), ...
                problem_instances(si,1), problem_instances(si,2), ...
                strrep(method_labels{m},',',';'), stat{m,si});
        fprintf(fid, '%s,%s,%s,%s,', numstr_(rt(m,si)), numstr_(ops(m,si)), ...
                numstr_(obj(m,si)), numstr_(gap(m,si)));
        fprintf(fid, '%s,%s,%s,%s\n', numstr_(flowlb(si)), numstr_(pprice(si)), ...
                numstr_(capcost(si)), numstr_(act(si)));
    end
end
fclose(fid);
fprintf('Saved a flat table to da_table_final_results.csv\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% ══════════════════════════ LOCAL FUNCTIONS ═══════════════════════════════
%% ═══════════════════════════════════════════════════════════════════════════

function out = tern(c, a, b)
    if c,  out = a;  else,  out = b;  end
end

%% ─────────────────────────────────────────────────────────────────────────
function S = build_instance(paths, nt, L, dt_min, eps_ent, eid, e_class, w_vals)
    P = numel(paths);
    t = linspace(0,1,nt)';
    [U,V] = ndgrid(t,t);
    supp = (V-U) >= dt_min;
    D = V-U;
    nC = numel(w_vals);
    Cc = cell(nC,1); Kc = cell(nC,1); CKc = cell(nC,1); CcInf = cell(nC,1);
    for c = 1:nC
        Ci = zeros(nt);  Ci(supp) = w_vals(c) ./ D(supp);
        Ki = zeros(nt);  Ki(supp) = exp(-Ci(supp)/eps_ent);
        %  NO per-class normalisation.  Dividing each class by its own max is
        %  NOT a gauge: it multiplies path p's kernel by 1/M_p with
        %  M_p = prod_l m_{e_l}, which depends on WHICH classes path p uses.
        %  The scalings that could absorb a constant -- u(t_0), v(t_T),
        %  W_v(t) -- are indexed by time and node, never by path, so M_p
        %  survives.  It is equivalent to adding eps*log(M_p) to each path's
        %  cost, i.e. reweighting the paths, and it changes the entropic
        %  optimizer.  Measured (P=10, n_t=50, eps=0.04):
        %      no normalisation      obj 1.40110031   peak 0.08812276
        %      per-class (was here)  obj 1.40432120   peak 0.08625430  (+0.23%)
        %      one global constant   obj 1.40110031   peak 0.08812276  (exact)
        %  A SINGLE global constant is a genuine gauge here only because every
        %  path has the same length L, so M_p = gmax^L is common to all paths.
        %  That would break if path lengths differed, so do not reintroduce it.
        %  Unnormalised magnitudes are safe in double precision at these
        %  parameters: kernel entries lie in [2.6e-4, 0.88], so an 11-edge
        %  product is >= 1e-39.  For eps well below 0.01 use log-domain
        %  updates (Remark 6) rather than rescaling.
        Cc{c} = Ci;  Kc{c} = Ki;  CKc{c} = Ci.*Ki;
        Ce = Ci;  Ce(~supp) = inf;  CcInf{c} = Ce;
    end
    ZInf = zeros(nt);  ZInf(~supp) = inf;

    reach = supp;                       % exact grid reachability over L edges
    for k = 2:L,  reach = (double(reach)*double(supp)) > 0;  end
    G = @(m,s) exp(-0.5*((t-m)/s).^2);
    mu0 = 0.9*G(0.10,0.05) + 0.6*G(0.22,0.05);
    muT = 0.8*G(0.76,0.05) + 0.7*G(0.88,0.05);
    mu0 = mu0 .* (sum(reach,2)>0);   mu0 = mu0/sum(mu0);
    muT = muT .* (sum(reach,1)'>0);  muT = muT/sum(muT);

    pcls = zeros(P,L);  pnod = zeros(P,L-1);
    for p = 1:P
        q = paths{p};
        for l = 1:L,    pcls(p,l) = e_class(full(eid(q(l),q(l+1))));  end
        for l = 1:L-1,  pnod(p,l) = q(l+1);  end
    end
    nodelist = unique(pnod(:));  nN = numel(nodelist);
    invmap = zeros(max(nodelist),1);  invmap(nodelist) = 1:nN;
    nidx = reshape(invmap(pnod), P, L-1);

    Sel = cell(L-1,1);
    for l = 1:L-1,  Sel{l} = sparse((1:P)', nidx(:,l), 1, P, nN);  end
    Grp = cell(L,1);
    for l = 1:L
        cls = unique(pcls(:,l))';  g = cell(numel(cls),1);
        for k = 1:numel(cls)
            g{k} = struct('c',cls(k),'idx',find(pcls(:,l)==cls(k))');
        end
        Grp{l} = g;
    end
    [rs,cs] = find(supp);
    S = struct('paths',{paths},'t',t,'nt',nt,'P',P,'L',L,'supp',supp, ...
               'Cc',{Cc},'Kc',{Kc},'CKc',{CKc},'CcInf',{CcInf},'ZInf',ZInf, ...
               'mu0',mu0,'muT',muT,'pcls',pcls,'nidx',nidx,'nN',nN, ...
               'nodelist',nodelist,'Sel',{Sel},'Grp',{Grp},'rs',rs,'cs',cs, ...
               'nnzE',numel(rs),'eps_ent',eps_ent,'w_vals',w_vals);
end

%% ─────────────────────────────────────────────────────────────────────────
function F = build_flow(S, adj, src, snk, nInner, nE, cap, eid, e_class, w_vals)
%BUILD_FLOW  Node-split time-expanded min-cost flow graph.
%   NODE SPLITTING is what turns the nodal capacity into an ARC capacity and
%   makes the whole thing a pure network:
%       (v,t)_in --[cap]--> (v,t)_out
%   every DAG edge u->v becomes (u,s)_out -> (v,t)_in for t-s >= Delta with
%   cost w_uv/(t-s);  supply mu0(t) at (src,t), demand muT(t) at (snk,t).
%   Size does NOT depend on P.
    nt = S.nt;  t = S.t;  rs = S.rs;  cs = S.cs;  m0 = S.nnzE;
    off_in  = nt;
    off_out = nt + nInner*nt;
    off_snk = nt + 2*nInner*nt;
    n = off_snk + nt;

    nsplit = nInner*nt;
    ntrans = nE*m0;
    tail = zeros(nsplit+ntrans,1);  head = tail;  cost = tail;  ub = tail;

    % splitting arcs carry the nodal capacity
    k = (1:nt)';
    for v = 1:nInner
        r = (v-1)*nt + (1:nt)';
        tail(r) = off_in  + (v-1)*nt + k;
        head(r) = off_out + (v-1)*nt + k;
        cost(r) = 0;
        ub(r)   = cap;
    end
    split = (1:nsplit)';

    q = nsplit;
    for u = 1:(nInner+2)
        for v = adj{u}
            w = w_vals(e_class(full(eid(u,v))));
            if u == src,  uu = rs;               else,  uu = off_out + (u-1)*nt + rs;  end
            if v == snk,  vv = off_snk + cs;     else,  vv = off_in  + (v-1)*nt + cs;  end
            r = q + (1:m0)';
            tail(r) = uu;  head(r) = vv;
            cost(r) = w ./ (t(cs) - t(rs));
            ub(r)   = inf;
            q = q + m0;
        end
    end
    tail = tail(1:q); head = head(1:q); cost = cost(1:q); ub = ub(1:q);

    b = zeros(n,1);
    b(1:nt) = S.mu0;
    b(off_snk+1:off_snk+nt) = -S.muT;
    F = struct('n',n,'m',q,'tail',tail,'head',head,'cost',cost,'cap',ub, ...
               'b',b,'split',split,'off_snk',off_snk);
end

%% ─────────────────────────────────────────────────────────────────────────
function [obj, x, pivots, status] = netsimplex(n, tailA, headA, cost, cap, b, ...
                                               maxiter, tol, time_lim)
%NETSIMPLEX  TEXTBOOK primal network simplex, Ahuja-Magnanti-Orlin Ch. 11.
%   min c'x  s.t.  A x = b (node-arc incidence),  0 <= x <= cap.
%   No shortcuts and no tuning:
%     * artificial-root basis for the initial spanning tree
%     * node potentials read off the tree
%     * entering arc by the Dantzig rule (full arc scan)
%     * cycle by walking both endpoints to their common ancestor
%     * leaving arc by the ratio test
%   O(m + n) per pivot.  The arc scan is one vectorised reduction (selecting
%   the most negative reduced cost IS a reduction, not a loop) -- that is a
%   MATLAB idiom, not an algorithmic change.
%
%   Validated against HiGHS: exact agreement on 12 random capacitated
%   min-cost flows, and 1.9e-16 relative difference on the node-split
%   time-expanded graph.
    m = numel(cost);
    ROOT = n+1;  N = n+1;
    BIG = 1 + 8*(max(abs(cost))+1)*max(1,sum(abs(b)));

    ii = (1:n)';
    a_tail = ii;              a_head = repmat(ROOT,n,1);
    neg = b < 0;
    a_tail(neg) = ROOT;       a_head(neg) = ii(neg);

    tail = [tailA(:); a_tail];
    head = [headA(:); a_head];
    c    = [cost(:);  BIG*ones(n,1)];
    u    = [cap(:);   inf(n,1)];
    M    = m + n;

    x = zeros(M,1);  x(m+1:M) = abs(b);
    parent = zeros(N,1);  pred = zeros(N,1);  depth = zeros(N,1);
    parent(1:n) = ROOT;   pred(1:n) = m + (1:n)';   depth(1:n) = 1;
    inT = false(M,1);  inT(m+1:M) = true;
    atU = false(M,1);

    t0 = tic;  pivots = 0;  status = 'maxiter';
    for it = 1:maxiter
        % ---- node potentials from the spanning tree.
        %  Processed one DEPTH LEVEL at a time: every node at depth d has its
        %  parent at depth d-1, already done, so a whole level updates in one
        %  vector operation.  Same arithmetic as the textbook node-by-node
        %  sweep, just not one MATLAB loop iteration per node.
        pi = zeros(N,1);
        maxd = max(depth);
        for d = 1:maxd
            vd = find(depth == d);
            if isempty(vd),  continue;  end
            a  = pred(vd);
            sg = 2*(tail(a) == vd) - 1;        % +1 child->parent, -1 otherwise
            pi(vd) = pi(parent(vd)) + sg .* c(a);
        end
        % ---- entering arc, Dantzig rule (vectorised full scan)
        rc = c - pi(tail) + pi(head);
        gain = -inf(M,1);
        s1 = ~inT & ~atU & (rc < -tol);   gain(s1) = -rc(s1);
        s2 = ~inT &  atU & (rc >  tol);   gain(s2) =  rc(s2);
        [g,e] = max(gain);
        if g <= tol,  status = 'optimal';  break;  end
        pivots = pivots + 1;

        if ~atU(e)
            i = tail(e);  j = head(e);  sgn_e =  1;  room_e = u(e) - x(e);
        else
            i = head(e);  j = tail(e);  sgn_e = -1;  room_e = x(e);
        end

        % ---- common ancestor
        mark = false(N,1);  v = i;
        while v ~= 0,  mark(v) = true;  v = parent(v);  end
        apex = j;
        while ~mark(apex),  apex = parent(apex);  end

        % ---- cycle arcs, signs and residual room
        arcs = zeros(2*N,1);  sgns = zeros(2*N,1);  rooms = zeros(2*N,1);
        nc = 1;  arcs(1) = e;  sgns(1) = sgn_e;  rooms(1) = room_e;
        v = j;                                   % j-branch: child -> parent
        while v ~= apex
            a = pred(v);
            s = tern(tail(a) == v, 1, -1);
            nc = nc+1;  arcs(nc) = a;  sgns(nc) = s;
            rooms(nc) = tern(s > 0, u(a)-x(a), x(a));
            v = parent(v);
        end
        v = i;                                   % i-branch: parent -> child
        while v ~= apex
            a = pred(v);
            s = tern(tail(a) == v, -1, 1);
            nc = nc+1;  arcs(nc) = a;  sgns(nc) = s;
            rooms(nc) = tern(s > 0, u(a)-x(a), x(a));
            v = parent(v);
        end
        arcs = arcs(1:nc);  sgns = sgns(1:nc);  rooms = rooms(1:nc);

        delta = min(rooms);
        if ~isfinite(delta),  status = 'unbounded';  break;  end

        % leaving arc: last blocking arc on the cycle
        blk = find(rooms <= delta + tol);
        leaving = arcs(blk(end));  leave_sgn = sgns(blk(end));

        x(arcs) = x(arcs) + sgns*delta;

        if leaving == e
            atU(e) = (sgn_e > 0);
        else
            inT(leaving) = false;  atU(leaving) = (leave_sgn > 0);
            inT(e) = true;         atU(e) = false;
            [parent, pred, depth] = rebuild(N, ROOT, inT, tail, head);
        end
        if toc(t0) > time_lim,  status = 'timelimit';  break;  end
    end
    obj = c(1:m)' * x(1:m);
    if strcmp(status,'optimal') && sum(x(m+1:M)) > 1e-7
        status = 'infeasible';
    end
    x = x(1:m);
end

function [parent, pred, depth] = rebuild(N, ROOT, inT, tail, head)
%REBUILD  Refresh parent / pred / depth by BFS over the basis arcs.
%   Breadth-first search, one LEVEL per loop iteration rather than one node,
%   over a CSR adjacency built by sorting.  Identical tree to a node-by-node
%   BFS; the vectorisation is a MATLAB idiom, not an algorithmic change.
%   Growing cell arrays per pivot (the obvious way) costs O(n) appends per
%   pivot and makes the whole solve hours long, which is an implementation
%   artefact and not a property of network simplex.
    ba = find(inT);
    s2 = [tail(ba); head(ba)];
    d2 = [head(ba); tail(ba)];
    a2 = [ba;       ba];
    [s2, o2] = sort(s2);
    d2 = d2(o2);  a2 = a2(o2);
    ptr = [0; cumsum(accumarray(s2, 1, [N 1]))];

    parent = zeros(N,1);  pred = zeros(N,1);  depth = zeros(N,1);
    vis = false(N,1);  vis(ROOT) = true;
    front = ROOT;  d = 0;
    while ~isempty(front)
        d = d + 1;
        st  = ptr(front) + 1;
        en  = ptr(front + 1);
        len = en - st + 1;
        len(len < 0) = 0;
        tot = sum(len);
        if tot == 0,  break;  end
        keepf = len > 0;
        stk = st(keepf);  enk = en(keepf);  frk = front(keepf);  lnk = len(keepf);
        ends   = cumsum(lnk);
        starts = ends - lnk + 1;
        inc = ones(tot,1);
        inc(starts) = stk - [0; enk(1:end-1)];
        idx = cumsum(inc);
        cand   = d2(idx);
        carc   = a2(idx);
        fowner = repelem(frk(:), lnk);
        k = ~vis(cand);
        cand = cand(k);  carc = carc(k);  fowner = fowner(k);
        if isempty(cand),  break;  end
        [cu, iu] = unique(cand, 'stable');
        parent(cu) = fowner(iu);
        pred(cu)   = carc(iu);
        depth(cu)  = d;
        vis(cu)    = true;
        front      = cu;
    end
end

%% ─────────────────────────────────────────────────────────────────────────
function LP = build_path_lp(S, cap)
%BUILD_PATH_LP  Path-restricted time-expanded LP (the paper's Problem 3).
%   Variables x_{p,l}(s,t).  Rows: departure | arrival | chain consistency |
%   node-wise capacity.  This is the smallest exact LP in which the interior
%   crossing-time marginals appear at all.
    nt = S.nt;  P = S.P;  L = S.L;  m = S.nnzE;
    rs = S.rs;  cs = S.cs;
    nV = P*L*m;
    base = @(p,l) ((p-1)*L + (l-1))*m;

    c = zeros(nV,1);
    dtv = S.t(cs) - S.t(rs);
    for p = 1:P
        for l = 1:L
            c(base(p,l)+(1:m)') = S.w_vals(S.pcls(p,l)) ./ dtv;
        end
    end
    nEq = 2*nt + P*(L-1)*nt;
    nz  = 2*P*m + 2*P*(L-1)*m;
    I = zeros(nz,1); J = zeros(nz,1); Vv = zeros(nz,1); q = 0;
    for p = 1:P
        j = base(p,1)+(1:m)';
        I(q+1:q+m)=rs;       J(q+1:q+m)=j; Vv(q+1:q+m)= 1; q=q+m;
        j = base(p,L)+(1:m)';
        I(q+1:q+m)=nt+cs;    J(q+1:q+m)=j; Vv(q+1:q+m)= 1; q=q+m;
    end
    off = 2*nt;
    for p = 1:P
        for l = 2:L
            r0 = off + ((p-1)*(L-1)+(l-2))*nt;
            j = base(p,l-1)+(1:m)';
            I(q+1:q+m)=r0+cs;  J(q+1:q+m)=j; Vv(q+1:q+m)= 1; q=q+m;
            j = base(p,l)+(1:m)';
            I(q+1:q+m)=r0+rs;  J(q+1:q+m)=j; Vv(q+1:q+m)=-1; q=q+m;
        end
    end
    Aeq = sparse(I(1:q),J(1:q),Vv(1:q),nEq,nV);
    beq = [S.mu0; S.muT; zeros(P*(L-1)*nt,1)];

    nz = P*(L-1)*m;  I = zeros(nz,1); J = zeros(nz,1); q = 0;
    for p = 1:P
        for l = 2:L
            k = S.nidx(p,l-1);
            j = base(p,l-1)+(1:m)';
            I(q+1:q+m) = (k-1)*nt+cs;  J(q+1:q+m) = j;  q = q+m;
        end
    end
    Aub = sparse(I(1:q),J(1:q),1,S.nN*nt,nV);
    LP = struct('c',c,'Aeq',Aeq,'beq',beq,'Aub',Aub, ...
                'bub',cap*ones(S.nN*nt,1),'nV',nV);
end

%% ─────────────────────────────────────────────────────────────────────────
function [verdict, detail] = feas_check(S, cap, do_lp, time_lim)
%FEAS_CHECK  (a) cheap per-stage necessary bound, (b) exact min-slack LP.
%   (a) is necessary but NOT sufficient -- it passed at 1.12 on an instance
%   that was genuinely infeasible, because mu0 is concentrated and each path's
%   causal window is narrow.  Only (b) settles it.
    nt = S.nt;  L = S.L;
    Fs = cell(L+1,1);  Fs{1} = (S.mu0>0)';
    for k = 1:L,  Fs{k+1} = (double(Fs{k})*double(S.supp)) > 0;  end
    Bs = cell(L+1,1);  Bs{1} = (S.muT>0);
    for k = 1:L,  Bs{k+1} = (double(S.supp)*double(Bs{k})) > 0;  end
    worst = inf;
    for k = 1:L-1
        reach = Fs{k+1}(:)' & Bs{L-k+1}(:)';
        worst = min(worst, numel(unique(S.nidx(:,k)))*sum(reach)*cap);
    end
    if ~do_lp
        if worst < 1
            verdict='infeasible';
            detail=sprintf('per-stage bound %.3f < 1 -> necessarily infeasible',worst);
        else
            verdict='unknown';
            detail=sprintf('per-stage bound %.3f >= 1 (necessary only; LP too large)',worst);
        end
        return
    end
    LP = build_path_lp(S, cap);
    ns = 2*nt;
    Sp = sparse(1:ns,1:ns,1,size(LP.Aeq,1),ns);
    A  = [LP.Aeq, Sp, -Sp];
    Au = [LP.Aub, sparse(size(LP.Aub,1), 2*ns)];
    cc = [zeros(LP.nV,1); ones(2*ns,1)];
    o = optimoptions('linprog','Algorithm','dual-simplex','Display','off', ...
                     'MaxTime',time_lim,'MaxIterations',3e6);
    [~, sl, fl] = linprog(cc, Au, LP.bub, A, LP.beq, zeros(LP.nV+2*ns,1), [], o);
    if fl ~= 1
        verdict='unknown';
        detail=sprintf('per-stage bound %.3f; slack LP inconclusive (flag %d)',worst,fl);
    elseif sl > 1e-8
        verdict='infeasible';
        detail=sprintf('per-stage bound %.3f but EXACT min slack %.4f -> only %.1f%% routable', ...
                       worst, sl, 100*(1-sl));
    else
        verdict='feasible';
        detail=sprintf('exact min boundary slack %.2e (per-stage bound %.3f)',sl,worst);
    end
end

%% ─────────────────────────────────────────────────────────────────────────
function [obj, lb, nCols, nRounds, flag, tlad] = colgen(S, cap, time_lim, add_per_round, glad)
%COLGEN  TEXTBOOK two-phase Dantzig-Wolfe.  Columns = (path, causal tuple).
%   Phase 1 minimises boundary slack at zero real cost until the restricted
%   master is feasible; phase 2 drops the slacks and minimises the true cost.
%   A big-M phase 1 is deliberately AVOIDED: with big-M the boundary duals
%   saturate at +/-M, the pricing signal is destroyed and the scheme cycles,
%   adding the same columns forever (observed, 4.2% above the optimum).
%
%   Reduced cost in MATLAB's linprog dual convention
%     Lag = f'x + lam_ineq'(Ax-b) + lam_eq'(Aeq x - beq),  lam_ineq >= 0
%   so   rc = c_col + lam_dep(t_0) + lam_arr(t_L) + sum_l lam_cap(v_l,t_l),
%   and pricing is an exact min-plus DP along the path, O(L n_t^2).
%
%   Lower bound: every column carries unit mass on exactly one departure row,
%   so the column weights sum to sum(mu0) = 1 and LB = z_master + min_p rc_p
%   is rigorous (verified: always below the true optimum, closing to 1.6e-16).
    nt = S.nt;  P = S.P;  nN = S.nN;
    t0 = tic;  flag = 0;  obj = NaN;  lb = -Inf;  nRounds = 0;
    nRows = 2*nt + nN*nt;
    tlad = nan(1, numel(glad));    % wall clock to reach each DW gap level
    colI = zeros(0,1);  colJ = zeros(0,1);  ccost = zeros(0,1);  nc = 0;
    keys = containers.Map('KeyType','char','ValueType','logical');
    o = optimoptions('linprog','Algorithm','dual-simplex','Display','off', ...
                     'MaxTime',max(1,time_lim/3),'MaxIterations',3e6);

    lam = zeros(nRows,1);
    for k = 1:5000                                    % phase 1
        nRounds = nRounds + 1;
        [rc, tt, tc] = price(S, lam, 1);
        cand = find(rc < -1e-12);
        if k == 1,  [~,ord] = sort(rc,'ascend');  cand = ord(1:min(add_per_round,P));  end
        if numel(cand) > add_per_round
            [~,ord] = sort(rc(cand),'ascend');  cand = cand(ord(1:add_per_round));
        end
        nAdd = 0;
        for p = cand(:)'
            if ~isfinite(tc(p)),  continue;  end
            key = sprintf('%d:%s', p, sprintf('%d,', tt(:,p)));
            if ~isKey(keys,key)
                keys(key) = true;  nc = nc+1;
                r = col_rows(S,p,tt(:,p));
                colI = [colI; r];  colJ = [colJ; nc*ones(numel(r),1)];  %#ok<AGROW>
                ccost(nc,1) = tc(p);  nAdd = nAdd + 1;
            end
        end
        cols = sparse(colI, colJ, 1, nRows, nc);
        A1 = [cols(1:2*nt,:), speye(2*nt), -speye(2*nt)];
        A2 = [cols(2*nt+1:end,:), sparse(nN*nt, 4*nt)];
        c1 = [zeros(nc,1); ones(4*nt,1)];
        [~, s1, f1, ~, lm] = linprog(c1, A2, cap*ones(nN*nt,1), A1, ...
                                     [S.mu0; S.muT], zeros(nc+4*nt,1), [], o);
        if f1 ~= 1,  flag = -1;  break;  end
        lam = [lm.eqlin; lm.ineqlin];
        if s1 <= 1e-10,  break;  end
        if nAdd == 0,  flag = -2;  break;  end
        if toc(t0) > time_lim,  flag = -1;  break;  end
    end
    if flag < 0,  nCols = nc;  return;  end

    for k = 1:5000                                    % phase 2
        nRounds = nRounds + 1;
        cols = sparse(colI, colJ, 1, nRows, nc);
        [~, o2, f2, ~, lm] = linprog(ccost, cols(2*nt+1:end,:), cap*ones(nN*nt,1), ...
                                     cols(1:2*nt,:), [S.mu0; S.muT], zeros(nc,1), [], o);
        if f2 ~= 1,  flag = -1;  break;  end
        obj = o2;  lam = [lm.eqlin; lm.ineqlin];
        [rc, tt, tc] = price(S, lam, 2);
        lb = obj + min(min(rc), 0);
        %  A cut-off run still has a rigorous bracket, so record WHEN each
        %  accuracy level was reached.  "time to 1% gap" is a direct number;
        %  ">300 s" is not.
        if lb > 0
            gnow = (obj-lb)/abs(lb);
            for q = 1:numel(glad)
                if isnan(tlad(q)) && gnow <= glad(q),  tlad(q) = toc(t0);  end
            end
        end
        cand = find(rc < -1e-10);
        if isempty(cand),  flag = 1;  break;  end
        if numel(cand) > add_per_round
            [~,ord] = sort(rc(cand),'ascend');  cand = cand(ord(1:add_per_round));
        end
        nAdd = 0;
        for p = cand(:)'
            if ~isfinite(tc(p)),  continue;  end
            key = sprintf('%d:%s', p, sprintf('%d,', tt(:,p)));
            if ~isKey(keys,key)
                keys(key) = true;  nc = nc+1;
                r = col_rows(S,p,tt(:,p));
                colI = [colI; r];  colJ = [colJ; nc*ones(numel(r),1)];  %#ok<AGROW>
                ccost(nc,1) = tc(p);  nAdd = nAdd + 1;
            end
        end
        if nAdd == 0,  flag = 1;  break;  end
        if toc(t0) > time_lim,  flag = -1;  break;  end
    end
    nCols = nc;
end

function [rc, tt, tc] = price(S, lam, phase)
%PRICE  One exact min-plus DP per path, O(L n_t^2) each.
    nt = S.nt;  P = S.P;  L = S.L;
    lam_dep = lam(1:nt);  lam_arr = lam(nt+1:2*nt);
    lam_cap = reshape(lam(2*nt+1:end), nt, S.nN);
    rc = zeros(P,1);  tc = zeros(P,1);  tt = zeros(L+1,P);
    lin = (1:nt)';
    for p = 1:P
        g = lam_dep;  tcv = zeros(nt,1);  arg = zeros(nt,L);
        for l = 1:L
            Ce = S.CcInf{S.pcls(p,l)};
            if phase == 1,  Ko = S.ZInf;  else,  Ko = Ce;  end
            [gv, ix] = min(g + Ko, [], 1);
            g = gv(:);  ix = ix(:);
            tcv = tcv(ix) + Ce(sub2ind([nt nt], ix, lin));
            arg(:,l) = ix;
            if l < L,  g = g + lam_cap(:, S.nidx(p,l));
            else,      g = g + lam_arr;  end
        end
        [rc(p), tL] = min(g);
        tt(L+1,p) = tL;
        for l = L:-1:1,  tt(l,p) = arg(tt(l+1,p), l);  end
        tc(p) = tcv(tL);
    end
end

function r = col_rows(S, p, tt)
    nt = S.nt;  L = S.L;
    r = zeros(L+1,1);
    r(1) = tt(1);  r(2) = nt + tt(L+1);
    for l = 1:L-1,  r(2+l) = 2*nt + (S.nidx(p,l)-1)*nt + tt(l+1);  end
end

%% ─────────────────────────────────────────────────────────────────────────
function R = graph_sinkhorn(S, cap, tol, maxIter, time_lim, damp, ladder)
%GRAPH_SINKHORN  Algorithm 1 of the paper, per-edge kernels.  OURS.
%
%   Cyclic (Gauss-Seidel) block updates per cycle:
%       u <- u .* (mu0 ./ m_0)^damp
%       v <- v .* (muT ./ m_T)^damp
%       W <- min{ W .* (r ./ m_v)^damp, 1 }
%   Two corrections relative to the original code, both required:
%     (i)  the capacity update is min{W.*(r./m),1}, NOT W.*min{1,r./m}.  The
%          latter is monotonically non-increasing, so a multiplier once pushed
%          down can never be released -- not the exact block maximiser that
%          the proof of Theorem 4 needs.
%     (ii) the updates are CYCLIC.  Updating u, v, w simultaneously from one
%          sweep (Jacobi) DIVERGES on a heterogeneous instance (multiplier
%          overflow -> NaN).  F depends only on (u,W) and B only on (v,W), so
%          recomputing both after the boundary updates makes every marginal
%          exactly consistent at the cost of one extra pass.
%
%   Three message passes per cycle, each O(L P n_t^2)  (Remark 7).  The chain
%   structure is PER PATH -- the union of the P paths is not a tree, and the
%   cross-path coupling is carried by the shared nodal multipliers.
%   Stopping rule max{E_0,E_T,V} <= tol  (Remark 6).
    nt = S.nt;  P = S.P;  L = S.L;  nN = S.nN;
    TINY = 1e-300;  t0 = tic;
    u = ones(nt,1);  v = ones(nt,1);  W = ones(nt,nN);
    R = struct('obj',NaN,'it',0,'flag',0,'E0',NaN,'ET',NaN,'V',NaN, ...
               'Vstall',false,'peak',NaN,'nAct',0,'nCell',nt*nN,'util',NaN);
    R.tl = nan(1,numel(ladder));
    Vh = nan(maxIter,1);
    mv = zeros(nt,nN);

    for it = 1:maxIter
        B  = bwd(S, v, W);                     % block 1: source
        s0 = sum(B{1},2);   m0 = u .* s0;
        u  = u .* (S.mu0 ./ max(m0,TINY)).^damp;

        F  = fwd(S, u, W);                     % block 2: sink
        sT = sum(F{L},2);   mT = v .* sT;
        v  = v .* (S.muT ./ max(mT,TINY)).^damp;

        % ── block 3: interior capacities, ONE STAGE AT A TIME ────────────
        %  Each path visits exactly one node per stage, so the nodes WITHIN a
        %  stage do not enter each other's marginals -- updating a whole stage
        %  at once is exact block coordinate ascent.  ACROSS stages it is not:
        %  F_l depends on W at stages < l and B_{l+1} on W at stages > l.
        %  Updating all interior nodes simultaneously is a Jacobi step across
        %  stages and OSCILLATES once the capacity actually binds.  Measured
        %  (P=10, n_t=50, eps=0.04, tol=1e-4):
        %      kappa   all-at-once        stage-wise
        %      0.70    281 cycles         29 cycles
        %      0.50    diverges           38 cycles
        %      0.35    diverges           56 cycles
        %      0.25    diverges           60 cycles
        %  Cost is preserved: B computed once at the start of the sweep is
        %  still valid for every stage > l (those W's are untouched), and F is
        %  extended incrementally carrying the W's just updated.  So one
        %  backward plus one forward pass per cycle -- Remark 7 still holds.
        B = bwd(S, v, W);
        if isfinite(cap)
            Mu = repmat(u,1,P);
            Fl = zeros(nt,P);
            g  = S.Grp{1};
            for k = 1:numel(g),  Fl(:,g{k}.idx) = S.Kc{g{k}.c}' * Mu(:,g{k}.idx);  end
            for l = 1:L-1
                Wl   = W(:, S.nidx(:,l)');
                mvl  = (Fl .* B{l+1} .* Wl) * S.Sel{l};
                cols = unique(S.nidx(:,l));
                W(:,cols) = min(W(:,cols) .* ...
                                (cap ./ max(mvl(:,cols),TINY)).^damp, 1);
                if l < L-1
                    Mu = W(:, S.nidx(:,l)') .* Fl;      % the JUST-updated W
                    Fl = zeros(nt,P);
                    g  = S.Grp{l+1};
                    for k = 1:numel(g)
                        Fl(:,g{k}.idx) = S.Kc{g{k}.c}' * Mu(:,g{k}.idx);
                    end
                end
            end
        end

        % ── residuals from a fresh consistent pass ───────────────────────
        B = bwd(S, v, W);  F = fwd(S, u, W);  Wc = wcols(S, W);
        mv(:) = 0;
        for l = 1:L-1
            mv = mv + (F{l} .* B{l+1} .* Wc{l}) * S.Sel{l};
        end
        m0 = u .* sum(B{1},2);
        mT = v .* sum(F{L},2);
        E0 = sum(abs(m0-S.mu0));
        ET = sum(abs(mT-S.muT));
        if isfinite(cap),  Vv = sum(sum(max(mv-cap,0)));  else,  Vv = 0;  end
        R.it = it;  R.E0 = E0;  R.ET = ET;  R.V = Vv;  Vh(it) = Vv;
        res = max([E0 ET Vv]);
        for q = 1:numel(ladder)
            if isnan(R.tl(q)) && res <= ladder(q),  R.tl(q) = toc(t0);  end
        end
        if res <= tol,          R.flag =  1;  break;  end
        if toc(t0) > time_lim,  R.flag = -1;  break;  end
        %  bail out as soon as the multipliers leave the representable range
        %  rather than grinding through the whole iteration budget
        if ~all(isfinite(u)) || ~all(isfinite(v)) || ~all(isfinite(W(:)))
            R.flag = -3;  break;
        end
    end
    % infeasible instance vs slow convergence
    if R.it >= 1000
        R.Vstall = max([R.E0 R.ET]) < 10*tol && R.V > tol && ...
                   abs(R.V - Vh(R.it-500)) <= 1e-3*max(R.V, eps);
    end

    B = bwd(S,v,W);  F = fwd(S,u,W);  Wc = wcols(S,W);
    o = 0;  mass = 0;
    for l = 1:L
        if l == 1,  Lf = repmat(u,1,P);  else,  Lf = Wc{l-1} .* F{l-1};  end
        if l == L,  Rg = repmat(v,1,P);  else,  Rg = Wc{l}   .* B{l+1};  end
        g = S.Grp{l};
        for k = 1:numel(g)
            ix = g{k}.idx;  c = g{k}.c;
            o = o + sum(sum(Lf(:,ix) .* (S.CKc{c} * Rg(:,ix))));
            if l == 1
                mass = mass + sum(sum(Lf(:,ix) .* (S.Kc{c} * Rg(:,ix))));
            end
        end
    end
    %  Numerical-breakdown guard.  If the iteration DIVERGES (which happens
    %  when the capacity binds harder than the scheme can handle at this eps),
    %  the multipliers run away, u reaches ~1e300, the next message pass
    %  overflows to Inf and Inf*0 = NaN.  Measured: at KAPPA=0.5 the residual
    %  stalls at E0 ~ 0.7 and the objective comes back NaN; at KAPPA=0.7 the
    %  objective is finite for every eps down to 0.01.  So this is overflow
    %  from divergence, NOT underflow from small eps -- report it as such.
    if ~isfinite(mass) || mass < 1e-200 || mass > 1e200 || ~isfinite(o)
        R.obj  = NaN;
        R.flag = -3;
    else
        R.obj  = o / mass;
    end
    R.peak = max(mv(:));
    if isfinite(cap)
        R.nAct = nnz(mv >= 0.99*cap);
        R.util = R.peak / cap;
    end
end

function Wc = wcols(S, W)
    Wc = cell(S.L-1,1);
    for l = 1:S.L-1,  Wc{l} = W(:, S.nidx(:,l)');  end
end

function B = bwd(S, v, W)
    nt=S.nt; P=S.P; L=S.L;  Wc = wcols(S,W);  B = cell(L,1);
    M = repmat(v,1,P);
    for l = L:-1:1
        nx = zeros(nt,P);  g = S.Grp{l};
        for k = 1:numel(g),  nx(:,g{k}.idx) = S.Kc{g{k}.c} * M(:,g{k}.idx);  end
        B{l} = nx;
        if l > 1,  M = Wc{l-1} .* nx;  end
    end
end

function F = fwd(S, u, W)
    nt=S.nt; P=S.P; L=S.L;  Wc = wcols(S,W);  F = cell(L,1);
    M = repmat(u,1,P);
    for l = 1:L
        nx = zeros(nt,P);  g = S.Grp{l};
        for k = 1:numel(g),  nx(:,g{k}.idx) = S.Kc{g{k}.c}' * M(:,g{k}.idx);  end
        F{l} = nx;
        if l < L,  M = Wc{l} .* nx;  end
    end
end

%% ─────────────────────────────────────────────────────────────────────────
function report(inst, labels, rt, obj, ops, gap, stat, nvar, narc, flowlb, ...
                pprice, act, cggap, milpdev, lad, epsT, epsB, TOL_LAD, ...
                EPS_SWEEP, TIME_LIM, TOL, MEM_LIM, capcost, cg_lad, nsproj, CG_LAD)
    nInst = size(inst,1);  nM = numel(labels);  Wc = 17;
    sep = repmat('-', 1, 33 + nInst*Wc);
    hdr = @(ttl) fprintf('\n\n%s\n%s\n%s\n%-33s', sep, ttl, sep, 'Method');

    hdr('RUNTIME (s)   [operation count in brackets]');
    for si=1:nInst, fprintf('%*s',Wc,sprintf('(%d,%d)',inst(si,1),inst(si,2))); end
    fprintf('\n%s\n', sep);
    for m = 1:nM
        fprintf('%-33s', labels{m});
        for si = 1:nInst
            ex = '';
            if m == 4 && ~isnan(cggap(si)),  ex = sprintf('(%.2f%% gap)', 100*cggap(si));  end
            if m == 3 && ~isnan(nsproj(si)), ex = sprintf('~%.1fh proj', nsproj(si)/3600); end
            fprintf('%*s', Wc, cellstr_(rt(m,si), ops(m,si), stat{m,si}, TIME_LIM, ex));
        end
        fprintf('\n');
    end
    fprintf('%s\n', sep);

    hdr('OBJECTIVE  <C,Pi>');
    for si=1:nInst, fprintf('%*s',Wc,sprintf('(%d,%d)',inst(si,1),inst(si,2))); end
    fprintf('\n%s\n', sep);
    for m = 1:nM
        fprintf('%-33s', labels{m});
        for si = 1:nInst
            if isnan(obj(m,si)), fprintf('%*s',Wc,'--');
            else,                fprintf('%*.6f',Wc,obj(m,si));  end
        end
        fprintf('\n');
    end
    fprintf('%s\n', sep);

    fprintf('\n%-33s','path LP #var');
    for si=1:nInst, fprintf('%*s',Wc,num_(nvar(si))); end
    fprintf('\n%-33s','flow model #arcs (P-indep.)');
    for si=1:nInst, fprintf('%*s',Wc,num_(narc(si))); end
    fprintf('\n%-33s','flow optimum (lower bound)');
    for si=1:nInst
        if isnan(flowlb(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*.6f',Wc,flowlb(si)); end
    end
    fprintf('\n%-33s','price of path restriction');
    for si=1:nInst
        if isnan(pprice(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*s',Wc,sprintf('%+.1f%%',100*pprice(si))); end
    end
    fprintf('\n%-33s','transport-cost change');
    for si=1:nInst
        if isnan(capcost(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*s',Wc,sprintf('%+.2f%%',100*capcost(si))); end
    end
    fprintf('\n%-33s','capacity ACTIVE cells');
    for si=1:nInst
        if isnan(act(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*s',Wc,sprintf('%.2f%%',100*act(si))); end
    end
    fprintf('\n%-33s','MILP cert. max |flow-int|');
    for si=1:nInst
        if isnan(milpdev(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*.2e',Wc,milpdev(si)); end
    end
    fprintf('\n%-33s','DW gap when cut off');
    for si=1:nInst
        if isnan(cggap(si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*s',Wc,sprintf('%.2f%%',100*cggap(si))); end
    end
    fprintf('\n%-33s','entropic bias (ours)');
    for si=1:nInst
        if isnan(gap(6,si)), fprintf('%*s',Wc,'--');
        else, fprintf('%*s',Wc,sprintf('%+.2f%%',100*gap(6,si))); end
    end
    fprintf('\n');

    fprintf('\n\n%s\nCOLUMN GENERATION: TIME (s) TO REACH A GIVEN DW OPTIMALITY GAP\n%s\n', sep, sep);
    for q = 1:numel(CG_LAD)
        fprintf('gap <= %-6.3g   ', CG_LAD(q));
        for si = 1:nInst
            vv = cg_lad(si,q);
            if isnan(vv), fprintf('%*s',Wc,'not reached');
            else,         fprintf('%*.2f',Wc,vv);  end
        end
        fprintf('\n');
    end
    fprintf('%s\n', sep);

    fprintf('\n\n%s\nGRAPHICAL SINKHORN: TIME (s) TO REACH max{E_0,E_T,V} <= tau\n%s\n', sep, sep);
    for q = 1:numel(TOL_LAD)
        fprintf('tau = %-8.0e   ', TOL_LAD(q));
        for si = 1:nInst
            vv = lad(si,q);
            if isnan(vv), fprintf('%*s',Wc,'not reached');
            else,         fprintf('%*.2f',Wc,vv);  end
        end
        fprintf('\n');
    end
    fprintf('%s\n', sep);

    if any(~isnan(epsB(:)))
        fprintf('\n\n%s\nACCURACY / TIME FRONTIER for graphical Sinkhorn\n%s\n', sep, sep);
        for si = 1:nInst
            if all(isnan(epsB(si,:))), continue; end
            fprintf('  (P,n_t) = (%d,%d)\n', inst(si,1), inst(si,2));
            fprintf('     %-10s','eps');
            for q=1:numel(EPS_SWEEP), fprintf('%12.3f',EPS_SWEEP(q)); end
            fprintf('\n     %-10s','time (s)');
            for q=1:numel(EPS_SWEEP), fprintf('%12.2f',epsT(si,q)); end
            fprintf('\n     %-10s','bias (%)');
            for q=1:numel(EPS_SWEEP), fprintf('%12.3f',100*epsB(si,q)); end
            fprintf('\n');
        end
        fprintf('%s\n', sep);
    end

    %% LaTeX
    tex = {'LP reference (industrial)$^{\ast}$', ...
           'Time-expanded MILP$^{\S}$', ...
           'Network simplex$^{\flat}$', ...
           'Column generation$^{\dagger}$', ...
           'Standard MM Sinkhorn', ...
           '\textbf{Graphical Sinkhorn (ours)}'};
    fprintf('\n\n%%%% --- LaTeX table ---\n');
    fprintf('\\begin{tabular}{l%s}\n\\toprule\n', repmat('r',1,nInst));
    fprintf('& \\multicolumn{%d}{c}{$(P,\\,n_t)$}\\\\\n\\cmidrule(lr){2-%d}\n', nInst, nInst+1);
    fprintf('Method');
    for si=1:nInst, fprintf(' & $(%d,\\;%d)$', inst(si,1), inst(si,2)); end
    fprintf(' \\\\\n\\midrule\n');
    for m = 1:nM
        fprintf('%s', tex{m});
        for si = 1:nInst
            fprintf(' & %s', texcell_(rt(m,si), ops(m,si), stat{m,si}, TIME_LIM, ...
                                      m, cggap(si), milpdev(si)));
        end
        fprintf(' \\\\\n');
        if m == 4,  fprintf('\\midrule\n');  end
    end
    fprintf('\\bottomrule\n\\end{tabular}\n');
    fprintf(['%%%%\n' ...
      '%%%% Baselines are TEXTBOOK implementations with no tuning; the proposed\n' ...
      '%%%% method exploits the chain structure of the path cost.  All MATLAB,\n' ...
      '%%%% same machine, shared tolerance %g on max{E_0,E_T,V}, shared %d s\n' ...
      '%%%% wall-clock limit applied to the proposed method as well.  Operation\n' ...
      '%%%% counts are in brackets and are language independent.  Commercial\n' ...
      '%%%% network-simplex codes are considerably faster than a textbook one;\n' ...
      '%%%% this comparison isolates algorithmic structure, not engineering.\n' ...
      '%%%% * industrial solver, ground truth for the entropic bias only.\n' ...
      '%%%% S not solved: the node-split time-expanded flow matrix is totally\n' ...
      '%%%%   unimodular, so with integer data the LP relaxation is already\n' ...
      '%%%%   integral (certificate: max |flow - integer| reported).\n' ...
      '%%%% b network simplex solves the AGGREGATED flow model (all DAG paths):\n' ...
      '%%%%   a relaxation of Problem 3, hence a lower bound; its size does not\n' ...
      '%%%%   depend on P.\n' ...
      '%%%% dagger two-phase Dantzig-Wolfe, columns = (path, causal time tuple).\n' ...

      '%%%% --- exceeds the %.0f GB construction limit.\n'], ...
      TOL, TIME_LIM, MEM_LIM/1e9);
end

function s = cellstr_(rt, op, st, TL, extra)
    if nargin < 5,  extra = '';  end
    switch st
        case 'mem',    s = '---';
        case 'n/a',    s = 'n/a';
        case 'infeas', s = 'infeas';
        case 'iter',   s = 'no conv';
        case 'skip',   s = 'not needed';
        case 'proj',   s = extra;      % calibrated projection, e.g. "~5.6h proj"
        case 'brkdn',  s = 'breakdown';
        case 'err',    s = 'err';
        case 'time'
            if isempty(extra),  s = sprintf('%.0fs (cut off)', rt);
            else,               s = sprintf('%.0fs %s', rt, extra);  end
        case 'ok'
            if ~isnan(op),  s = sprintf('%.2fs [%d]', rt, op);
            else,           s = sprintf('%.2fs', rt);  end
        otherwise,     s = '--';
    end
end

function s = texcell_(rt, op, st, TL, m, cgg, mdev)
    if m == 2 && strcmp(st,'ok')
        s = sprintf('int.\\ ($%.0e$)', mdev);  return
    end
    switch st
        case 'mem',    s = '---';
        case 'n/a',    s = 'n/a';
        case 'infeas', s = 'infeas.';
        case 'iter',   s = 'n.c.';
        case 'skip',   s = '--';
        case 'proj',   s = '---';
        case 'brkdn',  s = 'brkdn';
        case 'err',    s = 'err';
        case 'time'
            if m == 4 && ~isnan(cgg)
                s = sprintf('$>\\!%d$\\,s ($%.1f\\%%$)', TL, 100*cgg);
            else
                s = sprintf('$>\\!%d$\\,s', TL);
            end
        case 'ok'
            if ~isnan(op),  s = sprintf('$%.2f$\\,s [%d]', rt, op);
            else,           s = sprintf('$%.2f$\\,s', rt);  end
        otherwise,     s = '--';
    end
end

function s = numstr_(x)
    if isnan(x),  s = '';  else,  s = sprintf('%.6g', x);  end
end

function s = num_(n)
    if isnan(n),      s = '--';
    elseif n >= 1e9,  s = sprintf('%.1fB', n/1e9);
    elseif n >= 1e6,  s = sprintf('%.1fM', n/1e6);
    else,             s = sprintf('%.0fk', n/1e3);
    end
end

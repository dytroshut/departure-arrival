%% complex_graph_final_v10.m — final large-network visualization
%
%  Uses the same inverse travel-time Gibbs kernels and cyclic stage-wise
%  graphical Sinkhorn updates as da_table_final.m. The constrained and
%  unconstrained rows share one absolute color normalization. The convergence
%  figure reports the empirical dual-objective gap together with a
%  geometric upper bound on one logarithmic axis.

clear; clc; close all;
rng(42);

%% ── Network ───────────────────────────────────────────────────────────────
nStages = 10;
nWidth  = 10;
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

%% ── Edge list and heterogeneous edge weights ─────────────────────────────
Ledge = nStages + 1;
eu = zeros(1,0);  ev = zeros(1,0);
for uu = 1:nNodes
    for vv = adj{uu}
        eu(end+1) = uu;  ev(end+1) = vv;  %#ok<AGROW>
    end
end
nE  = numel(eu);
eid = sparse(eu, ev, 1:nE, nNodes, nNodes);

nWclass = 8;
w_ref   = 1/(Ledge^2);
w_vals  = w_ref * linspace(0.6,1.4,nWclass);
e_class = randi(nWclass,1,nE);

%% ── Generate 120 unique random paths ─────────────────────────────────────
nPaths = 120;
paths  = {};
seen   = containers.Map('KeyType','char','ValueType','logical');
ntry   = 0;
while numel(paths) < nPaths && ntry < 200000
    ntry = ntry + 1;
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
            seen(key)    = true;
            paths{end+1} = pth;  %#ok<AGROW>
        end
    end
end
nPaths = numel(paths);
fprintf('Generated %d unique paths  (%d attempts)\n', nPaths, ntry);

middles = unique([paths{:}]);
middles = middles(middles ~= src & middles ~= snk);
fprintf('Active interior nodes: %d / %d\n\n', numel(middles), nInner);

%% ── Time grid and feasible boundary marginals ─────────────────────────
nt = 50;
t  = linspace(0,1,nt)';
G  = @(m,s) exp(-0.5*((t-m)/s).^2);
mu0 = 0.9*G(0.10,0.05) + 0.6*G(0.22,0.05);
muT = 0.8*G(0.76,0.05) + 0.7*G(0.88,0.05);

%% ── Heterogeneous inverse travel-time Gibbs kernels ──────────────────────
eps_ent = 0.01;
dt_min  = 0.035;
[U,V]   = ndgrid(t,t);
D       = V-U;
kmask   = D >= dt_min;

Kclass  = cell(nWclass,1);
CKclass = cell(nWclass,1);
for cc = 1:nWclass
    Cc = zeros(nt);
    Cc(kmask) = w_vals(cc)./D(kmask);
    Kclass{cc} = zeros(nt);
    Kclass{cc}(kmask) = exp(-Cc(kmask)/eps_ent);
    CKclass{cc} = Cc.*Kclass{cc};
end

% Retain only boundary grid cells that can be joined by Ledge causal edges.
reach = kmask;
for ell = 2:Ledge
    reach = (double(reach)*double(kmask)) > 0;
end
mu0 = mu0 .* (sum(reach,2)>0);
muT = muT .* (sum(reach,1)'>0);
mu0 = mu0/sum(mu0);
muT = muT/sum(muT);

% Edge-weight class encountered by each path at every edge position.
pcls = zeros(nPaths,Ledge);
for pid = 1:nPaths
    pth = paths{pid};
    for ell = 1:Ledge
        pcls(pid,ell) = e_class(full(eid(pth(ell),pth(ell+1))));
    end
end

%% ── Correct cyclic graphical Sinkhorn solves ─────────────────────────────
KAPPA = 0.25;
opts = struct('eps0',1e-300,'damp',1.0,'maxIter',50000,'tol',1e-6);

% Solve the unconstrained problem first and use its peak nodal cell mass to
% define an active but feasible capacity for the constrained visualization.
fprintf('Solving UNCONSTRAINED (nt=%d, P=%d) ...\n',nt,nPaths);
tic;
[rib_free,m0_free,mT_free,nIt_free,diag_free] = ...
    runPathSinkhorn(Kclass,CKclass,pcls,mu0,muT, ...
        paths,middles,nt,opts,Inf,eps_ent);
t_free = toc;

peak_free = 0;
for vtx = middles
    peak_free = max(peak_free,max(rib_free(vtx)));
end
cap_tight = KAPPA*peak_free;
fprintf('  Done %.2f s (%d cycles), peak %.6f, cost %.6f\n', ...
    t_free,nIt_free,peak_free,diag_free.cost);
fprintf('  Capacity = %.2f x peak = %.6f per temporal cell\n\n', ...
    KAPPA,cap_tight);

fprintf('Solving CONSTRAINED (cap=%.6f, nt=%d, P=%d) ...\n', ...
    cap_tight,nt,nPaths);
tic;
[rib_cap,m0_cap,mT_cap,nIt_cap,diag_cap] = ...
    runPathSinkhorn(Kclass,CKclass,pcls,mu0,muT, ...
        paths,middles,nt,opts,cap_tight,eps_ent);
t_cap = toc;
cost_change = 100*(diag_cap.cost/diag_free.cost-1);
fprintf('  Done %.2f s (%d cycles), cost %.6f (%+.2f%% vs free)\n', ...
    t_cap,nIt_cap,diag_cap.cost,cost_change);
fprintf('  Capacity active at %d/%d node-time cells (%.2f%%)\n\n', ...
    diag_cap.nActive,diag_cap.nCell,100*diag_cap.nActive/diag_cap.nCell);

%% ── Node positions ────────────────────────────────────────────────────────
xy = zeros(nNodes, 2);
for s = 1:nStages
    for pp = 1:nWidth
        xy(node_of(s,pp),:) = [ s/(nStages+1),  (pp-1)/(nWidth-1) ];
    end
end
xy(src,:) = [0, 0.5];
xy(snk,:) = [1, 0.5];

%% ── Flux snapshots ────────────────────────────────────────────────────────
nSnap  = 6;
t_snap = linspace(0, 1, nSnap);

fluxCap = zeros(nt, nNodes);
for vtx = middles
    if isKey(rib_cap, vtx),  fluxCap(:,vtx) = rib_cap(vtx);  end
end
fluxCap(:,src) = m0_cap;   fluxCap(:,snk) = mT_cap;
snapCap = zeros(nSnap, nNodes);
for vi = 1:nNodes
    snapCap(:,vi) = max(0, interp1(t, fluxCap(:,vi), t_snap','linear',0));
end

fluxFree = zeros(nt, nNodes);
for vtx = middles
    if isKey(rib_free, vtx),  fluxFree(:,vtx) = rib_free(vtx);  end
end
fluxFree(:,src) = m0_free;  fluxFree(:,snk) = mT_free;
snapFree = zeros(nSnap, nNodes);
for vi = 1:nNodes
    snapFree(:,vi) = max(0, interp1(t, fluxFree(:,vi), t_snap','linear',0));
end

%% ── Edge geometry ─────────────────────────────────────────────────────────
Xe=[]; Ye=[]; Ze=[];
for kk = 1:nNodes
    for nb = adj{kk}(:)'
        Xe = [Xe,  xy(kk,1), xy(nb,1), NaN];  %#ok<AGROW>
        Ye = [Ye,  xy(kk,2), xy(nb,2), NaN];  %#ok<AGROW>
        Ze = [Ze,  0,        0,        NaN];  %#ok<AGROW>
    end
end

%% ── Visual parameters ─────────────────────────────────────────────────────
cSrc  = [0.06  0.28  0.92];
cSnk  = [0.92  0.04  0.04];
cEdge = [0.38  0.38  0.38];

thr       = 5e-4;
gainZ     = 2.5;
tubeW_i   = 0.015;
tubeW_ss  = 0.058;
alph_i    = 0.52;
alph_ss   = 0.92;

%% ── Global normalisations ─────────────────────────────────────────────────
maxInt_cap  = max(max(fluxCap(:, middles)))  + 1e-16;
maxSrc_cap  = max(m0_cap)  + 1e-16;
maxSnk_cap  = max(mT_cap)  + 1e-16;

maxInt_free = max(max(fluxFree(:, middles))) + 1e-16;
maxSrc_free = max(m0_free) + 1e-16;
maxSnk_free = max(mT_free) + 1e-16;

maxInt_global = max(maxInt_cap, maxInt_free);
maxSrc_global = max(maxSrc_cap, maxSrc_free);
maxSnk_global = max(maxSnk_cap, maxSnk_free);
% One absolute height reference for every node and both rows. This preserves
% the physical ratios between source, interior, and sink masses.
maxHeight_global = max([maxInt_global,maxSrc_global,maxSnk_global]);
capZ_global      = cap_tight / maxHeight_global * gainZ;

fprintf('Global colour scale:  interior max = %.4f  (cap=%.4f, free=%.4f)\n', ...
    maxInt_global, maxInt_cap, maxInt_free);

%% ── Figure 1: combined 2-row flux figure, shared colorbar ────────────────
drawCombinedFig(1, snapCap, snapFree, ...
    maxInt_global, maxHeight_global, ...
    Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
    thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ...
    cSrc, cSnk, cEdge, capZ_global, ...
    {sprintf('Constrained $(r{=}%.4f)$', cap_tight), ...
     'Unconstrained'});

%% ── Figure 2: standalone network topology ────────────────────────────────
drawTopologyFig(2, adj, xy, middles, src, snk, nNodes, nStages, nWidth, node_of);

%% ── Figure 3: dual convergence ──────────────────────────────────────────
% The dual objective is the diagnostic directly associated with the exact
% block-coordinate ascent argument in Theorem 5.
figure(3); clf;
set(gcf,'Color','w','Units','centimeters','Position',[5 5 11.5 8], ...
    'PaperPositionMode','auto');

ax = axes('Parent',gcf,'Position',[0.16 0.17 0.80 0.77]);
iters = (0:nIt_cap-1)';
plot_floor = 1e-14;

% For exact Gauss--Seidel block maximization, the dual objective should be
% nondecreasing, up to floating-point roundoff. Report any material decrease
% rather than concealing it in the plotted diagnostic.
Dhist = diag_cap.dual(:);
Dref  = max(Dhist);
dual_scale = max(1,abs(Dref));
dual_drop  = min(diff(Dhist));
if dual_drop < -1e-10*dual_scale
    warning(['The recorded dual objective is not monotone: ' ...
        'minimum cycle increment = %.3e.'],dual_drop);
end

% Dref is the largest computed dual value and estimates D^star.
dual_gap = Dref-Dhist;
if dual_gap(1) <= 0
    error('The initial estimated dual gap is not positive.');
end
dual_rel = dual_gap/dual_gap(1);
dual_rel(dual_rel < plot_floor) = plot_floor;

% Solid: computed dual gap. Dashed: empirical geometric upper bound.
hD = semilogy(ax,iters,dual_rel,'-', ...
    'Color',[0.92 0.04 0.04],'LineWidth',2.0);
hold(ax,'on');
[rhoD,ubD] = linearUpperBoundFromZero(iters,dual_rel);
if isfinite(rhoD) && rhoD < 1
    semilogy(ax,iters,ubD,'--', ...
        'Color',[0.55 0.00 0.00],'LineWidth',1.6, ...
        'HandleVisibility','off');
end

xlabel(ax,'Sinkhorn cycle $k$','Interpreter','latex','FontSize',13);
ylabel(ax,'Normalized dual gap','Interpreter','latex','FontSize',13);
legend(ax,hD,{sprintf('Dual gap, $\\rho_D=%.3f$',rhoD)}, ...
    'Interpreter','latex','FontSize',10.5,'Location','southwest','Box','on');
grid(ax,'on'); box(ax,'on');
set(ax,'TickLabelInterpreter','latex','FontSize',12,'GridAlpha',0.25, ...
    'XLim',[0,max(iters(end),1)],'YLim',[plot_floor 1]);

fprintf('  Empirical dual upper-bound rate: rho_D = %.4f\n',rhoD);

%% ── Scalability report ────────────────────────────────────────────────────
fprintf('\n── Scalability Report ──────────────────────────────\n');
fprintf('  Interior nodes : %d\n', nInner);
fprintf('  Paths          : %d\n', nPaths);
fprintf('  nt             : %d\n', nt);
fprintf('  Constrained    : %d iter,  %.2f s\n', nIt_cap,  t_cap);
fprintf('  Unconstrained  : %d iter,  %.2f s\n', nIt_free, t_free);
fprintf('  Free cost      : %.6f\n', diag_free.cost);
fprintf('  Constr. cost   : %.6f  (%+.2f%%)\n',diag_cap.cost,cost_change);
fprintf('────────────────────────────────────────────────────\n');


%% ======================================================================
%% Local functions
%% ======================================================================

function drawCombinedFig(figNo, snapFlux1, snapFlux2, ...
        maxInt, maxHeight, ...
        Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
        thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ...
        cSrc, cSnk, cEdge, capZ, rowLabels)
%DRAWCOMBINEDFIG  2×nSnap panel figure.  Panels are flush (gap ≈ 0) so the
%   figure fits compactly in a paper column.  A single colorbar on the right
%   is valid for both rows because maxInt is the GLOBAL maximum.

    figure(figNo); clf;
    set(gcf,'Color','w','Units','centimeters', ...
        'Position',[1, 2, 36, 16], ...
        'PaperPositionMode','auto','Renderer','painters');

    cmap = getVizMap();
    colormap(gcf, cmap);

    % ── Layout: panels flush, no gap ─────────────────────────────────────
    pad_l = 0.060;   % row labels
    pad_r = 0.090;   % colorbar
    pad_t = 0.11;    % column titles
    pad_b = 0.02;
    gap_w = 0;       % ← flush horizontally
    gap_h = 0.004;   % ← near-zero vertical (hairline between rows)

    w_p   = (1 - pad_l - pad_r - gap_w*(nSnap-1)) / nSnap;
    h_p   = (1 - pad_t - pad_b - gap_h) / 2;
    y_top = pad_b + h_p + gap_h;
    y_bot = pad_b;

    for row = 1:2
        if row == 1
            snapFlux = snapFlux1;
            y_row    = y_top;
        else
            snapFlux = snapFlux2;
            y_row    = y_bot;
        end

        for i = 1:nSnap
            x_pos = pad_l + (i-1)*w_p;   % gap_w = 0 → no offset
            ax = axes('Parent',gcf, ...
                'Position',[x_pos, y_row, w_p, h_p]);  %#ok<LAXES>
            hold(ax,'on');
            clim(ax,[0 1]);

            % Keep only the links as a quiet topology reference. The dense
            % gray node and capacity markers are omitted for visual clarity.
            plot3(ax, Xe, Ye, Ze, 'Color',cEdge,'LineWidth',0.70);

            for jj = 1:numel(middles)
                vi   = middles(jj);
                mass = snapFlux(i,vi);
                if mass >= thr*maxHeight
                    hv   = mass/maxHeight*gainZ;
                    c_vi = min(1,mass/maxInt);
                    drawTubeMapped(ax,xy(vi,1),xy(vi,2), ...
                        hv,tubeW_i,c_vi,alph_i);
                end
            end

            mass0 = snapFlux(i,src);
            if mass0 >= thr*maxHeight
                h0 = mass0/maxHeight*gainZ;
                drawTubeFixed(ax,xy(src,1),xy(src,2), ...
                    h0,tubeW_ss,cSrc,alph_ss);
            end

            massT = snapFlux(i,snk);
            if massT >= thr*maxHeight
                hT = massT/maxHeight*gainZ;
                drawTubeFixed(ax,xy(snk,1),xy(snk,2), ...
                    hT,tubeW_ss,cSnk,alph_ss);
            end

            % Column title on top row only — LARGER font
            if row == 1
                title(ax, sprintf('$t=%d/5$', i-1), ...
                    'Interpreter','latex','FontSize',12,'FontWeight','normal');
            end

            set(ax,'Visible','off');
            if row == 1, ax.Title.Visible = 'on'; end
            set(ax,'XLim',[-0.02 1.02],'YLim',[-0.03 1.03],'ZLim',[0 gainZ]);
            set(ax,'PlotBoxAspectRatio',[1.05 1.08 gainZ*0.68]);
            view(ax, 35, 28);
            hold(ax,'off');
        end

        % Row label — LARGER font
        lax = axes('Parent',gcf, ...
            'Position',[0.002, y_row, 0.055, h_p],'Visible','off');
        text(lax, 0.5, 0.5, rowLabels{row}, ...
            'Units','normalized','Rotation',90, ...
            'Interpreter','latex','FontSize',12, ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end

    % ── Single shared colorbar, LARGER labels ─────────────────────────────
    cb_x = 1 - pad_r + 0.020;
    cb_h = h_p * 2 + gap_h;
    cbax = axes('Parent',gcf, ...
        'Position',[cb_x, y_bot, 0.001, cb_h],'Visible','off');
    clim(cbax,[0 1]);
    colormap(cbax, cmap);
    cb = colorbar(cbax,'Position',[cb_x, y_bot, 0.022, cb_h]);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize             = 11;
    cb.Label.Interpreter    = 'latex';
    cb.Label.String         = '$m_k(t_i)/\!\max_{k,i}\,m_k(t_i)$';
    cb.Label.FontSize       = 12;
    cb.Ticks                = 0:0.2:1;
end

% ──────────────────────────────────────────────────────────────────────────

function drawTopologyFig(figNo, adj, xy, middles, src, snk, nNodes, nStages, nWidth, node_of)
%DRAWTOPOLOGYFIG  Standalone 2-D DAG topology figure for the paper.
%   Nodes coloured by stage (parula), source blue, sink red.
%   Edges are thin grey lines.  One colorbar shows the stage index.

    figure(figNo); clf;
    set(gcf,'Color','w','Units','centimeters', ...
        'Position',[42, 2, 18, 10], ...
        'PaperPositionMode','auto');

    ax = axes('Parent',gcf,'Color','w','Position',[0.04 0.04 0.82 0.92]);
    hold(ax,'on');

    % Edges — black
    for kk = 1:nNodes
        for nb = adj{kk}(:)'
            plot(ax, [xy(kk,1), xy(nb,1)], [xy(kk,2), xy(nb,2)], ...
                '-','Color',[0 0 0],'LineWidth',0.5);
        end
    end

    % Interior nodes coloured by stage
    all_x = zeros(1, nStages*nWidth);
    all_y = zeros(1, nStages*nWidth);
    all_s = zeros(1, nStages*nWidth);
    idx = 0;
    for s = 1:nStages
        for pp = 1:nWidth
            idx = idx + 1;
            vi = node_of(s,pp);
            all_x(idx) = xy(vi,1);
            all_y(idx) = xy(vi,2);
            all_s(idx) = s;
        end
    end
    scatter(ax, all_x, all_y, 28, all_s, 'filled', 'MarkerEdgeColor','none');
    colormap(ax, parula(nStages));
    clim(ax, [0.5, nStages + 0.5]);

    % Source v_0
    scatter(ax, xy(src,1), xy(src,2), 100, [0.06 0.28 0.92], ...
        'filled','MarkerEdgeColor',[0.02 0.12 0.55],'LineWidth',1.2);
    text(ax, xy(src,1) - 0.04, xy(src,2), '$v_0$', ...
        'Interpreter','latex','FontSize',13,'FontWeight','bold', ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color',[0.06 0.28 0.92]);

    % Sink v_T
    scatter(ax, xy(snk,1), xy(snk,2), 100, [0.92 0.04 0.04], ...
        'filled','MarkerEdgeColor',[0.55 0.02 0.02],'LineWidth',1.2);
    text(ax, xy(snk,1) + 0.04, xy(snk,2), '$v_{\mathcal{T}}$', ...
        'Interpreter','latex','FontSize',13,'FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'Color',[0.92 0.04 0.04]);

    % Stage x-tick labels
    for s = 1:nStages
        text(ax, s/(nStages+1), -0.09, sprintf('%d', s), ...
            'Interpreter','latex','FontSize',13, ...
            'HorizontalAlignment','center','Color',[0.35 0.35 0.35]);
    end
    text(ax, 0.5, -0.17, 'Stage', ...
        'Interpreter','latex','FontSize',12, ...
        'HorizontalAlignment','center','Color',[0.25 0.25 0.25]);

    set(ax,'Visible','off');
    set(ax,'XLim',[-0.10 1.10],'YLim',[-0.22 1.10]);

    hold(ax,'off');
end

% ──────────────────────────────────────────────────────────────────────────

function drawTubeMapped(ax, x0, y0, h, w, c_val, alph)
    if h <= 0, h = 1e-4; end
    cmap = getVizMap();
    nC   = size(cmap,1);
    idx  = max(1, min(nC, round(c_val*(nC-1))+1));
    col  = cmap(idx,:);
    verts = [x0-w, y0-w, 0; x0+w, y0-w, 0; x0+w, y0+w, 0; x0-w, y0+w, 0;
             x0-w, y0-w, h; x0+w, y0-w, h; x0+w, y0+w, h; x0-w, y0+w, h];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax,'Vertices',verts,'Faces',faces, ...
        'FaceColor',col,'FaceAlpha',alph, ...
        'EdgeColor',col*0.60,'EdgeAlpha',0.30,'LineWidth',0.4);
end

function drawTubeFixed(ax, x0, y0, h, w, col, alph)
    if h <= 0, h = 1e-4; end
    verts = [x0-w, y0-w, 0; x0+w, y0-w, 0; x0+w, y0+w, 0; x0-w, y0+w, 0;
             x0-w, y0-w, h; x0+w, y0-w, h; x0+w, y0+w, h; x0-w, y0+w, h];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax,'Vertices',verts,'Faces',faces, ...
        'FaceColor',col,'FaceAlpha',alph, ...
        'EdgeColor',col*0.60,'EdgeAlpha',0.45,'LineWidth',0.5);
end

function cmap = getVizMap()
    raw  = parula(256);
    cmap = raw(round(0.30*256)+1 : end, :);
end

function drawCapMarker(ax, x0, y0, z0, w)
    if z0 <= 0, return; end
    verts = [x0-w, y0-w, z0; x0+w, y0-w, z0;
             x0+w, y0+w, z0; x0-w, y0+w, z0];
    patch(ax,'Vertices',verts,'Faces',[1 2 3 4], ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.30, ...
        'EdgeColor',[0.45 0.45 0.45],'EdgeAlpha',0.55,'LineWidth',0.5);
end

function [ribbons,m0,mT,nIter,diag] = ...
        runPathSinkhorn(Kclass,CKclass,pcls,p,q,paths,middles,nt,opts,cap,eps_ent)
%RUNPATHSINKHORN Correct cyclic graphical Sinkhorn for the layered DAG.
%   Boundary blocks are updated first. Capacity blocks are then updated one
%   stage at a time using the latest preceding-stage multipliers. Nodes within
%   one stage are separable because every path visits exactly one such node.

    tiny    = opts.eps0;
    damp    = opts.damp;
    maxIter = opts.maxIter;
    tol     = opts.tol;

    nP      = numel(paths);
    nEdge   = size(pcls,2);
    nStage  = nEdge-1;
    nNode   = max(middles);

    stage_node = zeros(nP,nStage);
    for pid = 1:nP
        stage_node(pid,:) = paths{pid}(2:end-1);
    end

    Sel = cell(nStage,1);
    for ell = 1:nStage
        Sel{ell} = sparse((1:nP)',stage_node(:,ell),1,nP,nNode);
    end

    Grp = cell(nEdge,1);
    for ell = 1:nEdge
        cls = unique(pcls(:,ell))';
        gg  = cell(numel(cls),1);
        for kk = 1:numel(cls)
            gg{kk} = struct('c',cls(kk), ...
                'idx',find(pcls(:,ell)==cls(kk))');
        end
        Grp{ell} = gg;
    end

    if isscalar(cap)
        Rcap = cap*ones(nt,nNode);
    else
        Rcap = cap;
    end

    u = ones(nt,1);
    v = ones(nt,1);
    W = ones(nt,nNode);

    E0_hist   = nan(maxIter,1);
    ET_hist   = nan(maxIter,1);
    V_hist    = nan(maxIter,1);
    dual_hist = nan(maxIter,1);
    mv        = zeros(nt,nNode);
    nIter     = maxIter;

    for it = 1:maxIter
        % Source block using the current sink and capacity scalings.
        B  = backwardMessages(Kclass,Grp,v,W,stage_node);
        m0 = u .* sum(B{1},2);
        u  = u .* (p./max(m0,tiny)).^damp;

        % Sink block using the newly updated source scaling.
        F  = forwardMessages(Kclass,Grp,u,W,stage_node);
        mT = v .* sum(F{nEdge},2);
        v  = v .* (q./max(mT,tiny)).^damp;

        % Capacity blocks in forward stage order. The backward messages stay
        % valid because later-stage capacity blocks have not yet been changed.
        if isfinite(cap)
            B  = backwardMessages(Kclass,Grp,v,W,stage_node);
            Mu = repmat(u,1,nP);
            Fl = applyForwardKernel(Kclass,Grp{1},Mu,nt,nP);

            for ell = 1:nStage
                Wl   = W(:,stage_node(:,ell)');
                mvl  = (Fl.*B{ell+1}.*Wl)*Sel{ell};
                cols = unique(stage_node(:,ell));

                % Exact upper-bound marginal block update.
                W(:,cols) = min(W(:,cols).* ...
                    (Rcap(:,cols)./max(mvl(:,cols),tiny)).^damp,1);

                if ell < nStage
                    Mu = W(:,stage_node(:,ell)').*Fl;
                    Fl = applyForwardKernel( ...
                        Kclass,Grp{ell+1},Mu,nt,nP);
                end
            end
        end

        % Residuals evaluated from one fresh, internally consistent pass.
        B  = backwardMessages(Kclass,Grp,v,W,stage_node);
        F  = forwardMessages(Kclass,Grp,u,W,stage_node);
        Wc = weightColumns(W,stage_node);
        mv(:) = 0;
        for ell = 1:nStage
            mv = mv + (F{ell}.*B{ell+1}.*Wc{ell})*Sel{ell};
        end

        m0 = u .* sum(B{1},2);
        mT = v .* sum(F{nEdge},2);
        E0 = sum(abs(m0-p));
        ET = sum(abs(mT-q));
        if isfinite(cap)
            Vcap = sum(sum(max(mv(:,middles)-Rcap(:,middles),0)));
        else
            Vcap = 0;
        end

        E0_hist(it) = E0;
        ET_hist(it) = ET;
        V_hist(it)  = Vcap;

        % Concave dual objective in scaling variables. It is the appropriate
        % quantity for diagnosing the linear convergence of block-coordinate
        % ascent. The capacity term is absent in the unconstrained solve.
        dual_val = eps_ent*(dot(p,log(max(u,tiny))) + ...
            dot(q,log(max(v,tiny))) - sum(m0));
        if isfinite(cap)
            dual_val = dual_val + eps_ent*sum(sum( ...
                Rcap(:,middles).*log(max(W(:,middles),tiny))));
        end
        dual_hist(it) = dual_val;
        if max([E0,ET,Vcap]) <= tol
            nIter = it;
            break
        end

        if ~all(isfinite(u)) || ~all(isfinite(v)) || ...
                ~all(isfinite(W(:)))
            error('Graphical Sinkhorn encountered non-finite scalings.');
        end
    end

    E0_hist   = E0_hist(1:nIter);
    ET_hist   = ET_hist(1:nIter);
    V_hist    = V_hist(1:nIter);
    dual_hist = dual_hist(1:nIter);

    if max([E0,ET,Vcap]) > tol
        warning('Graphical Sinkhorn stopped at the iteration limit.');
    end

    stage_mass = zeros(nStage,1);
    for ell = 1:nStage
        stage_mass(ell) = sum(sum(mv(:,unique(stage_node(:,ell)))));
    end
    mass_ref = 0.5*(sum(m0)+sum(mT));
    cons_res = max(abs(stage_mass-mass_ref));
    fprintf('  Final residuals: E0=%.2e, ET=%.2e, V=%.2e, C=%.2e\n', ...
        E0,ET,Vcap,cons_res);

    ribbons = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        ribbons(vtx) = mv(:,vtx);
    end

    diag = struct;
    diag.E0      = E0_hist;
    diag.ET      = ET_hist;
    diag.V       = V_hist;
    diag.res     = max([E0_hist,ET_hist,V_hist],[],2);
    diag.dual    = dual_hist;
    diag.cost    = transportCost( ...
        Kclass,CKclass,Grp,u,v,W,stage_node);
    diag.cons    = cons_res;
    diag.nActive = 0;
    diag.nCell   = nt*numel(middles);
    if isfinite(cap)
        diag.nActive = nnz(mv(:,middles) >= 0.99*cap);
    end
end

function obj = transportCost(Kclass,CKclass,Grp,u,v,W,stage_node)
%TRANSPORTCOST Expected unregularized transport-cost component of the
%current path plans. The value is normalized by their total mass.
    nt     = numel(u);
    nP     = size(stage_node,1);
    nEdge  = numel(Grp);
    F      = forwardMessages(Kclass,Grp,u,W,stage_node);
    B      = backwardMessages(Kclass,Grp,v,W,stage_node);
    Wc     = weightColumns(W,stage_node);
    cost   = 0;
    mass   = 0;

    for ell = 1:nEdge
        if ell == 1
            left = repmat(u,1,nP);
        else
            left = Wc{ell-1}.*F{ell-1};
        end
        if ell == nEdge
            right = repmat(v,1,nP);
        else
            right = Wc{ell}.*B{ell+1};
        end

        gg = Grp{ell};
        for kk = 1:numel(gg)
            ix = gg{kk}.idx;
            cc = gg{kk}.c;
            cost = cost + sum(sum( ...
                left(:,ix).*(CKclass{cc}*right(:,ix))));
            if ell == 1
                mass = mass + sum(sum( ...
                    left(:,ix).*(Kclass{cc}*right(:,ix))));
            end
        end
    end
    obj = cost/max(mass,realmin);
end

function B = backwardMessages(Kclass,Grp,v,W,stage_node)
    nt     = numel(v);
    nP     = size(stage_node,1);
    nEdge  = numel(Grp);
    B      = cell(nEdge,1);
    M      = repmat(v,1,nP);

    for ell = nEdge:-1:1
        nx = zeros(nt,nP);
        gg = Grp{ell};
        for kk = 1:numel(gg)
            ix = gg{kk}.idx;
            nx(:,ix) = Kclass{gg{kk}.c}*M(:,ix);
        end
        B{ell} = nx;
        if ell > 1
            M = W(:,stage_node(:,ell-1)').*nx;
        end
    end
end

function F = forwardMessages(Kclass,Grp,u,W,stage_node)
    nt     = numel(u);
    nP     = size(stage_node,1);
    nEdge  = numel(Grp);
    F      = cell(nEdge,1);
    M      = repmat(u,1,nP);

    for ell = 1:nEdge
        nx = applyForwardKernel(Kclass,Grp{ell},M,nt,nP);
        F{ell} = nx;
        if ell < nEdge
            M = W(:,stage_node(:,ell)').*nx;
        end
    end
end

function nx = applyForwardKernel(Kclass,gg,M,nt,nP)
    nx = zeros(nt,nP);
    for kk = 1:numel(gg)
        ix = gg{kk}.idx;
        nx(:,ix) = Kclass{gg{kk}.c}'*M(:,ix);
    end
end

function Wc = weightColumns(W,stage_node)
    nStage = size(stage_node,2);
    Wc = cell(nStage,1);
    for ell = 1:nStage
        Wc{ell} = W(:,stage_node(:,ell)');
    end
end

function [rho,upper_y] = linearUpperBoundFromZero(k,y)
%LINEARUPPERBOUNDFROMZERO Tight empirical bound y_k <= rho^k, y_0=1.
    k = k(:);
    y = y(:);
    valid = k > 0 & isfinite(y) & y > 0;
    if ~any(valid)
        rho = NaN;
        upper_y = nan(size(y));
        return
    end
    rates = exp(log(y(valid))./k(valid));
    rho = max(rates);
    if ~isfinite(rho) || rho <= 0 || rho >= 1
        upper_y = nan(size(y));
        return
    end
    upper_y = rho.^k;
end


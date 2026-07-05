%% complex_graph.m  –  revised: shared colorbar + zero gaps + topology figure
%
%  Changes in this version
%  ────────────────────────
%  1. drawCombinedFig: gap_w = 0, gap_h ≈ 0 → panels flush for paper layout.
%  2. All text (titles, row labels, colorbar) made larger for readability.
%  3. Figure 2: standalone topology figure showing the DAG structure.

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

%% ── Time grid & marginals ─────────────────────────────────────────────────
nt  = 50;
t   = linspace(0,1,nt)';
G   = @(m,s) exp(-0.5*((t-m)/s).^2);
mu0 = 0.9*G(0.25,0.10) + 0.6*G(0.48,0.07);  mu0 = mu0/sum(mu0);
muT = 0.8*G(0.75,0.09) + 0.7*G(0.85,0.06);  muT = muT/sum(muT);

%% ── Kernel ────────────────────────────────────────────────────────────────
eps_ent = 0.04;
dt_min  = 0.035;
[U,V]   = ndgrid(t,t);
Kedge   = zeros(nt);
kmask   = (V - U) >= dt_min;
Kedge(kmask) = exp(-(V(kmask)-U(kmask)) / eps_ent);

%% ── Two capacity scenarios ────────────────────────────────────────────────
cap_tight = 0.005;
cap_free  = 1e4;

opts = struct('eps0',1e-8,'damp',0.5,'clip_min',1e-8,'clip_max',1e8,'maxIter',1000);

r_tight = containers.Map('KeyType','double','ValueType','any');
r_free  = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    r_tight(vtx) = cap_tight * ones(nt,1);
    r_free(vtx)  = cap_free  * ones(nt,1);
end

%% ── Solve ─────────────────────────────────────────────────────────────────
fprintf('Solving CONSTRAINED  (cap=%.3f, nt=%d, P=%d) ...\n', cap_tight, nt, nPaths);
tic;
[rib_cap,  m0_cap,  mT_cap,  nIt_cap,  bnd_hist_cap, cap_hist_cap, ~] = ...
    runPathSinkhorn(Kedge,mu0,muT,paths,middles,nt,opts,r_tight);
t_cap = toc;
fprintf('  Done %.2f s  (%d iter)\n\n', t_cap, nIt_cap);

fprintf('Solving UNCONSTRAINED (cap=%.0e, nt=%d, P=%d) ...\n', cap_free, nt, nPaths);
tic;
[rib_free, m0_free, mT_free, nIt_free, ~, ~, ~] = ...
    runPathSinkhorn(Kedge,mu0,muT,paths,middles,nt,opts,r_free);
t_free = toc;
fprintf('  Done %.2f s  (%d iter)\n\n', t_free, nIt_free);

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

thr       = 0.003;
gainZ     = 2.5;
tubeW_i   = 0.015;
tubeW_ss  = 0.058;
alph_i    = 0.52;
alph_ss   = 0.92;
ss_hscale = 0.55;

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
capZ_global   = cap_tight / maxInt_global * gainZ;

fprintf('Global colour scale:  interior max = %.4f  (cap=%.4f, free=%.4f)\n', ...
    maxInt_global, maxInt_cap, maxInt_free);

%% ── Figure 1: combined 2-row flux figure, shared colorbar ────────────────
drawCombinedFig(1, snapCap, snapFree, ...
    maxInt_global, maxSrc_global, maxSnk_global, ...
    Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
    thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ss_hscale, ...
    cSrc, cSnk, cEdge, capZ_global, ...
    {sprintf('Constraint $(r{=}%.3f)$', cap_tight), ...
     'Unconstraint'});

%% ── Figure 2: standalone network topology ────────────────────────────────
drawTopologyFig(2, adj, xy, middles, src, snk, nNodes, nStages, nWidth, node_of);

%% ── Figure 3: convergence ─────────────────────────────────────────────────
figure(3); clf;
set(gcf,'Color','w','Units','centimeters','Position',[5 5 11 8], ...
    'PaperPositionMode','auto');
iters = 1:nIt_cap;
semilogy(iters, bnd_hist_cap, '-', 'Color',[0.06 0.28 0.92],'LineWidth',2.5);
hold on;
semilogy(iters, cap_hist_cap, '--','Color',[0.92 0.04 0.04],'LineWidth',2.5);
xlabel('Iteration $k$','Interpreter','latex','FontSize',14);
ylabel('Residual','Interpreter','latex','FontSize',14);
lgd = legend({'Marginal residual','Capacity violation'}, ...
    'Interpreter','latex','FontSize',11,'Location','northeast');
lgd.Box = 'on';  lgd.ItemTokenSize = [14 8];
grid on; box on;
set(gca,'TickLabelInterpreter','latex','FontSize',14,'GridAlpha',0.3);

%% ── Scalability report ────────────────────────────────────────────────────
fprintf('\n── Scalability Report ──────────────────────────────\n');
fprintf('  Interior nodes : %d\n', nInner);
fprintf('  Paths          : %d\n', nPaths);
fprintf('  nt             : %d\n', nt);
fprintf('  Constrained    : %d iter,  %.2f s\n', nIt_cap,  t_cap);
fprintf('  Unconstrained  : %d iter,  %.2f s\n', nIt_free, t_free);
fprintf('────────────────────────────────────────────────────\n');


%% ======================================================================
%% Local functions
%% ======================================================================

function drawCombinedFig(figNo, snapFlux1, snapFlux2, ...
        maxInt, maxSrc, maxSnk, ...
        Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
        thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ss_hscale, ...
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
            showCap  = ~isempty(capZ);
        else
            snapFlux = snapFlux2;
            y_row    = y_bot;
            showCap  = false;
        end

        for i = 1:nSnap
            x_pos = pad_l + (i-1)*w_p;   % gap_w = 0 → no offset
            ax = axes('Parent',gcf, ...
                'Position',[x_pos, y_row, w_p, h_p]);  %#ok<LAXES>
            hold(ax,'on');
            clim(ax,[0 1]);

            plot3(ax, Xe, Ye, Ze, 'Color',cEdge,'LineWidth',0.85);
            scatter3(ax, xy(middles,1), xy(middles,2), zeros(numel(middles),1), ...
                5, [0.60 0.60 0.60],'filled','MarkerEdgeColor','none');
            scatter3(ax, xy([src;snk],1), xy([src;snk],2), [0;0], ...
                14, [0.20 0.20 0.20],'filled','MarkerEdgeColor','none');

            if showCap
                for jj = 1:numel(middles)
                    drawCapMarker(ax, xy(middles(jj),1), xy(middles(jj),2), ...
                        capZ, tubeW_i*1.6);
                end
            end

            for jj = 1:numel(middles)
                vi   = middles(jj);
                val  = snapFlux(i,vi) / maxInt;
                hv   = max(val, thr) * gainZ;
                c_vi = min(1, val);
                drawTubeMapped(ax, xy(vi,1), xy(vi,2), hv, tubeW_i, c_vi, alph_i);
            end

            h0 = max(snapFlux(i,src) / maxSrc, thr) * gainZ * ss_hscale;
            drawTubeFixed(ax, xy(src,1), xy(src,2), h0, tubeW_ss, cSrc, alph_ss);

            hT = max(snapFlux(i,snk) / maxSnk, thr) * gainZ * ss_hscale;
            drawTubeFixed(ax, xy(snk,1), xy(snk,2), hT, tubeW_ss, cSnk, alph_ss);

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

function [ribbons, m0, mT, nIter, bnd_hist, cap_hist, cons_hist] = ...
        runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)

    eps0    = opts.eps0;
    damp    = opts.damp;
    cmin    = opts.clip_min;
    cmax    = opts.clip_max;
    maxIter = opts.maxIter;
    KT      = Kedge';
    bnd_hist  = nan(maxIter, 1);
    cap_hist  = nan(maxIter, 1);
    cons_hist = nan(maxIter, 1);

    nP = numel(paths);
    u  = ones(nt,1);
    v  = ones(nt,1);
    w  = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles,  w(vtx) = ones(nt,1);  end

    A0 = zeros(nt,nP);  AT = zeros(nt,nP);
    Aw = containers.Map('KeyType','double','ValueType','any');
    Iw = containers.Map('KeyType','double','ValueType','any');
    Ow = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        Aw(vtx) = zeros(nt,nP);
        Iw(vtx) = zeros(nt,nP);
        Ow(vtx) = zeros(nt,nP);
    end

    nIter = maxIter;

    for it = 1:maxIter
        A0(:) = 0;  AT(:) = 0;
        for vtx = middles
            Aw(vtx) = zeros(nt,nP);
            Iw(vtx) = zeros(nt,nP);
            Ow(vtx) = zeros(nt,nP);
        end

        for pid = 1:nP
            Pth  = paths{pid};
            Pmid = Pth(2:end-1);
            L    = numel(Pmid);

            if L < 1
                A0(:,pid) = Kedge * v;
                AT(:,pid) = KT    * u;
                continue
            end

            B    = cell(L,1);
            B{L} = max(Kedge * v, eps0);
            for k = L-1:-1:1
                B{k} = max(Kedge * (w(Pmid(k+1)) .* B{k+1}), eps0);
            end

            F    = cell(L,1);
            F{1} = max(KT * u, eps0);
            for k = 2:L
                F{k} = max(KT * (w(Pmid(k-1)) .* F{k-1}), eps0);
            end

            A0(:,pid) = Kedge * (w(Pmid(1)) .* B{1});
            AT(:,pid) = KT    * (w(Pmid(L)) .* F{L});

            for k = 1:L
                vtx = Pmid(k);
                if k == 1
                    inc = (KT * u) .* B{k};
                else
                    inc = (KT * (w(Pmid(k-1)) .* F{k-1})) .* B{k};
                end
                if k == L
                    out = F{k} .* (Kedge * v);
                else
                    out = F{k} .* (Kedge * (w(Pmid(k+1)) .* B{k+1}));
                end
                tmp = Aw(vtx);  tmp(:,pid) = tmp(:,pid) + F{k}.*B{k};    Aw(vtx) = tmp;
                tmp = Iw(vtx);  tmp(:,pid) = tmp(:,pid) + max(inc,eps0);  Iw(vtx) = tmp;
                tmp = Ow(vtx);  tmp(:,pid) = tmp(:,pid) + max(out,eps0);  Ow(vtx) = tmp;
            end
        end

        sA0 = max(sum(A0,2), eps0);
        sAT = max(sum(AT,2), eps0);
        m0  = u .* sA0;
        mT  = v .* sAT;

        cons_res = 0;  cap_res = 0;
        mn = containers.Map('KeyType','double','ValueType','any');
        for vtx = middles
            mv       = w(vtx) .* max(sum(Aw(vtx),2), eps0);
            mn(vtx)  = mv;
            cons_res = cons_res + sum(abs(sum(Iw(vtx),2) - sum(Ow(vtx),2)));
            cap_res  = cap_res  + sum(max(mv - r_cap(vtx), 0));
        end

        u = min(max(u .* (p ./ max(m0,eps0)).^damp, cmin), cmax);
        v = min(max(v .* (q ./ max(mT,eps0)).^damp, cmin), cmax);
        for vtx = middles
            wv = w(vtx) .* min(1, r_cap(vtx) ./ max(mn(vtx),eps0)).^damp;
            w(vtx) = min(max(wv, cmin), cmax);
        end

        bnd_res       = sum(abs(m0-p)) + sum(abs(mT-q));
        bnd_hist(it)  = bnd_res;
        cap_hist(it)  = cap_res;
        cons_hist(it) = cons_res;
        if bnd_res < 2e-4 && cons_res < 2e-5 && cap_res < 2e-5
            nIter = it;  break
        end
    end
    bnd_hist  = bnd_hist(1:nIter);
    cap_hist  = cap_hist(1:nIter);
    cons_hist = cons_hist(1:nIter);
    fprintf('  Converged at iter %d  |  bnd=%.2e  cap=%.2e\n', nIter, bnd_res, cap_res);

    for kp = 1:12
        A0(:) = 0;  AT(:) = 0;
        for pid = 1:nP
            Pth=paths{pid}; Pmid=Pth(2:end-1); L=numel(Pmid);
            if L<1,  A0(:,pid)=Kedge*v;  AT(:,pid)=KT*u;  continue;  end
            B=cell(L,1); B{L}=max(Kedge*v,eps0);
            for ii=L-1:-1:1,  B{ii}=max(Kedge*(w(Pmid(ii+1)).*B{ii+1}),eps0);  end
            F=cell(L,1); F{1}=max(KT*u,eps0);
            for ii=2:L,     F{ii}=max(KT*(w(Pmid(ii-1)).*F{ii-1}),eps0);  end
            A0(:,pid)=Kedge*(w(Pmid(1)).*B{1});
            AT(:,pid)=KT   *(w(Pmid(L)).*F{L});
        end
        u = u .* (p ./ max(u .* max(sum(A0,2),eps0), eps0));
        AT(:) = 0;
        for pid = 1:nP
            Pth=paths{pid}; Pmid=Pth(2:end-1); L=numel(Pmid);
            if L<1,  AT(:,pid)=KT*u;  continue;  end
            Fp=KT*u;
            for ii=2:L,  Fp=KT*(w(Pmid(ii-1)).*Fp);  end
            AT(:,pid)=KT*(w(Pmid(L)).*Fp);
        end
        v = v .* (q ./ max(v .* max(sum(AT,2),eps0), eps0));
    end
    sA0=max(sum(A0,2),eps0);  sAT=max(sum(AT,2),eps0);
    m0=u.*sA0;  mT=v.*sAT;

    ribbons = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        ribbons(vtx) = sum(Aw(vtx) .* (w(vtx) * ones(1,nP)), 2);
    end
end






% %% complex_graph.m  –  revised: 3-D time-rising scatter
% %
% %  Large-scale Sinkhorn demo: 100-node series-parallel DAG, 120 paths.
% %
% %  Network : v0 (source) → 10 stages × 10 nodes → vT (sink)
% %            Each stage-s node connects to 3 random stage-(s+1) nodes.
% %            120 unique random paths sampled.
% %
% %  Solver  : path-wise Sinkhorn, linear cost C(s,t) = t-s
% %            dt_min = 0.035 enforced (causal, physical crossing times)
% %            eps_ent = 0.04, nt = 50
% %
% %  Figure  : single 3-D scatter — network in (x,y) plane at z=0,
% %            flux rises in z = t.  Each node leaves a vertical column
% %            of dots sized & coloured by m_k(t): the active wave
% %            travels left→right (source→sink) as t increases.
% %            · Light gray lines/dots = network structure at z=0
% %            · Hot-coloured dots     = interior node flux (all t)
% %            · Blue column           = source v0,  m0(t)
% %            · Red  column           = sink  vT,  mT(t)
% %            · z-axis = t,  x/y axes hidden
% 
% clear; clc; close all;
% rng(42);
% 
% %% ── Network ───────────────────────────────────────────────────────────────
% nStages = 10;
% nWidth  = 10;
% nInner  = nStages * nWidth;   % 100 interior nodes
% nNodes  = nInner + 2;
% 
% src     = nInner + 1;         % v0
% snk     = nInner + 2;         % vT
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
% %% ── Generate 120 unique random paths ─────────────────────────────────────
% nPaths = 120;
% paths  = {};
% seen   = containers.Map('KeyType','char','ValueType','logical');
% ntry   = 0;
% while numel(paths) < nPaths && ntry < 200000
%     ntry = ntry + 1;
%     kk   = src;  pth = kk;  ok = true;
%     while kk ~= snk
%         nb = adj{kk};
%         if isempty(nb),  ok = false;  break;  end
%         kk = nb(randi(numel(nb)));
%         pth(end+1) = kk;  %#ok<AGROW>
%     end
%     if ok
%         key = sprintf('%d,', pth);
%         if ~isKey(seen, key)
%             seen(key)    = true;
%             paths{end+1} = pth;  %#ok<AGROW>
%         end
%     end
% end
% nPaths = numel(paths);
% fprintf('Generated %d unique paths  (%d attempts)\n', nPaths, ntry);
% 
% middles = unique([paths{:}]);
% middles = middles(middles ~= src & middles ~= snk);
% fprintf('Active interior nodes: %d / %d\n\n', numel(middles), nInner);
% 
% %% ── Time grid & marginals ─────────────────────────────────────────────────
% nt  = 50;
% t   = linspace(0,1,nt)';
% G   = @(m,s) exp(-0.5*((t-m)/s).^2);
% mu0 = 0.9*G(0.25,0.10) + 0.6*G(0.48,0.07);  mu0 = mu0/sum(mu0);
% muT = 0.8*G(0.75,0.09) + 0.7*G(0.85,0.06);  muT = muT/sum(muT);
% 
% %% ── Kernel: C(s,t) = t-s, dt_min enforced ────────────────────────────────
% eps_ent = 0.04;
% dt_min  = 0.035;   % minimum travel time per edge (physical feasibility)
% [U,V]   = ndgrid(t,t);
% Kedge   = zeros(nt);
% kmask   = (V - U) >= dt_min;
% Kedge(kmask) = exp(-(V(kmask)-U(kmask)) / eps_ent);
% 
% %% ── Two capacity scenarios ────────────────────────────────────────────────
% cap_tight = 0.005;    % tight nodal capacity  (Figure 1)
% cap_free  = 1e4;     % effectively unconstrained  (Figure 2)
% 
% opts = struct('eps0',1e-8,'damp',0.5,'clip_min',1e-8,'clip_max',1e8,'maxIter',1000);
% 
% r_tight = containers.Map('KeyType','double','ValueType','any');
% r_free  = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles
%     r_tight(vtx) = cap_tight * ones(nt,1);
%     r_free(vtx)  = cap_free  * ones(nt,1);
% end
% 
% %% ── Solve: constrained ────────────────────────────────────────────────────
% fprintf('Solving CONSTRAINED  (cap = %.2f, nt=%d, %d paths) ...\n', ...
%     cap_tight, nt, nPaths);
% tic;
% [rib_cap,  m0_cap,  mT_cap,  nIt_cap,  bnd_hist_cap, cap_hist_cap, cons_hist_cap] = runPathSinkhorn(Kedge,mu0,muT,paths,middles,nt,opts,r_tight);
% t_cap = toc;
% fprintf('  Done %.2f s  (%d iter)\n\n', t_cap, nIt_cap);
% 
% %% ── Solve: unconstrained ──────────────────────────────────────────────────
% fprintf('Solving UNCONSTRAINED (cap = %.0e, nt=%d, %d paths) ...\n', ...
%     cap_free, nt, nPaths);
% tic;
% [rib_free, m0_free, mT_free, nIt_free, ~, ~, ~] = runPathSinkhorn(Kedge,mu0,muT,paths,middles,nt,opts,r_free);
% t_free = toc;
% fprintf('  Done %.2f s  (%d iter)\n\n', t_free, nIt_free);
% 
% %% ── Node (x,y) positions ──────────────────────────────────────────────────
% %   x = stage index (left → right),  y = position within stage (bottom → top)
% xy = zeros(nNodes, 2);
% for s = 1:nStages
%     for pp = 1:nWidth
%         xy(node_of(s,pp),:) = [ s/(nStages+1),  (pp-1)/(nWidth-1) ];
%     end
% end
% xy(src,:) = [0,   0.5];
% xy(snk,:) = [1,   0.5];
% 
% %% ── Flux matrices & snapshot interpolations ───────────────────────────────
% nSnap  = 6;
% t_snap = linspace(0, 1, nSnap);   % t = 0, 0.2, 0.4, 0.6, 0.8, 1.0
% 
% % ── Constrained ───────────────────────────────────────────────────────────
% fluxCap = zeros(nt, nNodes);
% for vtx = middles
%     if isKey(rib_cap, vtx),  fluxCap(:,vtx) = rib_cap(vtx);  end
% end
% fluxCap(:,src) = m0_cap;   fluxCap(:,snk) = mT_cap;
% snapCap = zeros(nSnap, nNodes);
% for vi = 1:nNodes
%     snapCap(:,vi) = max(0, interp1(t, fluxCap(:,vi), t_snap','linear',0));
% end
% 
% % ── Unconstrained ─────────────────────────────────────────────────────────
% fluxFree = zeros(nt, nNodes);
% for vtx = middles
%     if isKey(rib_free, vtx),  fluxFree(:,vtx) = rib_free(vtx);  end
% end
% fluxFree(:,src) = m0_free;  fluxFree(:,snk) = mT_free;
% snapFree = zeros(nSnap, nNodes);
% for vi = 1:nNodes
%     snapFree(:,vi) = max(0, interp1(t, fluxFree(:,vi), t_snap','linear',0));
% end
% 
% %% ── Pre-build edge geometry (NaN-separated, drawn once per panel) ─────────
% %   One plot3 call per panel instead of one per edge → much faster render.
% Xe=[]; Ye=[]; Ze=[];
% for kk = 1:nNodes
%     for nb = adj{kk}(:)'
%         Xe = [Xe,  xy(kk,1), xy(nb,1), NaN];  %#ok<AGROW>
%         Ye = [Ye,  xy(kk,2), xy(nb,2), NaN];  %#ok<AGROW>
%         Ze = [Ze,  0,        0,        NaN];  %#ok<AGROW>
%     end
% end
% 
% %% ── Visual parameters ─────────────────────────────────────────────────────
% cSrc  = [0.06  0.28  0.92];   % saturated blue  = initial distribution v0
% cSnk  = [0.92  0.04  0.04];   % saturated red   = target distribution  vT
% cEdge = [0.38  0.38  0.38];
% 
% thr       = 0.003;
% gainZ     = 2.5;
% tubeW_i   = 0.015;    % interior tube half-width
% tubeW_ss  = 0.058;    % source / sink tube half-width (wider → more visible)
% alph_i    = 0.52;     % interior: more transparent so graph reads through
% alph_ss   = 0.92;     % src/sink: nearly opaque so they stand out
% ss_hscale = 0.55;     % src/sink height scale (taller → more visible)
% 
% %% ── Per-scenario normalisations ───────────────────────────────────────────
% maxInt_cap  = max(max(fluxCap(:, middles)))  + 1e-16;
% maxSrc_cap  = max(m0_cap)  + 1e-16;
% maxSnk_cap  = max(mT_cap)  + 1e-16;
% capZ        = cap_tight / maxInt_cap * gainZ;  % capacity ceiling in z-units
% 
% maxInt_free = max(max(fluxFree(:, middles))) + 1e-16;
% maxSrc_free = max(m0_free) + 1e-16;
% maxSnk_free = max(mT_free) + 1e-16;
% 
% % %% ── Figure 1: WITH tight capacity constraint ──────────────────────────────
% % drawNetworkFig(1, snapCap,  maxInt_cap,  maxSrc_cap,  maxSnk_cap, ...
% %     Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
% %     thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ss_hscale, ...
% %     cSrc, cSnk, cEdge, capZ, ...
% %     sprintf('With capacity constraint  $(r = %.2f)$', cap_tight));
% % 
% % %% ── Figure 2: WITHOUT capacity constraint ─────────────────────────────────
% % drawNetworkFig(2, snapFree, maxInt_free, maxSrc_free, maxSnk_free, ...
% %     Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
% %     thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ss_hscale, ...
% %     cSrc, cSnk, cEdge, [], ...
% %     'Without capacity constraint');
% 
% %% ── Figure 3: All marginal violations — constrained case only ─────────────
% figure(3); clf;
% set(gcf,'Color','w','Units','centimeters','Position',[5 5 11 8], ...
%     'PaperPositionMode','auto');
% 
% iters = 1:nIt_cap;
% 
% semilogy(iters, bnd_hist_cap, '-', ...
%     'Color',[0.06 0.28 0.92],'LineWidth',2.5);
% hold on;
% 
% semilogy(iters, cap_hist_cap, '--', ...
%     'Color',[0.92 0.04 0.04],'LineWidth',2.5);
% 
% xlabel('Iteration $k$','Interpreter','latex','FontSize',14);
% ylabel('Residual','Interpreter','latex','FontSize',14);
% 
% lgd = legend({'Marginal residual', ...
%               'Capacity violation'}, ...
%     'Interpreter','latex', ...
%     'FontSize',11, ...
%     'Location','northeast');
% 
% lgd.Box = 'on';
% lgd.ItemTokenSize = [14 8];
% 
% grid on; box on;
% set(gca,'TickLabelInterpreter','latex','FontSize',14,'GridAlpha',0.3);
% 
% 
% %% ── Scalability report ────────────────────────────────────────────────────
% fprintf('\n── Scalability Report ──────────────────────────────\n');
% fprintf('  Interior nodes  : %d  (+ source + sink)\n', nInner);
% fprintf('  Paths           : %d\n', nPaths);
% fprintf('  Time grid nt    : %d\n', nt);
% fprintf('  Constrained     : %d iter,  %.2f s  (cap=%.2f)\n', nIt_cap,  t_cap,  cap_tight);
% fprintf('  Unconstrained   : %d iter,  %.2f s\n',              nIt_free, t_free);
% fprintf('────────────────────────────────────────────────────\n');
% 
% 
% %% ======================================================================
% %% Local functions
% %% ======================================================================
% 
% function drawNetworkFig(figNo, snapFlux, maxInt, maxSrc, maxSnk, ...
%         Xe, Ye, Ze, xy, middles, src, snk, nSnap, t_snap, ...
%         thr, gainZ, tubeW_i, tubeW_ss, alph_i, alph_ss, ss_hscale, ...
%         cSrc, cSnk, cEdge, capZ, figtitle)
% %DRAWNETWORKFIG  One-row 1×nSnap figure with manual tight layout.
% %   capZ = [] → no capacity markers (unconstrained figure).
% 
%     figure(figNo); clf;
%     set(gcf,'Color','w','Units','centimeters', ...
%         'Position',[1, 2 + (figNo-1)*17, 36, 8], ...
%         'PaperPositionMode','auto','Renderer','painters');
%     sgtitle(figtitle,'Interpreter','latex','FontSize',10);
% 
%     cmap = getVizMap();
%     colormap(gcf, cmap);
% 
%     % ── Manual tight layout (no tiledlayout gaps) ─────────────────────────
%     pad_l = 0.008;   pad_r = 0.080;   pad_t = 0.13;   pad_b = 0.02;
%     gap   = 0.001;   % pixel-thin gap between panels
%     w_p   = (1 - pad_l - pad_r - gap*(nSnap-1)) / nSnap;
%     h_p   = 1 - pad_t - pad_b;
% 
%     for i = 1:nSnap
%         x_pos = pad_l + (i-1) * (w_p + gap);
%         ax = axes('Parent', gcf, 'Position', [x_pos, pad_b, w_p, h_p]);  %#ok<LAXES>
%         hold(ax,'on');
%         clim(ax,[0 1]);
% 
%         % Edges + node anchors at z=0
%         plot3(ax, Xe, Ye, Ze, 'Color',cEdge,'LineWidth',0.85);
%         scatter3(ax, xy(middles,1), xy(middles,2), zeros(numel(middles),1), ...
%             5, [0.60 0.60 0.60],'filled','MarkerEdgeColor','none');
%         scatter3(ax, xy([src;snk],1), xy([src;snk],2), [0;0], ...
%             14, [0.20 0.20 0.20],'filled','MarkerEdgeColor','none');
% 
%         % Grey capacity ceiling (constrained only)
%         if ~isempty(capZ)
%             for jj = 1:numel(middles)
%                 vi = middles(jj);
%                 drawCapMarker(ax, xy(vi,1), xy(vi,2), capZ, tubeW_i*1.6);
%             end
%         end
% 
%         % Interior cubic tubes — colour = absolute flux / global peak
%         for jj = 1:numel(middles)
%             vi   = middles(jj);
%             hv   = max(snapFlux(i,vi) / maxInt * gainZ, thr * gainZ);
%             c_vi = min(1, snapFlux(i,vi) / maxInt);   % absolute normalisation
%             drawTubeMapped(ax, xy(vi,1), xy(vi,2), hv, tubeW_i, c_vi, alph_i);
%         end
% 
%         % Source v0 (blue)
%         h0 = max(snapFlux(i,src) / maxSrc, thr) * gainZ * ss_hscale;
%         drawTubeFixed(ax, xy(src,1), xy(src,2), h0, tubeW_ss, cSrc, alph_ss);
% 
%         % Sink vT (red)
%         hT = max(snapFlux(i,snk) / maxSnk, thr) * gainZ * ss_hscale;
%         drawTubeFixed(ax, xy(snk,1), xy(snk,2), hT, tubeW_ss, cSnk, alph_ss);
% 
%         title(ax, sprintf('$t=%d/5$', i-1), ...
%             'Interpreter','latex','FontSize',8,'FontWeight','normal');
%         set(ax,'Visible','off');
%         ax.Title.Visible = 'on';   % restore title hidden by Visible=off
%         set(ax,'XLim',[-0.07 1.07],'YLim',[-0.10 1.10],'ZLim',[0 gainZ*1.05]);
%         set(ax,'PlotBoxAspectRatio',[1.1 1 gainZ*0.85]);
%         view(ax, 35, 28);
%         hold(ax,'off');
%     end
% 
%     % Shared colorbar: standalone invisible axes so it doesn't shrink panels
%     cbax = axes('Parent', gcf, ...
%         'Position', [1-pad_r+0.018, pad_b, 0.001, h_p], 'Visible','off');
%     clim(cbax,[0 1]);
%     colormap(cbax, cmap);
%     cb = colorbar(cbax, 'Position', [1-pad_r+0.018, pad_b, 0.020, h_p]);
%     cb.TickLabelInterpreter = 'latex';
%     cb.FontSize             = 8;
%     cb.Label.Interpreter    = 'latex';
%     cb.Label.String         = 'Absolute flux $m_k(t_i)/\!\max_{k,i}\,m_k(t_i)$';
%     cb.Label.FontSize       = 9;
% end
% 
% % ──────────────────────────────────────────────────────────────────────────
% 
% function drawTubeMapped(ax, x0, y0, h, w, c_val, alph)
% %DRAWTUBEMAPPED  Square prism coloured via getVizMap() (parula 0.3–1).
%     if h <= 0, h = 1e-4; end
%     cmap = getVizMap();
%     nC   = size(cmap,1);
%     idx  = max(1, min(nC, round(c_val * (nC-1)) + 1));
%     col  = cmap(idx,:);
%     verts = [x0-w, y0-w, 0; x0+w, y0-w, 0; x0+w, y0+w, 0; x0-w, y0+w, 0;
%              x0-w, y0-w, h; x0+w, y0-w, h; x0+w, y0+w, h; x0-w, y0+w, h];
%     faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
%     patch(ax,'Vertices',verts,'Faces',faces, ...
%         'FaceColor',col,'FaceAlpha',alph, ...
%         'EdgeColor',col*0.60,'EdgeAlpha',0.30,'LineWidth',0.4);
% end
% 
% function drawTubeFixed(ax, x0, y0, h, w, col, alph)
% %DRAWTUBEFIXED  Square prism with a fixed RGB colour (source / sink).
%     if h <= 0, h = 1e-4; end
%     verts = [x0-w, y0-w, 0; x0+w, y0-w, 0; x0+w, y0+w, 0; x0-w, y0+w, 0;
%              x0-w, y0-w, h; x0+w, y0-w, h; x0+w, y0+w, h; x0-w, y0+w, h];
%     faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
%     patch(ax,'Vertices',verts,'Faces',faces, ...
%         'FaceColor',col,'FaceAlpha',alph, ...
%         'EdgeColor',col*0.60,'EdgeAlpha',0.45,'LineWidth',0.5);
% end
% 
% function cmap = getVizMap()
% %GETVIZMAP  Parula with the deep-blue bottom 30% removed.
% %   Returns a 180×3 colormap running from light-blue/cyan → green → yellow.
%     raw  = parula(256);
%     cmap = raw(round(0.30*256)+1 : end, :);
% end
% 
% function drawCapMarker(ax, x0, y0, z0, w)
% %DRAWCAPMARKER  Semi-transparent grey flat square at capacity ceiling z=z0.
%     if z0 <= 0, return; end
%     verts = [x0-w, y0-w, z0; x0+w, y0-w, z0;
%              x0+w, y0+w, z0; x0-w, y0+w, z0];
%     patch(ax,'Vertices',verts,'Faces',[1 2 3 4], ...
%         'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.30, ...
%         'EdgeColor',[0.45 0.45 0.45],'EdgeAlpha',0.55,'LineWidth',0.5);
% end
% 
% function [ribbons, m0, mT, nIter, bnd_hist, cap_hist, cons_hist] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)
% %RUNPATHSINKHORN  Path-wise Sinkhorn with nodal capacity multipliers.
% %   Returns ribbons, boundary marginals m0/mT, iteration count, and
% %   per-iteration residuals bnd_hist / cap_hist / cons_hist.
% 
%     eps0    = opts.eps0;
%     damp    = opts.damp;
%     cmin    = opts.clip_min;
%     cmax    = opts.clip_max;
%     maxIter = opts.maxIter;
%     KT      = Kedge';
%     bnd_hist  = nan(maxIter, 1);
%     cap_hist  = nan(maxIter, 1);
%     cons_hist = nan(maxIter, 1);
% 
%     nP = numel(paths);
%     u  = ones(nt,1);
%     v  = ones(nt,1);
%     w  = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles,  w(vtx) = ones(nt,1);  end
% 
%     A0 = zeros(nt,nP);  AT = zeros(nt,nP);
% 
%     Aw = containers.Map('KeyType','double','ValueType','any');
%     Iw = containers.Map('KeyType','double','ValueType','any');
%     Ow = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         Aw(vtx) = zeros(nt,nP);
%         Iw(vtx) = zeros(nt,nP);
%         Ow(vtx) = zeros(nt,nP);
%     end
% 
%     nIter = maxIter;
% 
%     for it = 1:maxIter
% 
%         A0(:) = 0;  AT(:) = 0;
%         for vtx = middles
%             Aw(vtx) = zeros(nt,nP);
%             Iw(vtx) = zeros(nt,nP);
%             Ow(vtx) = zeros(nt,nP);
%         end
% 
%         for pid = 1:nP
%             Pth  = paths{pid};
%             Pmid = Pth(2:end-1);
%             L    = numel(Pmid);
% 
%             if L < 1
%                 A0(:,pid) = Kedge * v;
%                 AT(:,pid) = KT    * u;
%                 continue
%             end
% 
%             % Backward messages
%             B    = cell(L,1);
%             B{L} = max(Kedge * v, eps0);
%             for k = L-1:-1:1
%                 B{k} = max(Kedge * (w(Pmid(k+1)) .* B{k+1}), eps0);
%             end
% 
%             % Forward messages
%             F    = cell(L,1);
%             F{1} = max(KT * u, eps0);
%             for k = 2:L
%                 F{k} = max(KT * (w(Pmid(k-1)) .* F{k-1}), eps0);
%             end
% 
%             A0(:,pid) = Kedge * (w(Pmid(1)) .* B{1});
%             AT(:,pid) = KT    * (w(Pmid(L)) .* F{L});
% 
%             for k = 1:L
%                 vtx = Pmid(k);
%                 if k == 1
%                     inc = (KT * u) .* B{k};
%                 else
%                     inc = (KT * (w(Pmid(k-1)) .* F{k-1})) .* B{k};
%                 end
%                 if k == L
%                     out = F{k} .* (Kedge * v);
%                 else
%                     out = F{k} .* (Kedge * (w(Pmid(k+1)) .* B{k+1}));
%                 end
%                 tmp = Aw(vtx);  tmp(:,pid) = tmp(:,pid) + F{k}.*B{k};    Aw(vtx) = tmp;
%                 tmp = Iw(vtx);  tmp(:,pid) = tmp(:,pid) + max(inc,eps0);  Iw(vtx) = tmp;
%                 tmp = Ow(vtx);  tmp(:,pid) = tmp(:,pid) + max(out,eps0);  Ow(vtx) = tmp;
%             end
%         end
% 
%         sA0 = max(sum(A0,2), eps0);
%         sAT = max(sum(AT,2), eps0);
%         m0  = u .* sA0;
%         mT  = v .* sAT;
% 
%         cons_res = 0;  cap_res = 0;
%         mn = containers.Map('KeyType','double','ValueType','any');
%         for vtx = middles
%             mv       = w(vtx) .* max(sum(Aw(vtx),2), eps0);
%             mn(vtx)  = mv;
%             cons_res = cons_res + sum(abs(sum(Iw(vtx),2) - sum(Ow(vtx),2)));
%             cap_res  = cap_res  + sum(max(mv - r_cap(vtx), 0));
%         end
% 
%         u = min(max(u .* (p  ./ max(m0,eps0)).^damp, cmin), cmax);
%         v = min(max(v .* (q  ./ max(mT,eps0)).^damp, cmin), cmax);
%         for vtx = middles
%             wv = w(vtx) .* min(1, r_cap(vtx) ./ max(mn(vtx),eps0)).^damp;
%             w(vtx) = min(max(wv, cmin), cmax);
%         end
% 
%         bnd_res       = sum(abs(m0-p)) + sum(abs(mT-q));
%         bnd_hist(it)  = bnd_res;
%         cap_hist(it)  = cap_res;
%         cons_hist(it) = cons_res;
%         if bnd_res < 2e-4 && cons_res < 2e-5 && cap_res < 2e-5
%             nIter = it;
%             break
%         end
%     end
%     bnd_hist  = bnd_hist(1:nIter);
%     cap_hist  = cap_hist(1:nIter);
%     cons_hist = cons_hist(1:nIter);
%     fprintf('  Converged at iter %d  |  bnd=%.2e  cap=%.2e\n', nIter, bnd_res, cap_res);
% 
%     %% Boundary polish (w frozen)
%     for kp = 1:12
%         A0(:) = 0;  AT(:) = 0;
%         for pid = 1:nP
%             Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
%             if L < 1,  A0(:,pid)=Kedge*v;  AT(:,pid)=KT*u;  continue;  end
%             B=cell(L,1); B{L}=max(Kedge*v,eps0);
%             for ii=L-1:-1:1,  B{ii}=max(Kedge*(w(Pmid(ii+1)).*B{ii+1}),eps0);  end
%             F=cell(L,1); F{1}=max(KT*u,eps0);
%             for ii=2:L,     F{ii}=max(KT*(w(Pmid(ii-1)).*F{ii-1}),eps0);  end
%             A0(:,pid)=Kedge*(w(Pmid(1)).*B{1});
%             AT(:,pid)=KT   *(w(Pmid(L)).*F{L});
%         end
%         u = u .* (p ./ max(u .* max(sum(A0,2),eps0), eps0));
% 
%         AT(:) = 0;
%         for pid = 1:nP
%             Pth=paths{pid}; Pmid=Pth(2:end-1); L=numel(Pmid);
%             if L<1,  AT(:,pid)=KT*u;  continue;  end
%             Fp=KT*u;
%             for ii=2:L,  Fp=KT*(w(Pmid(ii-1)).*Fp);  end
%             AT(:,pid)=KT*(w(Pmid(L)).*Fp);
%         end
%         v = v .* (q ./ max(v .* max(sum(AT,2),eps0), eps0));
%     end
%     sA0=max(sum(A0,2),eps0);  sAT=max(sum(AT,2),eps0);
%     m0=u.*sA0;  mT=v.*sAT;
% 
%     %% Collect node ribbons
%     ribbons = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         ribbons(vtx) = sum(Aw(vtx) .* (w(vtx) * ones(1,nP)), 2);
%     end
% end
% 
% 
% 

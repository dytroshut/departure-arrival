%% Path-graph: uniform vs time-varying capacity (dual-lobe 2-D plot)
%
%  Revised kernel:
%    - Cost  C(s,t) = t - s   (linear; prefers fast transit → concentrated lobes)
%    - No row-normalisation   (unnormalised Sinkhorn, standard formulation)
%    - eps_ent = 0.065        (spread ≈ eps × hops ≈ 0.065×6 ≈ 0.39, within [0,1])

clear; clc; close all;

%% ── Network topology ─────────────────────────────────────────────────────
nx = 4; ny = 3;
node_id = @(ix,iy) (iy-1)*nx + ix;

v0 = node_id(1,1);
vT = node_id(4,3);
J1 = node_id(2,2);
J2 = node_id(3,2);

pathA = [node_id(1,1), node_id(2,1), J1, J2, node_id(4,2), vT];
pathB = [node_id(1,1), node_id(1,2), J1, J2, node_id(3,3), vT];
pathC = [node_id(1,1), node_id(2,1), J1, J2, node_id(3,3), vT];
paths = {pathA, pathB, pathC};

nodes_on_paths = unique([pathA(:); pathB(:); pathC(:)])';
middles        = setdiff(nodes_on_paths, [v0 vT]);

% Display order (left to right)
seq_nodes = [v0, node_id(2,1), node_id(1,2), J1, J2, node_id(4,2), node_id(3,3), vT];
labels    = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','$v_{\mathcal T}$'};

% Order index for each interior node (used when building time-varying cap)
interior_seq = seq_nodes(2:end-1);
nodeIdx = containers.Map(num2cell(interior_seq), num2cell(1:numel(interior_seq)));

%% ── Time grid and marginals ──────────────────────────────────────────────
nt = 200;
t  = linspace(0,1,nt)';
Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);

p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07);  p = p/sum(p);
q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06);  q = q/sum(q);

%% ── Kernel: linear cost, no row-normalisation ────────────────────────────
%  C(s,t) = t - s  (zero for t<=s, linear for t>s)
%  Kedge(s,t) = exp(-C(s,t) / eps_ent)
%  With eps_ent small the mass hugs the diagonal → concentrated node profiles.
eps_ent = 0.065;

[U,V] = ndgrid(t,t);
mask   = (V > U);          % strictly causal

Kedge          = zeros(nt);
Kedge(mask)    = exp(-(V(mask)-U(mask)) / eps_ent);
% No row-normalisation: standard (unbalanced) Sinkhorn keeps the marginal
% constraints sharp via the u/v multipliers alone.

%% ── Solver options ───────────────────────────────────────────────────────
opts.eps0     = 1e-8;
opts.damp     = 0.4;
opts.clip_min = 1e-8;
opts.clip_max = 1e8;
opts.maxIter  = 4000;

%% ── Uniform capacity ─────────────────────────────────────────────────────
cap_scalar = 0.009;
r_capU = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    r_capU(vtx) = cap_scalar * ones(nt,1);
end

fprintf('Solving uniform-capacity case ...\n');
[ribbonsU, m0U, mTU] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_capU);

%% ── Time-varying capacity (sinusoidal, mean = cap_scalar) ────────────────
r_capV = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    k          = nodeIdx(vtx);
    phase      = (k-1) * (pi/2);
    modulation = 1 + 0.65 * sin(2*pi*t/0.60 + phase);
    profile    = cap_scalar * max(modulation, 0.05);
    profile    = profile * (cap_scalar * nt / sum(profile));   % fix mean
    r_capV(vtx) = max(profile, 5e-4);
end

fprintf('Solving time-varying-capacity case ...\n');
[ribbonsV, m0V, mTV] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_capV);

%% ── Global normalisation (both solutions on the same scale) ──────────────
epsn = 1e-16;
allVals = [m0U; mTU; m0V; mTV];
for vtx = middles
    allVals = [allVals; ribbonsU(vtx); ribbonsV(vtx); ...
               r_capU(vtx);            r_capV(vtx)];   %#ok<AGROW>
end
denG = max(allVals) + epsn;

%% ── Colours and layout ───────────────────────────────────────────────────
xAll = linspace(0,1,numel(seq_nodes));
xmap = containers.Map(num2cell(seq_nodes), num2cell(xAll));

cBlue    = [0.0000 0.4470 0.7410];
cArrival = [0.6350 0.0780 0.1840];
cUniform = [0.20   0.62   0.20  ];
cVarying = [0.88   0.28   0.10  ];
cRail    = [0.20   0.20   0.20  ];
cEdge    = [0.38   0.38   0.38  ];

gain    = 0.12;    % half-width of each lobe
arc_h   = 0.10;

%% ── Figure ───────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Color','w','Units','centimeters','Position',[2 2 20 10.5], ...
    'PaperPositionMode','auto','Renderer','painters');
hold on; box on; grid on;

% Vertical rails
for i = 1:numel(seq_nodes)
    x = xAll(i);
    plot([x x],[0 1],'Color',cRail,'LineWidth',1.4);
    plot(x, 0,'o','Color',cRail,'MarkerFaceColor',cRail,'MarkerSize',7);
end

% Source (blue, left)
drawLobe(xmap(v0),  m0U/denG, t, gain, cBlue,    0.35, 2.6, 'left');

% Interior nodes — uniform RIGHT, time-varying LEFT
for i = 2:(numel(seq_nodes)-1)
    vtx = seq_nodes(i);
    x   = xAll(i);

    if isKey(ribbonsU, vtx)
        mvU = ribbonsU(vtx) / denG;
        mvV = ribbonsV(vtx) / denG;
        cvU = r_capU(vtx)   / denG;
        cvV = r_capV(vtx)   / denG;
    else
        mvU = zeros(nt,1);  mvV = zeros(nt,1);
        cvU = zeros(nt,1);  cvV = zeros(nt,1);
    end

    % Uniform — right side + dashed cap
    drawLobe(x, mvU, t, gain, cUniform, 0.25, 2.0, 'right');
    plot(x + gain*cvU, t, '--', 'Color',cUniform, 'LineWidth',1.3);

    % Time-varying — left side + dashed cap
    drawLobe(x, mvV, t, gain, cVarying, 0.25, 2.0, 'left');
    plot(x - gain*cvV, t, '--', 'Color',cVarying, 'LineWidth',1.3);
end

% Sink (dark red, right)
drawLobe(xmap(vT), mTU/denG, t, gain, cArrival, 0.35, 2.6, 'right');

% Bottom arcs
for pid = 1:numel(paths)
    Pth = paths{pid};
    for k = 1:(numel(Pth)-1)
        drawArc(xmap(Pth(k)), xmap(Pth(k+1)), arc_h, cEdge, 1.2);
    end
end

% Node labels
for i = 1:numel(seq_nodes)
    col = [0 0 0];
    if i == 1,              col = cBlue;    end
    if i == numel(seq_nodes), col = cArrival; end
    text(xAll(i), -0.045, labels{i}, 'Interpreter','latex','Color',col, ...
        'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','top');
end

% Legend
hU = plot(nan,nan,'Color',cUniform,'LineWidth',2.2);
hV = plot(nan,nan,'Color',cVarying,'LineWidth',2.2);
legend([hU hV],{'Uniform capacity','Time-varying capacity'}, ...
    'Interpreter','latex','FontSize',12,'Location','northoutside', ...
    'Orientation','horizontal','Box','off');

xlim([-0.09 1.09]);  ylim([-0.13 1.03]);
set(gca,'XTick',[],'YTick',0:0.2:1,'TickLabelInterpreter','latex', ...
    'FontName','Times New Roman','FontSize',14,'LineWidth',0.8);
ylabel('$t$','Interpreter','latex','FontName','Times New Roman','FontSize',18);
hold off;

exportgraphics(gcf,'path_uniform_vs_time_varying_capacity.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(gcf,'path_uniform_vs_time_varying_capacity.png', ...
    'Resolution',600,'BackgroundColor','white');

%% ── Capacity profiles ────────────────────────────────────────────────────
figure(2); clf;
set(gcf,'Color','w','Units','centimeters','Position',[2 2 14 8], ...
    'PaperPositionMode','auto','Renderer','painters');
hold on;
cCap = lines(numel(interior_seq));
for i = 1:numel(interior_seq)
    vtx = interior_seq(i);
    plot(t, r_capV(vtx), 'Color',cCap(i,:), 'LineWidth',1.8);
end
plot(t, cap_scalar*ones(nt,1), 'k--', 'LineWidth',1.8);
lbs = arrayfun(@(k) labels{k+1}, 1:numel(interior_seq), 'UniformOutput',false);
lbs{end+1} = 'Uniform';
legend(lbs,'Interpreter','latex','Location','eastoutside','Box','off');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('Capacity $r_k(t)$','Interpreter','latex','FontSize',16);
set(gca,'FontName','Times New Roman','FontSize',13,'TickLabelInterpreter','latex');
grid on; box on; hold off;

exportgraphics(gcf,'path_time_varying_capacity_profiles.pdf', ...
    'ContentType','vector','BackgroundColor','white');

%% ======================================================================
%% Local functions
%% ======================================================================

function [ribbons, m0, mT] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)
%RUNPATHSINKHORN  Path-wise Sinkhorn with shared nodal capacity multipliers.

    eps0     = opts.eps0;
    damp     = opts.damp;
    clip_min = opts.clip_min;
    clip_max = opts.clip_max;
    maxIter  = opts.maxIter;

    u = ones(nt,1);
    v = ones(nt,1);
    w = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles,  w(vtx) = ones(nt,1);  end

    nP           = numel(paths);
    A0_paths     = zeros(nt, nP);
    AT_paths     = zeros(nt, nP);
    node_A_wo    = containers.Map('KeyType','double','ValueType','any');
    node_in_wo   = containers.Map('KeyType','double','ValueType','any');
    node_out_wo  = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        node_A_wo(vtx)   = zeros(nt, nP);
        node_in_wo(vtx)  = zeros(nt, nP);
        node_out_wo(vtx) = zeros(nt, nP);
    end

    %% Main iteration
    for it = 1:maxIter
        A0_paths(:) = 0;   AT_paths(:) = 0;
        for vtx = middles
            node_A_wo(vtx)  = 0*node_A_wo(vtx);
            node_in_wo(vtx) = 0*node_in_wo(vtx);
            node_out_wo(vtx)= 0*node_out_wo(vtx);
        end

        for pid = 1:nP
            Pth  = paths{pid};
            Pmid = Pth(2:end-1);
            L    = numel(Pmid);

            if L >= 1
                Bwo    = cell(L,1);
                Bwo{L} = max(Kedge * v, eps0);
                for k = L-1:-1:1
                    Bwo{k} = max(Kedge * (w(Pmid(k+1)).*Bwo{k+1}), eps0);
                end
                Fwo    = cell(L,1);
                Fwo{1} = max(Kedge' * u, eps0);
                for k = 2:L
                    Fwo{k} = max(Kedge' * (w(Pmid(k-1)).*Fwo{k-1}), eps0);
                end
                A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
                AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});

                for k = 1:L
                    vtx = Pmid(k);
                    if k == 1
                        incoming = (Kedge'*u) .* Bwo{k};
                    else
                        incoming = (Kedge'*(w(Pmid(k-1)).*Fwo{k-1})) .* Bwo{k};
                    end
                    if k == L
                        outgoing = Fwo{k} .* (Kedge*v);
                    else
                        outgoing = Fwo{k} .* (Kedge*(w(Pmid(k+1)).*Bwo{k+1}));
                    end
                    incoming = max(incoming, eps0);
                    outgoing = max(outgoing, eps0);

                    A = node_A_wo(vtx);  A(:,pid) = A(:,pid) + Fwo{k}.*Bwo{k};  node_A_wo(vtx) = A;
                    I = node_in_wo(vtx); I(:,pid) = I(:,pid) + incoming;          node_in_wo(vtx) = I;
                    O = node_out_wo(vtx);O(:,pid) = O(:,pid) + outgoing;          node_out_wo(vtx) = O;
                end
            else
                A0_paths(:,pid) = Kedge * v;
                AT_paths(:,pid) = Kedge' * u;
            end
        end

        sumA0 = max(sum(A0_paths,2), eps0);
        sumAT = max(sum(AT_paths,2), eps0);
        m0    = u .* sumA0;
        mT    = v .* sumAT;

        cons_res = 0;  cap_res = 0;
        m_node = containers.Map('KeyType','double','ValueType','any');
        for vtx = middles
            A_sum = max(sum(node_A_wo(vtx),2), eps0);
            mv    = w(vtx) .* A_sum;
            m_node(vtx) = mv;
            IN  = sum(node_in_wo(vtx),2);
            OUT = sum(node_out_wo(vtx),2);
            cons_res = cons_res + sum(abs(IN-OUT));
            cap_res  = cap_res  + sum(max(mv - r_cap(vtx), 0));
        end

        u = min(max(u .* (p  ./max(m0,eps0)).^damp, clip_min), clip_max);
        v = min(max(v .* (q  ./max(mT,eps0)).^damp, clip_min), clip_max);
        for vtx = middles
            mv = m_node(vtx);  rv = r_cap(vtx);
            wv = w(vtx) .* min(1, rv./max(mv,eps0)).^damp;
            w(vtx) = min(max(wv, clip_min), clip_max);
        end

        if sum(abs(m0-p)) < 2e-4 && sum(abs(mT-q)) < 2e-4 && cons_res < 2e-5 && cap_res < 2e-5
            break;
        end
    end
    fprintf('  Converged at iter %d  (bnd=%.2e  cap=%.2e)\n', ...
        it, sum(abs(m0-p))+sum(abs(mT-q)), cap_res);

    %% Boundary polish (w frozen)
    for kp = 1:12
        A0_paths(:) = 0;   AT_paths(:) = 0;
        for pid = 1:nP
            Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
            if L >= 1
                Bwo    = cell(L,1);  Bwo{L} = max(Kedge*v, eps0);
                for i = L-1:-1:1,  Bwo{i} = max(Kedge*(w(Pmid(i+1)).*Bwo{i+1}),eps0);  end
                Fwo    = cell(L,1);  Fwo{1} = max(Kedge'*u, eps0);
                for i = 2:L,        Fwo{i} = max(Kedge'*(w(Pmid(i-1)).*Fwo{i-1}),eps0); end
                A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
                AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
            else
                A0_paths(:,pid) = Kedge*v;  AT_paths(:,pid) = Kedge'*u;
            end
        end
        sumA0 = max(sum(A0_paths,2), eps0);
        u     = u .* (p ./ (u.*sumA0));

        for pid = 1:nP
            Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
            if L >= 1
                F1 = Kedge'*u;
                for i = 2:L,  F1 = Kedge'*(w(Pmid(i-1)).*F1);  end
                AT_paths(:,pid) = Kedge'*(w(Pmid(L)).*F1);
            else
                AT_paths(:,pid) = Kedge'*u;
            end
        end
        sumAT = max(sum(AT_paths,2), eps0);
        v     = v .* (q ./ (v.*sumAT));
    end

    sumA0 = max(sum(A0_paths,2), eps0);
    sumAT = max(sum(AT_paths,2), eps0);
    m0    = u .* sumA0;
    mT    = v .* sumAT;

    %% Collect node ribbons
    ribbons = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        A  = node_A_wo(vtx);
        Wv = w(vtx);
        ribbons(vtx) = sum(A .* (Wv.*ones(1,nP)), 2);
    end
end

%% ── Drawing helpers ────────────────────────────────────────────────────────

function drawLobe(x0, w, t, gain, col, alph, lw, side)
%DRAWLOBE  Filled density lobe to the left or right of x0.
    if strcmp(side,'left')
        x = x0 - gain*w(:);
    else
        x = x0 + gain*w(:);
    end
    fill([x0*ones(numel(t),1); flipud(x)], [t; flipud(t)], col, ...
        'FaceAlpha',alph,'EdgeColor','none');
    plot(x, t, 'Color',col, 'LineWidth',lw);
end

function drawArc(xL, xR, h, col, lw)
%DRAWARC  Quadratic Bezier arc curving below the t=0 axis.
    s  = linspace(0,1,80);
    xC = 0.5*(xL+xR);
    bx = (1-s).^2*xL + 2*(1-s).*s*xC + s.^2*xR;
    bz = 2*(-h)*(1-s).*s;
    plot(bx, bz, 'Color',col, 'LineWidth',lw);
end







% %% Path-graph: uniform vs time-varying capacity (dual-lobe 2-D plot)
% %
% %  Kernel and solver parameters are kept identical to the working original.
% %  The only additions are:
% %    1. A second Sinkhorn solve with sinusoidally modulated capacity.
% %    2. Interior nodes now show both solutions side-by-side
% %       (uniform → right, green;  time-varying → left, orange/red).
% 
% clear; clc; close all;
% 
% %% ── Network topology ─────────────────────────────────────────────────────
% nx = 4; ny = 3;
% node_id = @(ix,iy) (iy-1)*nx + ix;
% 
% v0 = node_id(1,1);
% vT = node_id(4,3);
% J1 = node_id(2,2);
% J2 = node_id(3,2);
% 
% pathA = [node_id(1,1), node_id(2,1), J1, J2, node_id(4,2), vT];
% pathB = [node_id(1,1), node_id(1,2), J1, J2, node_id(3,3), vT];
% pathC = [node_id(1,1), node_id(2,1), J1, J2, node_id(3,3), vT];
% paths = {pathA, pathB, pathC};
% 
% nodes_on_paths = unique([pathA(:); pathB(:); pathC(:)])';
% middles        = setdiff(nodes_on_paths, [v0 vT]);
% 
% % Display order (left to right)
% seq_nodes = [v0, node_id(2,1), node_id(1,2), J1, J2, node_id(4,2), node_id(3,3), vT];
% labels    = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','$v_{\mathcal T}$'};
% 
% % Order index for each interior node (used when building time-varying cap)
% interior_seq = seq_nodes(2:end-1);
% nodeIdx = containers.Map(num2cell(interior_seq), num2cell(1:numel(interior_seq)));
% 
% %% ── Time grid and marginals (unchanged from original) ────────────────────
% nt = 200;
% t  = linspace(0,1,nt)';
% Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);
% 
% p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07);  p = p/sum(p);
% q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06);  q = q/sum(q);
% 
% %% ── Kernel (unchanged from original: row-normalised, tau_gap smoothed) ───
% eps_ent    = 0.05;
% tau_gap    = 0.03;
% alpha_cost = 0.002;
% 
% [U,V] = ndgrid(t,t);
% mask  = (V >= U);
% Cedge = inf(nt);
% Cedge(mask) = 1./((V(mask)-U(mask)) + tau_gap);
% 
% Kedge = zeros(nt);
% Kedge(mask) = exp(-alpha_cost * Cedge(mask) / eps_ent);
% rowSum = sum(Kedge,2);  rowSum(rowSum==0) = 1;
% Kedge  = Kedge ./ rowSum;
% 
% %% ── Solver options (unchanged from original) ─────────────────────────────
% opts.eps0     = 1e-8;
% opts.damp     = 0.4;
% opts.clip_min = 1e-8;
% opts.clip_max = 1e8;
% opts.maxIter  = 4000;
% 
% %% ── Uniform capacity ─────────────────────────────────────────────────────
% cap_scalar = 0.009;
% r_capU = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles
%     r_capU(vtx) = cap_scalar * ones(nt,1);
% end
% 
% fprintf('Solving uniform-capacity case ...\n');
% [ribbonsU, m0U, mTU] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_capU);
% 
% %% ── Time-varying capacity (sinusoidal, mean = cap_scalar) ────────────────
% %  Build modulation relative to the uniform solution's ribbon at each node,
% %  then normalise so the mean capacity matches the uniform scalar.
% 
% r_capV = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles
%     k         = nodeIdx(vtx);                  % 1 … nMiddle
%     phase     = (k-1) * (pi/2);
%     modulation = 1 + 0.65 * sin(2*pi*t/0.60 + phase);
%     profile   = cap_scalar * max(modulation, 0.05);
%     % Normalise so mean(profile) == cap_scalar  (fair comparison)
%     profile   = profile * (cap_scalar * nt / sum(profile));
%     r_capV(vtx) = max(profile, 5e-4);
% end
% 
% fprintf('Solving time-varying-capacity case ...\n');
% [ribbonsV, m0V, mTV] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_capV);
% 
% %% ── Global normalisation (both solutions on the same scale) ──────────────
% epsn = 1e-16;
% allVals = [m0U; mTU; m0V; mTV];
% for vtx = middles
%     allVals = [allVals; ribbonsU(vtx); ribbonsV(vtx); ...
%                r_capU(vtx);            r_capV(vtx)];   %#ok<AGROW>
% end
% denG = max(allVals) + epsn;
% 
% %% ── Colours and layout ───────────────────────────────────────────────────
% xAll = linspace(0,1,numel(seq_nodes));
% xmap = containers.Map(num2cell(seq_nodes), num2cell(xAll));
% 
% cBlue    = [0.0000 0.4470 0.7410];
% cArrival = [0.6350 0.0780 0.1840];
% cUniform = [0.20   0.62   0.20  ];
% cVarying = [0.88   0.28   0.10  ];
% cRail    = [0.20   0.20   0.20  ];
% cEdge    = [0.38   0.38   0.38  ];
% 
% gain    = 0.12;    % half-width of each lobe
% arc_h   = 0.10;
% 
% %% ── Figure ───────────────────────────────────────────────────────────────
% figure(1); clf;
% set(gcf,'Color','w','Units','centimeters','Position',[2 2 20 10.5], ...
%     'PaperPositionMode','auto','Renderer','painters');
% hold on; box on; grid on;
% 
% % Vertical rails
% for i = 1:numel(seq_nodes)
%     x = xAll(i);
%     plot([x x],[0 1],'Color',cRail,'LineWidth',1.4);
%     plot(x, 0,'o','Color',cRail,'MarkerFaceColor',cRail,'MarkerSize',7);
% end
% 
% % Source (blue, left)
% drawLobe(xmap(v0),  m0U/denG, t, gain, cBlue,    0.35, 2.6, 'left');
% 
% % Interior nodes — uniform RIGHT, time-varying LEFT
% for i = 2:(numel(seq_nodes)-1)
%     vtx = seq_nodes(i);
%     x   = xAll(i);
% 
%     if isKey(ribbonsU, vtx)
%         mvU = ribbonsU(vtx) / denG;
%         mvV = ribbonsV(vtx) / denG;
%         cvU = r_capU(vtx)   / denG;
%         cvV = r_capV(vtx)   / denG;
%     else
%         mvU = zeros(nt,1);  mvV = zeros(nt,1);
%         cvU = zeros(nt,1);  cvV = zeros(nt,1);
%     end
% 
%     % Uniform — right side + dashed cap
%     drawLobe(x, mvU, t, gain, cUniform, 0.25, 2.0, 'right');
%     plot(x + gain*cvU, t, '--', 'Color',cUniform, 'LineWidth',1.3);
% 
%     % Time-varying — left side + dashed cap
%     drawLobe(x, mvV, t, gain, cVarying, 0.25, 2.0, 'left');
%     plot(x - gain*cvV, t, '--', 'Color',cVarying, 'LineWidth',1.3);
% end
% 
% % Sink (dark red, right)
% drawLobe(xmap(vT), mTU/denG, t, gain, cArrival, 0.35, 2.6, 'right');
% 
% % Bottom arcs
% for pid = 1:numel(paths)
%     Pth = paths{pid};
%     for k = 1:(numel(Pth)-1)
%         drawArc(xmap(Pth(k)), xmap(Pth(k+1)), arc_h, cEdge, 1.2);
%     end
% end
% 
% % Node labels
% for i = 1:numel(seq_nodes)
%     col = [0 0 0];
%     if i == 1,              col = cBlue;    end
%     if i == numel(seq_nodes), col = cArrival; end
%     text(xAll(i), -0.045, labels{i}, 'Interpreter','latex','Color',col, ...
%         'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','top');
% end
% 
% % Legend
% hU = plot(nan,nan,'Color',cUniform,'LineWidth',2.2);
% hV = plot(nan,nan,'Color',cVarying,'LineWidth',2.2);
% legend([hU hV],{'Uniform capacity','Time-varying capacity'}, ...
%     'Interpreter','latex','FontSize',12,'Location','northoutside', ...
%     'Orientation','horizontal','Box','off');
% 
% xlim([-0.09 1.09]);  ylim([-0.13 1.03]);
% set(gca,'XTick',[],'YTick',0:0.2:1,'TickLabelInterpreter','latex', ...
%     'FontName','Times New Roman','FontSize',14,'LineWidth',0.8);
% ylabel('$t$','Interpreter','latex','FontName','Times New Roman','FontSize',18);
% hold off;
% 
% exportgraphics(gcf,'path_uniform_vs_time_varying_capacity.pdf', ...
%     'ContentType','vector','BackgroundColor','white');
% exportgraphics(gcf,'path_uniform_vs_time_varying_capacity.png', ...
%     'Resolution',600,'BackgroundColor','white');
% 
% %% ── Capacity profiles ────────────────────────────────────────────────────
% figure(2); clf;
% set(gcf,'Color','w','Units','centimeters','Position',[2 2 14 8], ...
%     'PaperPositionMode','auto','Renderer','painters');
% hold on;
% cCap = lines(numel(interior_seq));
% for i = 1:numel(interior_seq)
%     vtx = interior_seq(i);
%     plot(t, r_capV(vtx), 'Color',cCap(i,:), 'LineWidth',1.8);
% end
% plot(t, cap_scalar*ones(nt,1), 'k--', 'LineWidth',1.8);
% lbs = arrayfun(@(k) labels{k+1}, 1:numel(interior_seq), 'UniformOutput',false);
% lbs{end+1} = 'Uniform';
% legend(lbs,'Interpreter','latex','Location','eastoutside','Box','off');
% xlabel('$t$','Interpreter','latex','FontSize',16);
% ylabel('Capacity $r_k(t)$','Interpreter','latex','FontSize',16);
% set(gca,'FontName','Times New Roman','FontSize',13,'TickLabelInterpreter','latex');
% grid on; box on; hold off;
% 
% exportgraphics(gcf,'path_time_varying_capacity_profiles.pdf', ...
%     'ContentType','vector','BackgroundColor','white');
% 
% %% ======================================================================
% %% Local functions
% %% ======================================================================
% 
% function [ribbons, m0, mT] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)
% %RUNPATHSINKHORN  Path-wise Sinkhorn with shared nodal capacity multipliers.
% %  Identical algorithm to the original working code, wrapped in a function
% %  so it can be called twice for the two capacity cases.
% 
%     eps0     = opts.eps0;
%     damp     = opts.damp;
%     clip_min = opts.clip_min;
%     clip_max = opts.clip_max;
%     maxIter  = opts.maxIter;
% 
%     u = ones(nt,1);
%     v = ones(nt,1);
%     w = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles,  w(vtx) = ones(nt,1);  end
% 
%     nP           = numel(paths);
%     A0_paths     = zeros(nt, nP);
%     AT_paths     = zeros(nt, nP);
%     node_A_wo    = containers.Map('KeyType','double','ValueType','any');
%     node_in_wo   = containers.Map('KeyType','double','ValueType','any');
%     node_out_wo  = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         node_A_wo(vtx)   = zeros(nt, nP);
%         node_in_wo(vtx)  = zeros(nt, nP);
%         node_out_wo(vtx) = zeros(nt, nP);
%     end
% 
%     %% Main iteration
%     for it = 1:maxIter
%         A0_paths(:) = 0;   AT_paths(:) = 0;
%         for vtx = middles
%             node_A_wo(vtx)  = 0*node_A_wo(vtx);
%             node_in_wo(vtx) = 0*node_in_wo(vtx);
%             node_out_wo(vtx)= 0*node_out_wo(vtx);
%         end
% 
%         for pid = 1:nP
%             Pth  = paths{pid};
%             Pmid = Pth(2:end-1);
%             L    = numel(Pmid);
% 
%             if L >= 1
%                 Bwo    = cell(L,1);
%                 Bwo{L} = max(Kedge * v, eps0);
%                 for k = L-1:-1:1
%                     Bwo{k} = max(Kedge * (w(Pmid(k+1)).*Bwo{k+1}), eps0);
%                 end
%                 Fwo    = cell(L,1);
%                 Fwo{1} = max(Kedge' * u, eps0);
%                 for k = 2:L
%                     Fwo{k} = max(Kedge' * (w(Pmid(k-1)).*Fwo{k-1}), eps0);
%                 end
%                 A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
%                 AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
% 
%                 for k = 1:L
%                     vtx = Pmid(k);
%                     if k == 1
%                         incoming = (Kedge'*u) .* Bwo{k};
%                     else
%                         incoming = (Kedge'*(w(Pmid(k-1)).*Fwo{k-1})) .* Bwo{k};
%                     end
%                     if k == L
%                         outgoing = Fwo{k} .* (Kedge*v);
%                     else
%                         outgoing = Fwo{k} .* (Kedge*(w(Pmid(k+1)).*Bwo{k+1}));
%                     end
%                     incoming = max(incoming, eps0);
%                     outgoing = max(outgoing, eps0);
% 
%                     A = node_A_wo(vtx);  A(:,pid) = A(:,pid) + Fwo{k}.*Bwo{k};  node_A_wo(vtx) = A;
%                     I = node_in_wo(vtx); I(:,pid) = I(:,pid) + incoming;          node_in_wo(vtx) = I;
%                     O = node_out_wo(vtx);O(:,pid) = O(:,pid) + outgoing;          node_out_wo(vtx) = O;
%                 end
%             else
%                 A0_paths(:,pid) = Kedge * v;
%                 AT_paths(:,pid) = Kedge' * u;
%             end
%         end
% 
%         sumA0 = max(sum(A0_paths,2), eps0);
%         sumAT = max(sum(AT_paths,2), eps0);
%         m0    = u .* sumA0;
%         mT    = v .* sumAT;
% 
%         cons_res = 0;  cap_res = 0;
%         m_node = containers.Map('KeyType','double','ValueType','any');
%         for vtx = middles
%             A_sum = max(sum(node_A_wo(vtx),2), eps0);
%             mv    = w(vtx) .* A_sum;
%             m_node(vtx) = mv;
%             IN  = sum(node_in_wo(vtx),2);
%             OUT = sum(node_out_wo(vtx),2);
%             cons_res = cons_res + sum(abs(IN-OUT));
%             cap_res  = cap_res  + sum(max(mv - r_cap(vtx), 0));
%         end
% 
%         u = min(max(u .* (p  ./max(m0,eps0)).^damp, clip_min), clip_max);
%         v = min(max(v .* (q  ./max(mT,eps0)).^damp, clip_min), clip_max);
%         for vtx = middles
%             mv = m_node(vtx);  rv = r_cap(vtx);
%             wv = w(vtx) .* min(1, rv./max(mv,eps0)).^damp;
%             w(vtx) = min(max(wv, clip_min), clip_max);
%         end
% 
%         if sum(abs(m0-p)) < 2e-4 && sum(abs(mT-q)) < 2e-4 && cons_res < 2e-5 && cap_res < 2e-5
%             break;
%         end
%     end
%     fprintf('  Converged at iter %d  (bnd=%.2e  cap=%.2e)\n', ...
%         it, sum(abs(m0-p))+sum(abs(mT-q)), cap_res);
% 
%     %% Boundary polish (w frozen)
%     for kp = 1:12
%         A0_paths(:) = 0;   AT_paths(:) = 0;
%         for pid = 1:nP
%             Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
%             if L >= 1
%                 Bwo    = cell(L,1);  Bwo{L} = max(Kedge*v, eps0);
%                 for i = L-1:-1:1,  Bwo{i} = max(Kedge*(w(Pmid(i+1)).*Bwo{i+1}),eps0);  end
%                 Fwo    = cell(L,1);  Fwo{1} = max(Kedge'*u, eps0);
%                 for i = 2:L,        Fwo{i} = max(Kedge'*(w(Pmid(i-1)).*Fwo{i-1}),eps0); end
%                 A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
%                 AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
%             else
%                 A0_paths(:,pid) = Kedge*v;  AT_paths(:,pid) = Kedge'*u;
%             end
%         end
%         sumA0 = max(sum(A0_paths,2), eps0);
%         u     = u .* (p ./ (u.*sumA0));
% 
%         for pid = 1:nP
%             Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
%             if L >= 1
%                 F1 = Kedge'*u;
%                 for i = 2:L,  F1 = Kedge'*(w(Pmid(i-1)).*F1);  end
%                 AT_paths(:,pid) = Kedge'*(w(Pmid(L)).*F1);
%             else
%                 AT_paths(:,pid) = Kedge'*u;
%             end
%         end
%         sumAT = max(sum(AT_paths,2), eps0);
%         v     = v .* (q ./ (v.*sumAT));
%     end
% 
%     sumA0 = max(sum(A0_paths,2), eps0);
%     sumAT = max(sum(AT_paths,2), eps0);
%     m0    = u .* sumA0;
%     mT    = v .* sumAT;
% 
%     %% Collect node ribbons
%     ribbons = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         A  = node_A_wo(vtx);
%         Wv = w(vtx);
%         ribbons(vtx) = sum(A .* (Wv.*ones(1,nP)), 2);
%     end
% end
% 
% %% ── Drawing helpers ────────────────────────────────────────────────────────
% 
% function drawLobe(x0, w, t, gain, col, alph, lw, side)
% %DRAWLOBE  Filled density lobe to the left or right of x0.
%     if strcmp(side,'left')
%         x = x0 - gain*w(:);
%     else
%         x = x0 + gain*w(:);
%     end
%     fill([x0*ones(numel(t),1); flipud(x)], [t; flipud(t)], col, ...
%         'FaceAlpha',alph,'EdgeColor','none');
%     plot(x, t, 'Color',col, 'LineWidth',lw);
% end
% 
% function drawArc(xL, xR, h, col, lw)
% %DRAWARC  Quadratic Bezier arc curving below the t=0 axis.
%     s  = linspace(0,1,80);
%     xC = 0.5*(xL+xR);
%     bx = (1-s).^2*xL + 2*(1-s).*s*xC + s.^2*xR;
%     bz = 2*(-h)*(1-s).*s;
%     plot(bx, bz, 'Color',col, 'LineWidth',lw);
% end
% 
% 


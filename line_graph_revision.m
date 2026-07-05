%% Line graph: uniform vs time-varying capacity
%
%  Single straight-line path  v0 → v1 → … → v6 → vT  (8 nodes, 7 edges).
%  Kernel and solver are identical to path_sinkhorn_2dplot.m.
%  Plot: uniform capacity → right lobes (green)
%        time-varying capacity → left lobes (orange/red)

clear; clc; close all;

%% ── Network: single line path ────────────────────────────────────────────
nNodes  = 8;
nMiddle = nNodes - 2;          % v1 … v6
nEdges  = nNodes - 1;
% labels  = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','$v_{\mathcal T}$'};
labels  = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','$v_7'};

v0      = 1;
vT      = nNodes;
middles = 2:(nNodes-1);        % [2 3 4 5 6 7]

pathLine = {1:nNodes};         % single path as a cell array (matches solver API)

%% ── Time grid & marginals (same as path_sinkhorn_2dplot) ─────────────────
nt = 200;
t  = linspace(0,1,nt)';
Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);

p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07);  p = p/sum(p);
q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06);  q = q/sum(q);

%% ── Kernel: linear transit-time cost  C(s,t) = t-s ──────────────────────
%  Old cost  1/(t-s+tau)  gave *small* cost for slow transit, so mass spread
%  across the full feasible window.  Linear cost penalises slow transit and
%  concentrates each node's profile near its natural crossing time.
%  eps_ent ≈ 0.065 sets mean hop time ≈ 0.065, giving 7 hops ≈ 0.455
%  which matches the source-to-sink travel distance (~0.80-0.35 = 0.45).

eps_ent = 0.065;

[U,V] = ndgrid(t,t);
causalMask = (V > U);                          % strict: no zero-time hops

Kedge = zeros(nt);
Kedge(causalMask) = exp(-(V(causalMask)-U(causalMask)) / eps_ent);
rowSum = sum(Kedge,2);  rowSum(rowSum==0) = 1;
Kedge  = Kedge ./ rowSum;

%% ── Solver options (identical to path_sinkhorn_2dplot) ───────────────────
opts.eps0     = 1e-8;
opts.damp     = 0.4;
opts.clip_min = 1e-8;
opts.clip_max = 1e8;
opts.maxIter  = 4000;

%% ── Uniform capacity ─────────────────────────────────────────────────────
cap_scalar = 0.01;
r_capU = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    r_capU(vtx) = cap_scalar * ones(nt,1);
end

fprintf('Solving uniform-capacity case ...\n');
[ribbonsU, m0U, mTU] = runPathSinkhorn(Kedge, p, q, pathLine, middles, nt, opts, r_capU);

%% ── Time-varying capacity (sinusoidal, mean = cap_scalar) ────────────────
r_capV = containers.Map('KeyType','double','ValueType','any');
for k = 1:nMiddle
    vtx   = middles(k);
    phase = (k-1) * (pi/2);
    modulation = 1 + 0.65 * sin(2*pi*t/0.60 + phase);
    profile    = cap_scalar * max(modulation, 0.05);
    profile    = profile * (cap_scalar * nt / sum(profile));  % normalise mean
    r_capV(vtx) = max(profile, 5e-4);
end

fprintf('Solving time-varying-capacity case ...\n');
[ribbonsV, m0V, mTV] = runPathSinkhorn(Kedge, p, q, pathLine, middles, nt, opts, r_capV);

%% ── Global normalisation ─────────────────────────────────────────────────
allVals = [m0U; mTU; m0V; mTV];
for vtx = middles
    allVals = [allVals; ribbonsU(vtx); ribbonsV(vtx); ...
               r_capU(vtx);            r_capV(vtx)];  %#ok<AGROW>
end
denG = max(allVals) + 1e-16;

%% ── Colours & layout ─────────────────────────────────────────────────────
xAll     = linspace(0,1,nNodes);
cBlue    = [0.0000 0.4470 0.7410];
cArrival = [0.6350 0.0780 0.1840];
cUniform = [0.20   0.62   0.20  ];
cVarying = [0.88   0.28   0.10  ];
cRail    = [0.20   0.20   0.20  ];
cEdge    = [0.38   0.38   0.38  ];

gain   = 0.09;
arc_h  = 0.075;
gainBd = 0.12;   % slightly wider lobe for source / sink

%% ── Figure ───────────────────────────────────────────────────────────────
figure(1); clf;
set(gcf,'Color','w','Units','centimeters','Position',[2 2 20 10.5], ...
    'PaperPositionMode','auto','Renderer','painters');
hold on; box off; grid on;

% Vertical rails and node dots
for i = 1:nNodes
    plot([xAll(i) xAll(i)],[0 1],'Color',cRail,'LineWidth',1.4);
    plot(xAll(i),0,'o','Color',cRail,'MarkerFaceColor',cRail,'MarkerSize',7);
end

% Source (blue, left)
drawLobe(xAll(1), m0U/denG, t, gainBd, cBlue, 0.35, 2.6, 'left');

% Visual t-offsets for v1 and v2 — shift their lobes upward for clarity
tOffsets = [0.08, 0.04, 0, 0, 0, 0];

% Interior nodes
for k = 1:nMiddle
    vtx   = middles(k);
    x     = xAll(k+1);
    tDraw = t + tOffsets(k);   % shifted only for k=1,2; zero elsewhere

    mvU = ribbonsU(vtx) / denG;
    mvV = ribbonsV(vtx) / denG;
    cvU = r_capU(vtx)   / denG;
    cvV = r_capV(vtx)   / denG;

    % Uniform → right + dashed cap
    drawLobe(x, mvU, tDraw, gain, cUniform, 0.25, 2.0, 'right');
    plot(x + gain*cvU, tDraw, '--', 'Color',cUniform, 'LineWidth',1.3);

    % Time-varying → left + dashed cap
    drawLobe(x, mvV, tDraw, gain, cVarying, 0.25, 2.0, 'left');
    plot(x - gain*cvV, tDraw, '--', 'Color',cVarying, 'LineWidth',1.3);
end

% Sink (dark red, right)
drawLobe(xAll(end), mTU/denG, t, gainBd, cArrival, 0.35, 2.6, 'right');

% Bottom arcs (single path = adjacent nodes)
for k = 1:nEdges
    drawArc(xAll(k), xAll(k+1), arc_h, cEdge, 1.2);
end

% Node labels
for i = 1:nNodes
    col = [0 0 0];
    if i == 1,      col = cBlue;    end
    if i == nNodes, col = cArrival; end
    text(xAll(i),-0.045,labels{i},'Interpreter','latex','Color',col, ...
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

% Remove bottom x-axis line; keep y-axis spine and tick labels
ax1 = gca;
ax1.XAxis.Axle.Visible = 'off';

exportgraphics(gcf,'line_uniform_vs_time_varying_capacity.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(gcf,'line_uniform_vs_time_varying_capacity.png', ...
    'Resolution',600,'BackgroundColor','white');

%% ── Capacity profiles figure ─────────────────────────────────────────────
figure(2); clf;
set(gcf,'Color','w','Units','centimeters','Position',[2 2 14 8], ...
    'PaperPositionMode','auto','Renderer','painters');
hold on;
cCap = lines(nMiddle);
for k = 1:nMiddle
    plot(t, r_capV(middles(k)), 'Color',cCap(k,:), 'LineWidth',1.8);
end
plot(t, cap_scalar*ones(nt,1), 'k--', 'LineWidth',1.8);
lbs = arrayfun(@(k) labels{k+1}, 1:nMiddle, 'UniformOutput',false);
lbs{end+1} = 'Uniform';
legend(lbs,'Interpreter','latex','Location','eastoutside','Box','off');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('Capacity $r_k(t)$','Interpreter','latex','FontSize',16);
set(gca,'FontName','Times New Roman','FontSize',13,'TickLabelInterpreter','latex');
grid on; box on; hold off;

exportgraphics(gcf,'line_time_varying_capacity_profiles.pdf', ...
    'ContentType','vector','BackgroundColor','white');

%% ======================================================================
%% Local functions (identical to path_sinkhorn_2dplot.m)
%% ======================================================================

function [ribbons, m0, mT] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)
    eps0     = opts.eps0;
    damp     = opts.damp;
    clip_min = opts.clip_min;
    clip_max = opts.clip_max;
    maxIter  = opts.maxIter;

    nP = numel(paths);
    u  = ones(nt,1);
    v  = ones(nt,1);
    w  = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles,  w(vtx) = ones(nt,1);  end

    A0_paths    = zeros(nt,nP);
    AT_paths    = zeros(nt,nP);
    node_A_wo   = containers.Map('KeyType','double','ValueType','any');
    node_in_wo  = containers.Map('KeyType','double','ValueType','any');
    node_out_wo = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        node_A_wo(vtx)   = zeros(nt,nP);
        node_in_wo(vtx)  = zeros(nt,nP);
        node_out_wo(vtx) = zeros(nt,nP);
    end

    for it = 1:maxIter
        A0_paths(:) = 0;  AT_paths(:) = 0;
        for vtx = middles
            node_A_wo(vtx)   = 0*node_A_wo(vtx);
            node_in_wo(vtx)  = 0*node_in_wo(vtx);
            node_out_wo(vtx) = 0*node_out_wo(vtx);
        end

        for pid = 1:nP
            Pth  = paths{pid};
            Pmid = Pth(2:end-1);
            L    = numel(Pmid);

            if L >= 1
                Bwo    = cell(L,1);
                Bwo{L} = max(Kedge*v, eps0);
                for k = L-1:-1:1
                    Bwo{k} = max(Kedge*(w(Pmid(k+1)).*Bwo{k+1}), eps0);
                end
                Fwo    = cell(L,1);
                Fwo{1} = max(Kedge'*u, eps0);
                for k = 2:L
                    Fwo{k} = max(Kedge'*(w(Pmid(k-1)).*Fwo{k-1}), eps0);
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
                    incoming = max(incoming,eps0);
                    outgoing = max(outgoing,eps0);
                    A = node_A_wo(vtx);  A(:,pid) = A(:,pid)+Fwo{k}.*Bwo{k};  node_A_wo(vtx) = A;
                    I = node_in_wo(vtx); I(:,pid) = I(:,pid)+incoming;          node_in_wo(vtx) = I;
                    O = node_out_wo(vtx);O(:,pid) = O(:,pid)+outgoing;          node_out_wo(vtx) = O;
                end
            else
                A0_paths(:,pid) = Kedge*v;
                AT_paths(:,pid) = Kedge'*u;
            end
        end

        sumA0 = max(sum(A0_paths,2),eps0);
        sumAT = max(sum(AT_paths,2),eps0);
        m0    = u.*sumA0;
        mT    = v.*sumAT;

        cons_res = 0;  cap_res = 0;
        m_node = containers.Map('KeyType','double','ValueType','any');
        for vtx = middles
            A_sum = max(sum(node_A_wo(vtx),2),eps0);
            mv    = w(vtx).*A_sum;
            m_node(vtx) = mv;
            IN  = sum(node_in_wo(vtx),2);
            OUT = sum(node_out_wo(vtx),2);
            cons_res = cons_res + sum(abs(IN-OUT));
            cap_res  = cap_res  + sum(max(mv-r_cap(vtx),0));
        end

        u = min(max(u.*(p ./max(m0,eps0)).^damp, clip_min), clip_max);
        v = min(max(v.*(q ./max(mT,eps0)).^damp, clip_min), clip_max);
        for vtx = middles
            mv = m_node(vtx);  rv = r_cap(vtx);
            w(vtx) = min(max(w(vtx).*min(1,rv./max(mv,eps0)).^damp, clip_min), clip_max);
        end

        if sum(abs(m0-p))<2e-4 && sum(abs(mT-q))<2e-4 && cons_res<2e-5 && cap_res<2e-5
            break;
        end
    end
    fprintf('  Converged at iter %d  (bnd=%.2e  cap=%.2e)\n', ...
        it, sum(abs(m0-p))+sum(abs(mT-q)), cap_res);

    % Boundary polish (w frozen)
    for kp = 1:12
        A0_paths(:) = 0;  AT_paths(:) = 0;
        for pid = 1:nP
            Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
            if L >= 1
                Bwo    = cell(L,1);  Bwo{L} = max(Kedge*v,eps0);
                for i = L-1:-1:1,  Bwo{i} = max(Kedge*(w(Pmid(i+1)).*Bwo{i+1}),eps0);  end
                Fwo    = cell(L,1);  Fwo{1} = max(Kedge'*u,eps0);
                for i = 2:L,        Fwo{i} = max(Kedge'*(w(Pmid(i-1)).*Fwo{i-1}),eps0); end
                A0_paths(:,pid) = Kedge *(w(Pmid(1)).*Bwo{1});
                AT_paths(:,pid) = Kedge'*(w(Pmid(L)).*Fwo{L});
            else
                A0_paths(:,pid) = Kedge*v;  AT_paths(:,pid) = Kedge'*u;
            end
        end
        sumA0 = max(sum(A0_paths,2),eps0);
        u     = u.*(p./(u.*sumA0));
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
        sumAT = max(sum(AT_paths,2),eps0);
        v     = v.*(q./(v.*sumAT));
    end

    sumA0 = max(sum(A0_paths,2),eps0);
    sumAT = max(sum(AT_paths,2),eps0);
    m0    = u.*sumA0;
    mT    = v.*sumAT;

    ribbons = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        A  = node_A_wo(vtx);
        Wv = w(vtx);
        ribbons(vtx) = sum(A.*(Wv.*ones(1,nP)), 2);
    end
end

function drawLobe(x0, w, t, gain, col, alph, lw, side)
    if strcmp(side,'left'),  x = x0 - gain*w(:);
    else,                    x = x0 + gain*w(:);
    end
    fill([x0*ones(numel(t),1); flipud(x)],[t; flipud(t)],col, ...
        'FaceAlpha',alph,'EdgeColor','none');
    plot(x, t, 'Color',col, 'LineWidth',lw);
end

function drawArc(xL, xR, h, col, lw)
    s  = linspace(0,1,80);
    xC = 0.5*(xL+xR);
    bx = (1-s).^2*xL + 2*(1-s).*s*xC + s.^2*xR;
    bz = 2*(-h)*(1-s).*s;
    plot(bx, bz, 'Color',col, 'LineWidth',lw);
end







% %% Line graph: uniform vs time-varying capacity
% %
% %  Single straight-line path  v0 → v1 → … → v6 → vT  (8 nodes, 7 edges).
% %  Kernel and solver are identical to path_sinkhorn_2dplot.m.
% %  Plot: uniform capacity → right lobes (green)
% %        time-varying capacity → left lobes (orange/red)
% 
% clear; clc; close all;
% 
% %% ── Network: single line path ────────────────────────────────────────────
% nNodes  = 8;
% nMiddle = nNodes - 2;          % v1 … v6
% nEdges  = nNodes - 1;
% labels  = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','$v_{\mathcal T}$'};
% 
% v0      = 1;
% vT      = nNodes;
% middles = 2:(nNodes-1);        % [2 3 4 5 6 7]
% 
% pathLine = {1:nNodes};         % single path as a cell array (matches solver API)
% 
% %% ── Time grid & marginals (same as path_sinkhorn_2dplot) ─────────────────
% nt = 200;
% t  = linspace(0,1,nt)';
% Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);
% 
% p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07);  p = p/sum(p);
% q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06);  q = q/sum(q);
% 
% %% ── Kernel (identical to path_sinkhorn_2dplot) ───────────────────────────
% eps_ent    = 0.05;
% tau_gap    = 0.03;
% alpha_cost = 0.002;
% 
% [U,V] = ndgrid(t,t);
% causalMask = (V >= U);
% Cedge = inf(nt);
% Cedge(causalMask) = 1./((V(causalMask)-U(causalMask)) + tau_gap);
% 
% Kedge = zeros(nt);
% Kedge(causalMask) = exp(-alpha_cost * Cedge(causalMask) / eps_ent);
% rowSum = sum(Kedge,2);  rowSum(rowSum==0) = 1;
% Kedge  = Kedge ./ rowSum;
% 
% %% ── Solver options (identical to path_sinkhorn_2dplot) ───────────────────
% opts.eps0     = 1e-8;
% opts.damp     = 0.4;
% opts.clip_min = 1e-8;
% opts.clip_max = 1e8;
% opts.maxIter  = 4000;
% 
% %% ── Uniform capacity ─────────────────────────────────────────────────────
% cap_scalar = 0.015;
% r_capU = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles
%     r_capU(vtx) = cap_scalar * ones(nt,1);
% end
% 
% fprintf('Solving uniform-capacity case ...\n');
% [ribbonsU, m0U, mTU] = runPathSinkhorn(Kedge, p, q, pathLine, middles, nt, opts, r_capU);
% 
% %% ── Time-varying capacity (sinusoidal, mean = cap_scalar) ────────────────
% r_capV = containers.Map('KeyType','double','ValueType','any');
% for k = 1:nMiddle
%     vtx   = middles(k);
%     phase = (k-1) * (pi/2);
%     modulation = 1 + 0.65 * sin(2*pi*t/0.60 + phase);
%     profile    = cap_scalar * max(modulation, 0.05);
%     profile    = profile * (cap_scalar * nt / sum(profile));  % normalise mean
%     r_capV(vtx) = max(profile, 5e-4);
% end
% 
% fprintf('Solving time-varying-capacity case ...\n');
% [ribbonsV, m0V, mTV] = runPathSinkhorn(Kedge, p, q, pathLine, middles, nt, opts, r_capV);
% 
% %% ── Global normalisation ─────────────────────────────────────────────────
% allVals = [m0U; mTU; m0V; mTV];
% for vtx = middles
%     allVals = [allVals; ribbonsU(vtx); ribbonsV(vtx); ...
%                r_capU(vtx);            r_capV(vtx)];  %#ok<AGROW>
% end
% denG = max(allVals) + 1e-16;
% 
% %% ── Colours & layout ─────────────────────────────────────────────────────
% xAll     = linspace(0,1,nNodes);
% cBlue    = [0.0000 0.4470 0.7410];
% cArrival = [0.6350 0.0780 0.1840];
% cUniform = [0.20   0.62   0.20  ];
% cVarying = [0.88   0.28   0.10  ];
% cRail    = [0.20   0.20   0.20  ];
% cEdge    = [0.38   0.38   0.38  ];
% 
% gain   = 0.09;
% arc_h  = 0.075;
% gainBd = 0.12;   % slightly wider lobe for source / sink
% 
% %% ── Figure ───────────────────────────────────────────────────────────────
% figure(1); clf;
% set(gcf,'Color','w','Units','centimeters','Position',[2 2 20 10.5], ...
%     'PaperPositionMode','auto','Renderer','painters');
% hold on; box on; grid on;
% 
% % Vertical rails and node dots
% for i = 1:nNodes
%     plot([xAll(i) xAll(i)],[0 1],'Color',cRail,'LineWidth',1.4);
%     plot(xAll(i),0,'o','Color',cRail,'MarkerFaceColor',cRail,'MarkerSize',7);
% end
% 
% % Source (blue, left)
% drawLobe(xAll(1), m0U/denG, t, gainBd, cBlue, 0.35, 2.6, 'left');
% 
% % Interior nodes
% for k = 1:nMiddle
%     vtx = middles(k);
%     x   = xAll(k+1);
% 
%     mvU = ribbonsU(vtx) / denG;
%     mvV = ribbonsV(vtx) / denG;
%     cvU = r_capU(vtx)   / denG;
%     cvV = r_capV(vtx)   / denG;
% 
%     % Uniform → right + dashed cap
%     drawLobe(x, mvU, t, gain, cUniform, 0.25, 2.0, 'right');
%     plot(x + gain*cvU, t, '--', 'Color',cUniform, 'LineWidth',1.3);
% 
%     % Time-varying → left + dashed cap
%     drawLobe(x, mvV, t, gain, cVarying, 0.25, 2.0, 'left');
%     plot(x - gain*cvV, t, '--', 'Color',cVarying, 'LineWidth',1.3);
% end
% 
% % Sink (dark red, right)
% drawLobe(xAll(end), mTU/denG, t, gainBd, cArrival, 0.35, 2.6, 'right');
% 
% % Bottom arcs (single path = adjacent nodes)
% for k = 1:nEdges
%     drawArc(xAll(k), xAll(k+1), arc_h, cEdge, 1.2);
% end
% 
% % Node labels
% for i = 1:nNodes
%     col = [0 0 0];
%     if i == 1,      col = cBlue;    end
%     if i == nNodes, col = cArrival; end
%     text(xAll(i),-0.045,labels{i},'Interpreter','latex','Color',col, ...
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
% exportgraphics(gcf,'line_uniform_vs_time_varying_capacity.pdf', ...
%     'ContentType','vector','BackgroundColor','white');
% exportgraphics(gcf,'line_uniform_vs_time_varying_capacity.png', ...
%     'Resolution',600,'BackgroundColor','white');
% 
% %% ── Capacity profiles figure ─────────────────────────────────────────────
% figure(2); clf;
% set(gcf,'Color','w','Units','centimeters','Position',[2 2 14 8], ...
%     'PaperPositionMode','auto','Renderer','painters');
% hold on;
% cCap = lines(nMiddle);
% for k = 1:nMiddle
%     plot(t, r_capV(middles(k)), 'Color',cCap(k,:), 'LineWidth',1.8);
% end
% plot(t, cap_scalar*ones(nt,1), 'k--', 'LineWidth',1.8);
% lbs = arrayfun(@(k) labels{k+1}, 1:nMiddle, 'UniformOutput',false);
% lbs{end+1} = 'Uniform';
% legend(lbs,'Interpreter','latex','Location','eastoutside','Box','off');
% xlabel('$t$','Interpreter','latex','FontSize',16);
% ylabel('Capacity $r_k(t)$','Interpreter','latex','FontSize',16);
% set(gca,'FontName','Times New Roman','FontSize',13,'TickLabelInterpreter','latex');
% grid on; box on; hold off;
% 
% exportgraphics(gcf,'line_time_varying_capacity_profiles.pdf', ...
%     'ContentType','vector','BackgroundColor','white');
% 
% %% ======================================================================
% %% Local functions (identical to path_sinkhorn_2dplot.m)
% %% ======================================================================
% 
% function [ribbons, m0, mT] = runPathSinkhorn(Kedge, p, q, paths, middles, nt, opts, r_cap)
%     eps0     = opts.eps0;
%     damp     = opts.damp;
%     clip_min = opts.clip_min;
%     clip_max = opts.clip_max;
%     maxIter  = opts.maxIter;
% 
%     nP = numel(paths);
%     u  = ones(nt,1);
%     v  = ones(nt,1);
%     w  = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles,  w(vtx) = ones(nt,1);  end
% 
%     A0_paths    = zeros(nt,nP);
%     AT_paths    = zeros(nt,nP);
%     node_A_wo   = containers.Map('KeyType','double','ValueType','any');
%     node_in_wo  = containers.Map('KeyType','double','ValueType','any');
%     node_out_wo = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         node_A_wo(vtx)   = zeros(nt,nP);
%         node_in_wo(vtx)  = zeros(nt,nP);
%         node_out_wo(vtx) = zeros(nt,nP);
%     end
% 
%     for it = 1:maxIter
%         A0_paths(:) = 0;  AT_paths(:) = 0;
%         for vtx = middles
%             node_A_wo(vtx)   = 0*node_A_wo(vtx);
%             node_in_wo(vtx)  = 0*node_in_wo(vtx);
%             node_out_wo(vtx) = 0*node_out_wo(vtx);
%         end
% 
%         for pid = 1:nP
%             Pth  = paths{pid};
%             Pmid = Pth(2:end-1);
%             L    = numel(Pmid);
% 
%             if L >= 1
%                 Bwo    = cell(L,1);
%                 Bwo{L} = max(Kedge*v, eps0);
%                 for k = L-1:-1:1
%                     Bwo{k} = max(Kedge*(w(Pmid(k+1)).*Bwo{k+1}), eps0);
%                 end
%                 Fwo    = cell(L,1);
%                 Fwo{1} = max(Kedge'*u, eps0);
%                 for k = 2:L
%                     Fwo{k} = max(Kedge'*(w(Pmid(k-1)).*Fwo{k-1}), eps0);
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
%                     incoming = max(incoming,eps0);
%                     outgoing = max(outgoing,eps0);
%                     A = node_A_wo(vtx);  A(:,pid) = A(:,pid)+Fwo{k}.*Bwo{k};  node_A_wo(vtx) = A;
%                     I = node_in_wo(vtx); I(:,pid) = I(:,pid)+incoming;          node_in_wo(vtx) = I;
%                     O = node_out_wo(vtx);O(:,pid) = O(:,pid)+outgoing;          node_out_wo(vtx) = O;
%                 end
%             else
%                 A0_paths(:,pid) = Kedge*v;
%                 AT_paths(:,pid) = Kedge'*u;
%             end
%         end
% 
%         sumA0 = max(sum(A0_paths,2),eps0);
%         sumAT = max(sum(AT_paths,2),eps0);
%         m0    = u.*sumA0;
%         mT    = v.*sumAT;
% 
%         cons_res = 0;  cap_res = 0;
%         m_node = containers.Map('KeyType','double','ValueType','any');
%         for vtx = middles
%             A_sum = max(sum(node_A_wo(vtx),2),eps0);
%             mv    = w(vtx).*A_sum;
%             m_node(vtx) = mv;
%             IN  = sum(node_in_wo(vtx),2);
%             OUT = sum(node_out_wo(vtx),2);
%             cons_res = cons_res + sum(abs(IN-OUT));
%             cap_res  = cap_res  + sum(max(mv-r_cap(vtx),0));
%         end
% 
%         u = min(max(u.*(p ./max(m0,eps0)).^damp, clip_min), clip_max);
%         v = min(max(v.*(q ./max(mT,eps0)).^damp, clip_min), clip_max);
%         for vtx = middles
%             mv = m_node(vtx);  rv = r_cap(vtx);
%             w(vtx) = min(max(w(vtx).*min(1,rv./max(mv,eps0)).^damp, clip_min), clip_max);
%         end
% 
%         if sum(abs(m0-p))<2e-4 && sum(abs(mT-q))<2e-4 && cons_res<2e-5 && cap_res<2e-5
%             break;
%         end
%     end
%     fprintf('  Converged at iter %d  (bnd=%.2e  cap=%.2e)\n', ...
%         it, sum(abs(m0-p))+sum(abs(mT-q)), cap_res);
% 
%     % Boundary polish (w frozen)
%     for kp = 1:12
%         A0_paths(:) = 0;  AT_paths(:) = 0;
%         for pid = 1:nP
%             Pth  = paths{pid};  Pmid = Pth(2:end-1);  L = numel(Pmid);
%             if L >= 1
%                 Bwo    = cell(L,1);  Bwo{L} = max(Kedge*v,eps0);
%                 for i = L-1:-1:1,  Bwo{i} = max(Kedge*(w(Pmid(i+1)).*Bwo{i+1}),eps0);  end
%                 Fwo    = cell(L,1);  Fwo{1} = max(Kedge'*u,eps0);
%                 for i = 2:L,        Fwo{i} = max(Kedge'*(w(Pmid(i-1)).*Fwo{i-1}),eps0); end
%                 A0_paths(:,pid) = Kedge *(w(Pmid(1)).*Bwo{1});
%                 AT_paths(:,pid) = Kedge'*(w(Pmid(L)).*Fwo{L});
%             else
%                 A0_paths(:,pid) = Kedge*v;  AT_paths(:,pid) = Kedge'*u;
%             end
%         end
%         sumA0 = max(sum(A0_paths,2),eps0);
%         u     = u.*(p./(u.*sumA0));
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
%         sumAT = max(sum(AT_paths,2),eps0);
%         v     = v.*(q./(v.*sumAT));
%     end
% 
%     sumA0 = max(sum(A0_paths,2),eps0);
%     sumAT = max(sum(AT_paths,2),eps0);
%     m0    = u.*sumA0;
%     mT    = v.*sumAT;
% 
%     ribbons = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         A  = node_A_wo(vtx);
%         Wv = w(vtx);
%         ribbons(vtx) = sum(A.*(Wv.*ones(1,nP)), 2);
%     end
% end
% 
% function drawLobe(x0, w, t, gain, col, alph, lw, side)
%     if strcmp(side,'left'),  x = x0 - gain*w(:);
%     else,                    x = x0 + gain*w(:);
%     end
%     fill([x0*ones(numel(t),1); flipud(x)],[t; flipud(t)],col, ...
%         'FaceAlpha',alph,'EdgeColor','none');
%     plot(x, t, 'Color',col, 'LineWidth',lw);
% end
% 
% function drawArc(xL, xR, h, col, lw)
%     s  = linspace(0,1,80);
%     xC = 0.5*(xL+xR);
%     bx = (1-s).^2*xL + 2*(1-s).*s*xC + s.^2*xR;
%     bz = 2*(-h)*(1-s).*s;
%     plot(bx, bz, 'Color',col, 'LineWidth',lw);
% end
%% 2D one-line chain with bottom arcs and one-sided distributions (fixed)
clear; clc; close all;

%% -------- Grid & admissible routes (for connectivity) --------------------
nx = 4; ny = 3;
node_id = @(ix,iy) (iy-1)*nx + ix;

v0  = node_id(1,1);
vT  = node_id(4,3);
J1  = node_id(2,2);
J2  = node_id(3,2);

% Paths exactly as specified
pathA = [ node_id(1,1), node_id(2,1), J1, J2, node_id(4,2), vT ];
pathB = [ node_id(1,1), node_id(1,2), J1, J2, node_id(3,3), vT ];
pathC = [ node_id(1,1), node_id(2,1), J1, J2, node_id(3,3), vT ];

paths = {pathA, pathB, pathC};
nodes_on_paths = unique([pathA(:); pathB(:); pathC(:)])';
middles = setdiff(nodes_on_paths, [v0 vT]);

%% -------- Time grid & boundary marginals ---------------------------------
nt = 200;  
t = linspace(0,1,nt)';

Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);
p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07); p = p/sum(p);
q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06); q = q/sum(q);

%% -------- Nearly-causal entropic kernel (row-normalized) -----------------
eps_ent   = 0.02;
tau_gap   = 0.03;
alpha_cost= 0.002;

[U,V] = ndgrid(t,t);
mask  = (V >= U);                           % allow diagonal
Cedge = inf(nt);  
Cedge(mask) = 1./((V(mask)-U(mask)) + tau_gap);
Kedge = zeros(nt); 
Kedge(mask) = exp(-alpha_cost * Cedge(mask) / eps_ent);
rowSum = sum(Kedge,2); rowSum(rowSum==0)=1;  
Kedge = Kedge ./ rowSum;                    % stochastic rows

%% -------- Uniform capacity at every interior node ------------------------
cap_scalar = 0.007;
r_cap = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    r_cap(vtx) = cap_scalar*ones(nt,1);
end

%% -------- Path-wise Sinkhorn with shared nodal multipliers ---------------
eps0     = 1e-8;
damp     = 0.4;
clip_min = 1e-8;
clip_max = 1e8;

u = ones(nt,1);
v = ones(nt,1);

w = containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    w(vtx) = ones(nt,1);
end

A0_paths = zeros(nt, numel(paths));
AT_paths = zeros(nt, numel(paths));

node_A_wo  = containers.Map('KeyType','double','ValueType','any');
node_in_wo = containers.Map('KeyType','double','ValueType','any');
node_out_wo= containers.Map('KeyType','double','ValueType','any');
for vtx = middles
    node_A_wo(vtx)  = zeros(nt, numel(paths));
    node_in_wo(vtx) = zeros(nt, numel(paths));
    node_out_wo(vtx)= zeros(nt, numel(paths));
end

histE0   = [];
histET   = [];
histCONS = [];
histCAP  = [];
maxIter  = 4000;

for it = 1:maxIter
    A0_paths(:) = 0;
    AT_paths(:) = 0;
    for vtx = middles
        node_A_wo(vtx)  = 0*node_A_wo(vtx);
        node_in_wo(vtx) = 0*node_in_wo(vtx);
        node_out_wo(vtx)= 0*node_out_wo(vtx);
    end

    for pid = 1:numel(paths)
        Pth  = paths{pid}; 
        Pmid = Pth(2:end-1); 
        L    = numel(Pmid);

        if L >= 1
            % Backward (without local w)
            Bwo      = cell(L,1);  
            Bwo{L}   = Kedge * v;  
            Bwo{L}   = max(Bwo{L},eps0);
            for k = L-1:-1:1
                Bwo{k} = Kedge * (w(Pmid(k+1)).*Bwo{k+1});
                Bwo{k} = max(Bwo{k},eps0);
            end

            % Forward (without local w)
            Fwo      = cell(L,1);  
            Fwo{1}   = Kedge' * u; 
            Fwo{1}   = max(Fwo{1},eps0);
            for k = 2:L
                Fwo{k} = Kedge' * (w(Pmid(k-1)).*Fwo{k-1});
                Fwo{k} = max(Fwo{k},eps0);
            end
        end

        % Source-side and sink-side messages
        if L >= 1
            A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
            AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
        else
            A0_paths(:,pid) = Kedge * v;
            AT_paths(:,pid) = Kedge' * u;
        end

        % Node crossing contributions (wo local w_v)
        if L >= 1
            for k = 1:L
                vtx = Pmid(k);

                % incoming at node time
                if k == 1
                    incoming = (Kedge'*u) .* Bwo{k};
                else
                    incoming = (Kedge'*(w(Pmid(k-1)).*Fwo{k-1})) .* Bwo{k};
                end

                % outgoing at node time
                if k == L
                    outgoing = Fwo{k} .* (Kedge*v);
                else
                    outgoing = Fwo{k} .* (Kedge*(w(Pmid(k+1)).*Bwo{k+1}));
                end

                incoming = max(incoming,eps0); 
                outgoing = max(outgoing,eps0);

                A = node_A_wo(vtx);  A(:,pid) = A(:,pid) + Fwo{k}.*Bwo{k}; node_A_wo(vtx) = A;
                I = node_in_wo(vtx); I(:,pid) = I(:,pid) + incoming;       node_in_wo(vtx) = I;
                O = node_out_wo(vtx);O(:,pid) = O(:,pid) + outgoing;       node_out_wo(vtx) = O;
            end
        end
    end

    % Aggregates
    sumA0 = sum(A0_paths,2); 
    sumAT = sum(AT_paths,2);
    sumA0 = max(sumA0,eps0);   
    sumAT = max(sumAT,eps0);
    m0    = u.*sumA0; 
    mT    = v.*sumAT;

    % Build node loads, conservation & cap diagnostics
    cons_res = 0; 
    cap_res  = 0; 
    m_node   = containers.Map('KeyType','double','ValueType','any');
    for vtx = middles
        A_sum = sum(node_A_wo(vtx),2); 
        A_sum = max(A_sum,eps0);
        mv    = w(vtx).*A_sum; 
        m_node(vtx) = mv;

        IN  = sum(node_in_wo(vtx),2); 
        OUT = sum(node_out_wo(vtx),2);

        cons_res = cons_res + sum(abs(IN-OUT));
        cap_res  = cap_res  + sum(max(mv - r_cap(vtx),0));
    end

    % Multiplicative projections
    u = u .* ((p./max(m0,eps0)).^damp); 
    u = min(max(u,clip_min),clip_max);

    v = v .* ((q./max(mT,eps0)).^damp); 
    v = min(max(v,clip_min),clip_max);

    for vtx = middles
        mv = m_node(vtx); 
        rv = r_cap(vtx);
        w(vtx) = w(vtx).*(min(1, rv./max(mv,eps0)).^damp);
        w(vtx) = min(max(w(vtx),clip_min),clip_max);
    end

    % Logs
    histE0(end+1)   = sum(abs(m0-p));
    histET(end+1)   = sum(abs(mT-q));
    histCONS(end+1) = cons_res;
    histCAP(end+1)  = cap_res;

    if histE0(end)<2e-4 && histET(end)<2e-4 && cons_res<2e-5 && cap_res<2e-5
        break; 
    end
end

% Boundary polish (fix p,q exactly; w frozen)
for k = 1:12
    A0_paths(:)=0; 
    AT_paths(:)=0;
    for pid = 1:numel(paths)
        Pth  = paths{pid}; 
        Pmid = Pth(2:end-1); 
        L    = numel(Pmid);
        if L >= 1
            Bwo    = cell(L,1);  
            Bwo{L} = Kedge * v;  
            Bwo{L} = max(Bwo{L},eps0);
            for i = L-1:-1:1
                Bwo{i} = Kedge * (w(Pmid(i+1)).*Bwo{i+1}); 
                Bwo{i} = max(Bwo{i},eps0);
            end
            Fwo    = cell(L,1);  
            Fwo{1} = Kedge' * u; 
            Fwo{1} = max(Fwo{1},eps0);
            for i = 2:L
                Fwo{i} = Kedge' * (w(Pmid(i-1)).*Fwo{i-1}); 
                Fwo{i} = max(Fwo{i},eps0);
            end
            A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
            AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
        else
            A0_paths(:,pid) = Kedge * v;  
            AT_paths(:,pid) = Kedge' * u;
        end
    end
    sumA0 = sum(A0_paths,2); 
    sumA0 = max(sumA0,eps0); 
    u     = u .* (p ./ (u.*sumA0));

    % refresh AT with new u
    for pid = 1:numel(paths)
        Pth  = paths{pid}; 
        Pmid = Pth(2:end-1); 
        L    = numel(Pmid);
        if L >= 1
            F1 = Kedge'*u; 
            for i = 2:L
                F1 = Kedge' * (w(Pmid(i-1)).*F1);
            end
            AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*F1);
        else
            AT_paths(:,pid) = Kedge' * u;
        end
    end
    sumAT = sum(AT_paths,2); 
    sumAT = max(sumAT,eps0); 
    v     = v .* (q ./ (v.*sumAT));
end

% Final m0, mT
sumA0 = sum(A0_paths,2); 
sumAT = sum(AT_paths,2);
sumA0 = max(sumA0,eps0);   
sumAT = max(sumAT,eps0);
m0    = u.*sumA0;  
mT    = v.*sumAT;

%% -------- Build node ribbons and normalize on one global scale -----------
ribbons = containers.Map('KeyType','double','ValueType','any');
all_vals = [];
epsn     = 1e-16;

for vtx = middles
    A  = node_A_wo(vtx); 
    Wv = w(vtx);
    mv = sum(A .* (Wv.*ones(1,numel(paths))), 2);
    ribbons(vtx) = mv; 
    all_vals     = [all_vals; mv]; %#ok<AGROW>
end

denG     = max([cap_scalar; m0; mT; all_vals]) + epsn;
rv_blue  = m0/denG;  
rv_red   = mT/denG;  
capn_node= cap_scalar/denG;

%% -------- 2D layout: 8 nodes in a row, bottom arcs, node dots ------------
seq_nodes = [ v0, node_id(2,1), node_id(1,2), J1, J2, node_id(4,2), node_id(3,3), vT ];
labels    = { '$v_0$', '$v_1$', '$v_2$', '$v_3$', '$v_4$', '$v_5$', '$v_6$', '$v_{\mathcal T}$' };

xAll = linspace(0, 1, numel(seq_nodes));
xmap = containers.Map(num2cell(seq_nodes), num2cell(xAll));

% colors
cBlue  = [0    0.4470 0.7410]; 
cRed   = [0.6350 0.0780 0.1840]; 
cGreen = [0.20 0.62   0.20]; 
cEdge  = [0.20 0.20   0.20];

gain  = 0.14;     % single global width for all lobes/slices
arc_h = 0.12;     % peak height of bottom arcs (curvature)

figure(1); clf; set(gcf,'Color','w'); hold on; box on; grid off;

% vertical rails
for i = 1:numel(seq_nodes)
    x = xAll(i);
    plot([x x],[0 1],'k','LineWidth',1.6);
end

% node dots at t=0
for i = 1:numel(seq_nodes)
    x = xAll(i);
    plot(x,0,'ko','MarkerFaceColor',[0.15 0.15 0.15],'MarkerSize',8);
end

% source lobe to LEFT (blue), using gain
x0 = xmap(v0);
draw_density_lobe_left2D(x0, rv_blue, t, gain, cBlue, 0.35);

% interior slices to RIGHT + dashed common cap (green), using SAME gain
for i = 2:(numel(seq_nodes)-1)
    vtx = seq_nodes(i);
    x   = xAll(i);
    if isKey(ribbons, vtx)
        mv = ribbons(vtx)/denG;
    else
        mv = zeros(nt,1);
    end
    draw_slice_right2D(x, mv, t, gain, cGreen, 0.30);
    plot(x + gain*capn_node*ones(size(t)), t, 'k--','LineWidth',1.4);
end

% sink lobe to RIGHT (red), using gain
xT = xmap(vT);
draw_density_lobe_right2D(xT, rv_red, t, gain, cRed, 0.35);

% bottom arcs for each path edge (anchored at node dots, curved downward)
for pid = 1:numel(paths)
    Pth = paths{pid};
    for k = 1:(numel(Pth)-1)
        uN = Pth(k); 
        vN = Pth(k+1);
        xL = xmap(uN); 
        xR = xmap(vN);
        draw_bottom_arc2D(xL, xR, arc_h, cEdge, 1.2);
    end
end

% labels under each node, shifted further down to avoid overlap
for i = 1:numel(seq_nodes)
    x   = xAll(i);
    col = [0 0 0];
    if i == 1
        col = cBlue;
    end
    if i == numel(seq_nodes)
        col = cRed;
    end
    text(x, -0.04, labels{i}, 'Interpreter','latex','Color',col, ...
        'FontSize',24, 'HorizontalAlignment','center','VerticalAlignment','top');
end

xlim([-0.1 1.15]); 
ylim([-0.15 1.02]);
set(gca,'XTick',[],'FontSize',16); 
ylabel('$t$','Interpreter','latex','FontSize',20);
grid on; box on; hold off;

%% ===================== Helpers (2D) ======================================
function draw_density_lobe_left2D(x0, width_vec, t, gain, col, alphaFace)
x_edge = x0 - gain*width_vec(:);
x_fill = [x0*ones(numel(t),1); x_edge(end:-1:1)];
z_fill = [t;                    t(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(x_edge, t, 'Color', col, 'LineWidth', 2.6);
end

function draw_density_lobe_right2D(x0, width_vec, t, gain, col, alphaFace)
x_edge = x0 + gain*width_vec(:);
x_fill = [x0*ones(numel(t),1); x_edge(end:-1:1)];
z_fill = [t;                    t(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(x_edge, t, 'Color', col, 'LineWidth', 2.6);
end

function draw_slice_right2D(xc, mv_norm, t, gain, col, alphaFace)
xR = xc + gain*mv_norm(:);
x_fill = [xc*ones(numel(t),1); xR(end:-1:1)];
z_fill = [t;                   t(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(xR, t, 'Color', col, 'LineWidth', 2.2);
end

function draw_bottom_arc2D(xL, xR, h, col, lw)
    % Quadratic Bezier between (xL,0) and (xR,0) with control point below
    % the axis, so the arc curves downward.
    xc = 0.5*(xL+xR);  
    tc = -h;                 % flip sign to make it "smile" downward
    s  = linspace(0,1,80);
    bx = (1-s).^2*xL + 2*(1-s).*s*xc + s.^2*xR;
    bz = (1-s).^2*0  + 2*(1-s).*s*tc + s.^2*0;
    plot(bx, bz, 'Color', col, 'LineWidth', lw);
end









% %% 2D one-line chain with bottom arcs and one-sided distributions (fixed)
% clear; clc; close all;
% 
% %% -------- Grid & admissible routes (for connectivity) --------------------
% nx = 4; ny = 3;
% node_id = @(ix,iy) (iy-1)*nx + ix;
% 
% v0  = node_id(1,1);
% vT  = node_id(4,3);
% J1  = node_id(2,2);
% J2  = node_id(3,2);
% 
% % Paths exactly as specified
% pathA = [ node_id(1,1), node_id(2,1), J1, J2, node_id(4,2), vT ];
% pathB = [ node_id(1,1), node_id(1,2), J1, J2, node_id(3,3), vT ];
% pathC = [ node_id(1,1), node_id(2,1), J1, J2, node_id(3,3), vT ];
% 
% paths = {pathA, pathB, pathC};
% nodes_on_paths = unique([pathA(:); pathB(:); pathC(:)])';
% middles = setdiff(nodes_on_paths, [v0 vT]);
% 
% %% -------- Time grid & boundary marginals ---------------------------------
% nt = 200;  t = linspace(0,1,nt)';
% Gaussian = @(z,m,s) exp(-0.5*((z-m)./s).^2);
% p = 0.9*Gaussian(t,0.25,0.10) + 0.6*Gaussian(t,0.48,0.07); p = p/sum(p);
% q = 0.8*Gaussian(t,0.75,0.09) + 0.7*Gaussian(t,0.85,0.06); q = q/sum(q);
% 
% %% -------- Nearly-causal entropic kernel (row-normalized) -----------------
% eps_ent   = 0.005;
% tau_gap   = 0.03;
% alpha_cost= 0.002;
% 
% [U,V] = ndgrid(t,t);
% mask  = (V >= U);                           % allow diagonal
% Cedge = inf(nt);  
% Cedge(mask) = 1./((V(mask)-U(mask)) + tau_gap);
% Kedge = zeros(nt); 
% Kedge(mask) = exp(-alpha_cost * Cedge(mask) / eps_ent);
% rowSum = sum(Kedge,2); rowSum(rowSum==0)=1;  
% Kedge = Kedge ./ rowSum;                    % stochastic rows
% 
% %% -------- Uniform capacity at every interior node ------------------------
% cap_scalar = 0.007;
% r_cap = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles, r_cap(vtx) = cap_scalar*ones(nt,1); end
% 
% %% -------- Path-wise Sinkhorn with shared nodal multipliers ----------------
% eps0=1e-8; damp=0.4; clip_min=1e-8; clip_max=1e8;
% u = ones(nt,1); v = ones(nt,1);
% w = containers.Map('KeyType','double','ValueType','any');
% for vtx = middles, w(vtx) = ones(nt,1); end
% 
% A0_paths = zeros(nt, numel(paths));
% AT_paths = zeros(nt, numel(paths));
% node_A_wo  = containers.Map('KeyType','double','ValueType','any');
% node_in_wo = containers.Map('KeyType','double','ValueType','any');
% node_out_wo= containers.Map('KeyType','double','ValueType','any');
% for vtx = middles
%     node_A_wo(vtx)  = zeros(nt, numel(paths));
%     node_in_wo(vtx) = zeros(nt, numel(paths));
%     node_out_wo(vtx)= zeros(nt, numel(paths));
% end
% 
% histE0=[]; histET=[]; histCONS=[]; histCAP=[];
% maxIter=4000;
% 
% for it=1:maxIter
%     A0_paths(:)=0; AT_paths(:)=0;
%     for vtx = middles
%         node_A_wo(vtx)  = 0*node_A_wo(vtx);
%         node_in_wo(vtx) = 0*node_in_wo(vtx);
%         node_out_wo(vtx)= 0*node_out_wo(vtx);
%     end
% 
%     for pid=1:numel(paths)
%         Pth = paths{pid}; 
%         Pmid = Pth(2:end-1); 
%         L = numel(Pmid);
% 
%         if L>=1
%             % Backward (without local w)
%             Bwo=cell(L,1);  
%             Bwo{L} = Kedge * v;  
%             Bwo{L}=max(Bwo{L},eps0);
%             for k=L-1:-1:1
%                 Bwo{k} = Kedge * (w(Pmid(k+1)).*Bwo{k+1});
%                 Bwo{k} = max(Bwo{k},eps0);
%             end
%             % Forward (without local w)
%             Fwo=cell(L,1);  
%             Fwo{1} = Kedge' * u; 
%             Fwo{1}=max(Fwo{1},eps0);
%             for k=2:L
%                 Fwo{k} = Kedge' * (w(Pmid(k-1)).*Fwo{k-1});
%                 Fwo{k} = max(Fwo{k},eps0);
%             end
%         end
% 
%         % Source-side and sink-side messages
%         if L>=1
%             A0_paths(:,pid) = Kedge  * (w(Pmid(1)).*Bwo{1});
%             AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
%         else
%             A0_paths(:,pid) = Kedge * v;
%             AT_paths(:,pid) = Kedge' * u;
%         end
% 
%         % Node crossing contributions (wo local w_v), with SAFE indexing
%         if L>=1
%             for k=1:L
%                 vtx=Pmid(k);
% 
%                 % incoming at node time
%                 if k==1
%                     incoming = (Kedge'*u) .* Bwo{k};
%                 else
%                     incoming = (Kedge'*(w(Pmid(k-1)).*Fwo{k-1})) .* Bwo{k};
%                 end
% 
%                 % outgoing at node time
%                 if k==L
%                     outgoing = Fwo{k} .* (Kedge*v);
%                 else
%                     outgoing = Fwo{k} .* (Kedge*(w(Pmid(k+1)).*Bwo{k+1}));
%                 end
% 
%                 incoming=max(incoming,eps0); 
%                 outgoing=max(outgoing,eps0);
% 
%                 A = node_A_wo(vtx);  A(:,pid)=A(:,pid)+Fwo{k}.*Bwo{k}; node_A_wo(vtx)=A;
%                 I = node_in_wo(vtx); I(:,pid)=I(:,pid)+incoming;       node_in_wo(vtx)=I;
%                 O = node_out_wo(vtx);O(:,pid)=O(:,pid)+outgoing;       node_out_wo(vtx)=O;
%             end
%         end
%     end
% 
%     % Aggregates
%     sumA0 = sum(A0_paths,2); sumAT = sum(AT_paths,2);
%     sumA0=max(sumA0,eps0);   sumAT=max(sumAT,eps0);
%     m0 = u.*sumA0; mT = v.*sumAT;
% 
%     % Build node loads, conservation & cap diagnostics
%     cons_res=0; cap_res=0; 
%     m_node = containers.Map('KeyType','double','ValueType','any');
%     for vtx = middles
%         A_sum = sum(node_A_wo(vtx),2); A_sum=max(A_sum,eps0);
%         mv = w(vtx).*A_sum; m_node(vtx)=mv;
%         IN = sum(node_in_wo(vtx),2); 
%         OUT= sum(node_out_wo(vtx),2);
%         cons_res = cons_res + sum(abs(IN-OUT));
%         cap_res  = cap_res  + sum(max(mv - r_cap(vtx),0));
%     end
% 
%     % Multiplicative projections
%     u = u .* ((p./max(m0,eps0)).^damp); u = min(max(u,clip_min),clip_max);
%     v = v .* ((q./max(mT,eps0)).^damp); v = min(max(v,clip_min),clip_max);
%     for vtx = middles
%         mv=m_node(vtx); rv=r_cap(vtx);
%         w(vtx) = w(vtx).*(min(1, rv./max(mv,eps0)).^damp);
%         w(vtx) = min(max(w(vtx),clip_min),clip_max);
%     end
% 
%     % Logs
%     histE0(end+1)=sum(abs(m0-p));
%     histET(end+1)=sum(abs(mT-q));
%     histCONS(end+1)=cons_res;
%     histCAP(end+1)=cap_res;
% 
%     if histE0(end)<2e-4 && histET(end)<2e-4 && cons_res<2e-5 && cap_res<2e-5
%         break; 
%     end
% end
% 
% % Boundary polish (fix p,q exactly; w frozen)
% for k=1:12
%     A0_paths(:)=0; AT_paths(:)=0;
%     for pid=1:numel(paths)
%         Pth=paths{pid}; Pmid=Pth(2:end-1); L=numel(Pmid);
%         if L>=1
%             Bwo=cell(L,1);  Bwo{L} = Kedge * v;  Bwo{L}=max(Bwo{L},eps0);
%             for i=L-1:-1:1, Bwo{i} = Kedge * (w(Pmid(i+1)).*Bwo{i+1}); Bwo{i}=max(Bwo{i},eps0); end
%             Fwo=cell(L,1);  Fwo{1} = Kedge' * u; Fwo{1}=max(Fwo{1},eps0);
%             for i=2:L,      Fwo{i} = Kedge' * (w(Pmid(i-1)).*Fwo{i-1}); Fwo{i}=max(Fwo{i},eps0); end
%             A0_paths(:,pid) = Kedge * (w(Pmid(1)).*Bwo{1});
%             AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*Fwo{L});
%         else
%             A0_paths(:,pid) = Kedge * v;  
%             AT_paths(:,pid) = Kedge' * u;
%         end
%     end
%     sumA0 = sum(A0_paths,2); sumA0=max(sumA0,eps0); u = u .* (p ./ (u.*sumA0));
% 
%     % refresh AT with new u
%     for pid=1:numel(paths)
%         Pth=paths{pid}; Pmid=Pth(2:end-1); L=numel(Pmid);
%         if L>=1
%             F1 = Kedge'*u; 
%             for i=2:L, F1 = Kedge' * (w(Pmid(i-1)).*F1); end
%             AT_paths(:,pid) = Kedge' * (w(Pmid(L)).*F1);
%         else
%             AT_paths(:,pid) = Kedge' * u;
%         end
%     end
%     sumAT = sum(AT_paths,2); sumAT=max(sumAT,eps0); 
%     v = v .* (q ./ (v.*sumAT));
% end
% 
% % Final m0, mT
% sumA0 = sum(A0_paths,2); sumAT = sum(AT_paths,2);
% sumA0=max(sumA0,eps0);   sumAT=max(sumAT,eps0);
% m0 = u.*sumA0;  mT = v.*sumAT;
% 
% %% -------- Build node ribbons and normalize on one global scale -----------
% ribbons = containers.Map('KeyType','double','ValueType','any');
% all_vals=[]; epsn=1e-16;
% for vtx = middles
%     A = node_A_wo(vtx); Wv = w(vtx);
%     mv = sum(A .* (Wv.*ones(1,numel(paths))), 2);
%     ribbons(vtx)=mv; all_vals=[all_vals; mv]; %#ok<AGROW>
% end
% denG  = max([cap_scalar; m0; mT; all_vals]) + epsn;
% rv_blue = m0/denG;   rv_red = mT/denG;   capn_node = cap_scalar/denG;
% 
% %% -------- 2D layout: 8 nodes in a row, bottom arcs, node dots ------------
% seq_nodes = [ v0, node_id(2,1), node_id(1,2), J1, J2, node_id(4,2), node_id(3,3), vT ];
% labels    = { '$v_0$', '$v_1$', '$v_2$', '$v_3$', '$v_4$', '$v_5$', '$v_6$', '$v_{\mathcal T}$' };
% 
% xAll = linspace(0, 1, numel(seq_nodes));
% xmap = containers.Map(num2cell(seq_nodes), num2cell(xAll));
% 
% % colors
% cBlue=[0 0.4470 0.7410]; cRed=[0.6350 0.0780 0.1840]; cGreen=[0.20 0.62 0.20]; cEdge=[0.20 0.20 0.20];
% gain = 0.18;     % global width for lobes/slices
% arc_h = 0.12;    % peak height of bottom arcs
% 
% figure(1); clf; set(gcf,'Color','w'); hold on; box on; grid off;
% 
% % vertical rails
% for i = 1:numel(seq_nodes)
%     x = xAll(i);
%     plot([x x],[0 1],'k','LineWidth',1.6);
% end
% 
% % node dots at t=0
% for i = 1:numel(seq_nodes)
%     x = xAll(i);
%     plot(x,0,'ko','MarkerFaceColor',[0.15 0.15 0.15],'MarkerSize',8);
% end
% 
% % source lobe to LEFT
% x0 = xmap(v0);
% draw_density_lobe_left2D(x0, rv_blue, t, gain, cBlue, 0.35);
% 
% % interior slices to RIGHT + dashed common cap
% for i = 2:(numel(seq_nodes)-1)
%     vtx = seq_nodes(i);
%     x = xAll(i);
%     if isKey(ribbons, vtx), mv = ribbons(vtx)/denG; else, mv = zeros(nt,1); end
%     draw_slice_right2D(x, mv, t, gain*0.55, cGreen, 0.30);
%     plot(x + gain*0.55*capn_node*ones(size(t)), t, 'k--','LineWidth',1.4);
% end
% 
% % sink lobe to RIGHT
% xT = xmap(vT);
% draw_density_lobe_right2D(xT, rv_red, t, gain, cRed, 0.35);
% 
% % bottom arcs for each path edge (anchored at node dots)
% for pid = 1:numel(paths)
%     Pth = paths{pid};
%     for k = 1:(numel(Pth)-1)
%         uN = Pth(k); vN = Pth(k+1);
%         xL = xmap(uN); xR = xmap(vN);
%         draw_bottom_arc2D(xL, xR, arc_h, cEdge, 1.2);
%     end
% end
% 
% % labels under each node
% for i = 1:numel(seq_nodes)
%     x = xAll(i);
%     col = [0 0 0];
%     if i==1, col = cBlue; end
%     if i==numel(seq_nodes), col = cRed; end
%     text(x, -0.055, labels{i}, 'Interpreter','latex','Color',col, ...
%         'FontSize',24, 'HorizontalAlignment','center','VerticalAlignment','top');
% end
% 
% xlim([-0.07 1.07]); ylim([-0.10 1.02]);
% set(gca,'XTick',[]); 
% ylabel('$t$','Interpreter','latex','FontSize',24);
% grid on; box on; hold off;
% 
% %% ===================== Helpers (2D) ======================================
% function draw_density_lobe_left2D(x0, width_vec, t, gain, col, alphaFace)
% x_edge = x0 - gain*width_vec(:);
% x_fill = [x0*ones(numel(t),1); x_edge(end:-1:1)];
% z_fill = [t;                    t(end:-1:1)];
% fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
% plot(x_edge, t, 'Color', col, 'LineWidth', 2.6);
% end
% 
% function draw_density_lobe_right2D(x0, width_vec, t, gain, col, alphaFace)
% x_edge = x0 + gain*width_vec(:);
% x_fill = [x0*ones(numel(t),1); x_edge(end:-1:1)];
% z_fill = [t;                    t(end:-1:1)];
% fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
% plot(x_edge, t, 'Color', col, 'LineWidth', 2.6);
% end
% 
% function draw_slice_right2D(xc, mv_norm, t, gain, col, alphaFace)
% xR = xc + gain*mv_norm(:);
% x_fill = [xc*ones(numel(t),1); xR(end:-1:1)];
% z_fill = [t;                   t(end:-1:1)];
% fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
% plot(xR, t, 'Color', col, 'LineWidth', 2.2);
% end
% 
% 
% function draw_bottom_arc2D(xL, xR, h, col, lw)
%     xc = 0.5*(xL+xR);  
%     tc = -h;                  % <--- flip the sign to reverse the arc
%     s  = linspace(0,1,80);
%     bx = (1-s).^2*xL + 2*(1-s).*s*xc + s.^2*xR;
%     bz = (1-s).^2*0  + 2*(1-s).*s*tc + s.^2*0;
%     plot(bx, bz, 'Color', col, 'LineWidth', lw);
% end
% 






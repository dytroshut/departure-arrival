%% Multi-time Dirac transport with 5 interior times, UNIFORM capacity, and plots
% Requires CVX.
clear; clc; close all;

%% ---------------------- Discrete supports -------------------------------
n1 = 80;                 % refined # points for t0
n2 = 80;                 % refined # points for t_T
t0 = (0.01:n1-0.01)'/n1;  % t0 in (0,1)
tf = (0.01:n2-0.01)'/n2;  % t_T in (0,1)

normalize = @(a) a/sum(a(:));
Gaussian  = @(x,m,s) exp(-(x-m).^2/(2*s^2));

sigma0 = 0.01; sigma1 = 0.02; sigma2 = 0.03; sigma3 = 0.04; sigma4 = 0.05;

p = 5.0*Gaussian(t0, .2, sigma3) + 1.5*Gaussian(t0, .3, sigma0) + 2.5*Gaussian(t0, .4, sigma2);
q = 3.0*Gaussian(tf, .8, sigma4) + 1.0*Gaussian(tf, .6, sigma2) + 2.0*Gaussian(tf, .9, sigma2);
p = normalize(p);  q = normalize(q);

%% ---------------------- Time grid (shared for all middles) --------------
nt = 100; 
T  = 1;
t  = T*(0+0.01:nt-0.01)'/nt;

%% ---------------------- Leg costs (additive) ----------------------------
BIG = 1e12;

% Leg 0: C0(i,k1) = 1/(t(k1)-t0(i))
[T0, T1] = ndgrid(t0, t);
C0 = 1./(T1 - T0); C0(~isfinite(C0) | C0<=0) = BIG;

% Leg 1..4: C(km,km+1) = 1/(t(km+1)-t(km))
[Tk1, Tk2] = ndgrid(t, t);
C1 = 1./(Tk2 - Tk1); C1(~isfinite(C1) | C1<=0) = BIG;
C2 = C1; C3 = C1; C4 = C1;

% Leg 5: C5(k5,j) = 1/(t_T(j)-t(k5))
[T5, TF] = ndgrid(t, tf);
C5 = 1./(TF - T5); C5(~isfinite(C5) | C5<=0) = BIG;

%% ---------------------- UNIFORM capacity per time & per rail ------------
% Same logic as the 1-rail Dirac example:
% capacity = cap * T / nt, with cap = 2, T = 1.
total_mass = 0.5*(sum(p) + sum(q));       % = 1 with our normalization
cap_frac   = 2;                           % "2" from your earlier example
cap0       = cap_frac * (total_mass/nt);  % per-time-slice bound (0.02)

% One uniform capacity profile, same for all 5 rails
cap_uni = cap0 * ones(nt,1);

cap1 = cap_uni;
cap2 = cap_uni;
cap3 = cap_uni;
cap4 = cap_uni;
cap5 = cap_uni;

%% ---------------------- CVX variables & problem -------------------------
cvx_begin
    % cvx_solver mosek
    cvx_precision best
    variables P0(n1,nt) P1(nt,nt) P2(nt,nt) P3(nt,nt) P4(nt,nt) P5(nt,n2)
    variables r1(nt) r2(nt) r3(nt) r4(nt) r5(nt)

    % Objective: sum of leg costs
    minimize( sum(sum(P0.*C0)) + sum(sum(P1.*C1)) + sum(sum(P2.*C2)) ...
            + sum(sum(P3.*C3)) + sum(sum(P4.*C4)) + sum(sum(P5.*C5)) )

    subject to
        % Marginal consistency along the chain
        sum(P0,2) == p;          sum(P0,1)' == r1;
        sum(P1,2) == r1;         sum(P1,1)' == r2;
        sum(P2,2) == r2;         sum(P2,1)' == r3;
        sum(P3,2) == r3;         sum(P3,1)' == r4;
        sum(P4,2) == r4;         sum(P4,1)' == r5;
        sum(P5,2) == r5;         sum(P5,1)' == q;

        % Nonnegativity
        P0 >= 0; P1 >= 0; P2 >= 0; P3 >= 0; P4 >= 0; P5 >= 0;
        r1 >= 0; r2 >= 0; r3 >= 0; r4 >= 0; r5 >= 0;

        % Uniform capacity bound on each middle rail
        r1 <= cap1;
        r2 <= cap2;
        r3 <= cap3;
        r4 <= cap4;
        r5 <= cap5;
cvx_end

%% ---------------------- Diagnostics -------------------------------------
fprintf('sum p_hat: %.6f | sum p: %.6f\n', sum(sum(P0,2)), sum(p));
fprintf('sum q_hat: %.6f | sum q: %.6f\n', sum(sum(P5,1)), sum(q));
viol = sum(max(0,r1-cap1))+sum(max(0,r2-cap2))+sum(max(0,r3-cap3))+sum(max(0,r4-cap4))+sum(max(0,r5-cap5));
fprintf('total capacity violation (should be 0): %.3e\n', viol);
fprintf('per-time capacity cap0 = %.4f\n', cap0);

%% ---------------------- Effective p–q coupling via chain ----------------
eps0 = 1e-12;
K0 = bsxfun(@rdivide, P0, p + eps0);   % n1 x nt
K1 = bsxfun(@rdivide, P1, r1 + eps0);  % nt x nt
K2 = bsxfun(@rdivide, P2, r2 + eps0);  % nt x nt
K3 = bsxfun(@rdivide, P3, r3 + eps0);  % nt x nt
K4 = bsxfun(@rdivide, P4, r4 + eps0);  % nt x nt
K5 = bsxfun(@rdivide, P5, r5 + eps0);  % nt x n2

Mchain     = K0 * K1 * K2 * K3 * K4 * K5;    % n1 x n2
pqCoupling = bsxfun(@times, Mchain, p);      % diag(p)*M
rk_sum     = r1 + r2 + r3 + r4 + r5;         % aggregate time marginal (not used in plot here)

%% ---------------------- Figure 1: 2D DA-style v0 · v1..v5 · v_T ---------
% Horizontal positions
x_src = 0.0;
x_snk = 1.0;
x_mid = linspace(0.15,0.85,5);   % five middle rails v_1,...,v_5

% Colors (match DA + capacity style)
cBlue   = [0 0.4470 0.7410];
cRed    = [0.6350 0.0780 0.1840];
greens  = [0.20 0.62 0.20;
           0.26 0.66 0.26;
           0.33 0.70 0.33;
           0.40 0.74 0.40;
           0.48 0.78 0.48];
cEdge   = [0.20 0.20 0.20];

% Gather rails and capacities for normalization
R_all   = [r1 r2 r3 r4 r5];                     % nt x 5
Cap_all = [cap1 cap2 cap3 cap4 cap5];           % nt x 5 (all identical now)
epsn    = 1e-16;
denG    = max([p; q; R_all(:); Cap_all(:)]) + epsn;

% Normalized profiles (same global scale)
rv_blue = p   / denG;         % n1 x 1
rv_red  = q   / denG;         % n2 x 1
rv_mid  = R_all   / denG;     % nt x 5
cap_mid = Cap_all / denG;     % nt x 5 (each column constant)

gain  = 0.18;   % global width for lobes/slices
arc_h = 0.10;   % bottom arc curvature

figure(1); clf; set(gcf,'Color','w'); hold on; box on; grid on;

% All node x-positions and labels
x_all  = [x_src, x_mid, x_snk];
labels = [{'$v_0$'}, ...
          arrayfun(@(k) sprintf('$v_{%d}$',k), 1:5, 'UniformOutput', false), ...
          {'$v_{\mathcal T}$'}];

% Vertical rails
for i = 1:numel(x_all)
    x = x_all(i);
    plot([x x],[0 1],'k','LineWidth',1.6);
end

% Node dots at t=0
for i = 1:numel(x_all)
    x = x_all(i);
    plot(x,0,'ko','MarkerFaceColor',[0.15 0.15 0.15],'MarkerSize',8);
end

% Source lobe at v0, using t0
draw_density_lobe_left2D(x_src, rv_blue, t0, gain, cBlue, 0.35);

% Middle rails v1..v5: each a right lobe with its uniform cap curve
for m = 1:5
    xm  = x_mid(m);
    col = greens(m,:);
    draw_slice_right2D(xm, rv_mid(:,m), t, gain, col, 0.30);
    % capacity: constant-per-time, so this is a straight vertical ribbon to the right
    plot(xm + gain*cap_mid(:,m), t, 'Color', 0.35*col, ...
         'LineWidth', 1.6, 'LineStyle', '--');
end

% Sink lobe at v_T, using tf
draw_density_lobe_right2D(x_snk, rv_red, tf, gain, cRed, 0.35);

% Optional bottom arcs between consecutive rails
for i = 1:(numel(x_all)-1)
    draw_bottom_arc2D(x_all(i), x_all(i+1), arc_h, cEdge, 1.1);
end

% Labels under each node, below arcs
for i = 1:numel(x_all)
    x   = x_all(i);
    col = [0 0 0];
    if i == 1,                     col = cBlue; end
    if i == numel(x_all),          col = cRed;  end
    text(x, -0.04, labels{i}, 'Interpreter','latex','Color',col, ...
        'FontSize',18, 'HorizontalAlignment','center','VerticalAlignment','top');
end

xlim([-0.10 1.15]);
ylim([-0.15 1.02]);
set(gca,'XTick',[],'FontSize',16);
ylabel('$t$','Interpreter','latex','FontSize',20);
grid on; box on; hold off;

%% ===================== 2D helper shapes ============================
function draw_density_lobe_left2D(x0, width_vec, tgrid, gain, col, alphaFace)
x_edge = x0 - gain*width_vec(:);
x_fill = [x0*ones(numel(tgrid),1); x_edge(end:-1:1)];
z_fill = [tgrid;                  tgrid(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(x_edge, tgrid, 'Color', col, 'LineWidth', 2.6);
end

function draw_density_lobe_right2D(x0, width_vec, tgrid, gain, col, alphaFace)
x_edge = x0 + gain*width_vec(:);
x_fill = [x0*ones(numel(tgrid),1); x_edge(end:-1:1)];
z_fill = [tgrid;                  tgrid(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(x_edge, tgrid, 'Color', col, 'LineWidth', 2.6);
end

function draw_slice_right2D(xc, mv_norm, tgrid, gain, col, alphaFace)
xR = xc + gain*mv_norm(:);
x_fill = [xc*ones(numel(tgrid),1); xR(end:-1:1)];
z_fill = [tgrid;                  tgrid(end:-1:1)];
fill(x_fill, z_fill, col, 'FaceAlpha',alphaFace, 'EdgeColor','none');
plot(xR, tgrid, 'Color', col, 'LineWidth', 2.2);
end

function draw_bottom_arc2D(xL, xR, h, col, lw)
xc = 0.5*(xL + xR);
tc = -h;
s  = linspace(0,1,80);
bx = (1-s).^2*xL + 2*(1-s).*s*xc + s.^2*xR;
bz = (1-s).^2*0  + 2*(1-s).*s*tc + s.^2*0;
plot(bx, bz, 'Color', col, 'LineWidth', lw);
end
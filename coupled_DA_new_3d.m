%% Taxi allocation (1D) • Dirac-style rails + DA surface + continuous time curve
clear all; clc; close all;

%% ---------------------- Supports ----------------------------------------
n1 = 120; n2 = 120;
t0 = (0.01:0.35*n1-0.01)'/n1;          % source support (0 .. 0.35)
tf = (0.65*n2:n2-0.01)'/n2;            % sink   support (0.65 .. 1.0)
nt0 = length(t0); ntf = length(tf);

normalize = @(a)a/sum(a,'all');
Gaussian  = @(x,m,s) exp(-(x-m).^2/(2*s^2));

% Fixed source–destination coupling (DAconstraint)
w0 = Gaussian(t0, .10, .03);
w  = Gaussian(t0, .20, .03);
w1 = Gaussian(t0, .30, .03);
w2 = Gaussian(tf, .80, .02);
w3 = Gaussian(tf, .90, .04);

W = 3.*w*w2' + 2.*w1*w3' + 5.*w*w3' + w2*w3';
DAconstraint = normalize(W);              % nt0 x ntf

% Spatial marginals
p = sum(DAconstraint,2);          % nt0 x 1
q = sum(DAconstraint,1)';         % ntf x 1

%% ---------------------- Time grid & cost --------------------------------
nt = 60;
t1 = (0.35*nt:0.65*nt-0.01)'/nt;         % crossing-time window in (0,1)
nt = length(t1);

% Per–time-slice capacity bound (inequality)
wt = 2.*ones(nt,1);
wt = normalize(wt);

% Cost: split over two legs, 1/(gap)
[X1, X2, Y] = ndgrid(t0, tf, t1);
C = 1./(Y - X1) + 1./(X2 - Y);
Cost = permute(C, [2,1,3]);              % (nt0 x ntf x nt)

%% ---------------------- CVX solve ---------------------------------------
cvx_begin
    % cvx_solver mosek
    cvx_precision best
    variables Pi(nt0,ntf,nt) rk(nt)
    minimize( sum(sum(sum(Cost .* Pi))) )
    subject to
        squeeze(sum(sum(Pi))) <= wt;      % time capacity
        Pi >= 0;                          % nonnegativity
        sum(Pi,3) == DAconstraint;        % fixed DA coupling
cvx_end

% Time marginal actually used
rk = squeeze(sum(sum(Pi)));               % nt x 1

% Couplings for plots
pqCoupling = sum(Pi,3);                   % nt0 x ntf
ptCoupling = squeeze(sum(Pi,2));          % nt0 x nt
qtCoupling = squeeze(sum(Pi,1));          % ntf x nt

%% ---------------------- Figure 1: source · middle · sink rails ----------
x_src = 0; x_mid = 0.5; x_snk = 1.0;
xx1 = t0; yy1 = tf; tt1 = t1;
xx2 = x_src*ones(size(xx1));
yy2 = x_snk*ones(size(yy1));
tt2 = x_mid*ones(size(tt1));

cBlue  = [0 0.4470 0.7410];
cRed   = [0.6350 0.0780 0.1840];
cGreen = [0.4660 0.6740 0.1880];

figure(1); clf; set(gcf,'Color','w'); hold on; box on; grid on;

% Filled ribbons under each rail
fill3([xx2;xx2],[xx1;flipud(xx1)],[zeros(nt0,1);flipud(p)], cBlue,  'FaceAlpha',0.25,'EdgeColor','none');
fill3([yy2;yy2],[yy1;flipud(yy1)],[zeros(ntf,1);flipud(q)], cRed,   'FaceAlpha',0.25,'EdgeColor','none');
fill3([tt2;tt2],[tt1;flipud(tt1)],[zeros(nt,1);flipud(rk)], cGreen, 'FaceAlpha',0.25,'EdgeColor','none');

% Rail outlines
plot3(xx2,xx1,p,  'LineWidth',2.5,'Color',cBlue);
plot3(yy2,yy1,q,  'LineWidth',2.5,'Color',cRed);
plot3(tt2,tt1,rk, 'LineWidth',2.5,'Color',cGreen);

% Capacity line on middle rail
plot3(tt2,tt1,wt,'LineWidth',2.5,'Color',[0.25 0.25 0.25],'LineStyle',':');

% Thin grey guides (optional)
[Ipt,Jpt,~] = find(ptCoupling > 2e-4);
for k = 1:length(Ipt)
    plot([x_src,x_mid],[t0(Ipt(k)),t1(Jpt(k))],'Color',[0.9 0.9 0.9],'LineWidth',2);
end
[Iqt,Jqt,~] = find(qtCoupling > 2e-4);
for k = 1:length(Iqt)
    plot([x_snk,x_mid],[tf(Iqt(k)),t1(Jqt(k))],'Color',[0.9 0.9 0.9],'LineWidth',2);
end

set(gca,'XTick',[x_src x_mid x_snk],'XTickLabel',{'source','middle','sink'},'FontSize',13);
xlim([-0.05 1.05]);
ylabel('time','FontSize',14);
zlabel('density','FontSize',14);
view([-27 28]); hold off;

%% ---------------------- Figure 2: EXACT x,y,z marginals + fancy colormap
list = find (Pi > 1e-6);
[I,J,KK] = ind2sub(size(Pi), list);
[Xs, Ys] = meshgrid(t0, tf);                % (ntf x nt0)

figure(2); clf; set(gcf,'Color','w'); hold on; box on; grid on;
set(gcf,'Renderer','opengl');

% Gains
ridgeGain = 3;       % height for p,q ridges
surfGain  = 30;      % DA surface height (bigger)
timeGain  = 1;       % x-scale for green time curve

% Compact ridges (trim tiny tails)
th_p = 1e-4*max(p);  th_q = 1e-4*max(q);
idxp = p > th_p;     idxq = q > th_q;

% Blue ridge along y = min(tf)
plot3(t0(idxp), min(tf)*ones(sum(idxp),1), ridgeGain*p(idxp), ...
      'LineWidth',3,'Color',cBlue);

% Red ridge along x = min(t0)
plot3(min(t0)*ones(sum(idxq),1), tf(idxq), ridgeGain*q(idxq), ...
      'LineWidth',3,'Color',cRed);

% DA coupling surface (larger, no edges, no shadows)
Zsurf = surfGain .* pqCoupling';            % (ntf x nt0)
S = surf(Xs, Ys, Zsurf, 'FaceAlpha',0.95, 'EdgeColor','none');
shading interp; lighting none;              % no shadowing

% Fancy, no-pink colormap
cm = makeFancyMap(256, 0.80, 1.20);
colormap(cm);
vmax = prctile(Zsurf(:), 99);
caxis([0, vmax]);

% Continuous green time curve over z in [0,1], zero outside active window
z_full  = linspace(0,1,1000)';                   % full vertical support
x_green = zeros(size(z_full));                   % zero outside window
mask    = (z_full >= t1(1)) & (z_full <= t1(end));
x_green(mask) = timeGain * interp1(t1, rk, z_full(mask), 'linear');
plot3(x_green, ones(size(z_full)), z_full, 'LineWidth',3, 'Color', cGreen);

% Time-colored dominant triplets (kept)
vals = Pi(list);
[~, ord] = sort(vals, 'descend');
M = min(300, numel(ord)); ord = ord(1:M);
x_pts = t0(I(ord)); y_pts = tf(J(ord)); z_pts = t1(KK(ord));
tt_norm = (z_pts - min(t1)) / (max(t1) - min(t1) + eps);
cmapT   = parula(256); idxT = max(1, min(256, round(1 + tt_norm*(256-1))));
Cpts    = cmapT(idxT, :);
w_pts   = vals(ord);
smin = 24; smax = 120; vmin = w_pts(end); vmaxv = w_pts(1);
if vmaxv > vmin
    sz = smin + (smax - smin)*(w_pts - vmin)/(vmaxv - vmin);
else
    sz = smin*ones(size(w_pts));
end
scatter3(x_pts, y_pts, z_pts, sz, Cpts, 'filled', ...
         'MarkerEdgeColor','k', 'LineWidth',0.25);

% Axes
xlabel('$t_0$','Interpreter','latex','FontSize',20);
ylabel('$t_{\mathcal T}$','Interpreter','latex','FontSize',20);
zlabel('$t_1$','Interpreter','latex','FontSize',20);
set(gca,'FontSize',16);
view([47 32]); axis tight; axis vis3d;
hold off;

%% ---------------------- Colormap helper ---------------------------------
function cm = makeFancyMap(n, lightness, gamma)
% makeFancyMap  Teal→Lime→Gold→Orange, no pinks, light low end
if nargin<1 || isempty(n),         n = 256;    end
if nargin<2 || isempty(lightness), lightness = 0.80; end
if nargin<3 || isempty(gamma),     gamma = 1.20; end
anchors = [ 0.05 0.25 0.55
            0.10 0.60 0.55
            0.35 0.78 0.30
            0.95 0.85 0.25
            0.98 0.55 0.05 ];
m  = size(anchors,1);
xi = linspace(0,1,n)'; 
xp = linspace(0,1,m);
cm = [interp1(xp, anchors(:,1), xi, 'pchip'), ...
      interp1(xp, anchors(:,2), xi, 'pchip'), ...
      interp1(xp, anchors(:,3), xi, 'pchip')];
w = lightness * (1 - xi).^gamma;          % weight of white
cm = cm .* (1 - w) + 1 .* w;
cm(end,:) = anchors(end,:);
cm = max(0, min(1, cm));
end
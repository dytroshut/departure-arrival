function S = example2_data(which, nt)
%EXAMPLE2_DATA Data for the two independent-DA experiments in Example 2.
% All discrete profiles are CELL MASSES, not continuous-time densities.
% Each edge has weight 1; the cost is 1/(arrival time - departure time).
if nargin < 2, nt = 200; end
validateattributes(nt, {'numeric'}, {'scalar','integer','>=',10});
S.t = linspace(0,1,nt)'; S.nt = nt; S.epsilon = 0.065;
G = @(m,s) exp(-0.5*((S.t-m)/s).^2);
p = 0.9*G(.25,.10) + 0.6*G(.48,.07);
q = 0.8*G(.75,.09) + 0.7*G(.85,.06);
switch lower(char(which))
    case 'line'
        S.paths = 1:8;
        S.middle = 2:7; S.displayNodes = 1:8;
        % Approved correction: the old .0100 sinusoidal case is infeasible.
        % Apply the same corrected mean to BOTH line-graph cases.
        S.capacity = .0105;
    case 'shared'
        % Node IDs follow the original (iy-1)*4+ix convention.
        % Display labels: v0,v1,v2,v3,v4,v5,v6,vT.
        S.paths = [1 2 6 7 8 12; 1 5 6 7 11 12; 1 2 6 7 11 12];
        S.middle = [2 5 6 7 8 11];
        S.displayNodes = [1 2 5 6 7 8 11 12];
        S.capacity = .009;
    otherwise
        error('example2:UnknownCase','Choose line or shared.');
end
S.name = lower(char(which));
S.source = S.paths(1,1); S.sink = S.paths(1,end);
S.nNodes = max(S.paths(:)); S.nEdges = size(S.paths,2)-1;
S.edgeWeight = 1;
S.labels = {'$v_0$','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$', ...
    '$v_{\mathcal T}$'};

% Every strict grid crossing advances at least one cell. Remove only
% boundary cells that cannot belong to a full source--sink time tuple.
S.removedDepartureMass = sum(p(end-S.nEdges+1:end))/sum(p);
S.removedArrivalMass = sum(q(1:S.nEdges))/sum(q);
p(end-S.nEdges+1:end) = 0; q(1:S.nEdges) = 0;
S.p = p/sum(p); S.q = q/sum(q);

[U,V] = ndgrid(S.t,S.t);
S.support = V > U;
S.cost = zeros(nt);
S.cost(S.support) = S.edgeWeight./(V(S.support)-U(S.support));
% Zero outside support is used only for cost contractions; the solver
% excludes those entries by setting the log kernel to -Inf.
S.uniform = inf(nt,S.nNodes); S.varying = S.uniform;
for k = 1:numel(S.middle)
    v = S.middle(k);
    S.uniform(:,v) = S.capacity;
    modulation = 1 + .65*sin(2*pi*S.t/.60 + (k-1)*pi/2);
    S.varying(:,v) = S.capacity*nt*modulation/sum(modulation);
end
end

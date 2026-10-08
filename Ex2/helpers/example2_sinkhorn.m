function R = example2_sinkhorn(S, capacity, opts)
%EXAMPLE2_SINKHORN Exact cyclic path-wise scaling in the log domain.
% Independent DA, a common source and sink, equal-length staged paths.
% Within a stage each path visits at most one node, so those node blocks
% are separable. Later-stage backward messages remain valid while the
% forward prefix is refreshed after each capacity block.
% No message floors, row normalization, damping, or boundary polishing.
if nargin < 3, opts = struct; end
opts = defaultOption(opts,'tol',1e-8);
opts = defaultOption(opts,'complementarityTol',1e-8);
opts = defaultOption(opts,'maxIter',100000);
opts = defaultOption(opts,'checkEvery',10);
opts = defaultOption(opts,'verbose',true);
opts = defaultOption(opts,'initial',[]);
nt = S.nt; paths = S.paths; [np,nn] = size(paths);
assert(all(paths(:,1)==S.source) && all(paths(:,end)==S.sink), ...
    'example2:Endpoints','All paths must share the prescribed endpoints.');
assert(isequal(size(capacity),[nt,S.nNodes]), ...
    'example2:CapacitySize','Capacity must be nt-by-nNodes.');
assert(all(all(isfinite(capacity(:,S.middle)) & capacity(:,S.middle)>0)), ...
    'example2:Capacity','Interior capacities must be finite and positive.');
% This implementation uses a global stage order, not arbitrary path orders.
for v = S.middle
    [~,pos] = find(paths==v);
    assert(isscalar(unique(pos)),'example2:NonStagedPaths', ...
        'Each interior node must occur in the same stage on every path.');
end
logK = -inf(nt); logK(S.support) = -S.cost(S.support)/S.epsilon;
logp = log(S.p); logq = log(S.q); logr = log(capacity);
logu = zeros(nt,1); logv = zeros(nt,1); logW = zeros(nt,S.nNodes);
logu(S.p==0) = -inf; logv(S.q==0) = -inf;
if ~isempty(opts.initial)
    logu = opts.initial.logu; logv = opts.initial.logv;
    logW = min(opts.initial.logW,0);
end
activeP = S.p>0; activeQ = S.q>0;
if np==1
    previous=cumsum(S.p);
    for pos=2:nn-1
        available=[0;previous(1:end-1)]; earliest=zeros(nt,1);
        for j=1:nt
            if j==1, before=0; else, before=earliest(j-1); end
            earliest(j)=min(available(j),before+capacity(j,paths(1,pos)));
        end
        previous=earliest;
    end
    shortfall=max(cumsum(S.q)-[0;previous(1:end-1)]);
    if shortfall>1e-12
        error('example2:InfeasibleLine', ...
            'The strict-causal line problem is infeasible (arrival shortfall %.8g).',shortfall);
    end
end
history = nan(opts.maxIter,6); started = tic; converged = false;
for it = 1:opts.maxIter
    % Exact source block with the current sink and interior multipliers.
    B = backward(logK,logv,logW,paths);
    A0 = logsum(B{1},2);
    if any(~isfinite(A0(activeP)))
        error('example2:SourceSupport','Positive supply has no causal route.');
    end
    logu(activeP) = logp(activeP)-A0(activeP);

    % Exact sink block, using the newly updated source multiplier.
    F = forward(logK,logu,logW,paths);
    AT = logsum(F{nn},2);
    if any(~isfinite(AT(activeQ)))
        error('example2:SinkSupport','Positive demand has no causal route.');
    end
    logv(activeQ) = logq(activeQ)-AT(activeQ);

    % Capacity blocks: logW <- min(log(r)-log(A_without_W),0).
    % Equivalently W <- min(W.*(r./m),1), allowing W to INCREASE.
    B = backward(logK,logv,logW,paths);
    prefix = repmat(logu,1,np);
    for pos = 2:nn-1
        incoming = logmultiply(logK',prefix);
        nodes = unique(paths(:,pos));
        for v = nodes'
            incident = paths(:,pos)==v;
            logA = logsum(incoming(:,incident)+B{pos}(:,incident),2);
            logW(:,v) = min(logr(:,v)-logA,0);
        end
        prefix = incoming+logW(:,paths(:,pos)');
    end

    % Check only after a complete cycle. Skipping intermediate diagnostics
    % changes neither the updates nor the final constrained problem.
    if it~=1 && mod(it,opts.checkEvery)~=0 && it~=opts.maxIter, continue; end
    % Evaluate every residual from the SAME final set of cycle scalings.
    B = backward(logK,logv,logW,paths);
    F = forward(logK,logu,logW,paths);
    [m0,mT,mnode,pathMass] = marginals(F,B,logu,logv,logW,paths,S.nNodes);
    E0 = norm(m0-S.p,1); ET = norm(mT-S.q,1);
    V = sum(sum(max(mnode(:,S.middle)-capacity(:,S.middle),0)));
    dual = S.epsilon*(sum(logu(activeP).*S.p(activeP)) + ...
        sum(logv(activeQ).*S.q(activeQ)) + ...
        sum(sum(logW(:,S.middle).*capacity(:,S.middle))) - sum(pathMass));
    comp = S.epsilon*sum(sum(abs(logW(:,S.middle).* ...
        (capacity(:,S.middle)-mnode(:,S.middle)))));
    history(it,:) = [it,E0,ET,V,dual,comp];
    if opts.verbose && (it==1 || mod(it,500)==0)
        fprintf('  cycle %d: boundary %.3e, capacity %.3e, complementarity %.3e\n', ...
            it,max(E0,ET),V,comp);
    end
    if max([E0,ET,V])<=opts.tol && ...
            comp<=opts.complementarityTol*max(1,abs(dual))
        converged = true; break;
    end
end
if ~converged
    error('example2:NotConverged', ...
        ['No convergence after %d cycles (E0 %.3g, ET %.3g, V %.3g, complementarity %.3g). ' ...
         'No figures should be reported from this solve.'],it,E0,ET,V,comp);
end

% Final adjacent-edge plans for independent consistency checks and cost.
edgePlans = cell(np,nn-1); totalCost = 0; chainError = 0;
sourceCheck = zeros(nt,1); sinkCheck = zeros(nt,1);
nodeCheck = zeros(nt,S.nNodes);
for pid = 1:np
    for pos = 1:nn-1
        if pos==1, localLeft=logu; else, localLeft=logW(:,paths(pid,pos)); end
        if pos==nn-1, localRight=logv; else, localRight=logW(:,paths(pid,pos+1)); end
        left = F{pos}(:,pid)+localLeft;
        right = B{pos+1}(:,pid)+localRight;
        edge = exp(logK+left+right');
        edgePlans{pid,pos} = edge;
        totalCost = totalCost+sum(sum(S.cost.*edge));
        if pos==1, sourceCheck=sourceCheck+sum(edge,2); end
        if pos==nn-1, sinkCheck=sinkCheck+sum(edge,1)'; end
        if pos<nn-1
            v=paths(pid,pos+1);
            nodeCheck(:,v)=nodeCheck(:,v)+sum(edge,1)';
        end
        if pos>1
            prev = edgePlans{pid,pos-1};
            chainError=max(chainError,norm(sum(prev,1)'-sum(edge,2),1));
        end
    end
end
freshError = max([norm(sourceCheck-m0,1),norm(sinkCheck-mT,1), ...
    max(sum(abs(nodeCheck(:,S.middle)-mnode(:,S.middle)),1))]);
assert(freshError<1e-9 && chainError<1e-9, ...
    'example2:InconsistentMessages','Returned profiles disagree with edge plans.');
assert(all(cellfun(@(x)all(x(~S.support)==0),edgePlans(:))), ...
    'example2:Causality','A returned edge plan has noncausal mass.');

% Entropy identity on the path tensor, without forming that tensor.
scalingTerm=sum(m0(activeP).*logu(activeP))+sum(mT(activeQ).*logv(activeQ))+ ...
    sum(sum(mnode(:,S.middle).*logW(:,S.middle)));
entropy=totalCost/S.epsilon-scalingTerm+sum(pathMass);
R.m0=m0; R.mT=mT; R.mnode=mnode; R.pathMass=pathMass;
R.edgePlans=edgePlans; R.cost=totalCost; R.entropy=entropy;
R.regularizedObjective=totalCost-S.epsilon*entropy;
R.dual=dual; R.objectiveDifference=R.regularizedObjective-dual;
% The difference is a certified gap only for an exactly feasible plan.
R.E0=E0; R.ET=ET; R.V=V; R.complementarity=comp;
R.cycles=it; R.seconds=toc(started); R.history=history(isfinite(history(:,1)),:);
R.chainError=chainError; R.freshMarginalError=freshError;
R.logu=logu; R.logv=logv; R.logW=logW; R.converged=converged;
if opts.verbose
    fprintf('  converged: %d cycles, %.2f s; transport cost %.10g\n',it,R.seconds,totalCost);
end
end

function B=backward(logK,logv,logW,paths)
[np,nn]=size(paths); nt=numel(logv);
B=cell(1,nn); B{nn}=zeros(nt,np);
for pos=nn-1:-1:1
    if pos==nn-1, s=repmat(logv,1,np); else, s=logW(:,paths(:,pos+1)'); end
    B{pos}=logmultiply(logK,s+B{pos+1});
end
end

function F=forward(logK,logu,logW,paths)
[np,nn]=size(paths); nt=numel(logu);
F=cell(1,nn); F{1}=zeros(nt,np);
for pos=2:nn
    if pos==2, s=repmat(logu,1,np); else, s=logW(:,paths(:,pos-1)'); end
    F{pos}=logmultiply(logK',s+F{pos-1});
end
end

function [m0,mT,mnode,pathMass]=marginals(F,B,logu,logv,logW,paths,nNodes)
[np,nn]=size(paths); nt=numel(logu);
m0=exp(logu+logsum(B{1},2)); mT=exp(logv+logsum(F{nn},2));
mnode=zeros(nt,nNodes);
for pos=2:nn-1
    for v=unique(paths(:,pos))'
        idx=paths(:,pos)==v;
        mnode(:,v)=exp(logW(:,v)+logsum(F{pos}(:,idx)+B{pos}(:,idx),2));
    end
end
pathMass=exp(logsum(repmat(logu,1,np)+B{1},1));
end

function Y=logmultiply(A,X)
Y=zeros(size(A,1),size(X,2));
for k=1:size(X,2), Y(:,k)=logsum(A+X(:,k)',2); end
end

function y=logsum(X,dim)
m=max(X,[],dim); safe=m; safe(~isfinite(safe))=0;
y=safe+log(sum(exp(X-safe),dim));
end

function s=defaultOption(s,key,value)
if ~isfield(s,key), s.(key)=value; end
end

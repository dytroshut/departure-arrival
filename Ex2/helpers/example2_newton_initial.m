function [initial,info] = example2_newton_initial(S,capacity,initial)
%EXAMPLE2_NEWTON_INITIAL Toolbox-free dual warm start for a stiff LINE case.
% This is initialization, not Algorithm 1. The final reported solve uses
% exact cyclic Sinkhorn updates with the same cost, capacities and epsilon.
% A projected, damped Newton step is accepted only if it improves the dual.
assert(size(S.paths,1)==1,'example2:NewtonLineOnly','This initializer is for a single path.');
nodes=S.paths; nt=S.nt; nn=numel(nodes);
if isempty(initial)
    initial=struct('logu',zeros(nt,1),'logv',zeros(nt,1), ...
        'logW',zeros(nt,S.nNodes));
    initial.logu(S.p==0)=-inf; initial.logv(S.q==0)=-inf;
end
K=-inf(nt); K(S.support)=-S.cost(S.support)/S.epsilon;
x=zeros(nt,nn); targets=zeros(nt,nn);
x(:,1)=initial.logu; x(:,nn)=initial.logv;
targets(:,1)=S.p; targets(:,nn)=S.q;
for j=2:nn-1
    x(:,j)=initial.logW(:,nodes(j)); targets(:,j)=capacity(:,nodes(j));
end
activeBoundary=isfinite(x(:,1)) & S.p>0;
[~,gauge]=max(S.p); started=tic; accepted=0; success=false;
% Normalize only the total starting mass; subsequent steps optimize duals.
[~,~,~,lm]=evaluate(x,K,targets);
x(activeBoundary,1)=x(activeBoundary,1)-lm;
for iteration=1:500
    % Exact blocks also restore very small endpoint masses that would be
    % lost by an active-set Hessian cutoff alone.
    x=cycle(x,K,targets);
    [value,m,edges,~,g]=evaluate(x,K,targets);
    assert(isfinite(value) && all(isfinite(m(:))) && ...
        all(isfinite(x(:,2:end-1)),'all'),'example2:NewtonNonfinite', ...
        'The initialization encountered nonfinite dual data.');
    err=max([norm(m(:,1)-S.p,1),norm(m(:,end)-S.q,1), ...
        sum(max(m(:,2:end-1)-targets(:,2:end-1),0),'all')]);
    comp=S.epsilon*sum(abs(x(:,2:end-1).* ...
        (targets(:,2:end-1)-m(:,2:end-1))),'all');
    if mod(iteration,20)==1
        fprintf('  Newton initialization %d: residual %.3e, complementarity %.3e\n', ...
            iteration,err,comp);
    end
    if err<=1e-10 && comp<=1e-9*max(1,abs(S.epsilon*value))
        success=true; break;
    end

    % Hessian of minus the dual is made of pairwise joint marginals.
    % Conditional transitions are computed from EDGE PLANS, not by
    % normalizing the Gibbs kernel or changing the transport cost.
    H=zeros(nt*nn); transitions=cell(1,nn-1);
    for i=1:nn-1
        transitions{i}=zeros(nt);
        nz=m(:,i)>0;
        transitions{i}(nz,:)=edges{i}(nz,:)./m(nz,i);
    end
    for i=1:nn
        rows=(i-1)*nt+(1:nt); joint=diag(m(:,i)); H(rows,rows)=joint;
        for j=i+1:nn
            joint=joint*transitions{j-1}; cols=(j-1)*nt+(1:nt);
            H(rows,cols)=joint; H(cols,rows)=joint';
        end
    end
    free=isfinite(x) & m>1e-13;
    free(:,2:end-1)=free(:,2:end-1) & ...
        (x(:,2:end-1)<-1e-8 | g(:,2:end-1)>0);
    free(gauge,1)=false; ids=find(free); gradient=g(ids);
    block=H(ids,ids); block=(block+block')/2;
    ridge=1e-12;
    for repair=1:8
        [factor,flag]=chol(block+ridge*eye(numel(ids)));
        if flag==0, break; end
        ridge=10*ridge;
    end
    direction=zeros(size(x));
    if flag==0
        direction(ids)=-(factor\(factor'\gradient));
    else
        direction(ids)=-gradient./max(m(ids),1e-8);
    end
    if sum(g(free).*direction(free))>=0
        direction(ids)=-gradient./max(m(ids),1e-8);
    end
    step=1; improved=false;
    for trial=1:40
        candidate=x+step*direction;
        candidate(:,2:end-1)=min(candidate(:,2:end-1),0);
        candidateValue=evaluate(candidate,K,targets);
        if isfinite(candidateValue) && candidateValue>=value
            x=candidate; improved=true; accepted=accepted+1; break;
        end
        step=step/2;
    end
    if ~improved
        % A fully exact cyclic update is a safe fallback at active-set kinks.
        x=cycle(x,K,targets);
    end
end
initial.logu=x(:,1); initial.logv=x(:,nn);
for j=2:nn-1, initial.logW(:,nodes(j))=x(:,j); end
info.iterations=iteration; info.seconds=toc(started);
info.acceptedSteps=accepted; info.converged=success;
info.residual=err; info.complementarity=comp;
fprintf('  initialization finished: %d Newton steps, %.2f s; residual %.3e\n', ...
    iteration,info.seconds,err);
% A warm start need not be optimal; only the following exact solve decides
% convergence and permits export. No figures are produced by this function.
end

function [value,m,edges,lm,g]=evaluate(x,K,targets)
[nt,nn]=size(x); F=zeros(nt,nn); B=F;
for i=2:nn, F(:,i)=lsum(K'+(F(:,i-1)+x(:,i-1))',2); end
for i=nn-1:-1:1, B(:,i)=lsum(K+(B(:,i+1)+x(:,i+1))',2); end
lm=lsum(F(:,nn)+x(:,nn),1); mass=exp(lm);
finite=isfinite(x); value=sum(targets(finite).*x(finite))-mass;
if nargout>1
    m=exp(F+B+x); edges=cell(1,nn-1); g=m-targets;
    for i=1:nn-1
        edges{i}=exp(K+(F(:,i)+x(:,i))+(B(:,i+1)+x(:,i+1))');
    end
end
end

function x=cycle(x,K,targets)
[nt,nn]=size(x); B=zeros(nt,nn);
for i=nn-1:-1:1, B(:,i)=lsum(K+(B(:,i+1)+x(:,i+1))',2); end
ap=targets(:,1)>0; x(:,1)=-inf;
x(ap,1)=log(targets(ap,1))-B(ap,1);
F=zeros(nt,nn);
for i=2:nn, F(:,i)=lsum(K'+(F(:,i-1)+x(:,i-1))',2); end
aq=targets(:,nn)>0; x(:,nn)=-inf;
x(aq,nn)=log(targets(aq,nn))-F(aq,nn);
for i=nn-1:-1:1, B(:,i)=lsum(K+(B(:,i+1)+x(:,i+1))',2); end
prefix=x(:,1);
for i=2:nn-1
    incoming=lsum(K'+prefix',2);
    x(:,i)=min(log(targets(:,i))-incoming-B(:,i),0);
    prefix=incoming+x(:,i);
end
end

function y=lsum(A,dim)
m=max(A,[],dim); safe=m; safe(~isfinite(safe))=0;
y=safe+log(sum(exp(A-safe),dim));
end

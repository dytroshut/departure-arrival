function verify_example2
%VERIFY_EXAMPLE2 Small exact regression tests; no extra toolboxes required.
repoRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'examples','example2_nodal_capacities','helpers'));
options=struct('verbose',false,'tol',1e-10,'complementarityTol',1e-10);

% Single intermediate node: compare with the analytically clipped Gibbs law.
S=example2_data('line',10);
S.paths=1:3; S.source=1; S.sink=3; S.middle=2; S.nNodes=3; S.nEdges=2;
S.p=zeros(S.nt,1); S.p(1)=1; S.q=zeros(S.nt,1); S.q(end)=1;
cap=inf(S.nt,3); cap(:,2)=.3;
R=example2_sinkhorn(S,cap,options);
idx=2:S.nt-1; c=1./S.t(idx)+1./(1-S.t(idx));
lo=-1000; hi=1000;
for k=1:100
    mid=(lo+hi)/2;
    if sum(min(exp(mid-c/S.epsilon),.3))>1, hi=mid; else, lo=mid; end
end
truth=zeros(S.nt,1); truth(idx)=min(exp((lo+hi)/2-c/S.epsilon),.3);
assert(norm(R.mnode(:,2)-truth,1)<1e-8,'example2:Test','Clipped Gibbs test failed.');
assert(abs(R.cost-sum(c.*truth(idx)))<1e-8,'example2:Test','Cost contraction failed.');
fprintf('PASS: capped one-node solution agrees with the exact Gibbs law.\n');

% The optional Newton initializer must agree with the same exact solution;
% the final checked result still comes from ordinary cyclic Sinkhorn.
N=example2_solve(S,cap,struct('epsilonSchedule',S.epsilon));
info=N.newtonInitialization;
assert(info.converged && norm(N.mnode(:,2)-truth,1)<1e-8, ...
    'example2:Test','The optional Newton initialization changed the solution.');
fprintf('PASS: Newton initialization followed by exact scaling preserves the solution.\n');

% Three paths, fixed endpoint cells: enumerate the FULL path tensors.
S=example2_data('shared',10);
S.p=zeros(S.nt,1); S.p(1)=1; S.q=zeros(S.nt,1); S.q(end)=1;
cap=inf(S.nt,S.nNodes); cap(:,S.middle)=2;
options.initial=struct('logu',zeros(S.nt,1),'logv',zeros(S.nt,1), ...
    'logW',-50*ones(S.nt,S.nNodes));
options.initial.logu(S.p==0)=-inf; options.initial.logv(S.q==0)=-inf;
R=example2_sinkhorn(S,cap,options);
tuples=nchoosek(2:S.nt-1,4); costs=zeros(size(tuples,1),1);
for k=1:size(tuples,1)
    times=S.t([1,tuples(k,:),S.nt]); costs(k)=sum(1./diff(times));
end
lw=-costs/S.epsilon; weights=exp(lw-max(lw)); weights=weights/sum(weights);
truth=zeros(S.nt,S.nNodes);
for pid=1:size(S.paths,1)
    for k=1:size(tuples,1)
        for pos=1:4
            v=S.paths(pid,pos+1);
            truth(tuples(k,pos),v)=truth(tuples(k,pos),v)+weights(k)/size(S.paths,1);
        end
    end
end
assert(max(sum(abs(R.mnode(:,S.middle)-truth(:,S.middle)),1))<1e-8, ...
    'example2:Test','Full-tensor marginal comparison failed.');
assert(abs(R.cost-sum(costs.*weights))<1e-8,'example2:Test','Full-tensor cost comparison failed.');
assert(max(abs(R.pathMass-1/3))<1e-8,'example2:Test','Path allocation is incorrect.');
assert(max(abs(R.logW(:,S.middle)),[],'all')<1e-10, ...
    'example2:Test','Previously active multipliers did not release.');
fprintf('PASS: shared-node marginals and costs agree with explicit full tensors.\n');
fprintf('PASS: capacity multipliers can increase back to one.\n');

% Active SHARED bounds: independently scale the explicitly enumerated
% full tensors, rather than reusing the graphical-message implementation.
% A less stiff epsilon keeps this independent regression test small.
% The full-sized paper experiments still use their prescribed .065.
S.epsilon=.5; cap(:,S.middle)=.30; options.initial=[];
R=example2_sinkhorn(S,cap,options);
np=size(S.paths,1); nc=size(tuples,1);
allTuples=repmat(tuples,np,1); logMass=repmat(-costs/S.epsilon,np,1);
membership=zeros(np*nc,4);
for pid=1:np
    membership((pid-1)*nc+(1:nc),:)=repmat(S.paths(pid,2:end-1),nc,1);
end
logW=zeros(S.nt,S.nNodes); done=false;
for cycle=1:100000
    update=-testLogSum(logMass); logMass=logMass+update;
    for pos=1:4
        for v=unique(S.paths(:,pos+1))'
            for j=1:S.nt
                selected=membership(:,pos)==v & allTuples(:,pos)==j;
                lm=testLogSum(logMass(selected));
                next=min(logW(j,v)+log(cap(j,v))-lm,0);
                logMass(selected)=logMass(selected)+next-logW(j,v);
                logW(j,v)=next;
            end
        end
    end
    masses=exp(logMass); truth=zeros(S.nt,S.nNodes);
    for pos=1:4
        for v=unique(S.paths(:,pos+1))'
            for j=1:S.nt
                selected=membership(:,pos)==v & allTuples(:,pos)==j;
                truth(j,v)=truth(j,v)+sum(masses(selected));
            end
        end
    end
    comp=S.epsilon*sum(abs(logW(:,S.middle).* ...
        (cap(:,S.middle)-truth(:,S.middle))),'all');
    if abs(sum(masses)-1)<1e-11 && comp<1e-10
        done=true; break;
    end
end
assert(done,'example2:Test','Independent full-tensor scaling did not converge.');
assert(max(sum(abs(R.mnode(:,S.middle)-truth(:,S.middle)),1))<1e-8, ...
    'example2:Test','Active shared-capacity marginals disagree with full tensors.');
assert(abs(R.cost-sum(repmat(costs,np,1).*masses))<1e-8, ...
    'example2:Test','Active shared-capacity cost disagrees with full tensors.');
assert(any(R.logW(:,S.middle)<-1e-6,'all'), ...
    'example2:Test','Shared-capacity test failed to activate a bound.');
assert(all(diff(R.history(:,5))>=-1e-8), ...
    'example2:Test','The cyclic dual values decreased.');
fprintf('PASS: active shared capacities agree with independent full-tensor scaling.\n');

% The infeasible original line settings must be rejected, not plotted.
S=example2_data('line'); oldCapacity=S.varying*(.0100/.0105);
rejected=false;
try
    example2_sinkhorn(S,oldCapacity,struct('verbose',false,'maxIter',1));
catch problem
    rejected=strcmp(problem.identifier,'example2:InfeasibleLine');
end
assert(rejected,'example2:Test','The infeasible original settings were not rejected.');
fprintf('PASS: original infeasible line capacities are rejected before scaling.\n');
fprintf('All Example 2 regression tests passed.\n');
end

function y=testLogSum(x)
if isempty(x), y=-inf; return; end
m=max(x); y=m+log(sum(exp(x-m)));
end

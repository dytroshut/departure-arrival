function R = example2_solve(S, capacity, runOptions)
%EXAMPLE2_SOLVE Initialize by epsilon continuation, then solve at .065.
% All reported costs, residuals, and profiles are from the FINAL epsilon.
% The higher-epsilon solves only provide an initial set of dual scalings.
% runOptions.epsilonSchedule overrides the continuation sequence.
% Set runOptions.useNewtonInitialization=false to disable that warm start.
schedule = S.epsilon*[16 4 1];
% The sinusoidal line bottleneck needs a coarser entropy-only warm start.
% This affects initialization only, not its final weights or epsilon.
predict = strcmp(S.name,'line') && any(any(abs(capacity(:,S.middle)-S.uniform(:,S.middle))>1e-12));
if predict
    schedule = S.epsilon*[256 64 16 4 1];
end
if nargin<3, runOptions=struct; end
if isfield(runOptions,'epsilonSchedule'), schedule=runOptions.epsilonSchedule; end
assert(~isempty(schedule) && all(schedule>0) && schedule(end)==S.epsilon, ...
    'example2:Schedule','The positive epsilon schedule must end at the target epsilon.');
initial=[]; previous=[]; warmCycles=0; totalSeconds=0; phaseCycles=zeros(size(schedule));
newtonInfo=[];
for k=1:numel(schedule)
    current=S; current.epsilon=schedule(k);
    options=struct('initial',initial,'checkEvery',10,'maxIter',100000);
    if isfield(runOptions,'maxIter'), options.maxIter=runOptions.maxIter; end
    if k<numel(schedule)
        options.tol=1e-5; options.complementarityTol=1e-5;
    else
        options.tol=1e-8; options.complementarityTol=1e-8;
    end
    fprintf('  initialization/solve phase %d: epsilon %.8g\n',k,current.epsilon);
    if predict && k==numel(schedule) && ...
            (~isfield(runOptions,'useNewtonInitialization') || runOptions.useNewtonInitialization)
        % The near-saturated varying line case is exceptionally stiff at
        % .065. Refine the starting dual, then run EXACT cyclic scaling.
        [initial,newtonInfo]=example2_newton_initial(current,capacity,initial);
        options.initial=initial;
        totalSeconds=totalSeconds+newtonInfo.seconds;
    end
    R=example2_sinkhorn(current,capacity,options);
    phaseCycles(k)=R.cycles; totalSeconds=totalSeconds+R.seconds;
    if isfield(runOptions,'checkpointDirectory')
        folder=runOptions.checkpointDirectory;
        if ~isfolder(folder), mkdir(folder); end
        save(fullfile(folder,sprintf('phase_%d.mat',k)),'S','capacity','R');
    end
    if k<numel(schedule)
        activeP=S.p>0; activeQ=S.q>0;
        % Fix the harmless source/sink gauge before extrapolating duals.
        shift=sum(S.p(activeP).*R.logu(activeP));
        u=R.logu-shift; v=R.logv+shift;
        potential.u=zeros(S.nt,1); potential.v=zeros(S.nt,1);
        potential.u(activeP)=schedule(k)*u(activeP);
        potential.v(activeQ)=schedule(k)*v(activeQ);
        potential.W=schedule(k)*R.logW;
        next=potential;
        if predict && ~isempty(previous) && k+1<numel(schedule)
            ratio=(schedule(k+1)-schedule(k))/(schedule(k)-schedule(k-1));
            next.u=potential.u+ratio*(potential.u-previous.u);
            next.v=potential.v+ratio*(potential.v-previous.v);
            next.W=potential.W+ratio*(potential.W-previous.W);
        end
        initial=struct('logu',-inf(S.nt,1),'logv',-inf(S.nt,1), ...
            'logW',min(next.W/schedule(k+1),0));
        initial.logu(activeP)=next.u(activeP)/schedule(k+1);
        initial.logv(activeQ)=next.v(activeQ)/schedule(k+1);
        previous=potential;
        warmCycles=warmCycles+R.cycles;
    end
end
R.epsilon=S.epsilon; R.initializationCycles=warmCycles;
R.epsilonSchedule=schedule; R.phaseCycles=phaseCycles;
R.totalCycles=sum(phaseCycles); R.totalSeconds=totalSeconds;
R.dualPredictionUsed=predict;
R.newtonInitialization=newtonInfo;
end

function example2_save(S,uniform,varying,outputDirectory)
%EXAMPLE2_SAVE Save checked numerical data, diagnostics, and paper figures.
assert(uniform.converged && varying.converged,'example2:UncheckedResults', ...
    'Both solves must converge before saving a comparison.');
assert(uniform.epsilon==S.epsilon && varying.epsilon==S.epsilon, ...
    'example2:WrongEpsilon','Only final-epsilon results may be reported.');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
save(fullfile(outputDirectory,[S.name '_results.mat']),'S','uniform','varying');
f=fopen(fullfile(outputDirectory,[S.name '_diagnostics.txt']),'w');
assert(f>=0,'example2:Output','Cannot open diagnostics file.');
cleaner=onCleanup(@()fclose(f));
fprintf(f,'Example 2: %s\nMATLAB: %s\n',S.name,version);
fprintf(f,'Grid points: %d; final epsilon: %.12g; edge weight: %.12g\n', ...
    S.nt,S.epsilon,S.edgeWeight);
fprintf(f,'Common mean cell capacity: %.12g\n',S.capacity);
fprintf(f,'Deleted impossible departure/arrival mass: %.12g / %.12g\n\n', ...
    S.removedDepartureMass,S.removedArrivalMass);
for k=1:2
    if k==1, R=uniform; label='Uniform'; else, R=varying; label='Time-varying'; end
    fprintf(f,'%s\n',label);
    fprintf(f,'Epsilon initialization schedule: %s\n',mat2str(R.epsilonSchedule));
    fprintf(f,'Cycles per phase: %s; final cycles: %d\n',mat2str(R.phaseCycles),R.cycles);
    fprintf(f,'Total solve seconds: %.6g\n',R.totalSeconds);
    if isfield(R,'newtonInitialization') && ~isempty(R.newtonInitialization)
        N=R.newtonInitialization;
        fprintf(f,'Newton initialization steps: %d; seconds: %.6g; residual: %.12g\n', ...
            N.iterations,N.seconds,N.residual);
    end
    fprintf(f,'Departure L1 error: %.12g\nArrival L1 error: %.12g\n',R.E0,R.ET);
    fprintf(f,'Total capacity excess: %.12g\nComplementarity residual: %.12g\n', ...
        R.V,R.complementarity);
    fprintf(f,'Edge conservation error: %.12g\nFresh marginal cross-check: %.12g\n', ...
        R.chainError,R.freshMarginalError);
    fprintf(f,'Transport cost: %.12g\nRegularized objective: %.12g\nDual: %.12g\n', ...
        R.cost,R.regularizedObjective,R.dual);
    fprintf(f,'Objective minus dual (approximately feasible plan): %.12g\n',R.objectiveDifference);
    fprintf(f,'Path masses: %s\n\n',mat2str(R.pathMass,12));
end
example2_plot(S,uniform,varying,fullfile(outputDirectory,'figures'));
end

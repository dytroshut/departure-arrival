function [S, uniform, varying] = example2_compare(which, outputDirectory)
%EXAMPLE2_COMPARE Solve and export both capacity cases without stale caches.
if nargin<2
    outputDirectory=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
end
S=example2_data(which);
fprintf('\nExample 2 %s: nt=%d, epsilon=%.5g, mean cell capacity=%.5g\n', ...
    S.name,S.nt,S.epsilon,S.capacity);
fprintf('Removed impossible boundary mass: departure %.3g; arrival %.3g\n', ...
    S.removedDepartureMass,S.removedArrivalMass);
fprintf('\nUniform capacity\n'); uniform=example2_solve(S,S.uniform);
fprintf('\nTime-varying capacity\n'); varying=example2_solve(S,S.varying);
example2_save(S,uniform,varying,outputDirectory);
end

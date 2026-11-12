% Re-test the pairs an exhaustive run rejected, with the per-denominator scaling bound (fork
% commit after 28350a9), then rebuild the classes and write a retested .mat and .json.
%
%   octave --no-gui retest_rejected.m <model> [<pairs.txt>]
% reads  reference_final/<model>_concordant_fixed.mat and reference_final/retest/pairs_<model>.txt
% writes reference_final/retest/<model>_concordant_fixed_retest.mat and reference_<model>_retest.json
args = argv(); m = args{1};
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat')); addpath(genpath(fullfile(getenv('UPSTREAM_ALGORITHM'), 'Code')));
bench = fullfile(here, '..');
load(fullfile(bench, 'reference_final', [m '_concordant_fixed.mat']));    % Results_balanced, class_with_balanced, B, CC, group
model = Results_balanced.MODEL_r{1}; nc = numel(model.complexes);
pairs_path = fullfile(bench, 'reference_final', 'retest', ['pairs_' m '.txt']);
if numel(args) > 1, pairs_path = args{2}; end
At = dlmread(pairs_path);
printf('%s: %d rejected pairs to re-test\n', m, size(At,1)); fflush(stdout);
t0 = tic;
CC2 = find_concordant_complexes(model, group, sparse(nc, nc), At, 1, size(At,1));
t_retest = toc(t0);
acc = arrayfun(@(k) full(CC2(At(k,1),At(k,2))) == 2, 1:size(At,1));
printf('re-accepted %d, residual %d, %.0f s\n', sum(acc), sum(~acc), t_retest); fflush(stdout);
for k = find(acc), CC(At(k,1),At(k,2)) = 2; end
% classes: connected components over 2 / -1, balanced excluded (as run_reference.m)
E = (CC == 2) | (CC == -1); E = E | E';
module_of = zeros(nc, 1); next_module = 0;
for i = 1:nc
  if module_of(i) ~= 0 || any(B == i), continue; end
  comp = i; frontier = i;
  while ~isempty(frontier)
    nb = find(any(E(frontier, :), 1)); new = setdiff(nb, comp); comp = [comp new]; frontier = new;
  end
  comp = setdiff(comp, B);
  if numel(comp) > 1, next_module = next_module + 1; module_of(comp) = next_module; else, module_of(comp) = -1; end
end
module_of(B) = 0;
class_with_balanced = {};
if ~isempty(B), class_with_balanced{end+1} = B(:); end
for q = 1:next_module, class_with_balanced{end+1} = find(module_of == q)(:); end
outdir = fullfile(bench, 'reference_final', 'retest');
save('-v7', fullfile(outdir, [m '_concordant_fixed_retest.mat']), 'Results_balanced', 'class_with_balanced', 'B', 'CC', 'group');
crn = jsondecode(fileread(fullfile(bench, 'unblocked', [m '_fixed.crn.json'])));
orig = jsondecode(fileread(fullfile(bench, 'reference_final', ['reference_' m '.json'])));
write_results_json(fullfile(outdir, ['reference_' m '_retest.json']), crn, model, module_of, orig.seconds + t_retest, ...
                   orig.seconds_balanced, orig.seconds_concordant + t_retest, orig.candidate_pairs, 'exhaustive+retest', ...
                   orig.pairs_tested, true, 0, linprog('nonoptimal'));
fid = fopen(fullfile(outdir, ['retest_' m '.txt']), 'w');
fprintf(fid, 'retested %d accepted %d residual %d seconds %.1f\n', size(At,1), sum(acc), sum(~acc), t_retest);
for k = find(~acc), fprintf(fid, 'residual %d %d\n', At(k,1)-1, At(k,2)-1); end
fclose(fid);
printf('modules %d, written retest outputs for %s\n', next_module, m); fflush(stdout);

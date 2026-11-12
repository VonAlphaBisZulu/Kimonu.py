% Unit test for the two concordance fixes, on pairs whose truth is known independently.
%
% The pairs below were adjudicated against the definition with flux vectors, outside both
% implementations: the first group is one concordance class (constant ratio, verified to a spread
% of 6e-11), the second group are pairs whose ratio demonstrably varies but stays inside a single
% round(.,2) bin near zero.
%
% Usage: octave --no-gui test_fixed_verdict.m <crn.json>

args = argv();
crn_path = args{1};
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat'));
addpath(genpath(fullfile(getenv('UPSTREAM_ALGORITHM'), 'Code')));

crn = jsondecode(fileread(crn_path));
nc = crn.n_complexes; nr = crn.n_reactions; nm = crn.n_metabolites;
model.A = triplets_to_sparse(crn.A, nc, nr);
model.Y = triplets_to_sparse(crn.Y, nm, nc);
model.S = triplets_to_sparse(crn.S, nm, nr);
model.lb = crn.lb(:); model.ub = crn.ub(:);
model.b = zeros(nm, 1); model.c = zeros(nr, 1);
model.complexes = 1:nc;

t0 = tic;
[B, group] = find_balanced_complexes(model);
printf('balanced %d in %.0f s\n', numel(B), toc(t0)); fflush(stdout);

% should be concordant: members of one verified class (1-based; the JSON indices are 0-based)
want_yes = [2636 2638; 2636 2640; 2638 2640; 2636 2649];
% should NOT be concordant: ratio varies, but within one round(.,2) bin near zero
want_no  = [2245 1626; 2404 3296; 957 3226; 1843 2369];

At = [want_yes; want_no];
CC = sparse(nc, nc);
CC = find_concordant_complexes(model, group, CC, At, 1, size(At,1));

ok = 0; bad = 0;
for k = 1:size(At,1)
  v = CC(At(k,1), At(k,2));
  expect = k <= size(want_yes,1);
  got = (v == 2);
  tag = 'FAIL';
  if got == expect, tag = 'ok'; ok = ok + 1; else bad = bad + 1; end
  printf('%-4s pair (%d,%d) expected %d got %d\n', tag, At(k,1), At(k,2), expect, got);
end
printf('\n%d ok, %d failed\n', ok, bad);

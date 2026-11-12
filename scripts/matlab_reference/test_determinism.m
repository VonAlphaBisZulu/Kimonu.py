% Run the reference's pair LP twice on a list of pairs and report whether verdicts agree.
% Usage: octave --no-gui test_determinism.m <crn.json> <pairs.txt>
args = argv(); here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat')); upstream = getenv('UPSTREAM_ALGORITHM'); if isempty(upstream), upstream = fullfile(here, '..', 'reference', 'Upstream_Algorithm'); end
addpath(genpath(fullfile(upstream, 'Code')));   % bash scripts/reference/fetch_upstream.sh
crn = jsondecode(fileread(args{1}));
nc = crn.n_complexes; nr = crn.n_reactions; nm = crn.n_metabolites;
model.A = triplets_to_sparse(crn.A, nc, nr); model.Y = triplets_to_sparse(crn.Y, nm, nc);
model.S = triplets_to_sparse(crn.S, nm, nr); model.lb = crn.lb(:); model.ub = crn.ub(:);
model.b = zeros(nm,1); model.c = zeros(nr,1); model.complexes = 1:nc;
[B, group] = find_balanced_complexes(model);
At = dlmread(args{2});
printf('%d pairs\n', size(At,1));
v = zeros(size(At,1), 2);
for rep = 1:2
  CC = sparse(nc, nc);
  CC = find_concordant_complexes(model, group, CC, At, 1, size(At,1));
  for k = 1:size(At,1), v(k,rep) = CC(At(k,1), At(k,2)); end
end
for k = 1:size(At,1)
  printf('pair (%d,%d): run1=%d run2=%d %s\n', At(k,1), At(k,2), v(k,1), v(k,2), ifelse(v(k,1)==v(k,2),'','<-- NONDETERMINISTIC'));
end
printf('\nconcordant in run1: %d  run2: %d  disagreements: %d\n', sum(v(:,1)==2), sum(v(:,2)==2), sum(v(:,1)~=v(:,2)));

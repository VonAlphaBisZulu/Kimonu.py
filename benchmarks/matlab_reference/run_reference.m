% Run the MATLAB reference implementation (ankueken/Upstream_Algorithm) under Octave.
%
% The reference's own .m files are used unmodified; only linprog and the two-argument form of
% round are supplied by compat/, because Octave has neither. The complex decomposition (A, Y) and
% the flux space come from the shared exported network, so this measures the reference's decision
% procedure against the other implementations on identical input.
%
% Every candidate pair is tested. The reference normally runs a sampling prefilter first
% (pre_concordant_sampling), which is stochastic and only ever *removes* candidates; testing all
% pairs gives the reference formulation's exact answer and makes the run deterministic.
%
% Usage:
%   octave --no-gui run_reference.m <crn.json> <out.json> [max_pairs]

args = argv();
if numel(args) < 2
  error('usage: run_reference.m <crn.json> <out.json> [max_pairs]');
end
crn_path = args{1};
out_path = args{2};
max_pairs = Inf;
if numel(args) >= 3
  max_pairs = str2double(args{3});
end

here = fileparts(mfilename('fullpath'));
if isempty(here)
  here = pwd();
end
addpath(fullfile(here, 'compat'));
upstream = getenv('UPSTREAM_ALGORITHM');
if isempty(upstream)
  upstream = fullfile(here, '..', '..', '..', 'Upstream_Algorithm');
end
if ~exist(fullfile(upstream, 'Code'), 'dir')
  error(['reference code not found at ' upstream '; set UPSTREAM_ALGORITHM to a clone of https://github.com/ankueken/Upstream_Algorithm']);
end
addpath(genpath(fullfile(upstream, 'Code')));

printf('reading %s\n', crn_path); fflush(stdout);
raw = fileread(crn_path);
crn = jsondecode(raw);

nc = crn.n_complexes;
nr = crn.n_reactions;
nm = crn.n_metabolites;

model.A = triplets_to_sparse(crn.A, nc, nr);
model.Y = triplets_to_sparse(crn.Y, nm, nc);
model.S = triplets_to_sparse(crn.S, nm, nr);
model.lb = crn.lb(:);
model.ub = crn.ub(:);
model.b = zeros(nm, 1);
model.c = zeros(nr, 1);
model.complexes = 1:nc;

printf('complexes %d, reactions %d, metabolites %d\n', nc, nr, nm); fflush(stdout);

%% balanced complexes -- reference function, unmodified
t0 = tic;
[B, group] = find_balanced_complexes(model);
t_balanced = toc(t0);
printf('balanced: %d complexes in %.1f s\n', numel(B), t_balanced); fflush(stdout);

%% candidate pairs
% CC encodes the relationships: 1 balanced, -1 trivially concordant, 2 concordant.
% This reproduces the bookkeeping of pre_concordant_sampling without its sampling step.
CC = sparse(nc, nc);
if ~isempty(B)
  CC(B, B) = 1;
end

species_degree_complexes = sum(model.Y ~= 0, 2);
idx = find(species_degree_complexes == 2);
for i = 1:numel(idx)
  tcc = find(model.Y(idx(i), :) ~= 0);
  if numel(tcc) == 2 && all(all(CC(tcc, tcc) ~= 1))
    CC(tcc, tcc) = [0 -1; -1 0];
  end
end

At = nchoosek(1:nc, 2);
CC_temp = CC;
if ~isempty(B)
  CC_temp(B, :) = 5;
  CC_temp(:, B) = 5;
end
[checked_rows, checked_cols] = find(triu(CC_temp, 1) ~= 0);
At = setdiff(At, [checked_rows checked_cols; checked_cols checked_rows], 'rows');
printf('candidate pairs: %d\n', size(At, 1)); fflush(stdout);

if size(At, 1) > max_pairs
  error('candidate pairs %d exceed max_pairs %d', size(At, 1), max_pairs);
end

%% concordance -- reference function, unmodified
t0 = tic;
CC = find_concordant_complexes(model, group, CC, At, 1, size(At, 1));
t_concordant = toc(t0);
printf('concordance done in %.1f s\n', t_concordant); fflush(stdout);

%% connected components over the concordance relation
% Edges are the non-trivial (2) and trivial (-1) concordances; value 1 marks balanced complexes,
% which are not concordance edges.
E = (CC == 2) | (CC == -1);
E = E | E';
module_of = zeros(nc, 1);
next_module = 0;
for i = 1:nc
  if module_of(i) ~= 0 || any(B == i)
    continue
  end
  comp = i;
  frontier = i;
  while ~isempty(frontier)
    nb = find(any(E(frontier, :), 1));
    new = setdiff(nb, comp);
    comp = [comp new];
    frontier = new;
  end
  comp = setdiff(comp, B);
  if numel(comp) > 1
    next_module = next_module + 1;
    module_of(comp) = next_module;
  else
    module_of(comp) = -1;
  end
end
module_of(B) = 0;

printf('modules: %d, balanced: %d\n', next_module, numel(B)); fflush(stdout);

%% write results in the format the comparison script expects
write_results_json(out_path, crn, model, module_of, t_balanced + t_concordant, ...
                   t_balanced, t_concordant, size(At, 1));
printf('written %s\n', out_path);

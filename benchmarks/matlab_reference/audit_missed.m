% For each pair the reference missed, report why: the sign class of both complexes and the LP
% outcome in the orientation the reference actually tests (numerator = lower index).
args = argv();
crn_path = args{1};
missed_path = args{2};

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat'));
upstream = getenv('UPSTREAM_ALGORITHM');
if isempty(upstream)
  upstream = fullfile(here, '..', '..', '..', 'Upstream_Algorithm');
end
if ~exist(fullfile(upstream, 'Code'), 'dir')
  error(['reference code not found at ' upstream '; set UPSTREAM_ALGORITHM to a clone of https://github.com/ankueken/Upstream_Algorithm']);
end
addpath(genpath(fullfile(upstream, 'Code')));

crn = jsondecode(fileread(crn_path));
nc = crn.n_complexes; nr = crn.n_reactions; nm = crn.n_metabolites;
model.A = triplets_to_sparse(crn.A, nc, nr);
model.S = triplets_to_sparse(crn.S, nm, nr);
model.lb = crn.lb(:); model.ub = crn.ub(:); model.b = zeros(nm,1); model.c = zeros(nr,1);
[B, group] = find_balanced_complexes(model);

missed = jsondecode(fileread(missed_path));
n = size(model.S, 2);
a = sparse([speye(n) -model.ub; -speye(n) model.lb]);
options = optimset('linprog'); options.Display = 'off';

printf('%-6s %-6s %-3s %-3s %-10s %-14s %-14s %s\n', 'num', 'den', 'gN', 'gD', 'branch', 'max', 'min', 'verdict');
counts = struct('infeasible', 0, 'bounds_differ', 0, 'accepted', 0);
for k = 1:size(missed, 1)
  inum = missed(k, 1); iden = missed(k, 2);
  gd = group(iden);
  if gd == 'N'
    branches = {'negative'};
  elseif gd == 'P'
    branches = {'positive'};
  else
    branches = {'positive', 'negative'};
  end
  Aeq = [model.S sparse(size(model.S,1),1); model.A(iden,:) 0];
  beq = [model.b; 1];
  for bi = 1:numel(branches)
    if strcmp(branches{bi}, 'positive'), tlb = 0; tub = 999; else tlb = -999; tub = 0; end
    lb = [-ones(n,1)*1e9; tlb]; ub = [ones(n,1)*1e9; tub];
    [x, ~, ef1] = linprog([-model.A(inum,:) 1], a, zeros(2*n,1), Aeq, beq, lb, ub, options);
    [y, ~, ef2] = linprog([ model.A(inum,:) 1], a, zeros(2*n,1), Aeq, beq, lb, ub, options);
    if ef1 == 1 && ef2 == 1
      mx = (model.A(inum,:)*x(1:end-1)) / (model.A(iden,:)*x(1:end-1));
      mn = (model.A(inum,:)*y(1:end-1)) / (model.A(iden,:)*y(1:end-1));
      if round(mx,2) == round(mn,2)
        verdict = 'accepted'; counts.accepted = counts.accepted + 1;
      else
        verdict = 'bounds differ'; counts.bounds_differ = counts.bounds_differ + 1;
      end
      printf('%-6d %-6d %-3s %-3s %-10s %-14.6g %-14.6g %s\n', inum, iden, group(inum), gd, branches{bi}, mx, mn, verdict);
    else
      counts.infeasible = counts.infeasible + 1;
      printf('%-6d %-6d %-3s %-3s %-10s %-14s %-14s LP infeasible (exitflags %d/%d)\n', ...
             inum, iden, group(inum), gd, branches{bi}, 'NaN', 'NaN', ef1, ef2);
    end
  end
end
printf('\nsummary: infeasible %d, bounds differ %d, accepted %d\n', counts.infeasible, counts.bounds_differ, counts.accepted);

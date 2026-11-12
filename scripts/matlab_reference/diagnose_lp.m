% Replicate the reference's Charnes-Cooper LPs for one ordered pair, printing what the reference
% only ever sees as a yes/no. Distinguishes "the LP failed" from "the bounds genuinely differ".
args = argv();
crn_path = args{1};
inum = str2double(args{2});   % numerator complex
iden = str2double(args{3});   % denominator complex

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat'));
upstream = getenv('UPSTREAM_ALGORITHM');
if isempty(upstream)
  upstream = fullfile(here, '..', 'reference', 'Upstream_Algorithm');
end
if ~exist(fullfile(upstream, 'Code'), 'dir')
  error(['reference code not found at ' upstream '; run bash scripts/reference/fetch_upstream.sh or set UPSTREAM_ALGORITHM to a clone of the fork at 8a6ed17']);
end
addpath(genpath(fullfile(upstream, 'Code')));

crn = jsondecode(fileread(crn_path));
nc = crn.n_complexes; nr = crn.n_reactions; nm = crn.n_metabolites;
model.A = triplets_to_sparse(crn.A, nc, nr);
model.S = triplets_to_sparse(crn.S, nm, nr);
model.lb = crn.lb(:); model.ub = crn.ub(:);
model.b = zeros(nm, 1);

n = size(model.S, 2);
a = sparse([speye(n) -model.ub; -speye(n) model.lb]);
Aeq = [model.S sparse(size(model.S,1),1); model.A(iden,:) 0];
beq = [model.b; 1];
options = optimset('linprog'); options.Display = 'off';

for branch = {'positive', 'negative'}
  if strcmp(branch{1}, 'positive')
    tlb = 0; tub = 999;
  else
    tlb = -999; tub = 0;
  end
  lb = [-ones(n,1)*1e9; tlb];
  ub = [ones(n,1)*1e9; tub];

  [x, ~, ef1] = linprog([-model.A(inum,:) 1], a, zeros(2*n,1), Aeq, beq, lb, ub, options);
  if ef1 == 1
    mx = (model.A(inum,:)*x(1:end-1)) / (model.A(iden,:)*x(1:end-1));
  else
    mx = NaN;
  end
  [x, ~, ef2] = linprog([model.A(inum,:) 1], a, zeros(2*n,1), Aeq, beq, lb, ub, options);
  if ef2 == 1
    mn = (model.A(inum,:)*x(1:end-1)) / (model.A(iden,:)*x(1:end-1));
  else
    mn = NaN;
  end
  printf('%-9s branch: exitflags %d/%d  max %s  min %s  round2 equal: %d\n', ...
         branch{1}, ef1, ef2, num2str(mx, '%.10g'), num2str(mn, '%.10g'), ...
         isequal(round(mx,2), round(mn,2)));
end

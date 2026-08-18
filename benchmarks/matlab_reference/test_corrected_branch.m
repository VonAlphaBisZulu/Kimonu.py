% Test a candidate correction to the reference's negative branch.
%
% The Charnes-Cooper transform sets t = 1/(denominator activity) and w = t*v. The box lb <= v <= ub
% becomes lb*t <= w <= ub*t only while t > 0; for t < 0 multiplying by t reverses the inequalities,
% so the correct constraint is ub*t <= w <= lb*t. The reference uses the t > 0 form in both
% branches:
%
%   a = sparse([eye(n) -model.ub; -eye(n) model.lb]);   %  w <= ub*t  and  w >= lb*t
%
% This script re-tests the pairs the reference missed, using the swapped form in the negative
% branch, and reports how many are then found concordant.
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
addpath(genpath(fullfile(upstream, 'Code')));

crn = jsondecode(fileread(crn_path));
nc = crn.n_complexes; nr = crn.n_reactions; nm = crn.n_metabolites;
model.A = triplets_to_sparse(crn.A, nc, nr);
model.S = triplets_to_sparse(crn.S, nm, nr);
model.lb = crn.lb(:); model.ub = crn.ub(:); model.b = zeros(nm,1); model.c = zeros(nr,1);
[B, group] = find_balanced_complexes(model);

missed = jsondecode(fileread(missed_path));
n = size(model.S, 2);
a_pos = sparse([speye(n) -model.ub; -speye(n) model.lb]);   % as in the reference
a_neg = sparse([speye(n) -model.lb; -speye(n) model.ub]);   % candidate correction for t < 0
options = optimset('linprog'); options.Display = 'off';

n_recovered = 0; n_still_missed = 0; n_neg_feasible = 0;
for k = 1:size(missed, 1)
  inum = missed(k, 1); iden = missed(k, 2);
  Aeq = [model.S sparse(nm, 1); model.A(iden, :) 0];
  beq = [model.b; 1];
  gd = group(iden);

  ok_p = false; ok_n = false; mx_p = NaN; mn_p = NaN; mx_n = NaN; mn_n = NaN;
  if gd ~= 'N'
    lb = [-ones(n,1)*1e9; 0]; ub = [ones(n,1)*1e9; 999];
    [x,~,e1] = linprog([-model.A(inum,:) 1], a_pos, zeros(2*n,1), Aeq, beq, lb, ub, options);
    [y,~,e2] = linprog([ model.A(inum,:) 1], a_pos, zeros(2*n,1), Aeq, beq, lb, ub, options);
    if e1 == 1 && e2 == 1
      mx_p = (model.A(inum,:)*x(1:end-1))/(model.A(iden,:)*x(1:end-1));
      mn_p = (model.A(inum,:)*y(1:end-1))/(model.A(iden,:)*y(1:end-1));
      ok_p = (round(mx_p,2) == round(mn_p,2));
    end
  end
  if gd ~= 'P'
    lb = [-ones(n,1)*1e9; -999]; ub = [ones(n,1)*1e9; 0];
    [x,~,e1] = linprog([-model.A(inum,:) 1], a_neg, zeros(2*n,1), Aeq, beq, lb, ub, options);
    [y,~,e2] = linprog([ model.A(inum,:) 1], a_neg, zeros(2*n,1), Aeq, beq, lb, ub, options);
    if e1 == 1 && e2 == 1
      n_neg_feasible = n_neg_feasible + 1;
      mx_n = (model.A(inum,:)*x(1:end-1))/(model.A(iden,:)*x(1:end-1));
      mn_n = (model.A(inum,:)*y(1:end-1))/(model.A(iden,:)*y(1:end-1));
      ok_n = (round(mx_n,2) == round(mn_n,2));
    end
  end

  if gd == 'P'
    found = ok_p;
  elseif gd == 'N'
    found = ok_n;
  else
    found = ok_p && ok_n && (round(mx_p,2) == round(mn_n,2));
  end
  if found
    n_recovered = n_recovered + 1;
  else
    n_still_missed = n_still_missed + 1;
    printf('still missed: %d/%d group %s/%s  p[%g %g] n[%g %g]\n', inum, iden, group(inum), gd, mx_p, mn_p, mx_n, mn_n);
  end
end
printf('\nwith the corrected negative branch: recovered %d, still missed %d (of %d)\n', ...
       n_recovered, n_still_missed, size(missed, 1));
printf('negative branch feasible in %d cases (was 0 before the correction)\n', n_neg_feasible);

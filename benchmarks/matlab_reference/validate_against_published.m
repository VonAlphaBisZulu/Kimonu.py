% Validate the Octave port against the reference's own published output.
%
% Loads one of the committed pre-split models, reruns find_balanced_complexes with the GLPK-backed
% linprog, and compares the resulting balanced set against the B stored in the committed
% Results/concordant/*.mat for the same model. Also reports the sign-class census, which bounds how
% much of the network is exposed to the negative-branch defect (see FINDINGS.md).
args = argv();
model_path = args{1};
published_path = args{2};

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat'));
upstream = getenv('UPSTREAM_ALGORITHM');
if isempty(upstream)
  upstream = fullfile(here, '..', '..', '..', 'Upstream_Algorithm');
end
addpath(genpath(fullfile(upstream, 'Code')));
warning('off', 'all');

function M = csc(s, nrows)
  jc = double(s.jc); ir = double(s.ir); dat = double(s.data);
  ncols = numel(jc) - 1;
  cols = zeros(numel(dat), 1);
  for j = 1:ncols
    cols(jc(j)+1 : jc(j+1)) = j;
  end
  M = sparse(ir + 1, cols, dat, nrows, ncols);
end

d = load(model_path);
fn = fieldnames(d);
mname = '';
for i = 1:numel(fn)
  if strncmp(fn{i}, 'model_elementary', 16), mname = fn{i}; end
end
m = d.(mname);
nr = numel(m.lb); nm = numel(m.b);
model.S = csc(m.S, nm);
model.A = csc(m.A, max(double(m.A.ir)) + 1);
model.Y = csc(m.Y, nm);
model.lb = m.lb(:); model.ub = m.ub(:); model.b = m.b(:); model.c = zeros(nr, 1);
nc = size(model.A, 1);
printf('model %s: %d complexes, %d reactions, %d metabolites\n', mname, nc, nr, nm);

t0 = tic;
[B, group] = find_balanced_complexes(model);
printf('rerun: %d balanced complexes in %.1f s\n', numel(B), toc(t0));
printf('sign classes: P %d, N %d, U %d, B %d\n', sum(group=='P'), sum(group=='N'), ...
       sum(group=='U'), sum(group=='B'));
printf('exposed denominators (N or U): %d of %d complexes (%.1f%%)\n', ...
       sum(group=='N') + sum(group=='U'), nc, ...
       100 * (sum(group=='N') + sum(group=='U')) / nc);

p = load(published_path);
if isfield(p, 'B')
  Bp = double(p.B(:));
  printf('published: %d balanced complexes\n', numel(Bp));
  printf('identical: %d  (only in rerun: %d, only in published: %d)\n', ...
         isequal(sort(double(B(:))), sort(Bp)), ...
         numel(setdiff(double(B(:)), Bp)), numel(setdiff(Bp, double(B(:)))));
else
  printf('published file has no B field; fields: %s\n', strjoin(fieldnames(p)', ', '));
end

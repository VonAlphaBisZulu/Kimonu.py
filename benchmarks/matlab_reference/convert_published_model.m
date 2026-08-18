% Convert a committed pre-split model from the reference repository into the shared crn.json
% format, and its committed result into the shared results format.
%
% The reference's models_with_elementary_steps/*.mat files are the elementary-step-split networks
% the published numbers were computed on, and they already carry A and Y. Converting them lets the
% other implementations be run on the paper's own inputs, with no preprocessing of our own in
% between and no need for the R stage that builds A and Y.
%
% Usage:
%   octave --no-gui convert_published_model.m <pre_balanced.mat> <out.crn.json> [<concordant.mat> <out.results.json>]

args = argv();
model_path = args{1};
crn_out = args{2};

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd(); end
addpath(fullfile(here, 'compat'));
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

function write_triplets(fid, name, M, last)
  [i, j, v] = find(M);
  fprintf(fid, '  "%s": [', name);
  parts = arrayfun(@(k) sprintf('[%d, %d, %.17g]', i(k), j(k), v(k)), 1:numel(v), ...
                   'UniformOutput', false);
  fprintf(fid, '%s', strjoin(parts, ', '));
  if last
    fprintf(fid, ']\n');
  else
    fprintf(fid, '],\n');
  end
end

d = load(model_path);
fn = fieldnames(d);
mname = '';
for i = 1:numel(fn)
  if strncmp(fn{i}, 'model_elementary', 16), mname = fn{i}; end
end
m = d.(mname);
nr = numel(m.lb); nm = numel(m.b);
S = csc(m.S, nm);
A = csc(m.A, max(double(m.A.ir)) + 1);
Y = csc(m.Y, nm);
nc = size(A, 1);
printf('%s: %d complexes, %d reactions, %d metabolites\n', mname, nc, nr, nm);

fid = fopen(crn_out, 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "n_complexes": %d,\n  "n_reactions": %d,\n  "n_metabolites": %d,\n', nc, nr, nm);
fprintf(fid, '  "reactions": [%s],\n', strjoin(arrayfun(@(k) sprintf('"r%d"', k), 1:nr, 'UniformOutput', false), ', '));
fprintf(fid, '  "metabolites": [%s],\n', strjoin(arrayfun(@(k) sprintf('"m%d"', k), 1:nm, 'UniformOutput', false), ', '));
fprintf(fid, '  "lb": [%s],\n', strjoin(arrayfun(@(k) sprintf('%.17g', m.lb(k)), 1:nr, 'UniformOutput', false), ', '));
fprintf(fid, '  "ub": [%s],\n', strjoin(arrayfun(@(k) sprintf('%.17g', m.ub(k)), 1:nr, 'UniformOutput', false), ', '));
write_triplets(fid, 'A', A, false);
write_triplets(fid, 'Y', Y, false);
write_triplets(fid, 'S', S, true);
fprintf(fid, '}\n');
fclose(fid);
printf('written %s\n', crn_out);

%% published result, if given
if numel(args) >= 4
  pub_path = args{3};
  res_out = args{4};
  p = load(pub_path);
  B = double(p.B(:));
  if isstruct(p.CC)
    CC = csc(p.CC, nc);
  else
    CC = sparse(p.CC);
  end
  printf('published: %d balanced, CC nnz %d\n', numel(B), nnz(CC));

  % Modules are the connected components of the concordance relation (values 2 and -1);
  % value 1 marks balanced complexes and is not a concordance edge.
  E = (CC == 2) | (CC == -1);
  E = E | E';
  module_of = zeros(nc, 1);
  next_module = 0;
  for i = 1:nc
    if module_of(i) ~= 0 || any(B == i), continue; end
    comp = i; frontier = i;
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
  printf('published modules (non-singleton): %d\n', next_module);

  crn.n_complexes = nc; crn.n_reactions = nr; crn.n_metabolites = nm;
  crn.metabolites = arrayfun(@(k) sprintf('m%d', k), 1:nm, 'UniformOutput', false);
  model.Y = Y;
  write_results_json(res_out, crn, model, module_of, 0, 0, 0, 0);
  printf('written %s\n', res_out);
end

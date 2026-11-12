function write_results_json(path, crn, model, module_of, seconds, t_balanced, t_concordant, n_pairs, mode, pairs_tested, completed, n_lp_fail, n_lp_nonopt)
% Write the reference result in the same shape as the Kimonu and COCOA runners.
%
% Complexes are keyed by composition rather than by name, matching the other two exporters, so the
% comparison never has to trust three separate naming schemes.

  fid = fopen(path, 'w');
  fprintf(fid, '{\n');
  fprintf(fid, '  "tool": "Upstream_Algorithm",\n');
  if nargin < 12, n_lp_fail = 0; end
  if nargin < 13, n_lp_nonopt = 0; end
  fprintf(fid, '  "seconds": %.6f,\n', seconds);
  % A non-optimal LP is mapped to nan by the reference and silently drops the pair or complex, so
  % the counts travel with the result. lp_nonoptimal includes the infeasible negative branch the
  % formulation produces on purpose; lp_degraded counts CPX_STAT_OPTIMAL_INFEAS, an optimal basis
  % carrying small infeasibilities after unscaling, which the 1e9 bounds in the pair LP provoke.
  fprintf(fid, '  "lp_nonoptimal": %d,\n', n_lp_nonopt);
  fprintf(fid, '  "lp_degraded": %d,\n', n_lp_fail);
  fprintf(fid, '  "seconds_balanced": %.6f,\n', t_balanced);
  fprintf(fid, '  "seconds_concordant": %.6f,\n', t_concordant);
  fprintf(fid, '  "candidate_pairs": %d,\n', n_pairs);
  if nargin >= 9
    fprintf(fid, '  "mode": "%s",\n', mode);
    fprintf(fid, '  "pairs_tested": %d,\n', pairs_tested);
    fprintf(fid, '  "completed": %s,\n', merge(completed, 'true', 'false'));
  end
  fprintf(fid, '  "n_reactions": %d,\n', crn.n_reactions);
  fprintf(fid, '  "n_metabolites": %d,\n', crn.n_metabolites);
  fprintf(fid, '  "complexes": [\n');

  nc = crn.n_complexes;
  for j = 1:nc
    rows = find(model.Y(:, j) ~= 0);
    parts = {};
    for k = 1:numel(rows)
      parts{end+1} = sprintf('"%s": %s', crn.metabolites{rows(k)}, ...
                             num2str(full(model.Y(rows(k), j)), '%.17g'));
    end
    sep = ',';
    if j == nc
      sep = '';
    end
    fprintf(fid, '    {"id": "complex_%d", "module": %d, "composition": {%s}}%s\n', ...
            j, module_of(j), strjoin(parts, ', '), sep);
  end
  fprintf(fid, '  ]\n}\n');
  fclose(fid);
end

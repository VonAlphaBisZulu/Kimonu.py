function write_results_json(path, crn, model, module_of, seconds, t_balanced, t_concordant, n_pairs)
% Write the reference result in the same shape as the Kimonu and COCOA runners.
%
% Complexes are keyed by composition rather than by name, matching the other two exporters, so the
% comparison never has to trust three separate naming schemes.

  fid = fopen(path, 'w');
  fprintf(fid, '{\n');
  fprintf(fid, '  "tool": "Upstream_Algorithm",\n');
  fprintf(fid, '  "seconds": %.6f,\n', seconds);
  fprintf(fid, '  "seconds_balanced": %.6f,\n', t_balanced);
  fprintf(fid, '  "seconds_concordant": %.6f,\n', t_concordant);
  fprintf(fid, '  "candidate_pairs": %d,\n', n_pairs);
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

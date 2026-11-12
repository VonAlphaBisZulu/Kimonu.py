function M = triplets_to_sparse(t, nrows, ncols)
% Build a sparse matrix from an n-by-3 array of [row, col, value] triplets (1-based).
  if isempty(t)
    M = sparse(nrows, ncols);
    return
  end
  M = sparse(t(:, 1), t(:, 2), t(:, 3), nrows, ncols);
end

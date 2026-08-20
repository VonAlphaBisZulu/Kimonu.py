function Z = null (A, opt)
% null with MATLAB's rational option, which Octave does not implement.
%
% F2C2 calls null(A, 'r') for a rational basis of the kernel; Octave's null only offers the
% singular-value form. The rational branch is read off the reduced row echelon form, giving the
% same basis MATLAB returns: one column per free variable, carrying a 1 in that variable's row.

  if nargin > 1 && ischar (opt) && numel (opt) > 0 && lower (opt(1)) == 'r'
    [R, pivots] = rref (A);
    n = columns (A);
    free = setdiff (1:n, pivots);
    Z = zeros (n, numel (free));
    for k = 1:numel (free)
      Z(free(k), k) = 1;
      Z(pivots, k) = -R(1:numel (pivots), free(k));
    end
    return
  end

  % default: singular-value form, matching Octave's own null
  if isempty (A)
    Z = [];
    return
  end
  [~, S, V] = svd (A, 0);
  if isempty (S)
    Z = V;
    return
  end
  s = diag (S);
  tol = max (size (A)) * s(1) * eps;
  rank_A = sum (s > tol);
  Z = V(:, rank_A+1:end);
end

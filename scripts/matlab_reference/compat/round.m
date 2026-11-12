function y = round(x, n)
% round with an optional digit count, as MATLAB's round has and Octave's built-in does not.
%
% The reference code calls round(x, 2) when deciding whether the maximum and minimum activity
% ratios of a complex pair coincide, so the two-argument form is load-bearing: without it the
% concordance test cannot run at all.

  if nargin < 2
    y = builtin('round', x);
  else
    scale = 10 ^ n;
    y = builtin('round', x * scale) / scale;
  end
end

function [x, fval, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options)
% Drop-in replacement for MATLAB's linprog, backed by Octave's bundled GLPK.
%
% The reference implementation calls linprog from the Optimization Toolbox. Providing the same
% interface on top of glpk lets the original .m files run unmodified and without a MATLAB licence.
% Only the calling convention used by the reference code is supported:
%
%   [x, fval, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options)
%
% minimising f'x subject to A*x <= b, Aeq*x = beq, lb <= x <= ub. exitflag is 1 on success, so
% that the reference code's `ExitFlag == 1` checks keep their meaning.

  % optimset(fun) queries a solver for its defaults by calling fun('defaults').
  if nargin == 1 && ischar(f) && strcmp(f, 'defaults')
    x = struct('Display', 'off', 'MaxIter', 1e5, 'TolFun', 1e-9);
    return
  end

  if nargin < 8, options = []; end
  if nargin < 7, ub = []; end
  if nargin < 6, lb = []; end
  if nargin < 5, beq = []; end
  if nargin < 4, Aeq = []; end
  if nargin < 3, b = []; end
  if nargin < 2, A = []; end

  f = full(f(:));
  n = numel(f);

  if isempty(lb), lb = -Inf(n, 1); else lb = full(lb(:)); end
  if isempty(ub), ub =  Inf(n, 1); else ub = full(ub(:)); end

  M = [];
  rhs = [];
  ctype = '';
  if ~isempty(A)
    M = [M; A];
    rhs = [rhs; full(b(:))];
    ctype = [ctype, repmat('U', 1, size(A, 1))];
  end
  if ~isempty(Aeq)
    M = [M; Aeq];
    rhs = [rhs; full(beq(:))];
    ctype = [ctype, repmat('S', 1, size(Aeq, 1))];
  end
  if isempty(M)
    % glpk needs at least one row; a vacuous one leaves the problem unchanged.
    M = sparse(1, n);
    rhs = 0;
    ctype = 'U';
  end

  vartype = repmat('C', 1, n);
  param.msglev = 0;

  [x, fval, errnum] = glpk(f, M, rhs, lb, ub, ctype, vartype, 1, param);

  % glpk returns 0 for "solved to optimality"; everything else is a failure as far as the
  % reference code is concerned, and it only ever tests for 1.
  if errnum == 0
    exitflag = 1;
  else
    exitflag = -errnum;
  end
end

function [x, fval, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options)
% Drop-in replacement for MATLAB's linprog, backed by a persistent CPLEX worker.
%
% Same calling convention the reference uses:
%   [x, fval, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, options)
% minimising f'x subject to A*x <= b, Aeq*x = beq, lb <= x <= ub; exitflag is 1 on success so the
% reference code's `== 1` checks keep their meaning.
%
% Consecutive calls in the reference differ only in the objective or in a single equality row, so
% the problem is held open in the worker and only the delta is sent. Several models are kept in
% slots because the pair walk alternates between two inequality blocks and two objective signs.

  persistent S

  if nargin == 1 && ischar(f) && strcmp(f, 'failures')
    if isempty(S), x = 0; else x = S.nfail; end
    return
  end

  if nargin == 1 && ischar(f) && strcmp(f, 'nonoptimal')
    if isempty(S), x = 0; else x = S.ninf; end
    return
  end

  if nargin == 1 && ischar(f) && strcmp(f, 'stats')
    if isempty(S)
      printf('no LP calls\n');
    else
      printf('modes NEW %d OBJ %d ROW %d OR %d | t_select %.1f s t_msg %.1f s t_io %.1f s\n', ...
             S.cnt.NEW, S.cnt.OBJ, S.cnt.ROW, S.cnt.OR, S.t_sel, S.t_msg, S.t_io);
      ks = S.sthist.keys();
      for ii = 1:numel(ks)
        printf('  status %d: %d\n', ks{ii}, S.sthist(ks{ii}));
      end
    end
    return
  end

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
  if isempty(b),   b   = zeros(0, 1); else b   = full(b(:));   end
  if isempty(beq), beq = zeros(0, 1); else beq = full(beq(:)); end
  mA = size(A, 1); mE = size(Aeq, 1);

  if isempty(S)
    S = struct('in', [], 'out', [], 'pid', [], 'slots', {{}}, 'tick', 0);
    here = fileparts(mfilename('fullpath'));
    srv = fullfile(here, 'lp_server.py');
    py = getenv('LP_SERVER_PYTHON');
    if isempty(py), py = 'python3'; end
    [S.in, S.out, S.pid] = popen2(py, {'-u', srv}, true);   % blocking pipes
    if S.pid < 0, error('linprog: could not start LP worker'); end
    S.slots = cell(1, 5);
    S.rv = [];
    S.cnt = struct('NEW',0,'OBJ',0,'ROW',0,'OR',0);
    S.nfail = 0; S.ninf = 0; S.sthist = containers.Map('KeyType','double','ValueType','double');
    S.t_sel=0; S.t_io=0; S.t_msg=0;
  end
  S.tick = S.tick + 1;
  t_sel0 = tic;
  % Row fingerprints: comparing Aeq against every slot with `~=` allocates a sparse logical each
  % time, which dominated selection. Two independent projections identify the changed rows in O(m).
  if isempty(S.rv) || size(S.rv, 1) ~= n
    rand('state', 42); S.rv = rand(n, 2);
  end
  if mE > 0, rh = Aeq * S.rv; else rh = zeros(0, 2); end
  if mA > 0, hA = [nnz(A), full(sum(nonzeros(A)))]; else hA = [0, 0]; end

  % LP_BRIDGE_NOCACHE forces a full rebuild every call, so a long cached run can be checked
  % against a stateless one over the same pairs.
  nocache = ~isempty(getenv('LP_BRIDGE_NOCACHE'));

  % Prefer a slot that needs only the one changed row; failing that claim a free slot, so the
  % four (inequality block x objective sign) combinations the pair walk alternates between each
  % keep their own basis instead of overwriting one another.
  best = 0; mode = 'NEW'; chrow = [];
  freeslot = 0; cand_obj = 0; cand_or = 0; cand_or_row = []; stale = 0;
  for k = 1:numel(S.slots)
    if nocache, break, end
    sl = S.slots{k};
    if isempty(sl)
      if freeslot == 0, freeslot = k; end
      continue
    end
    if sl.n ~= n || sl.mA ~= mA || sl.mE ~= mE || ~isequal(sl.hA, hA) || ...
       ~isequal(sl.lb, lb) || ~isequal(sl.ub, ub) || ~isequal(sl.b, b)
      % Reclaim only slots left over from an earlier stage. The pair walk cycles through its own
      % variants every pair, so a recently touched slot is still live and must not be evicted.
      if stale == 0 && (S.tick - sl.tick) > 16, stale = k; end
      continue
    end
    same_f = isequal(sl.f, f);
    if isequal(sl.rh, rh) && isequal(sl.beq, beq)
      if same_f, best = k; mode = 'OBJ'; break, end
      if cand_obj == 0, cand_obj = k; end
      continue
    end
    ch = find(any(rh ~= sl.rh, 2));
    bch = find(beq ~= sl.beq);
    if numel(ch) == 1 && (isempty(bch) || isequal(bch(:), ch(:)))
      if same_f
        best = k; mode = 'ROW'; chrow = ch; break
      elseif cand_or == 0
        cand_or = k; cand_or_row = ch;
      end
    end
  end

  if best == 0
    if freeslot > 0
      best = freeslot; mode = 'NEW';
    elseif stale > 0
      % Claim a slot holding an unrelated problem before mutating a matching one, so the four
      % (inequality block x objective sign) variants the pair walk cycles through each get a slot.
      best = stale; mode = 'NEW';
    elseif cand_obj > 0
      best = cand_obj; mode = 'OBJ';
    elseif cand_or > 0
      best = cand_or; mode = 'OR'; chrow = cand_or_row;
    else
      [~, best] = min(cellfun(@(s) ifelse_tick(s), S.slots));
      mode = 'NEW';
    end
  end

  S.t_sel = S.t_sel + toc(t_sel0);
  t_msg0 = tic;
  if ~strcmp(mode, 'NEW') && ~isequal(S.slots{best}.A, A)
    mode = 'NEW';
  end
  sid = sprintf('s%d', best);

  switch mode
    case 'NEW'
      d.n = n;
      d.lb = clampinf(lb)'; d.ub = clampinf(ub)';
      d.m_A = mA; d.b = b'; d.m_Aeq = mE; d.beq = beq';
      d.A = trip(A); d.Aeq = trip(Aeq);
      fi = find(f); d.obj = [fi(:) f(fi)(:)];
      msg = sprintf('NEW %s %s\n', sid, jsonencode(d));
    case 'OBJ'
      msg = sprintf('OBJ %s %s\n', sid, sparsetok(find(f), f(find(f))));
    case 'ROW'
      msg = sprintf('ROW %s %d %.17g %s\n', sid, chrow - 1, beq(chrow), rowtok(Aeq, chrow));
    case 'OR'
      msg = sprintf('OR %s %d %.17g %s | %s\n', sid, chrow - 1, beq(chrow), ...
                    rowtok(Aeq, chrow), sparsetok(find(f), f(find(f))));
  end

  S.t_msg = S.t_msg + toc(t_msg0);
  S.cnt.(mode) = S.cnt.(mode) + 1;
  t_io0 = tic;
  fputs(S.in, msg); fflush(S.in);
  hdr = fread(S.out, 2, 'double');
  if numel(hdr) < 2, error('linprog: LP worker died'); end
  x = fread(S.out, n, 'double');
  S.t_io = S.t_io + toc(t_io0);
  fval = hdr(2);
  st = hdr(1);
  if st == 1 || st == 101 || st == 102
    exitflag = 1;
  else
    exitflag = 0;
    % 2/3/4 are unbounded / infeasible / either: outcomes the reference's own formulation produces
    % on purpose. Only anything else is a solver failure worth reporting.
    if ~(st == 2 || st == 3 || st == 4)
      S.nfail = S.nfail + 1;
    end
    S.ninf = S.ninf + 1;
    if S.sthist.isKey(st), S.sthist(st) = S.sthist(st) + 1; else S.sthist(st) = 1; end
  end

  S.slots{best} = struct('n', n, 'mA', mA, 'mE', mE, 'A', A, 'b', b, 'rh', rh, 'hA', hA, ...
                         'beq', beq, 'lb', lb, 'ub', ub, 'f', f, 'tick', S.tick);
end

function varargout = linprog_stats()
end

function v = clampinf(v)
  v(v == Inf) = 1e30; v(v == -Inf) = -1e30;
end

function t = trip(M)
  if isempty(M), t = zeros(0, 3); return, end
  [i, j, v] = find(M);
  t = [i(:) j(:) v(:)];
end

function s = sparsetok(idx, val)
  s = sprintf('%d %s', numel(idx), sprintf('%d %.17g ', [idx(:)'; val(:)']));
end

function s = rowtok(Aeq, r)
  [~, j, v] = find(Aeq(r, :));
  s = sparsetok(j, v);
end

function t = ifelse_tick(s)
  if isempty(s), t = -1; else t = s.tick; end
end

"""Persistent CPLEX worker for the Octave linprog shim.

The reference calls linprog in tight loops where consecutive problems differ only in the objective
(find_balanced_complexes) or in one equality row (find_concordant_complexes). Octave's glpk
rebuilds and re-solves from scratch each call and throws the basis away. Keeping the model alive
lets each re-solve start from the previous optimal basis.

The pair walk alternates between two inequality blocks and flips the objective sign every call, so
several models are kept in numbered slots and the shim picks the one that needs the smallest edit.

Protocol: one text request line per LP, one binary response.
  NEW <slot> <json>                 build/replace a slot
  OBJ <slot> <k> i v ...            objective only (sparse, 1-based)
  ROW <slot> <r> <rhs> <k> i v ...  one equality row only (0-based row, 1-based cols)
  OR  <slot> <r> <rhs> <k> i v ... | <k2> i v ...   row and objective
Response: status(f64) fval(f64) x(n * f64), little-endian.
"""
import sys, json, struct, numpy as np, cplex

BIG = 1e30
INF = cplex.infinity


def _arr(v):
    """Octave's jsonencode collapses 1-element arrays to scalars and 1-row matrices to flat lists."""
    if v is None:
        return []
    if isinstance(v, (int, float)):
        return [v]
    return v


def _trip(v):
    v = _arr(v)
    if v and isinstance(v[0], (int, float)):
        return [v]
    return v


def _fin(v):
    return [INF if x >= BIG else (-INF if x <= -BIG else float(x)) for x in v]


class Slot:
    def __init__(self, d):
        c = cplex.Cplex()
        for s in (c.set_log_stream, c.set_results_stream, c.set_warning_stream, c.set_error_stream):
            s(None)
        c.parameters.threads.set(1)
        # The balanced test classifies at |value| < 5e-10, below default solver tolerances, and a
        # warm start from a neighbouring basis can otherwise stop ~1e-4 short of the true optimum.
        c.parameters.simplex.tolerances.optimality.set(1e-9)
        c.parameters.simplex.tolerances.feasibility.set(1e-9)
        n = int(d["n"])
        c.variables.add(obj=[0.0] * n, lb=_fin(_arr(d["lb"])), ub=_fin(_arr(d["ub"])))
        rows, senses, rhs = [], "", []
        for blk, sense, rhskey in (("A", "L", "b"), ("Aeq", "E", "beq")):
            m = int(d.get("m_" + blk, 0) or 0)
            if m == 0:
                continue
            byrow = [[[], []] for _ in range(m)]
            for t in _trip(d.get(blk)):
                i, j, v = int(t[0]) - 1, int(t[1]) - 1, float(t[2])
                byrow[i][0].append(j)
                byrow[i][1].append(v)
            rows += [cplex.SparsePair(idx, val) for idx, val in byrow]
            senses += sense * m
            rhs += [float(x) for x in _arr(d[rhskey])]
        if rows:
            c.linear_constraints.add(lin_expr=rows, senses=senses, rhs=rhs)
        c.objective.set_sense(c.objective.sense.minimize)
        self.c, self.n, self.prev = c, n, []
        self.set_obj(_trip(d.get("obj")))
        self.eq0 = int(d.get("m_A", 0) or 0)   # equality rows start after the inequality block
        last = int(d.get("m_Aeq", 0) or 0) - 1
        if last >= 0:
            self.prev_row = [int(t[1]) - 1 for t in _trip(d.get("Aeq")) if int(t[0]) - 1 == last]

    def set_obj(self, pairs):
        cur = [(int(i) - 1, float(v)) for i, v in pairs]
        new = {k for k, _ in cur}
        self.c.objective.set_linear([(k, 0.0) for k, _ in self.prev if k not in new] + cur)
        self.prev = cur

    def set_row(self, r, rhs, pairs):
        # set_linear_components writes the given coefficients but leaves the rest of the row
        # untouched, so the previous row's entries have to be zeroed explicitly or they persist.
        rr = self.eq0 + r
        idx = [int(i) - 1 for i, _ in pairs]
        val = [float(v) for _, v in pairs]
        new = set(idx)
        stale = [k for k in self.prev_row if k not in new]
        if stale:
            self.c.linear_constraints.set_linear_components(rr, [stale, [0.0] * len(stale)])
        self.c.linear_constraints.set_linear_components(rr, [idx, val])
        self.c.linear_constraints.set_rhs(rr, float(rhs))
        self.prev_row = idx

    def solve(self):
        # The raw CPLEX status travels back so the shim can tell an infeasible problem -- which the
        # reference creates deliberately for pairs that only normalise with t < 0 -- from a solver
        # failure. Both give exitflag ~= 1, as they did under glpk.
        try:
            self.c.solve()
            st = self.c.solution.get_status()
            if st in (1, 101, 102):
                return float(st), self.c.solution.get_objective_value(), self.c.solution.get_values()
            return float(st), 0.0, [0.0] * self.n
        except cplex.exceptions.errors.CplexSolverError as e:
            return float(getattr(e, "args", [0])[0] or 0) * 0 - 1.0, 0.0, [0.0] * self.n


def pairs_from(tok):
    k = int(tok[0])
    return [(tok[1 + 2 * i], tok[2 + 2 * i]) for i in range(k)], tok[1 + 2 * k:]


def main():
    slots, out = {}, sys.stdout.buffer
    rd = sys.stdin.readline   # not `for line in sys.stdin`: its read-ahead stalls the round trip
    while True:
        line = rd()
        if not line:
            return
        line = line.rstrip("\n")
        if not line:
            continue
        op, _, rest = line.partition(" ")
        if op == "BYE":
            return
        sid, _, rest = rest.partition(" ")
        if op == "NEW":
            slots[sid] = Slot(json.loads(rest))
        elif op == "OBJ":
            p, _ = pairs_from(rest.split())
            slots[sid].set_obj(p)
        elif op in ("ROW", "OR"):
            t = rest.split()
            r, rhs = int(t[0]), float(t[1])
            p, tail = pairs_from(t[2:])
            slots[sid].set_row(r, rhs, p)
            if op == "OR":
                assert tail and tail[0] == "|"
                po, _ = pairs_from(tail[1:])
                slots[sid].set_obj(po)
        else:
            continue
        st, fv, x = slots[sid].solve()
        out.write(struct.pack("<2d", st, fv))
        out.write(np.asarray(x, dtype="<f8").tobytes())
        out.flush()


if __name__ == "__main__":
    main()

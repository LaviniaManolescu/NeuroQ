"""Algoritmi clasici de fazare (arbore de acoperire cu paritate, căutare locală, recoacere simulată), adaptați pentru NeuroQ."""
import itertools
import numpy as np


def cost(x, a, d):
    """Vectorised pairwise cost (identical value to the notebook's pairwise_cost)."""
    x = np.asarray(x)
    diff = x[:, None] != x[None, :]
    return float(np.triu(np.where(diff, a, d), 1).sum())


def canon(x):
    """Complement-symmetry: return the copy that starts with 0 (SNP 1 first)."""
    x = np.asarray(x, dtype=int)
    return x if x[0] == 0 else 1 - x


# ---------------------------------------------------------------- 1. greedy
def greedy_tree(a, d):
    """Greedy maximum-spanning-tree phasing (Kruskal + parity union-find).

    Edges are taken in decreasing net strength |A-D| (ties: more votes first).
    An edge that links two different components is accepted and fixes their
    relative phase (same if A>D, opposite if A<D); an edge inside a component
    is a contradiction and is skipped. Exact when the evidence is consistent.
    Returns (x, n_components). Components are phased independently: relative
    phase between them is NOT determined by the data.
    """
    n = len(a)
    parent = list(range(n))
    par = [0] * n  # parity of i relative to parent[i]

    def find(i):
        if parent[i] == i:
            return i, 0
        r, p = find(parent[i])
        parent[i], par[i] = r, par[i] ^ p
        return r, par[i]

    edges = [(abs(a[i, j] - d[i, j]), a[i, j] + d[i, j], i, j)
             for i, j in itertools.combinations(range(n), 2) if a[i, j] + d[i, j] > 0]
    edges.sort(key=lambda e: (-e[0], -e[1]))
    for _, _, i, j in edges:
        (ri, pi), (rj, pj) = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            par[ri] = pi ^ pj ^ (0 if a[i, j] >= d[i, j] else 1)
    roots = [find(i) for i in range(n)]
    x = np.array([p for _, p in roots])
    return canon(x), len({r for r, _ in roots})


# ----------------------------------------------------------- 2. local search
def _weights(a, d):
    w = (a - d) + (a - d).T          # symmetric, zero diagonal
    return w


def local_search(x, a, d):
    """Steepest-descent 1-flip. Flipping i changes the cost by z_i * (W z)_i."""
    w = _weights(a, d)
    x = np.array(x, dtype=int)
    while True:
        z = 1 - 2 * x
        delta = z * (w @ z)
        i = int(np.argmin(delta))
        if delta[i] >= -1e-12:
            return canon(x)
        x[i] ^= 1


# ----------------------------------------------------- 3. simulated annealing
def anneal(a, d, steps=4000, t0=None, t1=0.05, seed=0, x0=None):
    """Single-spin-flip simulated annealing (geometric cooling). Returns best x seen."""
    rng = np.random.default_rng(seed)
    n = len(a)
    w = _weights(a, d)
    x = rng.integers(0, 2, n) if x0 is None else np.array(x0, dtype=int)
    z = 1 - 2 * x
    t0 = t0 if t0 is not None else max(1.0, np.abs(w).sum(axis=1).max() / 2)
    cur = cost(x, a, d)
    best, best_x = cur, x.copy()
    for s in range(steps):
        t = t0 * (t1 / t0) ** (s / max(1, steps - 1))
        i = rng.integers(n)
        delta = z[i] * (w[i] @ z)
        if delta <= 0 or rng.random() < np.exp(-delta / t):
            x[i] ^= 1
            z[i] = -z[i]
            cur += delta
            if cur < best - 1e-12:
                best, best_x = cur, x.copy()
    return canon(best_x)


# ------------------------------------------------- 4. exact (vectorised, n<=~24)
def exact_vectorized(a, d):
    """Brute force with numpy; fixes x_0 = 0 (complement symmetry), so 2**(n-1) states."""
    n = len(a)
    k = np.arange(2 ** (n - 1), dtype=np.int64)
    bits = [np.zeros_like(k)] + [(k >> i) & 1 for i in range(n - 1)]
    c = np.zeros(len(k))
    for i, j in itertools.combinations(range(n), 2):
        c += np.where(bits[i] != bits[j], a[i, j], d[i, j])
    best = int(np.argmin(c))
    return np.array([b[best] for b in bits]), float(c[best])


# ------------------------------------------------------------------ driver
def solve_classical(a, d, restarts=20, steps=4000, seed=0):
    """Greedy tree -> local search, plus SA restarts. Returns (x, cost, info)."""
    xg, ncomp = greedy_tree(a, d)
    cands = {'greedy': xg, 'greedy+LS': local_search(xg, a, d)}
    sa = [local_search(anneal(a, d, steps, seed=seed + r), a, d) for r in range(restarts)]
    cands['SA+LS'] = min(sa, key=lambda x: cost(x, a, d))
    scores = {k: cost(v, a, d) for k, v in cands.items()}
    name = min(scores, key=scores.get)
    return cands[name], scores[name], {'components': ncomp, 'scores': scores, 'winner': name}


# --------------------------------------------- optional: full-read MEC score
def mec_cost(x, reads):
    """Minimum Error Correction: per read, calls to flip to fit haplotype x or 1-x."""
    x = np.asarray(x)
    tot = 0
    for r in reads:
        m = r >= 0
        e = int(np.sum(r[m] != x[m]))
        tot += min(e, int(m.sum()) - e)
    return tot

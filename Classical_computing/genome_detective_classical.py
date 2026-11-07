import itertools
import numpy as np

# initial reads (rows = reads, columns = SNPs), -1 = missing
READS = np.array([
    [-1, -1, -1,  1,  0, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1,  0,  1],
    [-1, -1, -1, -1, -1, -1,  0,  1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1,  1,  0, -1],
    [-1, -1, -1, -1, -1, -1, -1,  1,  0,  1],
    [-1, -1, -1, -1,  1,  1, -1, -1, -1, -1],
    [ 0,  0,  1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  0,  0, -1, -1, -1, -1],
    [-1, -1,  1,  1,  0, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  1,  0,  1, -1],
    [-1, -1,  1,  1,  0, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  0,  0,  1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  1,  0, -1, -1],
    [-1, -1, -1, -1, -1, -1,  0,  1, -1, -1],
    [ 1,  1,  0, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  0,  1,  0, -1],
    [-1, -1, -1, -1,  0,  0,  1, -1, -1, -1],
    [-1,  1,  0, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1,  0,  1,  1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1,  0,  1,  0],
    [-1, -1, -1, -1, -1,  1,  0, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  0,  1, -1, -1],
    [-1, -1, -1, -1,  0,  0,  1, -1, -1, -1],
    [-1, -1, -1, -1, -1,  0,  1, -1, -1, -1]
])

N = READS.shape[1]      # no of SNPs
text = lambda x: ''.join(map(str, x))
canon = lambda x: np.asarray(x, int) if x[0] == 0 else 1 - np.asarray(x, int)

# Task A: build evidence (A, D) from reads
def build_evidence(reads):
    n = reads.shape[1]
    a = np.zeros((n, n)); d = np.zeros((n, n))
    for row in reads:
        for i, j in itertools.combinations(np.flatnonzero(row >= 0), 2):
            if row[i] == row[j]: a[i, j] += 1      # same
            else:                d[i, j] += 1      # opposite
    return a, d

# Task B: cost of a labeling x given evidence (A, D)
def pairwise_cost(x, a, d):
    x = np.asarray(x)
    return float(np.triu(np.where(x[:, None] != x[None, :], a, d), 1).sum())

# Task C: exact solution (brute-force, small n only)
def exact_reference(a, d):
    states = np.array(list(itertools.product([0, 1], repeat=len(a))))
    return states, np.array([pairwise_cost(x, a, d) for x in states])

def exact_vectorized(a, d):
    """x_0 = 0 fixed (complement symmetry): 2**(n-1) states, fine up to n~22."""
    n = len(a)
    k = np.arange(2 ** (n - 1), dtype=np.int64)
    bits = [np.zeros_like(k)] + [(k >> i) & 1 for i in range(n - 1)]
    c = np.zeros(len(k))
    for i, j in itertools.combinations(range(n), 2):
        c += np.where(bits[i] != bits[j], a[i, j], d[i, j])
    b = int(np.argmin(c))
    return np.array([v[b] for v in bits]), float(c[b])

# Task D - Ising check
def ising_from_evidence(a, d):
    return (d - a) / 2, float((a + d).sum() / 2)

def ising_energy(x, j, offset):
    z = 1 - 2 * np.asarray(x)
    return float(offset + sum(j[i, k] * z[i] * z[k]
                              for i, k in itertools.combinations(range(len(x)), 2)))

# Greedy spanning tree (Kruskal)
def greedy_tree(a, d):
    n = len(a)
    parent = list(range(n)); par = [0] * n

    def find(i):
        if parent[i] == i: return i, 0
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
    return canon(np.array([p for _, p in roots])), len({r for r, _ in roots})


def solve_classical(a, d, restarts=20, steps=4000, seed=0):
    xg, ncomp = greedy_tree(a, d)
    cands = {'greedy': xg}
    scores = {k: pairwise_cost(v, a, d) for k, v in cands.items()}
    win = min(scores, key=scores.get)
    return cands[win], scores[win], cands, scores, ncomp


if __name__ == '__main__':
    # Build check
    tiny = np.array([[0, 1, 1], [1, 0, 0], [0, 0, -1]])
    ta, td = build_evidence(tiny)
    assert pairwise_cost([0, 1, 1], ta, td) == 1 and pairwise_cost([0, 0, 1], ta, td) == 4

    A, D = build_evidence(READS)
    print('SAME votes A:\n', A.astype(int), '\n\nOPPOSITE votes D:\n', D.astype(int), '\n\n')

    states, costs = exact_reference(A, D)
    OPT = costs.min()
    print('\nExact minimum:', OPT, '| optimal labels:',
          [text(s) for s, c in zip(states, costs) if c == OPT])
    print('Uniform mean cost:', costs.mean())

    J, off = ising_from_evidence(A, D)
    assert all(np.isclose(ising_energy(x, J, off), c) for x, c in zip(states, costs))
    print('Ising energy == graph cost on all', 2 ** N, 'states: OK')

    x, c, cands, scores, ncomp = solve_classical(A, D)
    print('\n--- classical solvers ---')
    for k, v in cands.items():
        print(f'{k:10s} {text(v)}  cost={scores[k]}')
    print('connected components:', ncomp)
    print(f'\nRESULT: haplotype {text(x)} / {text(1 - x)}   cost={c}')
    assert np.isclose(c, OPT)
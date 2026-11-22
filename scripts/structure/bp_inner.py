"""D3: belief propagation on the inner test written as base-5 column additions with carries.
Variables: bits x_i (hex digit 6a + 9b per position) and carries c_j. Factor j (5-adic column j):
    sum_i w_ij x_i + c_j = t_j + 5 c_{j+1},   w_ij = j-th base-5 digit of w_i.
Messages to bits through each column are exact for the factor (convolution of the other bits' messages). Planted instances:
x* random, t = sum w_i x*_i mod 5^K. Measures (1) how many bits the converged marginals predict correctly, (2) whether BP-guided
decimation (fix the most biased bit, rerun) reaches a solution, compared with random guessing (50%)."""
import sys, random, math
import numpy as np
def digits5(v, K): return [(v // 5 ** j) % 5 for j in range(K)]
def run(w, t, K, iters=60, damp=0.5, fixed=None):
    n = len(w); W = np.array([digits5(x, K) for x in w])           # n x K
    T = digits5(t, K); C = 4 * n // 4 + 2                         # carry range bound
    Cmax = (4 * n) // 4 + 1
    # messages factor j -> bit i: prob(x_i = 1); bit -> factor: product of other factors
    mfb = np.full((K, n), 0.5); fixed = fixed or {}
    carry_msg = [np.ones(Cmax + 1) / (Cmax + 1) for _ in range(K + 1)]; carry_msg[0] = np.eye(Cmax + 1)[0]
    for it in range(iters):
        # bit beliefs
        lo = np.log(np.clip(mfb, 1e-12, 1)); l0 = np.log(np.clip(1 - mfb, 1e-12, 1))
        tot1 = lo.sum(0); tot0 = l0.sum(0)
        new = mfb.copy()
        for j in range(K):
            # bit -> factor j message: exclude factor j
            b1 = tot1 - lo[j]; b0 = tot0 - l0[j]
            p1 = 1 / (1 + np.exp(np.clip(b0 - b1, -50, 50)))
            for i, v in fixed.items(): p1[i] = float(v)
            # distribution of S = sum_i w_ij x_i (convolution), plus leave-one-out via forward/backward
            maxs = 4 * n + 1
            fwd = [np.zeros(maxs + 1)]; fwd[0][0] = 1
            for i in range(n):
                d = np.zeros(maxs + 1); a = W[i, j]
                d += fwd[-1] * (1 - p1[i]); d[a:] += fwd[-1][:maxs + 1 - a] * p1[i]; fwd.append(d)
            bwd = [None] * (n + 1); bwd[n] = np.zeros(maxs + 1); bwd[n][0] = 1
            for i in range(n - 1, -1, -1):
                d = np.zeros(maxs + 1); a = W[i, j]
                d += bwd[i + 1] * (1 - p1[i]); d[a:] += bwd[i + 1][:maxs + 1 - a] * p1[i]; bwd[i] = d
            # constraint: S + c_j = T_j + 5 c_{j+1}; carry messages give distributions of c_j (from below) and c_{j+1} (free above)
            cin = carry_msg[j]
            # admissible totals: value v = S + c_in with (v - T_j) % 5 == 0 and (v - T_j)//5 <= Cmax
            ok = np.zeros(2 * maxs + Cmax + 2)
            for v in range(ok.size):
                if (v - T[j]) % 5 == 0 and 0 <= (v - T[j]) // 5 <= Cmax: ok[v] = 1
            for i in range(n):
                other = np.convolve(fwd[i], bwd[i + 1])[: maxs + 1]
                other = np.convolve(other, cin)
                a = W[i, j]
                m0 = (other * ok[: other.size]).sum()
                sh = np.zeros_like(other); sh[a:] = other[: other.size - a] if a else other
                m1 = (sh * ok[: sh.size]).sum()
                new[j, i] = m1 / (m0 + m1) if m0 + m1 > 0 else 0.5
            # upward carry message (distribution of c_{j+1})
            full = np.convolve(np.convolve(fwd[n], np.ones(1)), cin)
            cm = np.zeros(Cmax + 1)
            for v in range(full.size):
                if (v - T[j]) % 5 == 0 and 0 <= (v - T[j]) // 5 <= Cmax: cm[(v - T[j]) // 5] += full[v]
            if cm.sum() > 0: carry_msg[j + 1] = 0.5 * carry_msg[j + 1] + 0.5 * cm / cm.sum()
        mfb = damp * mfb + (1 - damp) * new
    lo = np.log(np.clip(mfb, 1e-12, 1)).sum(0); l0 = np.log(np.clip(1 - mfb, 1e-12, 1)).sum(0)
    return 1 / (1 + np.exp(np.clip(l0 - lo, -50, 50)))
if __name__ == "__main__":
    npos, K, trials = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); M = 5 ** K
    inv = pow(16, -1, M); dig = [1] + [(pow(16, t, M) + pow(inv, t, M)) % M for t in range(1, npos)]
    w = []
    for v in dig: w += [6 * v % M, 9 * v % M]
    n = len(w); rnd = random.Random(1); acc = []; solved = 0
    for _ in range(trials):
        xs = [rnd.randrange(2) for _ in range(n)]; t = sum(a * b for a, b in zip(w, xs)) % M
        p = run(w, t, K)
        acc.append(np.mean([(p[i] > 0.5) == xs[i] for i in range(n)]))
        fixed = {}
        for step in range(n):                      # decimation
            p = run(w, t, K, iters=25, fixed=fixed)
            free = [i for i in range(n) if i not in fixed]
            i = max(free, key=lambda k: abs(p[k] - 0.5)); fixed[i] = int(p[i] > 0.5)
        x = [fixed[i] for i in range(n)]
        solved += (sum(a * b for a, b in zip(w, x)) % M == t)
    print(f"{npos} digit positions ({n} bits), K={K}, density {n/(K*math.log2(5)):.2f}: BP marginals predict {np.mean(acc)*100:.1f}% of planted bits "
          f"(random 50%); decimation solves {solved}/{trials}")

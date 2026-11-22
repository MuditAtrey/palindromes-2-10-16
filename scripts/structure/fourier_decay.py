"""B2: Fourier coefficients of the inner digit measure mu = (*)_tau uniform{a * w_tau : a in A} on Z/5^K, at every nonzero frequency.
Large |mu^(lambda)| = a frequency at which the inner sums are badly spread = a strong window test. Inner weights vs random weights."""
import sys, numpy as np
n, K = int(sys.argv[1]), int(sys.argv[2]); A = np.array([0, 6, 9, 15]); M = 5 ** K
inv16 = pow(16, -1, M)
w_in = [1] + [(pow(16, t, M) + pow(inv16, t, M)) % M for t in range(1, n)]
rng = np.random.default_rng(3); w_rd = [1] + [int(x) for x in rng.integers(1, M, n - 1)]
def spectrum(w):
    top = np.zeros(0); hist = np.zeros(6, dtype=np.int64); thr = [0.5, 0.3, 0.2, 0.1, 0.05, 0.02]
    for s in range(1, M, 1 << 20):
        lam = np.arange(s, min(M, s + (1 << 20)), dtype=np.float64)
        f = np.ones(lam.size, dtype=np.complex128)
        for wt in w:
            ph = (lam * (wt % M)) % M / M
            f *= np.exp(2j * np.pi * np.outer(ph, A) % (2j * np.pi) if False else 2j * np.pi * np.outer(ph, A)).mean(axis=1)
        a = np.abs(f)
        for i, t in enumerate(thr): hist[i] += int((a > t).sum())
        top = np.sort(np.concatenate([top, np.sort(a)[-10:]]))[-10:]
    return top[::-1], hist
for name, w in (("inner", w_in), ("random", w_rd)):
    top, hist = spectrum(w)
    print(f"{name:6s} n={n} K={K}: top |mu^| {np.round(top[:5], 4)} | frequencies with |mu^| > 0.5,0.3,0.2,0.1,0.05,0.02: {hist.tolist()}")

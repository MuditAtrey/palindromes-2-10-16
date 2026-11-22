# Theorem 1: triple palindromes are rarer than the trivial bound

Status: proof written 2026-11-21 in this project; every constant checked numerically (`scripts/proof/`). Needs an independent
careful reading before it is claimed anywhere. It proves rarity, not finiteness.

## Statement

Let T(X) be the number of integers N <= X that are palindromes in bases 2, 10 and 16. There are absolute constants C, c > 0 with

    T(X) <= C * X^(1/4) * 2^(-c sqrt(log_16 X))          (c = 0.04 works)

The trivial bound is T(X) <= #{base-2/16 palindromes <= X} ~ X^(1/4), so this is a (sub-power) saving over counting all candidates.

## Setting

Write B for the bit length of N and L = ceil(B/4) for its hex length.
- B = 4L: the hex digits lie in {0,6,9,15} = {6a + 9b : a, b in {0,1}} and form a palindrome (each hex digit is its own 4-bit
  reversal); about 2^L candidates.
- B = 4L - 1: the hex digits lie in {0,2,5,7} = {2a + 5b}, palindrome; about 2^L candidates.
- B = 4L - 2 and 4L - 3: at most about 2^(B/8) candidates (computed exactly for B <= 42: 1, 2, 2, 4, 4, 8, 8, 16, 16, 32;
  `scripts/proof/small_classes.py`). These contribute O(X^(1/8)) and are negligible. (A short structural proof of this count is still
  to be written; the computation covers the pattern.)

The argument below is for B = 4L; the B = 4L - 1 case differs only in the constants (given at the end).

## Step 1: one outer choice gives one 5-adic target

Fix d < L/2 and the d outermost digit pairs. Their contribution O is fixed, and N = O + I where I is the inner part,
I = sum_{j=d}^{L/2-1} e_j W_j, W_j = 16^j + 16^(L-1-j) (centre digit counted once if L is odd), with n = L/2 - d inner digits, and
0 <= I < 16^(L-d). All such N lie in an interval of length 16^(L-d) starting at O. If N has D decimal digits, then with
K = D - ceil(log_10 16^(L-d)) - 1 (so K >= 1.204 d - 3) the top K decimal digits of N take at most two values on this interval (one carry),
and D takes at most two values. A decimal palindrome has N = R (mod 10^K), R = the top K digits read backwards, so in particular
I = R - O (mod 5^K). Hence the number of triple palindromes extending this outer choice is at most

    sum over at most 4 targets t of  N_K(t),     N_K(t) = #{ e in A^n : sum_j e_j W_j = t (mod 5^K) },   A = {0,6,9,15}.

## Step 2: Fourier

With e(x) = exp(2 pi i x) and
    F(h) = prod_j (1/4) sum_{a in A} e(a h W_j / 5^K),
orthogonality gives N_K(t) = (4^n / 5^K) sum_h F(h) e(-h t / 5^K), so

    | N_K(t) - 4^n / 5^K |  <=  4^n * max_{h != 0 mod 5^K} |F(h)|.

Since a = 6a' + 9b', the factor is ((1 + e(6x))/2)((1 + e(9x))/2) with x = h W_j / 5^K, so

    |F(h)| = prod_j phi(x_j),   phi(x) = |cos(6 pi x) cos(9 pi x)|,   x_j = h W_j / 5^K mod 1.

phi(x) = 1 exactly when x is in (1/3)Z; away from those points it is smaller: with eps = 1/1000,
c1 := max{ phi(x) : dist(x, (1/3)Z) >= eps } = 0.999423 (computed).

## Step 3: the phases come from x16 and x16^(-1) dynamics

W_j = 4^(L-1) (4^s + 4^(-s)) with s = L - 1 - 2j, and 4^(L-1) is a unit mod 5^K, so after replacing h by h * 4^(L-1) (still nonzero
mod 5^K) the phases are x_j = y_j + z_j (mod 1) with

    y_j = { h 4^s / 5^K },   z_j = { h 4^(-s) / 5^K }   (4^(-s) taken mod 5^K).

Consecutive j change s by 2. Read the inner positions in order of increasing s (the order does not matter for runs): then
y' = {16 y} and 16 z' = z (mod 1), i.e. z' = (z + c)/16 with c in {0,...,15}. (In order of increasing j the roles of y and z are exchanged;
the argument below is symmetric in y and z.)

Let r be the exact 5-adic denominator: v_5(h) = K - r, 1 <= r <= K. Then y_j = A_j / 5^r with 5 not dividing A_j.

## Step 4: bad stretches are short (the key lemma)

Call position j bad if dist(x_j, (1/3)Z) < eps.

Lemma. If r >= 2 and eps <= 1/1000, every run of consecutive bad positions has length at most 0.581 r + 1.

Proof. Suppose j and j+1 are both bad: y_j + z_j = k/3 + m + e and y_{j+1} + z_{j+1} = k'/3 + m' + e' with |e|, |e'| < eps, integers k, k', m, m'.
Substituting y_{j+1} = 16 y_j - a (a an integer) and z_{j+1} = (z_j + c)/16 = (k/3 + m - y_j + e + c)/16 gives
(255/16) y_j = k'/3 + m' + a - (k/3 + m + c)/16 + e' - e/16, so y_j is within (16/255)(eps + eps/16) = eps/15 of a point of (1/765)Z.
Write y_j = q + delta, q in (1/765)Z, |delta| <= eps/15. Then 16 q is in (1/765)Z, and because 17 eps/15 < 1/765 for eps <= 1/1000 the point of
(1/765)Z nearest to y_{j+1} is {16 q}, so delta_{j+1} = 16 delta_j. Within a bad run the deviation therefore multiplies by 16 at every step
while staying below eps/15. It starts nonzero and at least 1/(153 * 5^r): y_j - q = (153 A_j - 5^(r-1) J)/(153 * 5^r) and the numerator is
nonzero because 5 does not divide 153 A_j while r >= 2. So a run of m + 1 consecutive bad positions needs 16^(m-1) / (153 * 5^r) <= eps / 15,
i.e. m <= 1 + log_16(153 * 5^r * eps / 15) = 0.5805 r - 0.654 for eps = 1/1000. QED.

(Numerical check, `scripts/proof/run_bound.py`, x4 version, moduli up to 5^24, 3,000 random and many adversarial h: longest run 0.69 of
the bound.)

For r = 1 every factor is small outright: x_j = a/5 with 5 not dividing a, and max phi(a/5) = 0.6545.

Consequence: for every h != 0 (mod 5^K), at least one position in every block of 0.581 K + 2 consecutive inner positions is good, so

    |F(h)| <= c1^floor(n / (0.581 K + 2)).

## Step 5: choosing the depth

Summing Step 1 over the 4^d outer choices (4^d * 4^n = 2^L):

    #{triple palindromes with B = 4L bits} <= 4 * 2^L * ( 5^(-K) + c1^floor(n / (0.581 K + 2)) ),   K >= 1.204 d - 3,  n = L/2 - d.

With kappa = log_2(1/c1) = 8.3e-4 and d = alpha sqrt(L): the first term saves 2^(-2.79 alpha sqrt(L)), the second
2^(-(kappa / (1.40 alpha)) sqrt(L)) (to leading order). Balancing gives alpha = sqrt(kappa / 3.9) ~ 0.0146 and a saving of 2^(-0.04 sqrt(L)).
Summing the geometric series over L <= log_16 X + 1 gives the theorem.

## Class B = 4L - 1 (digits {0,2,5,7} = 2a + 5b)

phi(x) = |cos(2 pi x) cos(5 pi x)|, equal to 1 only on Z. The same steps with bad set Z: (255/16) y_j = integer combination / 16 + errors,
so y_j is within eps/15 of (1/255)Z; separation from A/5^r is at least 1/(51 * 5^r); the nearest-point step needs eps < 1/289. With
eps = 1/300: c2 = 0.998409 (kappa = 2.3e-3), runs <= 0.581 r + 1, r = 1 factors <= 0.809. Same conclusion with a slightly better constant.

## What this does and does not say

- It is, as far as the literature search found, the first bound below the trivial count for numbers palindromic in two or three bases.
- It does not approach finiteness: finiteness needs the count per length to go below 1, and this argument saves only 2^(-0.04 sqrt(L)) from
  2^L. The obstruction to a power saving is exactly Step 4: for special h the bad runs can have length proportional to r, i.e. to K. A positive
  proportion of good positions for every h (bounded runs) would give a power saving; that is a statement about short segments of
  4^s mod 5^r for r much larger than log(n), equivalent to the "no unexpectedly short dual vectors" property measured in STATUS 39.
- Using the 2-adic half as well (only about 2^(0.8 d) of the 4^d outer choices pass it) would improve Step 5, but needs a rigorous count of
  those survivors, which is the Furstenberg-type statement that is not available for reversal-coupled sets.

## Adversarial checks (2026-11-21)

- The key lemma tested exactly as used (x16 dynamics, s in steps of 2, eps = 1/1000, both parities, moduli up to 5^28, 1,500 random and
  ~450 adversarial h per modulus including small h, h near 5^K and multiples of powers of 5): longest bad run 0.395 of the bound 0.581 r + 1
  (`scripts/proof/run_bound16.py`).
- Wording fixed: along increasing j the roles of y and z swap; the lemma is symmetric in them.

# Credits and references

This project stands on other people's work. Everything below was used, tested against, or cited in the documents, and every
link was checked on 2026-11-22 (DOIs resolved through Crossref, arXiv entries through the arXiv API).

## The question and earlier searches

| Who | Work | How it was used |
|---|---|---|
| Rick Regan | [In Search of Decimal/Binary/Hexadecimal Palindromes](https://www.exploringbinary.com/in-search-of-decimal-binary-hexadecimal-palindromes/), Exploring Binary, 2010 | Posed the question publicly and reported the first search; the reason this problem is known. |
| Matt Schoen | [mtschoen/multibase-palindrome-search](https://github.com/mtschoen/multibase-palindrome-search) (MIT license), 2026 | The previous record (no solution below 2^255). This engine follows the same line of work: it shares the program name `tripal` and the starting idea of building candidates from their base-2 and base-16 structure. It shares no source code (checked line by line). His open challenges (issues #1 to #6) are good starting points for newcomers. |
| Nev Delap | [nevdelap/pp](https://github.com/nevdelap/pp), 2026 | Search for primes that are palindromes in bases 2, 10 and 16, to 256 bits. Related work. |
| Eshed Shaham | [Finding Binary & Decimal Palindromes](https://ashdnazg.github.io/articles/22/Finding-Really-Big-Palindromes) | Search method behind most known base-2/10 double palindromes (OEIS A007632). Related work. |
| Charlton Harrison | [Binary/Decimal Palindromes](http://bach.dynet.com/palin/) | Early records for base-2/10 double palindromes. Related work. |
| Patrick De Geest | [Palindromic numbers beyond base 10](https://www.worldofnumbers.com/nobase10.htm), World!Of Numbers | Long-running collection of multi-base palindromes. Related work. |
| Project Euler | [Problem 36: Double-base palindromes](https://projecteuler.net/problem=36) | The classic small version of the question. |

## OEIS data

The On-Line Encyclopedia of Integer Sequences, [https://oeis.org](https://oeis.org), is licensed under
[CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

- [A007632](https://oeis.org/A007632), numbers palindromic in bases 2 and 10. The file `engine/data/oeis/b007632.txt` is a copy
  of its [b-file](https://oeis.org/A007632/b007632.txt) by Jeff Sponaugle (terms below 10^64; terms 1 to 175 from Eshed Shaham;
  terms 1 to 147 variously by Robert G. Wilson v, Charlton Harrison, Ilya Nikulshin and Andrey Astrelin; terms 184 to 197 and
  202 to 213 from Richard Jones, and 198 to 201 from Eshed Shaham, via Patrick De Geest's World!Of Numbers). It is the ground
  truth for the engine's self-test and the data that validated the random model. That file stays under CC BY-SA 4.0
  (see [engine/data/oeis/README.md](engine/data/oeis/README.md)).
- [A029731](https://oeis.org/A029731) (bases 10 and 16, Patrick De Geest; table by Jeff Sponaugle and Robert G. Wilson v) and
  [A340559](https://oeis.org/A340559) (bases 2 and 16, Glen Gilchrist; table by Chai Wah Wu): the two other pairs.
- [A146025](https://oeis.org/A146025), numbers written with only 0 and 1 in bases 2, 3, 4 and 5 (the "82000 problem"): cited as a
  problem of the same shape.

## Mathematics

| Reference | Used for |
|---|---|
| A. Bérczes, V. Ziegler, [On Simultaneous Palindromes](https://arxiv.org/abs/1403.0787), arXiv:1403.0787, 2014 | the theory of numbers palindromic in two bases |
| E. H. Goins, [Palindromes in Different Bases: A Conjecture of J. Ernest Wilkins](https://arxiv.org/abs/0909.5452), Integers 9, 2009 | known finiteness results for palindromes of equal digit length |
| H. G. Senge, E. G. Straus, [PV-numbers and sets of multiplicity](https://doi.org/10.1007/BF02018464), Periodica Math. Hungarica 3, 1973 | only finitely many numbers have few nonzero digits in two independent bases |
| C. L. Stewart, [On the representation of an integer in two different bases](https://doi.org/10.1515/crll.1980.319.63), J. reine angew. Math. 319, 1980 | the effective form of Senge-Straus (no sparse formula can give infinitely many solutions) |
| A. Cobham, [On the base-dependence of sets of numbers recognizable by finite automata](https://doi.org/10.1007/BF01746527), Math. Systems Theory 3, 1969 | why no finite automaton can read base 10 and base 16 together |
| S. A. Burrell, H. Yu, [Digit expansions of numbers in different bases](https://arxiv.org/abs/1905.00832), arXiv:1905.00832 | the closest known counting method; where it stops for palindromes (PROOF_RESEARCH 4b) |
| P. Shmerkin, [On Furstenberg's intersection conjecture...](https://doi.org/10.4007/annals.2019.189.2.1), Annals of Math. 189, 2019; M. Wu, [A proof of Furstenberg's conjecture on the intersections of xp- and xq-invariant sets](https://doi.org/10.4007/annals.2019.189.3.2), Annals of Math. 189, 2019 | the intersection bound that the search's survivor counts match |
| J. C. Lagarias, [Ternary expansions of powers of 2](https://doi.org/10.1112/jlms/jdn080), J. London Math. Soc. 79, 2009 ([arXiv](https://arxiv.org/abs/math/0512006)); M. A. Idowu, [Repunit Coordinates for the Erdős Ternary-Digit Problem](https://arxiv.org/abs/2610.03789), arXiv:2610.03789 | Erdős's ternary problem, the open problem of the same shape as this one (WALL.md) |
| A. Baker, [Linear forms in the logarithms of algebraic numbers](https://doi.org/10.1112/S0025579300003971), Mathematika 13, 1966 | the Baker route, tested and closed (WALL.md, section 1) |
| H. Weyl, [Über die Gleichverteilung von Zahlen mod. Eins](https://doi.org/10.1007/BF01475864), Math. Annalen 77, 1916 | equidistribution step in Theorem 2 |
| A. Tuxanidy, D. Panario, [Infinitude of palindromic almost-prime numbers](https://arxiv.org/abs/2307.16637), arXiv:2307.16637 | survey of exponential-sum methods for palindromes (PROOF_RESEARCH) |

## Algorithms

| Reference | Used for |
|---|---|
| A. K. Lenstra, H. W. Lenstra Jr., L. Lovász, [Factoring polynomials with rational coefficients](https://doi.org/10.1007/BF01457454), Math. Annalen 261, 1982 | LLL reduction, behind the dual-lattice certificate |
| C. P. Schnorr, M. Euchner, [Lattice basis reduction: improved practical algorithms and solving subset sum problems](https://doi.org/10.1007/BF01581144), Math. Programming 66, 1994 | BKZ and the enumeration used by the lattice test |
| B. H. Bloom, [Space/time trade-offs in hash coding with allowable errors](https://doi.org/10.1145/362686.362692), CACM 13, 1970 | the Bloom reach filter |
| D. Coppersmith, [Small solutions to polynomial equations, and low exponent RSA vulnerabilities](https://doi.org/10.1007/s001459900030), J. Cryptology 10, 1997 | tested and closed (WHAT_I_TRIED) |
| I. Dinur, O. Dunkelman, N. Keller, A. Shamir, [Efficient dissection of composite problems](https://doi.org/10.1007/978-3-642-32009-5_42), CRYPTO 2012 | tested and closed (WHAT_I_TRIED) |
| A. Montanaro, [Quantum walk speedup of backtracking algorithms](https://arxiv.org/abs/1509.02374), arXiv:1509.02374 | considered and closed (WHAT_I_TRIED) |

## Software

| Project | Used for |
|---|---|
| [FLINT](https://flintlib.org) | LLL and exact integer matrices in the dual-lattice certificate and lattice test |
| [GMP](https://gmplib.org) and [MPFR](https://www.mpfr.org) | arbitrary-precision arithmetic under FLINT and in the research C programs |
| [OpenTimestamps](https://opentimestamps.org) | the timestamped commitment (proof of priority) |
| [NumPy](https://numpy.org), [SciPy](https://scipy.org), [mpmath](https://mpmath.org) | the research scripts |

If your work is used here and is missing or misdescribed, please open an issue and it will be corrected.

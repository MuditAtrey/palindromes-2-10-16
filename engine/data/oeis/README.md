# OEIS data in this folder

`b007632.txt` is a copy of the b-file of OEIS [A007632](https://oeis.org/A007632), "Numbers that are palindromic in bases 2
and 10", from The On-Line Encyclopedia of Integer Sequences, [https://oeis.org](https://oeis.org).

- Original file: [https://oeis.org/A007632/b007632.txt](https://oeis.org/A007632/b007632.txt), by Jeff Sponaugle (all terms
  below 10^64; terms 1 to 175 from Eshed Shaham; terms 1 to 147 variously by Robert G. Wilson v, Charlton Harrison, Ilya
  Nikulshin and Andrey Astrelin; terms 184 to 197 and 202 to 213 from Richard Jones, and 198 to 201 from Eshed Shaham, via
  Patrick De Geest's World!Of Numbers).
- Change made here: three comment lines at the top saying where the file comes from and what it is used for.
- License: CC BY-SA 4.0, as for all OEIS content ([OEIS End-User License Agreement](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement)).
  This file, unlike the rest of the repository's data, is therefore share-alike: copies and adaptations must keep CC BY-SA 4.0.

The engine's self-test (`tripal selftest`) drops the base-16 condition and checks that the search reproduces exactly the terms
in this file below 2^48. `tools/gen_b007632.py` regenerates the smaller terms independently.

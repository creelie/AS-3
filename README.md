# Positive Black Hole Degeneracies for Dyons of Higher Torsion

Paper (`paper/main.tex`, `paper/main.pdf`) and the programs that check it.

- `code/consts.py`: every constant of the analytic proofs and the exact tables for P²/2 = 1, 2.
- `code/proofcheck.py`, `code/proofcheck2.py`: the coefficient lemmas checked against exact series.
- `code/psi.py`, `code/psif.py`, `code/dstar.py`: finite parts of the Fourier–Jacobi coefficients of 1/Φ₁₀ (data in `code/data/cF_240.pkl`).
- `code/numgen.py data/cF_240.pkl 240`, `code/numtors.py data/cF_240.pkl 240`: the doubling inequality and torsion positivity for general magnetic charge.
- `paper/figures/make_figs.py`: the figures.

Requires Python 3 with python-flint, sympy, numpy and matplotlib.

Author: Deep Bhattacharjee. Licence: see `LICENSE`.

# Programs for every magnetic charge

These programs check the finite arithmetic in "Sen's Positivity Conjecture at Every Torsion and Every Magnetic Charge". Every proof in the paper is analytic; the programs only re-check numbers.

- `expansion/`: exact coefficients of the finite part from the product formula (`check_exact.py`), the polar formula (`check_polar.py`) and the level-one Rademacher expansion with explicit tails (`check_rademacher.py`). Logs: `checks*.log`.
- `certs/`: ball-arithmetic certificates (outward rounding) for the error bound and the doubling inequality for every magnetic charge (`cert_rangeA.py`, `cert_rangeB.py`, `cert_doubling.py`), exact values and torsion sums (`cert_exact.py`), wall sums (`polar_sb.py`) and the expansion at sample charges (`rad_check.py`).

Run each script from its own folder with `python3 <script>`. Requires Python 3 with python-flint, mpmath and sympy.

# Dense direct solvers for equality constrained QPs

Code, data and results for my CPSC 542G midterm project, *Dense Direct Solvers for
Equality Constrained Quadratic Optimization*.

The problem is

    min_x  1/2 x'Hx - g'x   subject to   Cx = d

with `C` of full row rank. I compare four dense direct methods:

- **QR null space.** Full QR of `C'`, orthonormal null space basis `Z`, Cholesky of `Z'HZ`.
- **Variable reduction null space.** `Z = [I; -C2\C1]`, where the columns in `C2` are
  picked by LU with partial pivoting on `C'`. No QR and no SVD.
- **Schur complement.** Cholesky of `H`, then Cholesky of `S = C H^{-1} C'`.
- **LDL' of the KKT matrix.** Bunch-Kaufman LDL' of `[H C'; C 0]`. The signs of the
  pivots give the inertia, which should be `(n, p, 0)` when `Z'HZ` is positive definite.

The first three need `H` (or `Z'HZ`) to be positive definite. The LDL' solver doesn't, and
it reports the inertia so you can tell a minimizer from a saddle point.

## Layout

| Path | What's in it |
|---|---|
| `matlab/src/solve_qr_null.m` | QR null space solver |
| `matlab/src/solve_vr_null.m` | variable reduction null space solver |
| `matlab/src/solve_schur.m` | Schur complement solver |
| `matlab/src/solve_kkt_ldl.m` | LDL' of the KKT matrix with the inertia check |
| `matlab/src/ldl_inertia.m` | counts the inertia from the 1x1 and 2x2 blocks of `D` |
| `matlab/src/kkt_ldl_probe.m` | LDL' with a shift `H + alpha*I` or `-delta*I` in the constraint block, used for the regularization runs |
| `matlab/src/make_dense_case.m` | random dense problems with prescribed `cond(H)` and `cond(C)` |
| `matlab/src/make_indefinite_case.m` | problems with indefinite `H` where `Z'HZ` is either positive definite or has a known number of negative eigenvalues |
| `matlab/src/make_kriging_case.m` | ordinary kriging problem on the road elevation data |
| `matlab/src/orthogonal_benchmark.m` | the change of variables `x = Uy` that makes DTOC3 fully dense |
| `matlab/src/suite_measure.m` | timing and error measurement for one problem and a list of methods |
| `matlab/src/method_list.m` | maps method names to solvers |
| `matlab/run_*.m` | experiment drivers, see below |
| `matlab/data/` | stored problems (DTOC3 and AUG2DC from CUTEst, GEMAT1 least squares setups, the n = 80 conditioning problems, the road elevation sample) |
| `matlab/provenance/` | SIF files, export records, hashes and validation output for the stored problems |
| `python/` | scripts that built the stored problems and independent checks, see `python/README.md` |
| `results/` | `results.csv` and `environment.txt` for every experiment in the writeup |

## Running

Everything in `matlab/` runs in plain MATLAB without toolboxes. I ran it in MATLAB Online
R2026b. The problems are stored in `matlab/data/`, so Python and CUTEst aren't needed. From
the `matlab/` folder:

    run_literature('check'); run_kriging('check');   % small checks
    run_vr('check'); run_direct3('check'); run_indefinite('check');

    run_final_dense('dimension')      % QR and Schur, n = 1000 ... 12000, p/n = 0.5
    run_final_dense('fraction')       % QR and Schur, n = 12000, p/n = 0.1 ... 0.9
    run_final_dense('illcond')        % QR and Schur, cond(H) = 1e8, cond(C) = 1e6
    run_final_dense('conditioning')   % QR and Schur, n = 12000 conditioning sweep
    run_conditioning_tests()          % QR and Schur, n = 80 with 96 digit references

    run_vr('dimension'); run_vr('fraction'); run_vr('illcond');
    run_vr('conditioning'); run_vr('applications')

    run_direct3('pilot'); run_direct3('dimension'); run_direct3('conditioning')

    run_indefinite('full')            % indefinite H, inertia
    run_regularize('full')            % shifts in the constraint block and in H

    run_literature('full'); run_kriging('full')
    package_results()

Each driver writes a folder under `matlab/results/` with `results.csv`, `environment.txt`
and the solution vectors. The largest cases need about 9 GiB of RAM. The dense LDL' runs
stop at n = 8000 because the KKT matrix and its factor don't fit at n = 12000.

A timed solve is one full solver call. Generating or loading the problem, the reference
and the error computation are outside the timer.

## Sources

- M. Benzi, G. H. Golub, J. Liesen, *Numerical solution of saddle point problems*, Acta Numerica 14 (2005).
- N. I. M. Gould, M. E. Hribar, J. Nocedal, *On the solution of equality constrained quadratic programming problems arising in optimization*, SIAM J. Sci. Comput. 23 (2001).
- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer (2006), Chapter 16.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM (2002), Chapter 11.
- J. R. Bunch, L. Kaufman, *Some stable methods for calculating inertia and solving symmetric linear systems*, Math. Comp. 31 (1977).
- DTOC3 and AUG2DC from the CUTEst SIF collection, https://github.com/ralna/SIF, exported with PyCUTEst 1.8.2.
- GEMAT1 from the SuiteSparse Matrix Collection (HB/gemat1). The least squares setup follows J. Scott and M. Tůma, BIT Numer. Math. 62 (2022).
- Road elevation data from M. Kaul, 3D Road Network (North Jutland, Denmark), UCI Machine Learning Repository, CC BY 4.0.

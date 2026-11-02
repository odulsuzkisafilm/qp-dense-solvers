# Dense direct solvers for equality-constrained quadratic optimization

Code, data and results for the CPSC 542G midterm project
*Dense Direct Solvers for Equality-Constrained Quadratic Optimization:
a QR Null-Space Method versus a Schur-Complement Method* (Efe Cinar, November 2026).

The problem is

    min_x  1/2 x'Hx - g'x   subject to   Cx = d,

with `H` symmetric positive definite and `C` of full row rank. Two dense solvers are
compared: a null-space method built on a full QR of `C'` and a Schur-complement method
built on two Cholesky factorizations.

## Layout

| Path | Contents |
|---|---|
| `matlab/src/solve_qr_null.m`, `matlab/src/solve_schur.m` | the two solvers |
| `matlab/src/` | input checks, problem generators (`make_dense_case.m`, `make_kriging_case.m`), the orthogonal change of variables (`orthogonal_benchmark.m`), the timing harness (`suite_measure.m`) |
| `matlab/run_*.m` | experiment drivers: `run_final_dense` (size / fraction / illcond / conditioning sweeps), `run_literature` (DTOC3, GEMAT1), `run_kriging`, plus earlier `run_cutest_tests`, `run_conditioning_tests`, `run_dense_scaling` |
| `matlab/package_results.m`, `matlab/plot_suite_results.m` | zip results; summary tables |
| `matlab/data/*.mat` | stored problems: DTOC3 and AUG2DC exports from CUTEst, GEMAT1 least-squares setups, the n=80 conditioning problems, and the UCI road/elevation sample used for kriging |
| `matlab/provenance/` | SIF sources, export manifests, SHA-256 hashes, candidate-problem decisions, validation reports |
| `python/` | scripts that produced the `.mat` files and independent checks, see `python/README.md` |
| `results/` | raw `results.csv` and `environment.txt` of every experiment reported |

## Running the experiments

Everything in `matlab/` runs in plain MATLAB (tested in MATLAB Online R2026b); no
toolboxes, no Python and no CUTEst installation are needed because the problems are
stored in `matlab/data/`. From the `matlab/` folder:

    run_literature('check'); run_kriging('check');     % small checks with assertions
    run_final_dense('dimension');                       % n = 1000 ... 12000, p/n = 0.5
    run_final_dense('fraction');                        % n = 12000, p/n = 0.1 ... 0.9
    run_final_dense('illcond');                         % same sizes, cond(H)=1e8, cond(C)=1e6
    run_literature('full');  run_kriging('full');       % application problems
    run_final_dense('conditioning');                    % n = 12000 conditioning sweep
    run_conditioning_tests();                           % n = 80 sweep with stored references
    package_results();

Each driver writes a timestamped folder under `matlab/results/` with `results.csv`,
`environment.txt` and the solution vectors. The largest cases need about 9 GiB of RAM.
A timed solve is one complete solver call; generation, loading, references and error
evaluation are outside the timer.

The Python scripts are optional and document provenance; they need NumPy, SciPy and
mpmath, and `export_cutest.py` additionally needs PyCUTEst.

## Sources

- Methods: M. Benzi, G. H. Golub, J. Liesen, *Numerical solution of saddle point problems*, Acta Numerica 2005, Sections 5–6.
- DTOC3, AUG2DC: CUTEst / SIF collection, https://github.com/ralna/SIF (exported with PyCUTEst 1.8.2).
- GEMAT1: SuiteSparse Matrix Collection, HB/gemat1; equality-constrained least-squares setup after J. Scott and M. Tůma, BIT Numer. Math. 62 (2022).
- Elevation data: M. Kaul, 3D Road Network (North Jutland, Denmark), UCI Machine Learning Repository, CC BY 4.0.

All solver and experiment code is the author's own.

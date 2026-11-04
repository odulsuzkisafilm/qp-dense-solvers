# Python scripts

These scripts document how the stored problems in `matlab/data/` were produced and
provide independent checks. None of them is needed to run the MATLAB experiments.
They require NumPy, SciPy and mpmath. `export_cutest.py` also needs PyCUTEst with
a working CUTEst installation and the SIF files in `matlab/provenance/`.

| Script | Purpose |
|---|---|
| `reproduce/export_cutest.py NAME '{"PARAM": value}'` | export a CUTEst problem (DTOC3, AUG2DC) to a `.mat` file with a sparse double precision KKT reference |
| `reproduce/build_conditioning.py` | build the n = 80 conditioning problems and their 96 digit references |
| `reproduce/build_literature.py` | build the GEMAT1 least squares problems from the SuiteSparse matrix |
| `reproduce/check_lse_reference.py` | recompute the GEMAT1 reference with an augmented system |
| `reproduce/validate_kriging_python.py` | run the QR and Schur solvers in SciPy on small kriging problems and compare with a direct KKT solve |
| `reproduce/validate_new_methods.py` | NumPy versions of the LDL inertia count, the indefinite problem generator and the variable reduction solver, run on small cases and on the n = 4000 indefinite problems |
| `precision_check/check_precision.py` | rerun the ten failed Schur cases in mpmath at 64 and 96 digits |
| `make_deck_data.py RESULTS_DIR` | turn `results.csv` files into the plot data used by the slides and the writeup |

Outputs of the first five scripts go to `matlab/data/` and `matlab/provenance/`.
`check_precision.py` writes `precision_results.csv` next to itself.

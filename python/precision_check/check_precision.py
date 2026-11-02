from pathlib import Path
import csv, json, hashlib, time
import numpy as np
from scipy.io import loadmat
import scipy.linalg as la
import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "matlab" / "data"
RESULTS = ROOT / "results" / "conditioning_20261029_040845_691" / "results.csv"
with RESULTS.open() as fh:
    failures = [r for r in csv.DictReader(fh) if r["method"] == "Schur" and r["status"] == "failed"]
assert len(failures) == 10
hashes = json.loads((ROOT / "matlab" / "provenance" / "SHA256.json").read_text())


def cholsolve(L, B):
    n = L.rows
    Y = mp.matrix(n, B.cols)
    X = mp.matrix(n, B.cols)
    for j in range(B.cols):
        for i in range(n):
            Y[i, j] = (B[i, j] - mp.fsum(L[i, k] * Y[k, j] for k in range(i))) / L[i, i]
        for i in range(n - 1, -1, -1):
            X[i, j] = (Y[i, j] - mp.fsum(L[k, i] * X[k, j] for k in range(i + 1, n))) / L[i, i]
    return X


def mp_schur(Hd, Cd, gd, dd, dps):
    with mp.workdps(dps):
        H = mp.matrix(Hd.tolist())
        C = mp.matrix(Cd.tolist())
        g = mp.matrix(gd.tolist())
        d = mp.matrix(dd.tolist())
        L = mp.cholesky(H)
        u = cholsolve(L, g)
        W = cholsolve(L, C.T)
        S = C * W
        S = (S + S.T) / 2
        LS = mp.cholesky(S)
        lam = cholsolve(LS, C * u - d)
        x = u - W * lam
        ev = mp.eigsy(S, eigvals_only=True)
        return x, lam, S, ev[0], ev[-1]


def ref_kkt(Hd, Cd, gd, dd):
    with mp.workdps(110):
        p = len(dd)
        K = mp.matrix(np.block([[Hd, Cd.T], [Cd, np.zeros((p, p))]]).tolist())
        b = mp.matrix(np.r_[gd, dd].tolist())
        v = mp.lu_solve(K, b)
        return v, K, b


rows = []
for f in failures:
    begin = time.time()
    name = f["case_id"]
    path = DATA / (name + ".mat")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == hashes["data/" + path.name]
    P = loadmat(path, squeeze_me=True)
    H = P["H"]
    C = P["C"]
    g = P["g"]
    d = P["d"]
    n = len(g)
    L = la.cholesky(H, lower=True)

    def ds(B):
        return la.solve_triangular(L.T, la.solve_triangular(L, B, lower=True), lower=False)

    Sd = C @ ds(C.T)
    Sd = (Sd + Sd.T) / 2
    try:
        la.cholesky(Sd, lower=True)
        double_status = "success"
    except la.LinAlgError:
        double_status = "failed"
    vref, K, b = ref_kkt(H, C, g, d)
    x64, l64, S64, m64, M64 = mp_schur(H, C, g, d, 64)
    x96, l96, S96, m96, M96 = mp_schur(H, C, g, d, 96)
    with mp.workdps(110):
        xr = vref[:n]
        lr = vref[n:]
        v = mp.matrix(list(x96) + list(l96))
        sd = mp.matrix(Sd.tolist())
        emin_double = mp.eigsy(sd, eigvals_only=True)[0]

        def rel(a, b):
            return mp.norm(a - b) / mp.norm(b)

        vals = dict(
            x_error_64=rel(x64, xr),
            lambda_error_64=rel(l64, lr),
            x_error_96=rel(x96, xr),
            lambda_error_96=rel(l96, lr),
            x_precision_agreement=rel(x64, x96),
            lambda_precision_agreement=rel(l64, l96),
            S_min_eigenvalue_96=m96,
            S_condition_96=M96 / m96,
            rounded_double_S_min_eigenvalue=emin_double,
            S_relative_Frobenius_error=rel(sd, S96),
            KKT_residual_96=mp.norm(K * v - b) / (mp.norm(K) * mp.norm(v) + mp.norm(b)),
        )
        assert (
            m96 > 0
            and vals["x_error_96"] < mp.mpf("1e-60")
            and vals["lambda_error_96"] < mp.mpf("1e-60")
        )
    row = dict(
        case_id=name,
        kC_target=float(P["kC_target"]),
        input_sha256=digest,
        MATLAB_double_status=f["error_id"],
        SciPy_double_status=double_status,
        mp64_status="success",
        mp96_status="success",
        **{k: mp.nstr(v, 18) for k, v in vals.items()},
        elapsed_seconds=time.time() - begin
    )
    rows.append(row)
    with (OUT / "precision_results.csv").open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (OUT / "precision_results.json").write_text(json.dumps(rows, indent=2))
    print(
        name,
        double_status,
        "mp64/96 OK",
        "min eigen double S",
        row["rounded_double_S_min_eigenvalue"],
        "x96",
        row["x_error_96"],
        flush=True,
    )
print("Completed", len(rows), "cases", flush=True)

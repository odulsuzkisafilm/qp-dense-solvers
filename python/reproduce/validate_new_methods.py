import json
from pathlib import Path
import numpy as np
import scipy.linalg as la


def inertia(D):
    m = D.shape[0]
    pos = neg = zer = 0
    i = 0
    while i < m:
        if i < m - 1 and D[i + 1, i] != 0:
            e = np.linalg.eigvalsh(D[i : i + 2, i : i + 2])
            pos += int((e > 0).sum())
            neg += int((e < 0).sum())
            zer += int((e == 0).sum())
            i += 2
        else:
            if D[i, i] > 0:
                pos += 1
            elif D[i, i] < 0:
                neg += 1
            else:
                zer += 1
            i += 1
    return pos, neg, zer


def solve_ldl(H, C, g, d):
    n, p = H.shape[0], C.shape[0]
    K = np.block([[H, C.T], [C, np.zeros((p, p))]])
    L, D, perm = la.ldl(K, lower=True)
    ine = inertia(D)
    if ine != (n, p, 0):
        raise ValueError(f"inertia {ine}")
    v = la.solve(K, np.r_[g, d], assume_a="sym")
    return v[:n], v[n:], ine


def solve_vr(H, C, g, d):
    n, p = H.shape[0], C.shape[0]
    P, L, U = la.lu(C.T)
    perm = np.argmax(P, axis=0)
    i2, i1 = perm[:p], perm[p:]
    L1, L2 = L[:p], L[p:]
    Y = -la.solve_triangular(L1, L2.T, lower=True, trans="T", unit_diagonal=True)
    x0 = np.zeros(n)
    x0[i2] = la.solve_triangular(L1, la.solve_triangular(U, d, trans="T"), lower=True, trans="T", unit_diagonal=True)
    HZ = H[:, i1] + H[:, i2] @ Y
    G = HZ[i1] + Y.T @ HZ[i2]
    r = g - H @ x0
    z = la.cho_solve(la.cho_factor((G + G.T) / 2), r[i1] + Y.T @ r[i2])
    x = x0.copy()
    x[i1] += z
    x[i2] += Y @ z
    r = g - H @ x
    lam = la.solve_triangular(U, la.solve_triangular(L1, r[i2], lower=True, unit_diagonal=True))
    return x, lam, float(la.norm(Y))


def indefinite_case(n, p, rng, kind, kH=100, kC=10):
    Q, _ = la.qr(rng.standard_normal((n, n)))
    V, Z = Q[:, :p], Q[:, p:]
    q = n - p
    W, _ = la.qr(rng.standard_normal((p, p)))
    C = (W * np.logspace(0, -np.log10(kC), p)) @ V.T
    a = np.logspace(0, -np.log10(kH), p)
    a[::2] *= -1
    b = np.logspace(0, -np.log10(kH), q)
    m = 0
    if kind == "reduced_indef":
        m = max(1, round(q / 10))
        b[-m:] *= -1
    E = 0.1 * rng.standard_normal((p, q))
    H = V @ (a[:, None] * V.T) + Z @ (b[:, None] * Z.T) + V @ E @ Z.T + Z @ E.T @ V.T
    H = (H + H.T) / 2
    x = rng.standard_normal(n)
    lam = rng.standard_normal(p)
    return H, C, H @ x + C.T @ lam, C @ x, x, lam, (n - m, p + m, 0)


report = []
rng = np.random.default_rng(11)
for kind in ["reduced_pd", "reduced_indef"]:
    H, C, g, d, xr, lr, expected = indefinite_case(200, 100, rng, kind)
    e = np.linalg.eigvalsh(np.block([[H, C.T], [C, np.zeros((100, 100))]]))
    true = (int((e > 0).sum()), int((e < 0).sum()), int((abs(e) < 1e-10).sum()))
    entry = dict(kind=kind, expected_inertia=expected, eigenvalue_inertia=true)
    try:
        x, lam, ine = solve_ldl(H, C, g, d)
        entry["ldl"] = dict(inertia=ine, x_error=float(la.norm(x - xr) / la.norm(xr)))
    except ValueError as err:
        entry["ldl"] = dict(error=str(err))
    try:
        x, lam, ynorm = solve_vr(H, C, g, d)
        entry["vr"] = dict(y_norm=ynorm, x_error=float(la.norm(x - xr) / la.norm(xr)))
    except la.LinAlgError as err:
        entry["vr"] = dict(error=str(err))
    try:
        la.cholesky(H)
        entry["cholesky_H"] = "succeeded"
    except la.LinAlgError:
        entry["cholesky_H"] = "failed"
    report.append(entry)

for kind in ["reduced_pd", "reduced_indef"]:
    for seed in [11, 22, 33]:
        H, C, g, d, xr, lr, expected = indefinite_case(4000, 2000, np.random.default_rng(seed), kind)
        entry = dict(kind=kind, n=4000, seed=seed)
        try:
            x, lam, ynorm = solve_vr(H, C, g, d)
            entry["vr"] = dict(y_norm=ynorm, x_error=float(la.norm(x - xr) / la.norm(xr)))
        except la.LinAlgError as err:
            entry["vr"] = dict(error=str(err))
        report.append(entry)

n, p = 300, 30
U, _ = la.qr(rng.standard_normal((n, n)))
H = (U * np.logspace(0, -4, n)) @ U.T
H = (H + H.T) / 2
V, _ = la.qr(rng.standard_normal((n, p)), mode="economic")
W, _ = la.qr(rng.standard_normal((p, p)))
C = (W * np.logspace(0, -1, p)) @ V.T
xr = rng.standard_normal(n)
lr = rng.standard_normal(p)
g = H @ xr + C.T @ lr
d = C @ xr
x, lam, ine = solve_ldl(H, C, g, d)
report.append(dict(problem="dense_kH1e4", method="ldl", inertia=ine, x_error=float(la.norm(x - xr) / la.norm(xr))))
x, lam, ynorm = solve_vr(H, C, g, d)
report.append(dict(problem="dense_kH1e4", method="vr", y_norm=ynorm, x_error=float(la.norm(x - xr) / la.norm(xr)), lambda_error=float(la.norm(lam - lr) / la.norm(lr)), feasibility=float(la.norm(C @ x - d))))

root = Path(__file__).resolve().parents[2] / "matlab"
(root / "provenance" / "new_methods_python_validation.json").write_text(json.dumps(report, indent=2))
for r in report:
    print(r)

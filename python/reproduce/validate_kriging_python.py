from pathlib import Path
import numpy as np, json
from scipy.io import loadmat
from scipy.linalg import qr, cho_factor, cho_solve, solve

root = Path(__file__).resolve().parents[2] / "matlab"
D = loadmat(root / "data" / "road_elevation.mat")
Xall = D["xy_km"]
report = []
for n in [40, 200]:
    ids = np.random.default_rng(11).choice(len(Xall), n + 1, replace=False)
    X = Xall[ids[:-1]]
    s = Xall[ids[-1]]
    H = np.exp(-np.sum((X[:, None] - X[None, :]) ** 2, axis=2) / 200) + 1e-3 * np.eye(n)
    C = np.ones((1, n))
    g = np.exp(-np.sum((X - s) ** 2, axis=1) / 200)
    Q, R = qr(C.T, mode="full")
    Z = Q[:, 1:]
    Q1 = Q[:, :1]
    r = R[:1]
    x0 = (Q1 @ solve(r.T, np.ones(1))).ravel()
    G = Z.T @ H @ Z
    xq = x0 + Z @ cho_solve(cho_factor((G + G.T) / 2), Z.T @ (g - H @ x0))
    lq = solve(r, Q1.T @ (g - H @ xq))
    ch = cho_factor(H)
    u = cho_solve(ch, g)
    W = cho_solve(ch, C.T)
    ls = solve(C @ W, C @ u - 1)
    xs = u - W @ ls
    K = np.block([[H, C.T], [C, np.zeros((1, 1))]])
    b = np.r_[g, 1]
    ref = solve(K, b)
    errors = [float(np.linalg.norm(x - ref[:-1]) / np.linalg.norm(ref[:-1])) for x in [xq, xs]]
    assert max(errors) < 1e-8
    assert abs(xq.sum() - 1) < 1e-10
    assert abs(xs.sum() - 1) < 1e-10
    report.append(
        dict(n=n, relative_disagreement_with_double_KKT=errors, cond_H=float(np.linalg.cond(H)))
    )
(root / "provenance" / "kriging_python_validation.json").write_text(
    json.dumps(
        dict(
            note="Independent small SciPy mirror using real packaged data; not MATLAB execution, not timing data. NumPy seeds do not select MATLAB-identical rows.",
            cases=report,
        ),
        indent=2,
    )
)
print(report)

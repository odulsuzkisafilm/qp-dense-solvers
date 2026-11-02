import os, sys, json, hashlib
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.io import savemat
import pycutest

root = Path(__file__).resolve().parents[2] / "matlab"
out = root
name = sys.argv[1]
par = json.loads(sys.argv[2])
tag = name + "_" + "_".join(f"{k}{v}" for k, v in par.items())
p = pycutest.import_problem(name, sifParams=par, sifOptions=["-force"], drop_fixed_variables=True)
n, pn = p.n, p.m
assert 0 < pn < n and np.all(p.is_eq_cons) and np.all(p.is_linear_cons)
assert np.all(p.bl <= -1e19) and np.all(p.bu >= 1e19)
z = np.zeros(n)
H = p.isphess(z).tocsc()
c0, C = p.scons(z, gradient=True)
C = C.tocsc()
g = -p.grad(z)
d = p.cl - c0
f0 = float(p.obj(z))
assert sp.linalg.norm(H - H.T) < 1e-14

assert (H - sp.diags(H.diagonal())).nnz == 0 and H.diagonal().min() > 0
rng = np.random.default_rng(19)
err = []
for k in range(2):
    x = rng.standard_normal(n)
    f = float(p.obj(x))
    grad = p.grad(x)
    c = p.cons(x)
    fq = 0.5 * x @ (H @ x) - g @ x + f0
    err.append(
        max(
            abs(f - fq) / (1 + abs(f)),
            np.linalg.norm(grad - (H @ x - g)) / (1 + np.linalg.norm(grad)),
            np.linalg.norm(c - (C @ x + c0)) / (1 + np.linalg.norm(c)),
        )
    )
assert max(err) < 1e-12
K = sp.bmat([[H, C.T], [C, None]], format="csc")
b = np.r_[g, d]
v = spsolve(K, b)
r = K @ v - b
eta = np.linalg.norm(r) / (sp.linalg.norm(K) * np.linalg.norm(v) + np.linalg.norm(b))

W = sp.diags(1 / H.diagonal()) @ C.T
ls = spsolve(C @ W, C @ (g / H.diagonal()) - d)
xs = g / H.diagonal() - W @ ls
agreement = np.linalg.norm(v[:n] - xs) / max(np.linalg.norm(v[:n]), np.finfo(float).tiny)
source = "https://raw.githubusercontent.com/ralna/SIF/master/" + name + ".SIF"
sha = hashlib.sha256((root / "provenance" / f"{name}.SIF").read_bytes()).hexdigest()
savemat(
    out / "data" / f"{tag}.mat",
    dict(
        H=H,
        C=C,
        g=g[:, None],
        d=d[:, None],
        x_ref=v[:n, None],
        lambda_ref=v[n:, None],
        f_constant=f0,
        ref_kind="double_sparse_KKT",
        ref_eta_fro=eta,
        ref_crosscheck=agreement,
        case_id=tag,
        family=name,
        category="CUTEst_original",
        seed=0,
        kH_target=float(H.diagonal().max() / H.diagonal().min()),
        kC_target=np.nan,
        kH_measured=float(H.diagonal().max() / H.diagonal().min()),
        kC_measured=np.nan,
        source_url=source,
        sif_sha256=sha,
        n=n,
        p=pn,
    ),
    do_compression=True,
)
record = dict(
    case_id=tag,
    n=n,
    p=pn,
    parameters=par,
    n_fixed_eliminated=int(p.n_fixed),
    objective_constant=f0,
    assembly_errors=err,
    reference_eta_fro=eta,
    reference_crosscheck=agreement,
    source_url=source,
    sif_sha256=sha,
)
(out / "provenance" / f"{tag}.json").write_text(json.dumps(record, indent=2))
print(json.dumps(record), flush=True)

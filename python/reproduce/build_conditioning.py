from pathlib import Path
import numpy as np, scipy.linalg as la
from scipy.io import savemat, loadmat
import mpmath as mp
import json, time

root = Path(__file__).resolve().parents[2] / "matlab"
out = root


def reference(H, C, g, d, dps):
    with mp.workdps(dps):
        n = len(g)
        p = len(d)
        K = mp.matrix(np.block([[H, C.T], [C, np.zeros((p, p))]]).tolist())
        b = mp.matrix(np.r_[g, d].tolist())
        v = mp.lu_solve(K, b)
        eta = mp.norm(K * v - b) / (mp.norm(K) * mp.norm(v) + mp.norm(b))
        return [mp.mpf(t) for t in v], float(eta)


records = []
for family, targets in [("H", [0, 4, 8, 10, 12]), ("C", [0, 4, 6, 8, 10, 12])]:
    for seed in [11, 22, 33, 44, 55]:
        rng = np.random.default_rng(seed)
        n = 80
        p = 20
        U, _ = la.qr(rng.standard_normal((n, n)))
        V, _ = la.qr(rng.standard_normal((n, p)), mode="economic")
        W, _ = la.qr(rng.standard_normal((p, p)))
        x = rng.standard_normal(n)
        lam = rng.standard_normal(p)
        for exponent in targets:
            kh = 10.0**exponent if family == "H" else 100.0
            kc = 10.0**exponent if family == "C" else 10.0
            H = (U * np.geomspace(1, 1 / kh, n)) @ U.T
            H = (H + H.T) / 2
            C = (W * np.geomspace(1, 1 / kc, p)) @ V.T
            g = H @ x + C.T @ lam
            d = C @ x
            tag = f"COND_{family}_log{exponent:02d}_seed{seed}"
            begin = time.time()
            v64, eta = reference(H, C, g, d, 64)
            v96, eta96 = reference(H, C, g, d, 96)
            with mp.workdps(100):
                delta = float(mp.norm(mp.matrix(v64) - mp.matrix(v96)) / mp.norm(mp.matrix(v96)))
            assert delta < 1e-30, (tag, delta)
            v = np.array([float(t) for t in v96])
            savemat(
                out / "data" / f"{tag}.mat",
                dict(
                    H=H,
                    C=C,
                    g=g[:, None],
                    d=d[:, None],
                    x_ref=v[:n, None],
                    lambda_ref=v[n:, None],
                    f_constant=0,
                    ref_kind="stored_input_KKT_96digits_checked64",
                    ref_precision_agreement=delta,
                    ref_eta_mp=eta96,
                    case_id=tag,
                    family=f"COND_{family}",
                    category="synthetic_conditioning",
                    seed=seed,
                    kH_target=kh,
                    kC_target=kc,
                    kH_measured=np.linalg.cond(H),
                    kC_measured=np.linalg.cond(C),
                    x_manufactured=x[:, None],
                    lambda_manufactured=lam[:, None],
                ),
                do_compression=True,
            )
            records.append(
                dict(
                    case_id=tag, precision_agreement=delta, eta96=eta96, seconds=time.time() - begin
                )
            )
            (out / "provenance/conditioning_reference_validation.json").write_text(
                json.dumps(records, indent=2)
            )
            print(tag, "agreement", delta, "seconds", round(time.time() - begin, 2), flush=True)

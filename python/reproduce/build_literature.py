from pathlib import Path
import numpy as np, scipy.sparse as sp, json, hashlib
from scipy.io import loadmat, savemat
from scipy.linalg import cholesky
from scipy.sparse.linalg import spsolve

root = Path(__file__).resolve().parents[2] / "matlab"
records = []
for name, path in [
    ("gemat1", "provenance/gemat1_source.mat"),
    ("deter3", "provenance/deter3_source.mat"),
]:
    P = loadmat(root / path, struct_as_record=False, squeeze_me=True)["Problem"]
    B = P.A.T.tocsr()
    B.eliminate_zeros()
    counts = np.diff(B.indptr)
    order = np.lexsort((np.arange(B.shape[0]), -counts))
    sel = order[:20]
    rest = np.sort(order[20:])
    keep = np.asarray(abs(B[rest]).sum(axis=0)).ravel() > 0
    B = B[:, keep]
    scale = 1 / np.sqrt(np.asarray(B.power(2).sum(axis=0)).ravel())
    B = (B @ sp.diags(scale)).tocsr()
    A = B[rest]
    n = B.shape[1]
    H = (A.T @ A).tocsc()
    H = (H + H.T) * 0.5
    record = dict(
        name=name,
        n=n,
        selected_constraint_rows_1based=(sel + 1).tolist(),
        removed_null_columns=int((~keep).sum()),
        source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        rule="Inspired by Scott-Tuma: 20 densest numerical rows removed from A, descending nnz ties by original index; normalize full extended columns; b=d=1. Tie policy and explicit-zero handling may differ from paper. Not a claimed exact reproduction.",
    )
    print("Checking", name, n, flush=True)
    try:
        L = cholesky(H.toarray(), lower=True, overwrite_a=True, check_finite=False)
        del L
        record["stored_H_cholesky"] = "passed"
    except Exception as e:
        record["stored_H_cholesky"] = str(e)
        records.append(record)
        continue
    for p in [5, 20]:
        C = B[sel[:p]].tocsc()
        sv = np.linalg.svd(C.toarray(), compute_uv=False)
        assert sv[-1] > sv[0] * 1e-12
        g = np.asarray(A.T @ np.ones(A.shape[0])).ravel()
        d = np.ones(p)
        K = sp.bmat([[H, C.T], [C, None]], format="csc")
        v = spsolve(K, np.r_[g, d])
        assert np.isfinite(v).all()
        eta = np.linalg.norm(K @ v - np.r_[g, d]) / (
            sp.linalg.norm(K) * np.linalg.norm(v) + np.linalg.norm(np.r_[g, d])
        )
        tag = f"LSE_{name}_p{p}"
        savemat(
            root / "data" / f"{tag}.mat",
            dict(
                H=H,
                C=C,
                g=g[:, None],
                d=d[:, None],
                x_ref=v[:n, None],
                lambda_ref=v[n:, None],
                f_constant=0.5 * A.shape[0],
                ref_kind="double_sparse_normal_KKT_not_certified",
                case_id=tag,
                family="LSE_" + name,
                category="paper_inspired_application_LSE",
                seed=0,
                kH_target=np.nan,
                kC_target=np.nan,
                kH_measured=np.nan,
                kC_measured=sv[0] / sv[-1],
                A_ls=A,
                b_ls=np.ones((A.shape[0], 1)),
                source_rows=(sel[:p] + 1)[:, None],
            ),
            do_compression=True,
        )
        record[f"p{p}"] = dict(
            rank_C=p,
            cond_C=float(sv[0] / sv[-1]),
            reference_residual=float(eta),
            H_density=H.nnz / n**2,
            C_density=C.nnz / (p * n),
        )
    records.append(record)
    (root / "provenance/literature_lse_audit.json").write_text(json.dumps(records, indent=2))
    print(record, flush=True)
(root / "provenance/literature_lse_audit.json").write_text(json.dumps(records, indent=2))

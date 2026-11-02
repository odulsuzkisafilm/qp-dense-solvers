from scipy.io import loadmat, savemat
import scipy.sparse as sp, numpy as np
from scipy.sparse.linalg import spsolve
from pathlib import Path

root = Path(__file__).resolve().parents[2] / "matlab"
for path in (root / "data").glob("LSE*.mat"):
    D = loadmat(path)
    A = D["A_ls"]
    C = D["C"]
    m, n = A.shape
    p = C.shape[0]
    b = D["b_ls"].ravel()
    d = D["d"].ravel()

    K = sp.bmat([[sp.eye(m), A, None], [A.T, None, -C.T], [None, C, None]], format="csc")
    v = spsolve(K, np.r_[b, np.zeros(n), d])
    x = v[m : m + n]
    l = v[m + n :]
    delta = np.linalg.norm(x - D["x_ref"].ravel()) / np.linalg.norm(x)
    print(path.name, "normal vs augmented", delta, flush=True)
    D["normal_vs_augmented_discrepancy"] = delta
    D["x_ref"] = x[:, None]
    D["lambda_ref"] = l[:, None]
    D["ref_kind"] = "double_augmented_LS_not_certified"
    D = {k: v for k, v in D.items() if not k.startswith("__")}
    savemat(path, D, do_compression=True)

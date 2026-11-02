import sys, glob, os, json, hashlib
import numpy as np, pandas as pd

root = sys.argv[1]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
prov = {}


def newest(kind):
    f = sorted(glob.glob(os.path.join(root, kind + "_*")))
    return f[-1] if f else None


def load(kind):
    d = newest(kind)
    if d is None:
        return None
    T = pd.read_csv(os.path.join(d, "results.csv"))
    prov[kind] = {
        "folder": os.path.basename(d),
        "sha256": hashlib.sha256(open(os.path.join(d, "results.csv"), "rb").read()).hexdigest(),
        "rows": int(len(T)),
        "status": T.status.value_counts().to_dict(),
    }
    return T


def seed_medians(T, keys):
    ok = T[T.status == "ok"]
    c = (
        ok.groupby(keys + ["seed", "method"])
        .agg(
            t=("seconds", "median"),
            xerr=("x_reference_error", "median"),
            lerr=("lambda_reference_error", "median"),
            kkt=("kkt_fro_scaled", "median"),
        )
        .reset_index()
    )
    return (
        c.groupby(keys + ["method"])
        .agg(
            y=("t", "median"),
            lo=("t", "min"),
            hi=("t", "max"),
            seeds=("seed", "nunique"),
            xerr=("xerr", "median"),
            lerr=("lerr", "median"),
            kkt=("kkt", "median"),
        )
        .reset_index()
    )


def write_plot(S, xcol, name):
    for m, tag in [("QR", "qr"), ("Schur", "schur")]:
        A = S[S.method == m].sort_values(xcol)
        pd.DataFrame({"x": A[xcol], "y": A.y, "minus": A.y - A.lo, "plus": A.hi - A.y}).to_csv(
            os.path.join(out, f"{name}_{tag}.csv"), index=False
        )


nums = {}

T = load("final_dense_dimension")
if T is not None:
    S = seed_medians(T, ["n"])
    write_plot(S, "n", "v2_size")
    q = S[S.method == "QR"].sort_values("n")
    s = S[S.method == "Schur"].sort_values("n")
    nums["SizeQRTop"] = f"{q.y.iloc[-1]:.1f}"
    nums["SizeSchurTop"] = f"{s.y.iloc[-1]:.2f}"
    nums["SizeNTop"] = f"{int(q.n.iloc[-1])}"
    nums["SizeRatioTop"] = f"{q.y.iloc[-1]/s.y.iloc[-1]:.2f}"
    nums["SizeRatioLow"] = f"{q.y.iloc[0]/s.y.iloc[0]:.2f}"
    for m, A in [("QR", q), ("Schur", s)]:
        sl = np.polyfit(np.log(A.n), np.log(A.y), 1)[0]
        dy = A[A.n <= 8000]
        sl2 = np.polyfit(np.log(dy.n), np.log(dy.y), 1)[0]
        nums[f"Slope{m}"] = f"{sl:.2f}"
        nums[f"SlopeDyadic{m}"] = f"{sl2:.2f}"
        r = dy.y.values[1:] / dy.y.values[:-1]
        nums[f"Doubling{m}"] = ", ".join(f"{v:.1f}" for v in r)
        nums[f"XerrSize{m}"] = f"{A.xerr.median():.0e}"
    print(S.to_string())

T = load("final_dense_illcond")
if T is not None:
    S = seed_medians(T, ["n"])
    write_plot(S, "n", "v2_illcond")
    for m in ["QR", "Schur"]:
        A = S[S.method == m].sort_values("n")
        nums[f"Ill{m}XLow"] = f"{A.xerr.iloc[0]:.0e}"
        nums[f"Ill{m}XTop"] = f"{A.xerr.iloc[-1]:.0e}"
        nums[f"Ill{m}LTop"] = f"{100*A.lerr.iloc[-1]:.1f}"
        nums[f"Ill{m}KKTTop"] = f"{A.kkt.iloc[-1]:.0e}"
        nums[f"Ill{m}TTop"] = f"{A.y.iloc[-1]:.1f}"
    ok = T[T.status == "ok"].iloc[0]
    nums["IllKH"] = f"{ok.kH_target:.0e}"
    nums["IllKC"] = f"{ok.kC_target:.0e}"
    print(S.to_string())
T = load("final_dense_dimension")
if T is not None:
    S = seed_medians(T, ["n"])
    for m in ["QR", "Schur"]:
        A = S[S.method == m].sort_values("n")
        nums[f"Well{m}XLow"] = f"{A.xerr.iloc[0]:.0e}"
        nums[f"Well{m}XTop"] = f"{A.xerr.iloc[-1]:.0e}"

T = load("final_dense_fraction")
if T is not None:
    S = seed_medians(T, ["ratio"])
    write_plot(S, "ratio", "v2_fraction")
    q = S[S.method == "QR"].sort_values("ratio")
    s = S[S.method == "Schur"].sort_values("ratio")
    nums["FracN"] = f"{int(T.n.iloc[0])}"
    nums["FracRatioLow"] = f"{q.y.iloc[0]/s.y.iloc[0]:.2f}"
    nums["FracRatioHigh"] = f"{q.y.iloc[-1]/s.y.iloc[-1]:.2f}"
    nums["FracYmax"] = f"{np.ceil(max(q.hi.max(), s.hi.max())):.0f}"
    print(S.to_string())

T = load("literature_full")
if T is not None:
    T = T.copy()
    T["kind"] = T.case_id.str.replace(r"_seed\d+$", "", regex=True)
    S = seed_medians(T, ["kind", "n", "p"])
    S.to_csv(os.path.join(out, "v2_literature_summary.csv"), index=False)
    print(S.to_string())

    def row(kind, m):
        return S[(S.kind == kind) & (S.method == m)].iloc[0]

    for kind, tag in [
        ("DTOC3_N4001", "DtocNative"),
        ("DTOC3_N4001_orthogonal", "DtocOrth"),
        ("LSE_gemat1_p5", "GematFive"),
        ("LSE_gemat1_p20", "GematTwenty"),
    ]:
        for m in ["QR", "Schur"]:
            r = row(kind, m)
            nums[f"{tag}{m}T"] = f"{r.y:.2f}" if r.y >= 1 else f"{r.y:.3f}"
            nums[f"{tag}{m}Lo"] = f"{r.lo:.2f}" if r.lo >= 1 else f"{r.lo:.3f}"
            nums[f"{tag}{m}Hi"] = f"{r.hi:.2f}" if r.hi >= 1 else f"{r.hi:.3f}"
            nums[f"{tag}{m}X"] = f"{r.xerr:.1e}"

T = load("kriging_full")
if T is not None:
    S = seed_medians(T, ["n"])
    write_plot(S, "n", "v2_kriging")
    q = S[S.method == "QR"].sort_values("n")
    s = S[S.method == "Schur"].sort_values("n")
    nums["KrigQRTop"] = f"{q.y.iloc[-1]:.1f}"
    nums["KrigSchurTop"] = f"{s.y.iloc[-1]:.2f}"
    nums["KrigRatioTop"] = f"{q.y.iloc[-1]/s.y.iloc[-1]:.1f}"
    nums["KrigNTop"] = f"{int(q.n.iloc[-1])}"
    ok = T[T.status == "ok"]
    nums["KrigWsum"] = f"{ok.weight_sum_error.max():.0e}"
    nums["KrigKKT"] = f"{ok.kkt_fro_scaled.median():.0e}"
    try:
        import scipy.io as sio

        d = newest("kriging_full")
        nt = int(q.n.iloc[-1])
        xq = sio.loadmat(os.path.join(d, "solutions", f"KRIGING_n{nt}_seed11_QR_r1.mat"))[
            "x"
        ].ravel()
        xs = sio.loadmat(os.path.join(d, "solutions", f"KRIGING_n{nt}_seed11_Schur_r1.mat"))[
            "x"
        ].ravel()
        nums["KrigAgree"] = f"{np.linalg.norm(xq-xs)/np.linalg.norm(xq):.0e}"
    except Exception as e:
        print("agreement skipped:", e)
    print(S.to_string())


def texnum(v):
    if "e" in v:
        a, b = v.split("e")
        return f"{a}\\times10^{{{int(b)}}}"
    return v


with open(os.path.join(out, "numbers_v2.tex"), "w") as f:
    for k, v in nums.items():
        f.write(f"\\newcommand{{\\{k}}}{{{texnum(v)}}}\n")
json.dump(prov, open(os.path.join(out, "provenance_v2.json"), "w"), indent=1)
print("\n".join(f"{k} = {v}" for k, v in nums.items()))

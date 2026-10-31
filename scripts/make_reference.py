#!/usr/bin/env python3
"""Generate the synthetic example datasets and the reference values used by the unit tests.

All data produced here are SYNTHETIC (simulated with a fixed random seed); they do not
come from any real experiment.

Reference values are computed independently of the TypeScript code with SciPy
(scipy.stats.ttest_ind, f_oneway, tukey_hsd, dunnett, studentized_range,
scipy.optimize.curve_fit) and plain NumPy for the ddCt arithmetic.

Usage:
    pip install numpy scipy
    python3 scripts/make_reference.py

Outputs:
    examples/*.txt|csv|tsv|json   synthetic example datasets (loaded by the app)
    tests/fixtures/reference.json reference values for the vitest suite
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import scipy
from scipy import integrate, optimize, stats

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "examples"
FIX = ROOT / "tests" / "fixtures"
EX.mkdir(exist_ok=True)
FIX.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261031)


def clean(o):
    """Convert numpy types and NaN to JSON-friendly values (NaN -> None)."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    return o


ref: dict = {"generator": "scripts/make_reference.py", "scipy": scipy.__version__, "numpy": np.__version__}

# ---------------------------------------------------------------- distributions
ref["dist"] = {
    "t_two_sided": [
        {"t": t, "df": df, "p": 2 * stats.t.sf(abs(t), df)}
        for t in (0.3, 1.0, 2.1, 3.7, 6.5) for df in (2, 5.5, 13, 40)
    ],
    "f_upper": [
        {"f": f, "df1": d1, "df2": d2, "p": stats.f.sf(f, d1, d2)}
        for f in (0.5, 2.0, 4.5, 12.0) for d1, d2 in ((1, 6), (2, 9), (3, 20), (5, 44))
    ],
    "t_quantile": [{"p": 0.975, "df": df, "q": stats.t.ppf(0.975, df)} for df in (1, 2, 3, 5, 8, 14, 30, 100)],
    "tukey_cdf": [
        {"q": q, "k": k, "df": df, "cdf": stats.studentized_range.cdf(q, k, df)}
        for q in (1.5, 3.0, 4.2, 6.0) for k in (3, 4, 6) for df in (6, 12, 30)
    ],
}

# ---------------------------------------------------------------- tests on random data
welch_cases = []
for n1, n2, shift, s2 in ((4, 4, 1.2, 1.0), (3, 6, 0.4, 2.5), (8, 5, -2.0, 0.5)):
    a = rng.normal(0, 1, n1)
    b = rng.normal(shift, s2, n2)
    r = stats.ttest_ind(a, b, equal_var=False)
    welch_cases.append({"a": a, "b": b, "t": r.statistic, "df": r.df, "p": r.pvalue})
ref["welch"] = welch_cases


def dunnett_exact(c: float, lam: np.ndarray, df: float) -> float:
    """P(max|T_j| <= c), product-correlation multivariate t, by nested adaptive quadrature."""
    sc = np.sqrt(1 - lam ** 2)

    def inner(s):
        f = lambda z: np.prod(stats.norm.cdf((c * s - lam * z) / sc) - stats.norm.cdf((-c * s - lam * z) / sc)) * stats.norm.pdf(z)
        return integrate.quad(f, -9, 9, epsabs=1e-13, epsrel=1e-12, limit=200)[0]

    dens = lambda s: stats.chi.pdf(s * math.sqrt(df), df) * math.sqrt(df)
    lo, hi = stats.chi.ppf([1e-15, 1 - 1e-15], df) / math.sqrt(df)
    return integrate.quad(lambda s: dens(s) * inner(s), lo, hi, epsabs=1e-13, epsrel=1e-11, limit=200)[0]


multi_cases = []
for sizes, shifts in (((5, 4, 6, 3), (0, 0.8, 1.6, -0.5)), ((4, 4, 4), (0, 1.5, 0.3)), ((6, 5, 5, 4, 6), (0, 0.2, 1.1, 2.0, -0.7))):
    groups = [rng.normal(m, 0.8, n) for n, m in zip(sizes, shifts)]
    an = stats.f_oneway(*groups)
    tk = stats.tukey_hsd(*groups)
    pairs = []
    k = len(groups)
    N = sum(sizes)
    mse = sum(((g - g.mean()) ** 2).sum() for g in groups) / (N - k)
    for i in range(k):
        for j in range(i + 1, k):
            se = math.sqrt(mse / 2 * (1 / sizes[i] + 1 / sizes[j]))
            pairs.append({"i": i, "j": j, "diff": groups[i].mean() - groups[j].mean(), "q": abs(tk.statistic[i, j]) / se, "p": tk.pvalue[i, j]})
    dn = stats.dunnett(*groups[1:], control=groups[0], rng=np.random.default_rng(1))
    lam = np.array([math.sqrt(n / (n + sizes[0])) for n in sizes[1:]])
    exact = [1 - dunnett_exact(abs(t), lam, N - k) for t in dn.statistic]
    multi_cases.append({
        "groups": groups,
        "anova": {"f": an.statistic, "p": an.pvalue, "df1": k - 1, "df2": N - k},
        "tukey": pairs,
        "dunnett": [{"t": t, "p_scipy": p, "p_exact": pe} for t, p, pe in zip(dn.statistic, dn.pvalue, exact)],
    })
ref["multi"] = multi_cases

# ---------------------------------------------------------------- qPCR synthetic data
GENES = {"GAPDH": 18.0, "ACTB": 17.2, "IL6": 27.0, "BCL2": 24.1}
EFFECT = {"IL6": {"Ctrl": 0, "TrtA": -2.0, "TrtB": -0.6}, "BCL2": {"Ctrl": 0, "TrtA": 1.1, "TrtB": 1.6}}
EFF = {"GAPDH": 98.0, "ACTB": 101.0, "IL6": 94.0, "BCL2": 91.0}
ROWS16 = "ABCDEFGHIJKLMNOP"

wells = []  # (pos, name, ct or None)
samples = [f"{g}{i}" for g in ("Ctrl", "TrtA", "TrtB") for i in range(1, 5)]
idx = 0
for gi, gene in enumerate(GENES):
    for s in samples:
        grp = re.sub(r"\d+$", "", s)
        bio = rng.normal(0, 0.25)
        base = GENES[gene] + EFFECT.get(gene, {}).get(grp, 0) + bio
        for rep in range(3):
            ct = base + rng.normal(0, 0.07)
            if s == "TrtB2" and gene == "IL6" and rep == 2:
                ct += 1.1  # technical outlier for the QC demo
            pos = f"{ROWS16[idx // 24]}{idx % 24 + 1}"
            idx += 1
            if s == "TrtA3" and gene == "BCL2" and rep == 1:
                wells.append((pos, f"{s}_{gene}", None))  # no amplification
            else:
                wells.append((pos, f"{s}_{gene}", round(ct, 2)))
    pos = f"{ROWS16[idx // 24]}{idx % 24 + 1}"
    idx += 1
    wells.append((pos, f"NTC_{gene}", 36.84 if gene == "GAPDH" else None))

# sample-level loading offset shared by all genes of one sample
offsets = {s: rng.normal(0, 0.3) for s in samples}
wells = [(p, n, None if ct is None else round(ct + offsets.get(n.split("_")[0], 0), 2)) for p, n, ct in wells]

lc_lines = ["Experiment: SYNTHETIC example generated by labcalc  Selected Filter: SYBR Green I / HRM Dye (465-510)",
            "Include\tColor\tPos\tName\tCp\tConcentration\tStandard\tStatus"]
for p, n, ct in wells:
    lc_lines.append(f"True\t255\t{p}\t{n}\t{'' if ct is None else f'{ct:.2f}'}\t\t0\t")
(EX / "qpcr_lc480_synthetic.txt").write_text("\r\n".join(lc_lines) + "\r\n", encoding="utf-8")

# Two-group tidy CSV (two reference genes)
csv_rows = ["Well,Sample,Group,Gene,Ct"]
two = []
w = 0
for gene, base in (("RPLP0", 19.4), ("HPRT1", 23.0), ("MYC", 25.5)):
    for grp, eff in (("Control", 0.0), ("Knockdown", 1.4)):
        for i in range(1, 5):
            s = f"{grp[:2].upper()}-{i}"
            mu = base + (eff if gene == "MYC" else 0) + rng.normal(0, 0.3)
            for rep in range(3):
                ct = round(mu + rng.normal(0, 0.06), 2)
                pos = f"{'ABCDEFGH'[w // 12]}{w % 12 + 1}"
                w += 1
                ctxt = "Undetermined" if (s == "KN-4" and gene == "HPRT1" and rep == 0) else f"{ct:.2f}"
                csv_rows.append(f"{pos},{s},{grp},{gene},{ctxt}")
                two.append((pos, s, grp, gene, None if ctxt == "Undetermined" else ct))
for gene in ("RPLP0", "HPRT1", "MYC"):
    pos = f"{'ABCDEFGH'[w // 12]}{w % 12 + 1}"
    w += 1
    csv_rows.append(f"{pos},NTC,,{gene},Undetermined")
(EX / "qpcr_two_groups_synthetic.csv").write_text("\n".join(csv_rows) + "\n", encoding="utf-8")


def ddct_reference(rows, refs, targets, groups_of, control, eff=None, multi="tukey"):
    """Independent ddCt implementation. rows: (sample, gene, ct) with NTC and missing Ct already removed."""
    reps = {}
    for s, g, ct in rows:
        reps.setdefault((s, g), []).append(ct)
    mean_ct = {k: float(np.mean(v)) for k, v in reps.items()}
    samples_ = list(dict.fromkeys(s for s, _, _ in rows))
    per = []
    for s in samples_:
        if any((s, r) not in mean_ct for r in refs):
            continue
        ref_ct = float(np.mean([mean_ct[(s, r)] for r in refs]))
        for t in targets:
            if (s, t) not in mean_ct:
                continue
            d = {"sample": s, "group": groups_of[s], "gene": t, "ctTarget": mean_ct[(s, t)], "refCt": ref_ct, "dCt": mean_ct[(s, t)] - ref_ct}
            if eff:
                d["dCqE"] = math.log2(1 + eff[t] / 100) * mean_ct[(s, t)] - np.mean([math.log2(1 + eff[r] / 100) * mean_ct[(s, r)] for r in refs])
            per.append(d)
    out_groups, out_stats = [], []
    order = [control] + [g for g in dict.fromkeys(groups_of[s] for s in samples_) if g != control]
    for t in targets:
        rows_t = [d for d in per if d["gene"] == t]
        cm = np.mean([d["dCt"] for d in rows_t if d["group"] == control])
        for d in rows_t:
            d["ddCt"] = d["dCt"] - cm
            d["fc"] = 2 ** -d["ddCt"]
        if eff:
            # Pfaffl ratio written in its original form: E_t^dCt_t / geomean_r(E_r^dCt_r), dCt = control mean - sample
            ctrl = [d["sample"] for d in rows_t if d["group"] == control]
            mc = {g: np.mean([mean_ct[(s, g)] for s in ctrl]) for g in [t] + refs}
            for d in rows_t:
                s = d["sample"]
                num = (1 + eff[t] / 100) ** (mc[t] - mean_ct[(s, t)])
                den = np.prod([(1 + eff[r] / 100) ** (mc[r] - mean_ct[(s, r)]) for r in refs]) ** (1 / len(refs))
                d["ratio"] = num / den
        for g in order:
            v = [d for d in rows_t if d["group"] == g]
            if not v:
                continue
            arr = lambda key: np.array([d[key] for d in v])
            entry = {"gene": t, "group": g, "n": len(v)}
            for key in ("dCt", "ddCt", "fc") + (("ratio",) if eff else ()):
                a = arr(key)
                entry[key] = {"mean": a.mean(), "sd": a.std(ddof=1), "sem": a.std(ddof=1) / math.sqrt(len(a))}
            entry["fcOfMean"] = 2 ** -arr("ddCt").mean()
            out_groups.append(entry)
        key = "dCqE" if eff else "dCt"
        data = [np.array([d[key] for d in rows_t if d["group"] == g]) for g in order]
        if len(order) == 2:
            r = stats.ttest_ind(data[1], data[0], equal_var=False)
            out_stats.append({"gene": t, "test": "welch", "a": order[1], "b": order[0], "diff": data[1].mean() - data[0].mean(), "stat": r.statistic, "df": r.df, "p": r.pvalue})
        else:
            an = stats.f_oneway(*data)
            out_stats.append({"gene": t, "test": "anova", "stat": an.statistic, "p": an.pvalue})
            if multi == "tukey":
                tk = stats.tukey_hsd(*data)
                for i in range(len(order)):
                    for j in range(i + 1, len(order)):
                        out_stats.append({"gene": t, "test": "tukey", "a": order[j], "b": order[i], "diff": data[j].mean() - data[i].mean(), "p": tk.pvalue[i, j]})
            else:
                dn = stats.dunnett(*data[1:], control=data[0], rng=np.random.default_rng(1))
                k, N = len(data), sum(len(x) for x in data)
                lam = np.array([math.sqrt(len(x) / (len(x) + len(data[0]))) for x in data[1:]])
                for j in range(1, len(order)):
                    pe = 1 - dunnett_exact(abs(dn.statistic[j - 1]), lam, N - k)
                    out_stats.append({"gene": t, "test": "dunnett", "a": order[j], "b": order[0], "diff": data[j].mean() - data[0].mean(), "stat": dn.statistic[j - 1], "p_scipy": dn.pvalue[j - 1], "p_exact": pe})
    return {"perSample": per, "groups": out_groups, "stats": out_stats}


NTC_RE = re.compile(r"^(ntc|h2o|water|nrt|neg)", re.I)
lc_rows = []
groups_of = {}
for p, n, ct in wells:
    s, g = n.rsplit("_", 1)
    if NTC_RE.match(s) or ct is None:
        continue
    lc_rows.append((s, g, ct))
    groups_of[s] = re.sub(r"\d+$", "", s)
ref["qpcr_lc480"] = {
    "file": "examples/qpcr_lc480_synthetic.txt",
    "refs": ["GAPDH", "ACTB"],
    "efficiency": EFF,
    "control": "Ctrl",
    "livak_tukey": ddct_reference(lc_rows, ["GAPDH", "ACTB"], ["IL6", "BCL2"], groups_of, "Ctrl"),
    "livak_dunnett": ddct_reference(lc_rows, ["GAPDH", "ACTB"], ["IL6", "BCL2"], groups_of, "Ctrl", multi="dunnett"),
    "pfaffl_tukey": ddct_reference(lc_rows, ["GAPDH", "ACTB"], ["IL6", "BCL2"], groups_of, "Ctrl", eff=EFF),
    "n_wells": len(wells),
}
two_rows = [(s, g, ct) for _, s, grp, g, ct in two if ct is not None]
ref["qpcr_two_groups"] = {
    "file": "examples/qpcr_two_groups_synthetic.csv",
    "refs": ["RPLP0", "HPRT1"],
    "control": "Control",
    "livak": ddct_reference(two_rows, ["RPLP0", "HPRT1"], ["MYC"], {s: grp for _, s, grp, _, _ in two}, "Control"),
}

# ---------------------------------------------------------------- viability synthetic plate
conc = [100 / 3 ** k for k in range(9)]  # uM, columns 3..11
truth = {"Cmpd-A": (5, 100, 2.0, 1.2), "Cmpd-B": (15, 100, 0.15, 0.9), "Cmpd-C": (8, 100, 300.0, 1.0)}
row_cmpd = ["Cmpd-A"] * 3 + ["Cmpd-B"] * 3 + ["Cmpd-C"] * 2
layout = []
A = np.zeros((8, 12))
R = np.zeros((8, 12))
for r in range(8):
    for c in range(12):
        pos = f"{'ABCDEFGH'[r]}{c + 1}"
        bg = 0.04 + rng.normal(0, 0.004)
        R[r, c] = bg
        if c == 0:
            layout.append({"pos": pos, "role": "blank"})
            signal = 0.0
        elif c == 1:
            layout.append({"pos": pos, "role": "control"})
            signal = 1.10 * rng.normal(1, 0.05)
        elif c <= 10:
            name = row_cmpd[r]
            x = conc[c - 2]
            b, t, ic, h = truth[name]
            v = b + (t - b) / (1 + (x / ic) ** h)
            layout.append({"pos": pos, "role": "compound", "compound": name, "conc": x})
            signal = 1.10 * v / 100 * rng.normal(1, 0.05)
        else:
            A[r, c] = np.nan
            R[r, c] = np.nan
            continue
        A[r, c] = bg + 0.05 + signal + rng.normal(0, 0.006)


def matrix_text(M):
    lines = ["\t" + "\t".join(str(i) for i in range(1, 13))]
    for r in range(8):
        lines.append("ABCDEFGH"[r] + "\t" + "\t".join("" if np.isnan(v) else f"{v:.4f}" for v in M[r]))
    return "\n".join(lines) + "\n"


(EX / "viability_570nm_synthetic.tsv").write_text(matrix_text(A), encoding="utf-8")
(EX / "viability_650nm_reference_synthetic.tsv").write_text(matrix_text(R), encoding="utf-8")
(EX / "viability_layout_synthetic.json").write_text(json.dumps({
    "format": "labcalc-plate-layout", "version": 1, "unit": "µM",
    "note": "SYNTHETIC example layout", "wells": layout}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")

# Reference: use the values as written to the files (4 decimals)
A4 = np.round(A, 4)
R4 = np.round(R, 4)
Cc = A4 - R4
lay = {w["pos"]: w for w in layout}
cell = lambda pos: Cc["ABCDEFGH".index(pos[0]), int(pos[1:]) - 1]
blank = np.array([cell(p) for p, w in lay.items() if w["role"] == "blank"])
ctrl = np.array([cell(p) for p, w in lay.items() if w["role"] == "control"]) - blank.mean()
viab = {p: 100 * (cell(p) - blank.mean()) / ctrl.mean() for p, w in lay.items()}


def f4(x, b, t, c, h):
    return b + (t - b) / (1 + 10 ** (h * (x - c)))


def fit_ref(xc, y, fix_bottom=None, fix_top=None):
    """Least-squares 4PL fit with SciPy (Levenberg-Marquardt), best of several starts.
    Parameters are reported in canonical form (hill > 0, top = low-concentration plateau)."""
    x = np.log10(xc)
    free = [fix_bottom is None, fix_top is None, True, True]

    def full(q):
        it = iter(q)
        return [next(it) if free[0] else fix_bottom, next(it) if free[1] else fix_top, next(it), next(it)]

    model = lambda x_, *q: f4(x_, *full(q))
    best = None
    for h0 in (1, 0.5, 2, 4):
        for s in (1, -1):
            for c0 in (np.median(x), x.min(), x.max()):
                p0 = [v for v, f in zip([y.min(), y.max(), c0, s * h0], free) if f]
                try:
                    p, cov = optimize.curve_fit(model, x, y, p0=p0, method="lm", maxfev=200000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
                except RuntimeError:
                    continue
                ssr = float(((y - model(x, *p)) ** 2).sum())
                if best is None or ssr < best[2]:
                    best = (p, cov, ssr)
    p, cov, ssr = best
    k = len(p)
    n, df = len(y), len(y) - k
    se_free = np.sqrt(np.diag(cov)) if df > 0 else np.full(k, np.nan)
    it = iter(se_free)
    se = [next(it) if f else 0.0 for f in free]
    b, t, c, h = full(p)
    if h < 0:  # canonical form: same curve with hill > 0
        b, t, h = t, b, -h
        se[0], se[1] = se[1], se[0]
    tq = stats.t.ppf(0.975, df) if df > 0 else np.nan
    absic = np.nan
    if min(b, t) < 50 < max(b, t):
        ratio = (t - b) / (50 - b) - 1
        if ratio > 0:
            absic = 10 ** (c + math.log10(ratio) / h)
    return {"bottom": b, "top": t, "logIC50": c, "hill": h, "se": se, "ic50": 10 ** c,
            "ic50CI": [10 ** (c - tq * se[2]), 10 ** (c + tq * se[2])], "r2": 1 - ssr / ((y - y.mean()) ** 2).sum(),
            "ssr": ssr, "absIC50": absic, "df": df}


fits = {}
for name in truth:
    pts = [(w["conc"], viab[p]) for p, w in lay.items() if w.get("compound") == name]
    xc = np.array([q[0] for q in pts])
    y = np.array([q[1] for q in pts])
    if name == "Cmpd-C":
        # Only the upper part of this curve is observed, so the free 4PL is not identifiable.
        fits[name + " (bottom fixed at 0)"] = fit_ref(xc, y, fix_bottom=0.0)
    else:
        fits[name] = fit_ref(xc, y)
        fits[name + " (top fixed at 100)"] = fit_ref(xc, y, fix_top=100.0)
ref["viability"] = {
    "layout": "examples/viability_layout_synthetic.json",
    "abs": "examples/viability_570nm_synthetic.tsv",
    "ref": "examples/viability_650nm_reference_synthetic.tsv",
    "blankMean": blank.mean(), "controlMean": ctrl.mean(),
    "viability": viab,
    "fits": fits,
}

# Extra standalone fit cases: increasing curve, and too few points for a CI.
xc = np.array([0.01, 0.01, 0.1, 0.1, 1, 1, 10, 10, 100, 100])
y_inc = f4(np.log10(xc), 10, 90, 0.3, -1.5) + rng.normal(0, 3, xc.size)
xc2 = np.array([0.1, 1, 10, 100, 1000])
y2 = f4(np.log10(xc2), 0, 100, 1.0, 1.0) + np.array([1.5, -2.0, 2.2, -1.0, 0.8])
ref["fit_cases"] = [
    {"conc": xc, "y": y_inc, **fit_ref(xc, y_inc)},
    {"conc": xc2, "y": y2, **fit_ref(xc2, y2)},
    {"conc": xc2, "y": y2, "fixBottom": 0.0, "fixTop": 100.0, **fit_ref(xc2, y2, fix_bottom=0.0, fix_top=100.0)},
]

(FIX / "reference.json").write_text(json.dumps(clean(ref), indent=1) + "\n", encoding="utf-8")
print("wrote", FIX / "reference.json")

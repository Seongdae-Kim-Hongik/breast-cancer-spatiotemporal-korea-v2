"""89g — 88g 민감도 결과 종합(89c 정신질환판 일반화, 2026-10-21). 사용: python3 01_스크립트/89g_sens_combine.py <TAG> <VER>
산출 results_<VER>/S41_retention_by_cov.csv, S41_noise_rubin.csv, S41_nco_detail.csv, S41_evalues.csv, S41_leaveout_ranges.csv, S41_summary.md"""
import os, sys, glob, math, numpy as np, pandas as pd
from decimal import Decimal, ROUND_HALF_UP
TAG, VER = sys.argv[1], sys.argv[2]; SUF = f"_{VER}" if VER else ""
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(ROOT, "04_분석", TAG, f"results{SUF}")
P = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/S41_part_*.csv"))], ignore_index=True)
S7 = pd.read_csv(f"{R}/S7_multiverse_summary.csv"); ROB = S7.loc[S7.robust, "cov"].tolist()
M3 = pd.read_csv(f"{R}/S3_main_IRR.csv").set_index("cov"); sgn = np.sign(np.log(M3.IRR))
alt = P[~P.analysis.str.startswith(("G2_noise", "H1_nco", "bw_twoway"))]
rows = []
for v in ROB:
    a = alt[alt["cov"] == v]; ok = (a.cred & (np.sign(np.log(a.IRR)) == sgn[v])).sum()
    rows.append(dict(cov=v, kept=int(ok), n_specs=a.analysis.nunique(), IRR_main=M3.IRR[v], IRR_min=a.IRR.min(), IRR_max=a.IRR.max()))
ret = pd.DataFrame(rows); ret.to_csv(f"{R}/S41_retention_by_cov.csv", index=False)
# leave-out ranges
lo = []
for pre, lab in (("F1_drop_sido", "province"), ("F2_drop_year", "year")):
    a = P[P.analysis.str.startswith(pre)]
    for v in ROB:
        b = a[a["cov"] == v]; lo.append(dict(set=lab, cov=v, n=b.analysis.nunique(), n_cred_same_sign=int((b.cred & (np.sign(np.log(b.IRR)) == sgn[v])).sum())))
    k = a.groupby("analysis").first(); lo.append(dict(set=lab, cov="__robust_kept", n=len(k), n_cred_same_sign=f"{k.robust_kept.min()}-{k.robust_kept.max()} (median {k.robust_kept.median():g})", hotJ=f"{k.hot_jaccard_vs_main.min():.3f}-{k.hot_jaccard_vs_main.max():.3f}"))
pd.DataFrame(lo).to_csv(f"{R}/S41_leaveout_ranges.csv", index=False)
# Rubin pooling over noise replicates (log scale, SE from CrI width)
N = P[P.analysis.str.startswith("G2_noise")]; rr = []
for v in sorted(set(N["cov"])):
    b = N[N["cov"] == v]; q = np.log(b.IRR); se = (np.log(b.hi) - np.log(b.lo)) / (2 * 1.96); m = len(b)
    qb = q.mean(); W = (se ** 2).mean(); B = q.var(ddof=1); T = W + (1 + 1 / m) * B
    rr.append(dict(cov=v, m=m, IRR_pooled=round(math.exp(qb), 3), lo=round(math.exp(qb - 1.96 * math.sqrt(T)), 3), hi=round(math.exp(qb + 1.96 * math.sqrt(T)), 3), n_cred=int(b.cred.sum())))
nz = pd.DataFrame(rr); nz.to_csv(f"{R}/S41_noise_rubin.csv", index=False)
# comparison outcomes
nc = P[P.analysis.str.startswith("H1_nco") & P["cov"].isin(ROB)][["analysis", "cov", "IRR", "lo", "hi", "cred"]]; nc.to_csv(f"{R}/S41_nco_detail.csv", index=False)
# E-values (VanderWeele & Ding; RR scale, per 1 SD), from unrounded? main CSV is 3-decimal; use as given
def ev(r):
    r = 1 / r if r < 1 else r; return r + math.sqrt(r * (r - 1))
def r2(x): return float(Decimal(str(x)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
E = []
for v in ROB:
    r, l, h = M3.IRR[v], M3.lo[v], M3.hi[v]; lim = l if r > 1 else h
    E.append(dict(cov=v, IRR=r, E_point=r2(ev(r)), E_limit=r2(ev(lim)) if (lim > 1) == (r > 1) and lim != 1 else 1.0))
pd.DataFrame(E).to_csv(f"{R}/S41_evalues.csv", index=False)
oos = pd.read_csv(f"{R}/S41_oos2024.csv") if os.path.exists(f"{R}/S41_oos2024.csv") else None
bw = P[P.analysis == "bw_twoway"]; bw = bw[bw["cov"].str.contains("__")][["cov", "IRR", "lo", "hi", "cred"]]
k = P.groupby("analysis").first()[["robust_kept", "n_robust", "hot_jaccard_vs_main", "WAIC", "note"]]
with open(f"{R}/S41_summary.md", "w") as f:
    f.write(f"# {TAG} {VER} 민감도 묶음 요약(88g, {len(k)}개 적합)\n\n## 강건 공변량 유지(대안 분석 {alt.analysis.nunique()}개)\n{ret.to_markdown(index=False)}\n\n## 지역 내·지역 간(이원 중심화)\n{bw.to_markdown(index=False)}\n\n")
    f.write(f"## 대체값 교란 Rubin 결합\n{nz[nz['cov'].isin(ROB)].to_markdown(index=False)}\n\n## 비교 결과\n{nc.to_markdown(index=False)}\n\n## E-value\n{pd.DataFrame(E).to_markdown(index=False)}\n\n")
    f.write(f"## 2024 조건부 검증\n{oos.to_markdown(index=False) if oos is not None else '-'}\n\n## 분석별\n{k.to_markdown()}\n")
print(open(f"{R}/S41_summary.md").read()[:6000])

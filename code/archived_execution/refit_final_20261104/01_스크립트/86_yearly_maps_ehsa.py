#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 86 — 연도별 평활 SIR 지도(F45/F45b) + EHSA(Esri 17분류, T24) + EB SIR 연도별 Moran(T24b) + 지속성(Kendall W·인구보존 영모형, T34/T34c)
#  사용: python3 86_yearly_maps_ehsa.py <TAG> <VER>     VER = v1(=results/, figs/) | std | v2 ...(=results_<VER>/, figs_<VER>/)
#  원형: FBD 595_salm_yearly_persistence.R(9/15 Kendall W 버그 수정판) · 411b_ehsa_esri.R/578(Esri 표준) — 같은 정의를 파이썬으로 옮김. INLA 미사용.
#  입력: 04_분석/<TAG>/panel_<TAG><SUF>.csv(cases, E = 그 버전의 기대건수) · results<SUF>/nb_<TAG>.graph(81 이 만든 queen+섬 최근접 그래프, final.shp 순서)
#  정의
#   - 평활 SIR = 연도마다 Marshall 전역 경험적 베이즈(spdep::EBest estmm와 같은 식, 기대건수 E를 '인구' 자리에) — 관측/기대의 비(SIR)이지 조율이 아님.
#   - EHSA: 연도별 EB SIR 의 Getis-Ord Gi*(자기 포함 행표준화 W, 같은 그래프), |z|>1.96 = hot/cold 연도; Mann-Kendall = 연도 vs z 의 Kendall τ(P<0.05).
#           Esri 분류(411b와 동일): 마지막 해 hot/cold 여부 → ≥90% 연도 조건(ceil(0.9·T)) → Intensifying/Diminishing/Persistent, New, Consecutive, Oscillating, Sporadic; 마지막 해 아님 + ≥90% = Historical.
#   - 지속성: 연도별 원 SIR 순위(1=최고, 동률 평균)의 Kendall W(동률 보정; 미보정 병기). 영모형 (i) 연도 내 라벨 순열 (ii) Poisson(E × 그 해 전국 SIR) — 시군구 인구규모·연도 수준 보존,
#            영모형 순위도 관측과 같은 SIR(모의건수/E)·동률 처리(9/15 수정). 시군구별 'SIR>1 연수'의 영모형 기대·P·BH-FDR.
#  산출: results<SUF>/S11_ehsa.csv(시군구별) · S11_ehsa_summary.csv(17분류 수) · S11b_moran_yearly_EB.csv · S12_persistence.csv(Kendall W) · S12b_persistence_district.csv · S12c_yearly_EB_SIR.csv
#        figs<SUF>/<TAG>_F45_yearly_EB_SIR.png · _F45b_yearly_relSIR.png · _F24_EHSA.png · _F46_persistence.png · <TAG>_F45_F24_F46_captions.json
import os, sys, json, math, unicodedata
import numpy as np, pandas as pd, geopandas as gpd
from scipy.stats import kendalltau, rankdata
import warnings; warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*matmul.*")   # macOS Accelerate 가짜 경고
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
TAG = sys.argv[1]; VER = sys.argv[2] if len(sys.argv) > 2 else "v1"; SUF = "" if VER in ("", "v1") else f"_{VER}"
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, "04_분석", TAG)
RES = os.path.join(OUT, f"results{SUF}"); FIG = os.path.join(OUT, f"figs{SUF}"); os.makedirs(FIG, exist_ok=True)
SHP = "${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
NFC = lambda s: unicodedata.normalize("NFC", str(s)).strip()
mpl.rcParams.update({"font.family": ["Helvetica Neue", "Arial", "Apple SD Gothic Neo"], "font.size": 9, "pdf.fonttype": 42})
B = 2000; rng = np.random.default_rng(20261019)
log = lambda *a: print(f"[86 {TAG}{SUF}]", *a, flush=True)
# ---- 자료 ----
D = pd.read_csv(os.path.join(OUT, f"panel_{TAG}{SUF}.csv"), encoding="utf-8-sig"); D["fbd_region"] = D.fbd_region.map(NFC)
G = gpd.read_file(SHP); G["region"] = G.region.map(NFC); n = len(G); regs = list(G.region)
lines = open(os.path.join(RES, f"nb_{TAG}.graph")).read().split("\n"); assert int(lines[0]) == n
NB = [[] for _ in range(n)]
for ln in lines[1:]:
    t = ln.split()
    if len(t) >= 2: i = int(t[0]) - 1; NB[i] = [int(x) - 1 for x in t[2:]]; assert len(NB[i]) == int(t[1])
assert all(i in NB[j] for i in range(n) for j in NB[i]), "graph not symmetric"
YRS = sorted(D.year.unique()); T = len(YRS)
obs = np.array([r in set(D.fbd_region) for r in regs]); idx = np.where(obs)[0]; m = len(idx)   # 미관측 시군구는 모든 계산에서 제외(부분그래프)
pos = {g: k for k, g in enumerate(idx)}
O = np.full((T, m), np.nan); E = np.full((T, m), np.nan)
for (y, r), g in D.groupby(["year", "fbd_region"]):
    O[YRS.index(y), pos[regs.index(r)]] = g.cases.sum(); E[YRS.index(y), pos[regs.index(r)]] = g.E.sum()
assert not np.isnan(O).any() and not np.isnan(E).any(), "panel not complete over observed districts × years"
nbo = [[pos[j] for j in NB[i] if j in pos] for i in idx]
# ---- (a) 연도별 EB 평활 SIR (Marshall global, spdep EBest 식) ----
def ebest(o, e):
    r = o / e; b = o.sum() / e.sum(); s2 = (e * (r - b) ** 2).sum() / e.sum(); a = max(s2 - b / e.mean(), 0.0)
    return b + a * (r - b) / (a + b / e)
RAW = O / E; EB = np.vstack([ebest(O[t], E[t]) for t in range(T)]); NAT = O.sum(1) / E.sum(1)
yr = pd.DataFrame([{"fbd_region": regs[idx[k]], "year": YRS[t], "cases": O[t, k], "E": round(E[t, k], 3), "rawSIR": round(RAW[t, k], 4), "ebSIR": round(EB[t, k], 4),
                    "relSIR": round(EB[t, k] / NAT[t], 4)} for t in range(T) for k in range(m)])
yr.to_csv(os.path.join(RES, "S12c_yearly_EB_SIR.csv"), index=False, encoding="utf-8-sig")
BR = [0, 0.67, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 99]; LAB8 = ["<0.67", "0.67-0.8", "0.8-0.9", "0.9-1", "1-1.1", "1.1-1.25", "1.25-1.5", ">=1.5"]
PAL8 = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]
G["sido"] = G.sigungu_cd.astype(str).str[:2]; SIDO = G.dissolve("sido")
def mainland_limits():   # 울릉·독도 같은 먼 섬 때문에 지도가 작아지지 않게 경계는 전체 범위 그대로 두되 여백만 줄임
    b = G.total_bounds; return (b[0], b[2]), (b[1], b[3])
XL, YL = mainland_limits()
def draw(ax, vals, title, br=BR, pal=PAL8, outline=None, ocol="black"):
    gg = G.copy(); gg["v"] = np.nan; gg.loc[idx, "v"] = vals
    gg.plot(ax=ax, column="v", cmap=ListedColormap(pal), norm=BoundaryNorm(br, len(pal)), edgecolor="#b0b0b0", linewidth=0.08, missing_kwds={"color": "#e6e6e6"})
    SIDO.boundary.plot(ax=ax, color="#555555", linewidth=0.25)
    if outline is not None:
        for flag, col in outline:
            sel = np.zeros(n, bool); sel[idx[flag]] = True
            if sel.any(): G[sel].boundary.plot(ax=ax, color=col, linewidth=0.5)
    ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_axis_off(); ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
def grid14(mat, fname, legtitle):
    nc = 5; nr = math.ceil(T / nc); fig, axs = plt.subplots(nr, nc, figsize=(13, 3.0 * nr))
    for a in axs.flat: a.set_axis_off()
    for t in range(T): draw(axs.flat[t], mat[t], str(YRS[t]))
    h = [Patch(facecolor=c, edgecolor="#777777", label=l) for c, l in zip(PAL8[::-1], LAB8[::-1])]
    axs.flat[-1].legend(handles=h, title=legtitle, loc="center", frameon=False, fontsize=8, title_fontsize=8.5)
    plt.tight_layout(); fig.savefig(os.path.join(FIG, fname), dpi=250, facecolor="white"); plt.close(fig)
grid14(EB, f"{TAG}_F45_yearly_EB_SIR.png", "EB-smoothed SIR"); grid14(EB / NAT[:, None], f"{TAG}_F45b_yearly_relSIR.png", "EB-smoothed SIR\nrelative to national\nSIR of the year")
log("F45 done; national SIR by year", dict(zip(YRS, NAT.round(3))))
# ---- (b) EHSA: Gi* (include.self, 행표준화) + MK ----
def gistar(x):   # spdep::localG(x, include.self(nb) W 스타일) 와 같은 식: Gi* = (Σ w x − x̄ Σw)/(s √((n Σw² − (Σw)²)/(n−1))), x̄·s는 전체(n) 기준
    nn = len(x); xb = x.mean(); s = math.sqrt((x ** 2).sum() / nn - xb ** 2); z = np.zeros(nn)
    for i in range(nn):
        nbr = [i] + nbo[i]; k = len(nbr); w = 1.0 / k   # 자기 포함 행표준화: Σw=1, Σw²=1/k
        z[i] = (x[nbr].sum() * w - xb) / (s * math.sqrt((nn * (k * w * w) - 1.0) / (nn - 1)))
    return z
Z = np.vstack([gistar(EB[t]) for t in range(T)])   # T × m
yrs = np.array(YRS, float); thr = math.ceil(0.9 * T)
def mk(zv):
    tau, p = kendalltau(yrs, zv); return (np.nan, np.nan) if np.isnan(tau) else (tau, p)
def last_run(b):
    c = 0
    for v in b[::-1]:
        if v: c += 1
        else: break
    return c
def esri(zv):
    hot = zv > 1.96; cold = zv < -1.96; nh, nc = hot.sum(), cold.sum(); tau, p = mk(zv); up = (p < 0.05) and tau > 0; dn = (p < 0.05) and tau < 0
    if hot[-1]:
        if nh >= thr: return "Intensifying hot" if up else ("Diminishing hot" if dn else "Persistent hot")
        if nh == 1: return "New hot"
        if last_run(hot) >= 2 and last_run(hot) == nh: return "Consecutive hot"
        return "Oscillating hot" if cold.any() else "Sporadic hot"
    if cold[-1]:
        if nc >= thr: return "Intensifying cold" if dn else ("Diminishing cold" if up else "Persistent cold")
        if nc == 1: return "New cold"
        if last_run(cold) >= 2 and last_run(cold) == nc: return "Consecutive cold"
        return "Oscillating cold" if hot.any() else "Sporadic cold"
    if nh >= thr: return "Historical hot"
    if nc >= thr: return "Historical cold"
    return "No pattern"
CATS = [f"{a} {b}" for b in ("hot", "cold") for a in ("New", "Consecutive", "Intensifying", "Persistent", "Diminishing", "Sporadic", "Oscillating", "Historical")] + ["No pattern"]
cat = [esri(Z[:, k]) for k in range(m)]; mkr = [mk(Z[:, k]) for k in range(m)]
eh = pd.DataFrame({"fbd_region": [regs[i] for i in idx], "sigungu_cd": G.sigungu_cd.values[idx], "EHSA": cat, "n_hot_years": (Z > 1.96).sum(0), "n_cold_years": (Z < -1.96).sum(0),
                   "giz_first": Z[0].round(2), "giz_final": Z[-1].round(2), "mk_tau": [round(a, 3) for a, _ in mkr], "mk_p": [round(b, 4) for _, b in mkr],
                   "ebSIR_mean": EB.mean(0).round(4), "rawSIR_pooled": (O.sum(0) / E.sum(0)).round(4)})
eh.to_csv(os.path.join(RES, "S11_ehsa.csv"), index=False, encoding="utf-8-sig")
es = pd.DataFrame({"category": CATS, "n": [cat.count(c) for c in CATS]}); es["regions_top"] = [";".join(eh[eh.EHSA == c].sort_values("ebSIR_mean", ascending=("cold" in c)).fbd_region.head(8)) for c in CATS]
es.to_csv(os.path.join(RES, "S11_ehsa_summary.csv"), index=False, encoding="utf-8-sig")
log("EHSA", {r.category: r.n for r in es.itertuples() if r.n})
# Moran I (EB SIR, 행표준화 W 자기 미포함, 999 순열) + Kendall 추세
WR = np.zeros((m, m))
for i in range(m):
    if nbo[i]: WR[i, nbo[i]] = 1.0 / len(nbo[i])
S0 = WR.sum()
def moranM(X):   # X: (k, m) 여러 벡터 한 번에
    Xc = X - X.mean(1, keepdims=True); return (m / S0) * np.einsum("ij,ij->i", Xc, Xc @ WR.T) / (Xc ** 2).sum(1)
mr = []
for t in range(T):
    I = moranM(EB[t][None])[0]; P = np.array([rng.permutation(EB[t]) for _ in range(999)]); Ip = moranM(P)
    mr.append({"year": YRS[t], "moran_I": round(I, 4), "p_perm": round((np.sum(Ip >= I) + 1) / 1000, 4)})
mr = pd.DataFrame(mr); kt = kendalltau(mr.year, mr.moran_I); mr["trend_tau"] = round(kt[0], 3); mr["trend_p"] = round(kt[1], 4)
mr.to_csv(os.path.join(RES, "S11b_moran_yearly_EB.csv"), index=False)
# EHSA 지도
EPAL = {"New hot": "#fdd0a2", "Consecutive hot": "#fdae6b", "Intensifying hot": "#a50f15", "Persistent hot": "#de2d26", "Diminishing hot": "#fc9272", "Sporadic hot": "#fee0d2", "Oscillating hot": "#d6a5c9", "Historical hot": "#bcbddc",
        "New cold": "#c6dbef", "Consecutive cold": "#9ecae1", "Intensifying cold": "#08519c", "Persistent cold": "#3182bd", "Diminishing cold": "#6baed6", "Sporadic cold": "#deebf7", "Oscillating cold": "#a1d99b", "Historical cold": "#c7e9c0", "No pattern": "#f7f7f7"}
fig, ax = plt.subplots(figsize=(7.5, 7.5)); gg = G.copy(); gg["c"] = "unobserved"; gg.loc[idx, "c"] = cat
gg.plot(ax=ax, color=[EPAL.get(c, "#e6e6e6") for c in gg.c], edgecolor="#9a9a9a", linewidth=0.1); SIDO.boundary.plot(ax=ax, color="#555555", linewidth=0.3)
ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_axis_off()
present = [c for c in CATS if c in cat]
ax.legend(handles=[Patch(facecolor=EPAL[c], edgecolor="#777777", label=f"{c} ({cat.count(c)})") for c in present], loc="lower right", frameon=False, fontsize=8, title="EHSA category (n)", title_fontsize=8.5)
plt.tight_layout(); fig.savefig(os.path.join(FIG, f"{TAG}_F24_EHSA.png"), dpi=250, facecolor="white"); plt.close(fig)
# ---- (c) 지속성: Kendall W(SIR 순위) + 인구보존 Poisson 영모형 ----
def ranks(M): return np.vstack([rankdata(-r, method="average") for r in M])   # 행=연도, 1=최고
def kW(R, corrected=True):
    mm, nn = R.shape; S = ((R.sum(0) - R.sum(0).mean()) ** 2).sum(); den = mm ** 2 * (nn ** 3 - nn)
    if corrected:
        Tt = 0.0
        for r in R: _, c = np.unique(r, return_counts=True); Tt += (c ** 3 - c).sum()
        den -= mm * Tt
    return 12 * S / den
RK = ranks(RAW); EMU = E * NAT[:, None]
post = [y for y in YRS if y >= 2020]
SETS = {"2011-2024": [True] * T, "2011-2019": [y <= 2019 for y in YRS], "2020-2024": [y >= 2020 for y in YRS]}
for wn, (a, b) in {"W1_2011-15": (2011, 2015), "W2_2016-19": (2016, 2019), "W3_2020-24": (2020, 2024)}.items(): SETS[wn] = [a <= y <= b for y in YRS]
SIMS = [rng.poisson(EMU) for _ in range(B)]   # 공용 모의자료(B × T × m)
kw = []
for nm, sel in SETS.items():
    sel = np.array(sel); R = RK[sel]; Wv = kW(R)
    Wp = np.array([kW(np.vstack([rng.permutation(r) for r in R])) for _ in range(B)])
    Wq = np.array([kW(ranks((S_ / E)[sel])) for S_ in SIMS])
    kw.append({"set": nm, "n_years": int(sel.sum()), "kendall_W": round(Wv, 4), "kendall_W_uncorrected": round(kW(R, False), 4),
               "perm_p": round((np.sum(Wp >= Wv) + 1) / (B + 1), 4), "perm_null_97.5": round(np.quantile(Wp, 0.975), 4), "perm_z": round((Wv - Wp.mean()) / Wp.std(ddof=1), 1),
               "poisson_p": round((np.sum(Wq >= Wv) + 1) / (B + 1), 4), "poisson_null_97.5": round(np.quantile(Wq, 0.975), 4), "poisson_z": round((Wv - Wq.mean()) / Wq.std(ddof=1), 1), "B": B})
kw = pd.DataFrame(kw); kw.to_csv(os.path.join(RES, "S12_persistence.csv"), index=False)
gt1 = (RAW > 1) & (O > 0); q80 = np.array([np.quantile(RAW[t][O[t] > 0], 0.8) for t in range(T)]); top = (RAW >= q80[:, None]) & (O > 0)
simG = np.array([((S_ > 0) & (S_ > E)).sum(0) for S_ in SIMS])   # B × m
obsG = gt1.sum(0); pv = (np.sum(simG >= obsG, 0) + 1) / (B + 1)
def bh(p):
    o = np.argsort(p); r = np.empty_like(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; r[o] = np.minimum(q, 1); return r
pst = pd.DataFrame({"fbd_region": [regs[i] for i in idx], "cases": O.sum(0), "E": E.sum(0).round(1), "SIR_pooled": (O.sum(0) / E.sum(0)).round(4), "years_SIR_gt1": obsG,
                    "exp_years_SIR_gt1": simG.mean(0).round(2), "p_years_SIR_gt1": pv.round(4), "fdr_years_SIR_gt1": bh(pv).round(4), "years_top_quintile": top.sum(0), "mean_rank": RK.mean(0).round(1)})
pst = pst.sort_values(["p_years_SIR_gt1", "years_SIR_gt1"], ascending=[True, False]); pst.to_csv(os.path.join(RES, "S12b_persistence_district.csv"), index=False, encoding="utf-8-sig")
PB = [-0.5, 0.5, 2.5, 4.5, 6.5, 9.5, 13.5, 14.5]; PL = ["0", "1-2", "3-4", "5-6", "7-9", "10-13", "14"]; PP = ["#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026", "#67000d"]
fig, ax = plt.subplots(figsize=(7.5, 7.5)); vals = pst.set_index("fbd_region").loc[[regs[i] for i in idx], "years_SIR_gt1"].values
sig = pst.set_index("fbd_region").loc[[regs[i] for i in idx], "fdr_years_SIR_gt1"].values < 0.05
draw(ax, vals, "", br=PB, pal=PP, outline=[(sig, "black")])
ax.legend(handles=[Patch(facecolor=c, edgecolor="#777777", label=l) for c, l in zip(PP[::-1], PL[::-1])] + [Patch(facecolor="none", edgecolor="black", label="FDR < 0.05 vs population-\npreserving Poisson null")],
          title=f"Years with SIR > 1 (of {T})", loc="lower right", frameon=False, fontsize=8, title_fontsize=8.5)
plt.tight_layout(); fig.savefig(os.path.join(FIG, f"{TAG}_F46_persistence.png"), dpi=250, facecolor="white"); plt.close(fig)
Etxt = "indirect age-sex standardized expected counts (E_std)" if SUF == "_std" else "expected counts from the national crude rate of each year (not age-standardized)"
cap = {"F45": {"file": f"{TAG}_F45_yearly_EB_SIR.png", "caption": f"Yearly empirical-Bayes smoothed standardized incidence ratio (SIR = observed / expected; {Etxt}) by district, {YRS[0]}-{YRS[-1]}. Marshall global EB smoothing within each year; common log-symmetric classes. This is a ratio to expected counts, not a crude rate."},
       "F45b": {"file": f"{TAG}_F45b_yearly_relSIR.png", "caption": "Same as F45 divided by the national SIR of the same year (removes the national level of each year)."},
       "F24": {"file": f"{TAG}_F24_EHSA.png", "caption": f"Emerging hot spot analysis (Esri classification): yearly Getis-Ord Gi* of the EB-smoothed SIR (self-inclusive row-standardized queen weights with island links; |z| > 1.96), Mann-Kendall trend of Gi* z over {T} years; persistent/intensifying/diminishing/historical require at least {thr} of {T} years."},
       "F46": {"file": f"{TAG}_F46_persistence.png", "caption": f"Number of years (of {T}) with raw SIR above 1 by district; black outline = more years than expected under a Poisson null preserving district expected counts and the national level of each year ({B} draws, BH-FDR < 0.05)."}}
json.dump(cap, open(os.path.join(FIG, f"{TAG}_F45_F24_F46_captions.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
log("Kendall W", kw[["set", "kendall_W", "perm_p", "poisson_p", "poisson_z"]].to_dict("records")); log("persistent FDR<0.05:", int(sig.sum()), "| Moran EB", mr.moran_I.min(), "-", mr.moran_I.max(), "tau", mr.trend_tau[0])
log("DONE")

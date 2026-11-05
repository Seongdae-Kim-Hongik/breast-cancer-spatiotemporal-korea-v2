"""88m prep — 출산·혼인 변수 장기 시차 노출(2026-10-21, 유방암 /333·정독 반영). 사용: python3 01_스크립트/88m_long_lag_prep.py <TAG> <VER>
 산출 results_<VER>/S45_panel_longlag.csv: <v>__lag10(t-10 원값 → z), <v>__base0004(2000-2004 지역 평균 → z). 결측 지역(경계 변경)은 가장 가까운 관측 연도 값."""
import os, sys, numpy as np, pandas as pd
TAG, VER = sys.argv[1], sys.argv[2]; SUF = f"_{VER}" if VER else ""
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); B = os.path.join(ROOT, "04_분석", TAG); R = os.path.join(B, f"results{SUF}")
P = pd.read_csv(f"{B}/panel_{TAG}{SUF}.csv"); META = pd.read_csv(f"{B}/panel_META_{TAG}{SUF}.csv").set_index("variable")
V = ["total_fertility_rate", "mean_age_first_marriage_wife"]
M = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet"), columns=["year", "fbd_region"] + [META.loc[v, "source_column"] for v in V])
z = lambda s: (s - s.mean()) / s.std(ddof=0); X = P[["year", "fbd_region"]].copy(); rep = {}
for v in V:
    c = META.loc[v, "source_column"]; w = M.pivot_table(index="fbd_region", columns="year", values=c)
    def at(r, y):
        if r not in w.index: return np.nan
        s = w.loc[r].dropna();
        if not len(s): return np.nan
        return s.iloc[(abs(s.index - y)).argmin()]
    X[f"{v}__lag10"] = [at(r, y - 10) for r, y in zip(P.fbd_region, P.year)]
    base = {r: (w.loc[r, [y for y in range(2000, 2005) if y in w.columns]].mean() if r in w.index and w.loc[r, [y for y in range(2000, 2005) if y in w.columns]].notna().any() else at(r, 2002)) for r in P.fbd_region.unique()}
    X[f"{v}__base0004"] = P.fbd_region.map(base)
    for k in (f"{v}__lag10", f"{v}__base0004"):
        rep[k] = int(X[k].isna().sum()); X[k] = z(X[k].fillna(X[k].mean()))
X.to_csv(f"{R}/S45_panel_longlag.csv", index=False); print("missing filled:", rep)

"""89i — 진료 유병(BYM2 창별 RR) 대 국가암등록 시군구 발생률 검증(2026-10-21, 유방암 /333 반영). 사용: python3 01_스크립트/89i_incidence_validation.py <TAG> <VER> <암등록 접두어>
 예) python3 01_스크립트/89i_incidence_validation.py BRCA std 암_유방
 입력: 마스터 canc__<접두어>__연령표준화발생률·발생자수(5년 시기, year=종료연도), results_<VER>/S4b_window_spatialRR.csv(창별 BYM2 RR·핫스팟)
 비교: 2009-2013 발생 ↔ 창 W1(2011-15), 2014-2018 발생 ↔ 창 W2(2016-19). 지표 = Spearman, 발생 상·하위 20% 대 핫·콜드스팟 일치, 수도권 대 기타 비.
 산출: results_<VER>/S44_incidence_validation.csv, S44_incidence_pairs.csv"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
TAG, VER, PFX = sys.argv[1], sys.argv[2], sys.argv[3]; SUF = f"_{VER}" if VER else ""
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(ROOT, "04_분석", TAG, f"results{SUF}")
M = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet"), columns=["year", "fbd_region", f"canc__{PFX}__연령표준화발생률", f"canc__{PFX}__발생자수"])
M.columns = ["year", "fbd_region", "asir", "inc_n"]
W = pd.read_csv(f"{R}/S4b_window_spatialRR.csv"); print(W.columns.tolist())
wcol = [c for c in W.columns if c.lower() in ("window",)][0]; rrc = [c for c in W.columns if c.upper().startswith("RR")][0]
clc = [c for c in W.columns if "class" in c.lower()]
pairs = {2013: [w for w in W[wcol].unique() if "2011" in str(w)][0], 2018: [w for w in W[wcol].unique() if "2016" in str(w)][0]}
S6 = pd.read_csv(f"{R}/S6_spatial_RR.csv")[["fbd_region", "RR", "class", "SIR_obs"]]
out, allp = [], []
cap = lambda r: r.str.startswith(("서울", "경기", "인천"))
for yr, w in pairs.items():
    a = M[M.year == yr].dropna(subset=["asir"]); b = W[W[wcol] == w][["fbd_region", rrc] + clc]
    p = a.merge(b, on="fbd_region"); p["window"] = w; p["inc_period"] = f"{yr-4}-{yr}"; allp.append(p)
    rho, prho = spearmanr(p.asir, p[rrc]); r, pr = pearsonr(np.log(p.asir), np.log(p[rrc]))
    q_hi, q_lo = p.asir.quantile(.8), p.asir.quantile(.2)
    row = dict(inc_period=f"{yr-4}-{yr}", window=w, n=len(p), spearman=round(rho, 3), p_spearman=prho, pearson_log=round(r, 3),
               asir_ratio_capital_vs_other=round(p[cap(p.fbd_region)].asir.mean() / p[~cap(p.fbd_region)].asir.mean(), 3),
               rr_ratio_capital_vs_other=round(p[cap(p.fbd_region)][rrc].mean() / p[~cap(p.fbd_region)][rrc].mean(), 3))
    if clc:
        c = clc[0]; row["hot_in_top20_inc"] = f"{int(((p[c]=='hotspot') & (p.asir>=q_hi)).sum())}/{int((p[c]=='hotspot').sum())}"
        row["cold_in_bottom20_inc"] = f"{int(((p[c]=='coldspot') & (p.asir<=q_lo)).sum())}/{int((p[c]=='coldspot').sum())}"
        row["hot_mean_asir_vs_cold"] = round(p[p[c]=='hotspot'].asir.mean() / p[p[c]=='coldspot'].asir.mean(), 3)
    out.append(row)
# 통합 RR ↔ 두 시기 평균 발생률
# Pooled concordance requires both registry periods; do not silently average one period.
registry = M[M.year.isin([2013, 2018])].pivot(index="fbd_region", columns="year", values="asir")
registry[registry.isna().any(axis=1)].to_csv(f"{R}/S44_registry_missing_periods.csv")
a = registry.dropna(subset=[2013,2018]).mean(axis=1).rename("asir_mean").reset_index()
p = a.merge(S6, on="fbd_region", validate="one_to_one")
assert p.asir_mean.gt(0).all() and p.RR.gt(0).all()
rho, prho = spearmanr(p.asir_mean, p.RR); out.append(dict(inc_period="2009-2018 mean (both periods observed)", window="pooled 2011-2024", n=len(p), spearman=round(rho, 3), p_spearman=prho,
    pearson_log=round(pearsonr(np.log(p.asir_mean), np.log(p.RR))[0], 3), hot_mean_asir_vs_cold=round(p[p["class"]=="hotspot"].asir_mean.mean() / p[p["class"]=="coldspot"].asir_mean.mean(), 3),
    asir_ratio_capital_vs_other=round(p[cap(p.fbd_region)].asir_mean.mean() / p[~cap(p.fbd_region)].asir_mean.mean(), 3), rr_ratio_capital_vs_other=round(p[cap(p.fbd_region)].RR.mean() / p[~cap(p.fbd_region)].RR.mean(), 3)))
o = pd.DataFrame(out); o.to_csv(f"{R}/S44_incidence_validation.csv", index=False); pd.concat(allp).to_csv(f"{R}/S44_incidence_pairs.csv", index=False)
print(o.to_string(index=False))

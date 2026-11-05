"""88g prep — 공간 트랙 질환 공통 민감도 묶음 입력 패널(88f_mental_sens_r3_prep.py 일반화, 2026-10-21 유방암 2호부터).
사용(데이터셋 루트): python3 01_스크립트/88g_sens_prep.py <TAG> <VER> [비교결과 TAG,TAG]   예) python3 01_스크립트/88g_sens_prep.py BRCA std GCA,INFECT
산출 04_분석/<TAG>/results_<VER>/S41_panel_*.csv (year, fbd_region 키 + 해당 분석 열). 원 패널·주모형 산출은 읽기만.
 sido: 시도 코드(sigungu_cd 앞 두 자리), sido_year
 lag1/lag2: 가설 공변량 z 를 지역 내 1·2년 지연(첫 해 행은 결측 → 해당 분석에서 그 연도 제외)
 ma3: KCHS 공변량(열 이름 kchs__ 원천) 원값 지역 내 3년 중심 이동평균(양 끝 2년 평균) → 전체 z
 extreme: 가설 공변량 원값이 연도별 사분위범위 5배 밖 → 결측 → 지역 내 선형보간 → z (A형간염 R7 규칙)
 splus: Spatial+ — 공변량 z 를 그래프 라플라시안 저주파 고유벡터 k개(k=round(0.1n))에 회귀한 잔차 → z (지역 자료용 변형)
 noise_r01..r20: 마스터 원천 열에서 결측이던(=대체된) 셀에 지역 내 관측 표준편차 크기 정규잡음 → z (대체 셀이 없으면 만들지 않음)
 nco_<TAG2>: 비교 결과 질환 cases·E_std 를 이 질환 키에 붙임
"""
import os, sys, json, numpy as np, pandas as pd
TAG, VER = sys.argv[1], sys.argv[2]; NCO = [t for t in (sys.argv[3].split(",") if len(sys.argv) > 3 else []) if t and t != TAG]
SUF = f"_{VER}" if VER else ""
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); B = os.path.join(ROOT, "04_분석", TAG); R = os.path.join(B, f"results{SUF}")
P = pd.read_csv(f"{B}/panel_{TAG}{SUF}.csv"); S = json.load(open(f"{B}/config_{TAG}{SUF}.json"))
best = open(f"{B}/HM_best_{TAG}{SUF}.txt").read().split("best_set=")[1].split()[0]; HYP = S["sets_after_vif"][best]; FORCED = S["forced"]; ALLV = list(dict.fromkeys(HYP + FORCED))
META = pd.read_csv(f"{B}/panel_META_{TAG}{SUF}.csv").set_index("variable")
P = P.sort_values(["fbd_region", "year"]).reset_index(drop=True)
z = lambda s: (s - s.mean()) / s.std(ddof=0)
key = P[["year", "fbd_region"]].copy()
s6 = pd.read_csv(f"{R}/S6_spatial_RR.csv", dtype={"sigungu_cd": str})
# sido
cd = dict(zip(s6.fbd_region, s6.sigungu_cd.str[:2]))
X = key.copy(); X["sido"] = P.fbd_region.map(cd); assert X.sido.notna().all(); X["sido_year"] = X.sido + "_" + P.year.astype(str); X.to_csv(f"{R}/S41_panel_sido.csv", index=False)
# lags
for L in (1, 2):
    X = key.copy()
    for v in HYP: X[v] = P.groupby("fbd_region")[v].shift(L)
    X.to_csv(f"{R}/S41_panel_lag{L}.csv", index=False)
# ma3 (KCHS 원천 공변량)
src = {v: META.loc[v, "source_column"] for v in ALLV if v in META.index and isinstance(META.loc[v, "source_column"], str)}
KCHS = [v for v in HYP if str(src.get(v, "")).startswith("kchs__") and f"{v}_raw" in P.columns]
X = key.copy()
for v in HYP: X[v] = P[v]
for v in KCHS: X[v] = z(P.groupby("fbd_region")[f"{v}_raw"].transform(lambda s: s.rolling(3, center=True, min_periods=2).mean()))
X.to_csv(f"{R}/S41_panel_ma3.csv", index=False)
# extreme (R7)
X = key.copy(); nflag = {}
for v in HYP:
    if f"{v}_raw" not in P.columns: X[v] = P[v]; continue
    r = P[f"{v}_raw"].astype(float).copy(); q = P.groupby("year")[f"{v}_raw"]; q1, q3 = q.transform(lambda s: s.quantile(.25)), q.transform(lambda s: s.quantile(.75)); iqr = q3 - q1
    bad = ((r < q1 - 5 * iqr) | (r > q3 + 5 * iqr)) & (iqr > 0); nflag[v] = int(bad.sum()); r[bad] = np.nan
    r = r.groupby(P.fbd_region).transform(lambda s: s.interpolate(limit_direction="both")); X[v] = z(r.fillna(r.mean()))
X.to_csv(f"{R}/S41_panel_extreme.csv", index=False); json.dump(nflag, open(f"{R}/S41_extreme_flags.json", "w"), ensure_ascii=False, indent=1)
# Component-aware Spatial+ is a separately verified frozen input, not regenerated here.
K = 23
# noise replicates on imputed cells
M = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet")); M = M[(M.year >= 2011) & (M.year <= 2024)]
MM = P[["year", "fbd_region"]].merge(M[["year", "fbd_region"] + sorted(set(c for c in src.values() if c in M.columns))], on=["year", "fbd_region"], how="left")
flag = {v: MM[c].isna().values | ((MM[c] == 0).values & (P[f"{v}_raw"].values != 0)) for v, c in src.items() if c in MM.columns and f"{v}_raw" in P.columns}
nimp = {v: int(f.sum()) for v, f in flag.items() if f.sum()}
json.dump(nimp, open(f"{R}/S41_imputed_cells.json", "w"), ensure_ascii=False, indent=1)
if nimp:
    rng = np.random.default_rng(20261021)
    for r_ in range(1, 21):
        X = key.copy()
        for v in ALLV:
            if f"{v}_raw" not in P.columns: X[v] = P[v]; continue
            raw = P[f"{v}_raw"].values.astype(float).copy()
            if v in flag and flag[v].any():
                sd = pd.Series(np.where(flag[v], np.nan, raw)).groupby(P.fbd_region.values).transform("std").fillna(pd.Series(raw).std() * 0.5).values
                raw[flag[v]] = raw[flag[v]] + rng.normal(0, sd[flag[v]])
            X[v] = z(pd.Series(raw))
        X.to_csv(f"{R}/S41_panel_noise_r{r_:02d}.csv", index=False)
# comparison outcomes
for t2 in NCO:
    g = pd.read_csv(os.path.join(ROOT, "04_분석", t2, f"panel_{t2}.csv"), usecols=["year", "fbd_region", "cases"])
    e = pd.read_csv(os.path.join(ROOT, "02_패널", f"{t2}_Estd.csv"), usecols=["year", "fbd_region", "E_std"])
    X = key.merge(g, on=["year", "fbd_region"], how="left").merge(e, on=["year", "fbd_region"], how="left"); X.columns = ["year", "fbd_region", f"cases_{t2}", f"E_{t2}"]
    print(t2, "missing", int(X.isna().any(axis=1).sum())); X.to_csv(f"{R}/S41_panel_nco_{t2}.csv", index=False)
print(TAG, VER, "HYP", len(HYP), "FORCED", len(FORCED), "KCHS", KCHS, "Spatial+ k", K, "imputed", nimp, "extreme", {k: v for k, v in nflag.items() if v}, "NCO", NCO)

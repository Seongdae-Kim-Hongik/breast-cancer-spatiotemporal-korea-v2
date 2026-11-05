#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 87 — 창별 결정요인 격자(F42b 방식: 변수 행 × 창 열, IRR·95% CrI, pooled 열 없음, 라벨 "English (한글)"). 원형 FBD 592_ehec_window_grid.py.
#  사용: python3 87_window_grid.py <TAG> <VER>   VER = v1(=results/, figs/) | std | v2 ...
#  입력: results<SUF>/S4_window_IRR.csv(81 S4) · S4_window_fit.csv(창별 n·사례) · S4_interaction.csv(있으면 † = 창×변수 상호작용 credible, 기준창 W1) · S3_main_IRR.csv(행 순서)
#  산출: figs<SUF>/<TAG>_F42b_window_grid.png · results<SUF>/S4c_window_grid_table.csv(그림과 같은 수치, 긴 형식)
#  색: credible(95% CrI가 1 제외)만 빨강(IRR>1)/파랑(IRR<1) 칸, 나머지 흰 칸 회색 글자. 가설 변수 먼저, 강제통제는 구분선 아래.
import os, sys, unicodedata
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
TAG = sys.argv[1]; VER = sys.argv[2] if len(sys.argv) > 2 else "v1"; SUF = "" if VER in ("", "v1") else f"_{VER}"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "04_분석", TAG)
RES = os.path.join(OUT, f"results{SUF}"); FIG = os.path.join(OUT, f"figs{SUF}"); os.makedirs(FIG, exist_ok=True)
mpl.rcParams.update({"font.family": ["Helvetica Neue", "Arial", "Apple SD Gothic Neo", "AppleGothic"], "font.size": 9, "pdf.fonttype": 42})
NFC = lambda s: unicodedata.normalize("NFC", str(s)).strip()
S = pd.read_csv(os.path.join(RES, "S4_window_IRR.csv")); F = pd.read_csv(os.path.join(RES, "S4_window_fit.csv")); M = pd.read_csv(os.path.join(RES, "S3_main_IRR.csv"))
ixp = os.path.join(RES, "S4_interaction.csv"); IX = pd.read_csv(ixp) if os.path.exists(ixp) else pd.DataFrame(columns=["cov", "window", "cred"])
ixc = {(r.cov, r.window): bool(r.cred) for r in IX.itertuples()}
hyp = [c for c in M["cov"] if M.set_index("cov").loc[c, "role"] == "hypothesis"]; frc = [c for c in M["cov"] if c not in hyp]
order = hyp + frc; cols = list(F.window)
lab = {r.cov: f"{NFC(r.label_en)} ({NFC(r.label_ko)})" for r in S.drop_duplicates("cov").itertuples()}
# 범위 표기 괄호 중첩 방지: 영문 라벨 속 "(SR %)" 같은 괄호는 그대로 두고 한글만 바깥 괄호
RED, BLUE, RL, BL = "#c0392b", "#1f4e79", "#f4d0cb", "#d3dff0"
nr, nc = len(order), len(cols)
fig, ax = plt.subplots(figsize=(4.2 + 2.0 * nc, 0.42 * nr + 1.3))
rows = []
for i, v in enumerate(order):
    for j, c in enumerate(cols):
        r = S[(S.window == c) & (S["cov"] == v)]
        if not len(r): continue
        r = r.iloc[0]; cred = bool(r.cred); posi = r.IRR > 1
        ax.add_patch(Rectangle((j, nr - 1 - i), 1, 1, facecolor=(RL if posi else BL) if cred else "white", edgecolor=(RED if posi else BLUE) if cred else "#cccccc", linewidth=1.4 if cred else 0.5))
        dp = 3 if max(abs(r.lo - 1), abs(r.hi - 1)) < 0.1 else 2   # 칸 안 세 숫자는 같은 자릿수(1 근처 효과는 3자리)
        ax.text(j + 0.5, nr - 1 - i + 0.5, f"{r.IRR:.{dp}f} ({r.lo:.{dp}f}-{r.hi:.{dp}f})",
                ha="center", va="center", fontsize=8.2, color=((RED if posi else BLUE) if cred else "#555555"), fontweight="bold" if cred else "normal")
        dag = ixc.get((v, c), False)
        if dag: ax.text(j + 0.97, nr - 1 - i + 0.93, "†", ha="right", va="top", fontsize=8.5, color="#333333")
        rows.append({"cov": v, "label": lab[v], "role": "hypothesis" if v in hyp else "forced", "window": c, "IRR": r.IRR, "lo": r.lo, "hi": r.hi, "cred": cred, "interaction_vs_W1_cred": dag})
if frc and hyp: ax.axhline(nr - len(hyp), color="#333333", linewidth=0.8)
ax.set_xlim(0, nc); ax.set_ylim(0, nr); ax.set_xticks(np.arange(nc) + 0.5)
ax.set_xticklabels([f"{r.years}\n(n = {int(r.cases):,})" for r in F.itertuples()], fontsize=9); ax.xaxis.tick_top()
ax.set_yticks(np.arange(nr) + 0.5); ax.set_yticklabels([lab[v] for v in order][::-1], fontsize=8.6)
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(length=0); plt.tight_layout()
fn = os.path.join(FIG, f"{TAG}_F42b_window_grid.png"); fig.savefig(fn, dpi=220, facecolor="white"); plt.close(fig)
pd.DataFrame(rows).to_csv(os.path.join(RES, "S4c_window_grid_table.csv"), index=False, encoding="utf-8-sig")
nd = pd.DataFrame(rows)
print(f"[87 {TAG}{SUF}] saved {os.path.basename(fn)} | credible hyp per window:", nd[(nd.role == "hypothesis")].groupby("window").cred.sum().to_dict(), "| interaction † :", int(nd.interaction_vs_W1_cred.sum()))

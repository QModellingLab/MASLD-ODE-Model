# Equilibrate-then-Inhibit (Two-Window) Control Analysis

**目的**：回應 npj Systems Biology and Applications Reviewer 2 的質疑——
「所有 active state 從 0 開始，未平衡就給藥，NAFL/NASH 早治療優勢是否為
未平衡的人工假象？」

## 協定摘要

- **Window 1（t=0–300h，undrugged）**：與原稿 Fig 5 / Supplementary Fig S2
  完全相同的模擬（NAFL/NASH 各自從 pyDESeq2 ratio 起始）。已數值驗證此
  終點與真正數學平衡態（t≈3000h）對三個 core output 及其上游驅動者而言
  差距 0.0000%——即稿件 Methods 已聲稱的「300h 足以達到 steady state or
  stable saturating plateau」為真。
- **Window 2（Window 1 終點為新起始點，再模擬 300h）**：此時才引入
  silymarin 8 靶點抑制（ki = 1.0/0.7/0.5/0.3）。

介入點仍由同一個 NAFL 或 NASH transcriptomic snapshot 驅動（節點層級
inactive+active 質量守恆，等於原始 ratio），不涉及外推到未觀察的晚期
疾病階段。

## 檔案結構

```
analyses/equilibrate_then_inhibit/
├── equilibrate_then_inhibit.py   ← 分析腳本（獨立可執行）
├── README.md                      ← 本檔
└── outputs/                       ← 執行後自動產生
    ├── TwoWindow_silymarin_summary.xlsx
    │     - Summary_ki0.3：三個 core output 的 NAFL/NASH ratio
    │     - Full_dose_response：完整 dose-response 表
    │     - Window1_convergence_check：Window-1 終點 vs 真正平衡態的收斂診斷
    └── two_window_results.json
```

## 執行方式

腳本會自動在下列路徑搜尋 `ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py`
與 `ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py`（依序嘗試，找到第一個
存在的就用）：

1. `<repo_root>/data/models/`
2. `<repo_root>/data/models_pydeseq2_mean/`
3. `analyses/equilibrate_then_inhibit/models_pydeseq2_mean/`
4. `analyses/equilibrate_then_inhibit/`（腳本同層）

若你的模型檔案不在 `data/` 底下，最簡單的做法是把兩個模型 .py 檔複製一份
到 `analyses/equilibrate_then_inhibit/` 這個資料夾內（第 4 個路徑會抓到）。

```powershell
cd F:\MASLD\ODEModeling\MASLD-ODE-Model\analyses\equilibrate_then_inhibit
conda activate deseq2_env
python equilibrate_then_inhibit.py
```

## 主要結果（GSE126848，ki=0.3）

| Output | NAFL reduction (%) | NASH reduction (%) | NAFL/NASH ratio |
|---|---|---|---|
| P_Hepatocyte_injury | 2.37 | 1.35 | 1.76 |
| P_Cell_death | 1.01 | 0.50 | 2.03 |
| P_Inflammation | 0.88 | 0.44 | 2.00 |

方向保留（均 >1，對照原 Table 4 的 2.10–3.03×），絕對量級大幅萎縮
（9–34% → 0.4–2.4%）——與 receptor-reserve 機制（Buchwald 2020，稿件已引用）
一致，非模型建構的必然假象。

## 建議 response letter 措辭

見對話紀錄（2026-09-02/03）。核心論點：
1. 未外推到未觀察的晚期疾病（節點總量守恆 + t=300h 已收斂）
2. 方向保留、量級萎縮 → 支持將 claim 由 "therapeutic window in established
   disease" 收斂為 "stage-dependent transient/early-response window"

Alternative Multi-gene Node Aggregation Rules (R1-3)
=======================================================
用途：回應 npj SBA Reviewer 1, Major Comment 3——
「多基因節點（包含大型呼吸鏈複合體與 AMPK）用不加權算術平均代表，但這些
成員基因的功能與化學計量比不同。應該釐清這個對應關係，並測試不同的
聚合規則是否影響結果。」

------------------------------------------------------------
方法
------------------------------------------------------------
1. 用你原本 GSE126848 的 raw counts + mapped.xlsx（跟
   patient_bootstrap_resampling 用的同一組檔案），跑「一次」（不重抽樣）
   完整的 PyDESeq2 分析，取得每個基因真實的 log2FC 點估計值。

2. 對 16 個多基因節點（adipoR, Akt, CxIV, AP_1, GSK_3, IL_1, IRS_1_2,
   p38, JNK1_2, ChREBP, CxI, NF_kB, PI3K, AMPK, CxII, CxIII），分別用
   4 種聚合方法重新計算節點初始 ratio：
     mean            — manuscript 原本方法（線性空間算術平均，即 v10 Mean）
     geomean         — 幾何平均（log2FC 空間先平均、再轉換回線性空間）
     median          — 中位數
     max_abs_log2fc  — 取 |log2FC| 最大的那個成員基因當代表值（「最具
                        差異表現的代表基因」，是比 KEGG 任意排序的
                        「First」更有原則性的替代方案）

3. 用 4 種聚合結果分別重跑 manuscript 同一套 Table 4 指標（NASH/NAFL
   %reduction、advantage ratio），比較是否穩健。

單基因節點（42 個）在 4 種方法下完全不變（沒有聚合選擇的問題），只有
16 個多基因節點的初始值會隨方法改變。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
alt_node_aggregation/
├── altagg_common.py         <- 共用函式（含 DESeq2 pipeline 與 4 種聚合方法）
├── node_gene_map.py         <- 58 節點 -> 基因清單對照表（與 bootstrap 相同）
├── run_alt_aggregation.py   <- 主程式
├── requirements.txt
├── models/
│   ├── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
│   └── ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py   <- 已內建
└── data/                     <- ***需自行放入以下兩個檔案***
    ├── GSE126848_Gene_counts_raw.txt
    └── GSE126848_Gene_counts_mapped.xlsx

這兩個資料檔案跟 patient_bootstrap_resampling 用的完全一樣，可以直接複製
那邊 data/ 資料夾裡的檔案過來，不需要重新準備。

------------------------------------------------------------
執行方式
------------------------------------------------------------
    conda activate deseq2_env
    python run_alt_aggregation.py

約 1-2 分鐘跑完（主要時間花在跑一次 DESeq2，ODE 模擬部分很快）。

------------------------------------------------------------
已知結果（雲端測試，供比對）—— 對你有利的穩健性結果
------------------------------------------------------------
全部 4 種聚合方法、全部 3 個核心 output，NAFL/NASH advantage ratio 全部
維持 >1：

  Method            P_Cell_death   P_Hepatocyte_injury   P_Inflammation
  mean（原本方法）    2.458          2.105                  3.089
  geomean            2.504          2.105                  2.995
  median             2.406          2.105                  2.874
  max_abs_log2fc     2.171          2.104                  3.535

  P_Cell_death           全部 4 方法皆 >1，range=[2.171, 2.504]
  P_Hepatocyte_injury    全部 4 方法皆 >1，range=[2.104, 2.105]（幾乎不變！）
  P_Inflammation         全部 4 方法皆 >1，range=[2.874, 3.535]

值得注意：P_Hepatocyte_injury 的 ratio 在 4 種方法下幾乎完全不變
（2.104-2.105），因為它的兩個主要驅動者 CASP3、CASP7 都是「單基因節點」，
完全不受聚合方法選擇影響——這點可以直接說明為什麼這個 output 對這項
特定的方法學選擇特別穩健。

------------------------------------------------------------
⚠️ 一個小的技術性差異，需要在文字中誠實說明
------------------------------------------------------------
這裡用標準設定（refit_cooks=True）重新跑的「mean」基準結果
（NASH%=5.568，manuscript Table 4 是 5.669），跟 manuscript 正式報告的
數字有約 1.8% 的微小落差。這是因為 DESeq2 的細部參數設定（例如
refit_cooks 開關）不同所致，**不影響「跨聚合方法比較」這個分析本身的
有效性**，因為 4 種方法都在同一套自洽的 pipeline、同一次 DESeq2 fit
下比較。如果要在 Response to Reviewers 中引用，建議說明這裡報告的是
「相對於同一 pipeline 內 mean 方法的變化幅度」，而非要取代 manuscript
正式的 Table 4 數字。

------------------------------------------------------------
建議在 Response to Reviewers 中的寫法
------------------------------------------------------------
- 直接、正面回應：測試了 4 種聚合規則（含 manuscript 原本的算術平均、
  幾何平均、中位數、以及以最大差異表現基因為代表的替代方案），核心的
  NAFL/NASH advantage 方向在所有測試方法下完全穩健維持
- 特別說明：P_Hepatocyte_injury 因為其兩個主要驅動者（CASP3、CASP7）
  都是單基因節點，完全不受聚合方法選擇的影響，ratio 幾乎維持不變
  （2.104-2.105），這是這個 output 相對於聚合方法選擇特別穩健的具體證據
- 誠實補充：其餘兩個 output（P_Cell_death、P_Inflammation）的 ratio
  絕對量級會隨聚合方法在 2.17-2.50、2.87-3.54 之間變動，但方向本身
  不受影響

Cohort-level Non-independence Tests (R1-8b / R3-6)
====================================================
用途：回應 npj SBA
  - Reviewer 1, Major Comment 8(b)：separation index 與 intervention advantage
    是用同一組模擬 AUC 算出的，可能有數學耦合；10% threshold 是 post-hoc。
  - Reviewer 3, Major Comment 6：主要檢定是 18 個「非獨立」觀測值（同一 cohort
    的 3 個 output 共用同一個模擬/轉錄體），應該用考慮 cohort 內非獨立性的方法
    （cohort-level permutation 或 mixed-effects model）重新檢定。

兩件事完全不需要重跑 DESeq2 或 ODE 模擬，資料/腳本都已內建，執行只需幾秒鐘。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
cohort_permutation_test/
├── table4_directional_test.py          <- 檢定 1：Table 4 方向性結論
├── requirements.txt
└── si_advantage_relationship/          <- 檢定 2：SI 預測介入優勢的關聯
    ├── permutation_core18.py               (3 核心 output x 6 cohort, N=18)
    ├── permutation_SI_advantage.py         (11 個 disease-positive output x 6 cohort)
    └── Table_AllOutputs_SeparationIndex_vs_AUCratio.xlsx   (已算好的資料，不需重算)

------------------------------------------------------------
檢定 1：table4_directional_test.py
------------------------------------------------------------
問題：manuscript 的基本方向性主張——「NAFL/NASH 早期介入優勢比值系統性 >1」
（Table 4，3 核心 output × 6 cohort）——如果把每個 cohort 內的 3 個 output
當作一個不可分割的觀測單位（因為它們共用同一個轉錄體/模擬），統計上還站得住嗎？

方法：
  - 對 6 個 cohort 做 exact sign-flip permutation（2^6=64 種排列組合，
    每次同時翻轉一整個 cohort 的 3 個 output，而非個別翻轉）
  - Wilcoxon signed-rank（6 個 cohort-level 平均 log-ratio vs 0）
  - Mixed-effects model（cohort 為 random intercept）作為參考，但因僅 6 個
    cluster、變異數估計在邊界，此模型的 P 值不應引用

執行：
    cd cohort_permutation_test
    python table4_directional_test.py

已知結果（雲端跑過）：
  Exact cohort-level sign-flip permutation: one-sided P = 0.0156（mean 與
    median 統計量皆同，這是 n=6 cohort 下可能達到的最小 P 值，因為 6 個
    cohort 的 log-ratio 總和皆為正）
  Wilcoxon signed-rank (n=6): P = 0.0156
  → 結論：方向性主張穩健，即使用最保守的「整個 cohort 視為一個觀測值」
    的檢定方式依然顯著。

------------------------------------------------------------
檢定 2：si_advantage_relationship/ 內兩支腳本
------------------------------------------------------------
問題：manuscript 宣稱「separation index (SI) 能預測 intervention advantage」
（Fig. 8, Spearman r=0.47 P=0.048；Mann-Whitney P=0.0059）。這個檢定用的
naive Spearman/Mann-Whitney 把 18 或 63 個 cohort-output 組合當作獨立觀測，
但同一 cohort 的多個 output 共用同一組模擬，並不獨立。

方法：Cohort-block permutation——只在同一 cohort 內部打亂 output 標籤，
保留 cohort 間的結構，重新估計 P 值。

執行：
    cd cohort_permutation_test\si_advantage_relationship
    python permutation_core18.py
    python permutation_SI_advantage.py

已知結果（雲端跑過，seed=42, N_PERM=20000）：
  [核心 N=18 (3 output x 6 cohort)]
    naive Spearman r=0.4716 (asymptotic P=0.048，與 manuscript 一致)
    global permutation P=0.024
    *cohort-block permutation P=0.209（不顯著）*
    naive Mann-Whitney (asymptotic) P=0.0059（與 manuscript 一致）
    *cohort-block permutation P=0.074（邊緣/不顯著）*

  [延伸 N=63 (11 output x 6 cohort)]
    naive Spearman r=0.15 (P=0.239，本來就不顯著)
    cohort-block permutation P=0.187
    naive Mann-Whitney asymptotic P=0.072
    cohort-block permutation P=0.073

→ 結論（重要，需要你決定如何處理）：一旦正確考慮 cohort 內非獨立性，
  「SI 預測介入優勢」這個較強的機制性主張（Results/Fig.8 的核心論點）
  **無法通過 cohort-block permutation 檢定**（P 從 0.0059/0.048 上升到
  0.07–0.21）。但檢定 1 的基本方向性主張（NAFL/NASH ratio 系統性 >1）
  依然穩健。

  建議：
  (a) 在 Response to Reviewers 誠實同時報告 naive 與 cohort-block
      permutation 兩種結果
  (b) 考慮把 manuscript 中 "predicts" 這類強機制性用詞，針對 SI-advantage
      關係軟化為 hedged 說法，並在 Limitations 承認此關聯在校正
      非獨立性後不穩健
  (c) 基本方向性主張（Table 4-based）可以維持現有措辭，因為統計上站得住

------------------------------------------------------------
注意
------------------------------------------------------------
- 兩支 permutation 腳本用 N_PERM=20000，執行時間約數秒到十幾秒。
- si_advantage_relationship 資料夾內的 xlsx 是根據你上傳的 6 個 cohort
  ratio 檔案 + 主模型 ODE 腳本產生，已與 manuscript Table 4 數字逐一核對
  一致（見對話紀錄），不需重新產生；如果你更新了任何一個 cohort 的
  ratio 檔案，需要另外重跑 quantify_all_outputs_separation_vs_intervention.py
  （此檔未包含在此包裡，如需要請告訴我另外打包）。

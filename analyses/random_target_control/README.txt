Random-target / Matched-control Analysis (R1-7 / R2-4 / R3-5)
================================================================
用途：回應 npj SBA
  - Reviewer 1, Major Comment 7：silymarin 8 個靶點是否只是「剛好」會降低
    injury outputs（因為多個靶點本來就是 output-proximal 節點）；建議做
    random target 的 background distribution 並報告 FDR/percentile 顯著性。
  - Reviewer 2, Major Comment 4：預測的反應是否 silymarin 特異，還是同時
    抑制任意 8 個 output-proximal 或高連結度節點就會有類似效果？
  - Reviewer 3, Major Comment 5：4/5 主要 sensitivity driver 同時也是
    silymarin 靶點，可能只是反映「研究熱點重疊」；建議對「既非靶點、也非
    top driver」的 matched control 節點集做相同介入，證明 window 是
    藥物靶點特異的。

------------------------------------------------------------
方法
------------------------------------------------------------
使用完全真實的 hsa04932 拓樸、真實動力學參數（Vmax=ksp=n=kcat=2.0）、
真實的 GSE126848 主 cohort 初始條件（直接讀取模型檔內建的
GENE_METADATA[node]['initial_ratio']，與 manuscript Table 2/S1、Table 4
完全一致）——不對網路做任何 scramble（那是 edge_scramble_control 在做的
事）。唯一改變的變數是：**被抑制（ki=0.3）的是哪 8 個節點**。

兩個 null model：

1. random_target：從全部 58 個分子節點中隨機抽 8 個做抑制，重複多次，
   建立「任意 8 節點」的 null distribution。回應 R1-7 / R2-4。

2. matched_control：從 58 個節點中先剔除
   - 8 個真實 silymarin 靶點：CASP3, CASP7, CASP8, CYP2E1, IL_8, TNFa,
     NF_kB, TGF_b1
   - Fig.6 報告的 top sensitivity driver：CASP7, CASP3, Cytc, Bax, TNFa,
     FasL, IL_8, TGF_b1
   （聯集共 11 個節點排除，剩 47 個節點的 pool），再從中隨機抽 8 個做抑制。
   直接回應 R3-5「排除 targets 和 top drivers 的 matched control」的要求。

因為未給藥（no-drug）的 baseline 軌跡不受「抑制哪些節點」影響，程式只在
每次批次執行開始時算一次 baseline、重複利用，所以每個 replicate 只需要
2 次新的 ODE 求解（NASH-drug、NAFL-drug），速度比 patient_bootstrap 或
edge_scramble 快很多（實測約 1 秒/replicate）。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
random_target_control/
├── common.py              <- 共用函式（模型載入、Y0、通用 ki 抑制注入等）
├── step1_run_batch.py     <- 主要執行腳本
├── step2_analyze.py       <- 統計分析（差值/比值檢定、joint test）
├── requirements.txt
├── models/
│   ├── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
│   └── ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py   <- 已內建
└── outputs/                <- 執行後結果會寫在這裡

------------------------------------------------------------
環境設定
------------------------------------------------------------
    conda activate deseq2_env
    pip install -r requirements.txt

（deseq2_env 應該已經有 numpy/pandas/scipy；若沒有 openpyxl 才需另外裝。
本分析完全不需要 pydeseq2，因為不做任何 DESeq2 重跑。）

------------------------------------------------------------
執行方式
------------------------------------------------------------
1. 先算一次「真實 8 靶點」的參考結果（幾秒鐘完成，會自動核對 manuscript
   Table 4 數字是否吻合作為 sanity check）：

       cd path\to\random_target_control
       python step1_run_batch.py true

   應該會印出（與 manuscript Table 4 GSE126848 完全一致）：
       P_Cell_death           ratio=2.44
       P_Hepatocyte_injury    ratio=2.10
       P_Inflammation         ratio=3.03

2. 跑 random_target null（建議總數 2000，約 30-40 分鐘；可一次跑完或分批）：

       python step1_run_batch.py random_target 0 2000

   分批範例（支援中斷續跑，已完成的 replicate 會自動跳過）：
       python step1_run_batch.py random_target 0 500
       python step1_run_batch.py random_target 500 1000
       python step1_run_batch.py random_target 1000 1500
       python step1_run_batch.py random_target 1500 2000

3. 跑 matched_control null（同樣建議 2000）：

       python step1_run_batch.py matched_control 0 2000

4. 兩個 null 都跑完（或先跑一部分想看目前結果）後：

       python step2_analyze.py

   會印出：
   - 兩個 null model 各自的邊際檢定（diff 與 ratio 兩種指標）：
     true value、null median、方向一致比例、P(>=true) 及其 95% CI
   - Joint test（三個 output 同時達到/超過真實效應量的比例）
   - Min-diff-across-3 統計量及其 P 值
   並輸出 outputs/random_target_control_summary.xlsx 彙總表。

------------------------------------------------------------
如何解讀（結合已完成的 edge_scramble_control 一起看）
------------------------------------------------------------
- 若 random_target 的邊際 P(>=true) 偏大（例如 >0.05），代表單看某一個
  output，隨機抽 8 個節點做抑制，有不小機率達到跟真實 silymarin 一樣強的
  效應——這符合 R1-7/R2-4 的疑慮：部分效應量可能只是「抑制夠多節點、
  夠靠近 output」的一般現象，不是 silymarin 特異的。
- 若 joint test（三個 output 同時達到真實效應量）的 P 值遠小於邊際 P
  值，代表雖然單軸效應量不算稀有，但「三軸同時、協調一致」地達到這個量級
  是 silymarin 這組特定靶點組合特有的——這與 edge_scramble_control 的
  發現模式一致（單軸邊際不夠強，但 joint pattern 顯著），可以在 Response
  to Reviewers 中兩個控制分析一起呈現，形成一致的論述。
- matched_control（排除已知 driver）通常會比 random_target 更嚴格：如果
  排除已知 driver 後，null median 明顯下降、真實效應量的 P(>=true) 更小，
  代表 silymarin 之所以有效，部分原因確實來自於它剛好命中了高 sensitivity
  的節點（CASP3/7、TNFa、CXCL8、TGF_b1 等）——這點需要在 Discussion
  誠實承認（manuscript 本來就已經在 Discussion 提到這個 ascertainment
  bias 的可能性），而不是說 silymarin 的選擇完全是任意的。

------------------------------------------------------------
注意事項
------------------------------------------------------------
- 隨機抽樣使用固定 seed（BASE_SEED 依 mode 不同，定義在
  step1_run_batch.py 開頭），同一台機器重跑會得到完全相同的抽樣序列，
  具重現性。
- 若某個隨機抽到的 8 節點組合導致 ODE 積分數值不穩定（NaN/Inf），該筆
  record 會被標記 nan_flag=True 並在 step2_analyze.py 自動排除，不會
  中斷整體流程。
- N=2000 是建議值，足以得到穩定的 P(>=true) 估計（即使 P 很小，
  1/2000=0.0005 的解析度足夠寫在 manuscript 裡；若要更精細，可以加大
  N，程式支援分批續跑，直接把 <end> 改大再執行一次即可）。

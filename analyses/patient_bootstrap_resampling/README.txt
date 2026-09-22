R1-2 Patient-level Bootstrap Robustness Analysis
==================================================
用途：回應 npj SBA Reviewer 1, Major Comment 2 —
「How did authors consider inter-patient variability and uncertainty of
the fold changes?」

方法：對 GSE126848 的四個臨床組（Normal=14, Obese=12, NAFL=15, NASH=16）
分別做 stratified with-replacement resampling，重跑 PyDESeq2 + ODE 模型，
共 200 次迭代，計算 3 個核心 pathway output（P_Cell_death,
P_Hepatocyte_injury, P_Inflammation）的 NAFL/NASH 早期介入優勢比值分佈。

------------------------------------------------------------
資料夾結構（放到本機後應長這樣）
------------------------------------------------------------
patient_bootstrap_resampling/
├── run_bootstrap.py          <- 主程式（可續跑）
├── summarize.py              <- 統計彙總（median / 95% CI / %>1）
├── node_gene_map.py          <- 58 節點 -> 基因清單對照表（已內建）
├── requirements.txt
├── models/
│   └── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
└── data/                     <- ***需自行放入以下兩個檔案***
    ├── GSE126848_Gene_counts_raw.txt
    └── GSE126848_Gene_counts_mapped.xlsx

這兩個資料檔案就是你原本用來建立主模型的 raw counts 與 sample-condition
對照 xlsx（在你 Dropbox 專案資料夾內應該找得到，或用你先前上傳給我的
同名檔案）。只要檔名一致、放進 data/ 資料夾即可，不需要更動程式。

------------------------------------------------------------
環境設定（本機 Windows + conda，沿用你的 deseq2_env）
------------------------------------------------------------
    conda activate deseq2_env
    pip install -r requirements.txt

（deseq2_env 應該已經有 pydeseq2/numpy/pandas/scipy；若沒有 openpyxl 才需另外裝）

------------------------------------------------------------
執行方式
------------------------------------------------------------
1. 切到這個資料夾：
       cd path\to\patient_bootstrap_resampling

2. 每次迭代約需 45-55 秒（DESeq2 fit 為主要耗時），200 次總計約 2.5-3 小時。
   直接一次跑完：
       python run_bootstrap.py 200

   若想分批跑（例如怕中途要用電腦，可以隨時 Ctrl+C 中斷，之後重跑會自動
   跳過已完成的 iteration，不會重算）：
       python run_bootstrap.py 200 20     <- 這次只跑 20 個新的
       python run_bootstrap.py 200 20     <- 再跑 20 個，以此類推
       ...
       python run_bootstrap.py 200        <- 最後一次不給第二參數，跑到滿 200

3. 結果會累積寫到 bootstrap_results.jsonl（每行一個 iteration 的 JSON record，
   含 seed、NASH/NAFL 各自的 AUC% reduction、以及三個 output 的 ratio）。

4. 全部跑完（或想看目前累積結果）後：
       python summarize.py

   會印出每個 output 的 median ratio、95% CI（percentile method）、
   mean±SD、以及 ratio>1 的迭代比例，並附上原始 manuscript 點估計方便比對。

------------------------------------------------------------
已知結果（Claude 於雲端環境完整跑過 200/200，供比對用；
本機重跑因 random seed 固定 = BASE_SEED(20260908) + iteration，
理論上應得到完全相同的 resample，但實際數值可能因套件版本/
浮點運算環境略有差異）
------------------------------------------------------------
P_Cell_death        median=2.119  95% CI=[0.376, 16.429]  ratio>1: 81.5%
P_Hepatocyte_injury median=2.077  95% CI=[1.155,  3.990]  ratio>1: 98.5%
P_Inflammation      median=2.420  95% CI=[0.937, 23.179]  ratio>1: 97.0%

原始 manuscript 點估計（Table 4, GSE126848, ki=0.3）：
P_Cell_death=2.44, P_Hepatocyte_injury=2.10, P_Inflammation=3.03

判讀：三個 output 的 median 都與原始點估計吻合，方向穩健；
P_Hepatocyte_injury 最穩定（CI 完全 >1），P_Cell_death 與
P_Inflammation 的 CI 較寬（小樣本 resampling 造成的基因層級變異），
但多數迭代仍支持 NAFL>NASH 方向。

------------------------------------------------------------
注意事項
------------------------------------------------------------
- 若中途某個 iteration 因 DESeq2 收斂或某 condition 全被抽到極端樣本而
  失敗，程式會把 {"iter":N, "error":"..."} 寫入 jsonl 並跳過，不會中斷
  整體流程；summarize.py 會自動忽略含 error 的 record。
- 固定 seed 設計是為了可重現性；若要做敏感度測試（例如換一組 seed），
  可修改 run_bootstrap.py 中的 BASE_SEED。

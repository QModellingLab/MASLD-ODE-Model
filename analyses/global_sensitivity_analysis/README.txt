Global Sensitivity Analysis + Topology Comparison (R1-6)
==========================================================
用途：回應 npj SBA Reviewer 1, Major Comment 6——
「這些 sensitivity 結果跟純拓樸排序（這些基因本來就是 output 的直接上游
activator）有什麼不同？需要 global sensitivity analysis，並報告這些基因
在所有基因中的排名與顯著性/信賴區間。」

------------------------------------------------------------
方法：為什麼選 Morris 而非 Sobol
------------------------------------------------------------
Manuscript master notes 原本預估 global sensitivity 要花 1-2 天，那是假設
用 Sobol variance-based indices（需要 (2k+2)×N 次模擬，58 個參數、N=1024
時要跑 12 萬次以上）。這裡改用 **Morris elementary-effects screening**
（SALib），只需要 R_TRAJ×(58+1) 次模擬，同樣能：
  (a) 同時變動全部 58 個節點的初始條件（而非 manuscript 原本 local
      sensitivity 的「一次只動一個節點 ±1%」）
  (b) 探索更寬的範圍（節點真實 NASH baseline 值的 0.3-3.0 倍，而非 ±1%）
  (c) 透過 SALib 內建的 bootstrap resampling 給出 mu*（整體重要性）的
      95% 信賴區間，直接回應 reviewer 要的「p-value/confidence」

三個核心 output（P_Cell_death、P_Hepatocyte_injury、P_Inflammation）各自
獨立計算 Morris mu*、mu、sigma、信賴區間，並報告 CASP3/CASP7/TNFa/FasL/
IL_8/TGF_b1/Cytc/Bax 這些「有興趣的」節點在全部 58 個節點中排第幾名。

第二支腳本額外做「與純拓樸排序比較」：從 Supplementary Table S2 的 106
條交互作用建立真實網路，計算每個節點的 out-degree（下游連結數）以及到
每個 output 的最短路徑距離（hop 數），再用 Spearman 相關係數（附
p-value）檢驗 Morris 排名是否可以單純用「這個基因連得多」或「這個基因
離 output 很近」來解釋。

------------------------------------------------------------
需要跑多久？可以 GPU 加速嗎？
------------------------------------------------------------
**不能用 GPU 加速。** 瓶頸是 scipy 的 odeint（LSODA 求解器），這是一個
140 維、中度剛性（stiff）的小型 ODE 系統，用的是序列型 Fortran 常微分
方程求解器，沒有現成的 GPU 實作可以直接套用；要上 GPU 需要把整個模型
改寫成向量化的 JAX/PyTorch ODE 求解器，工程量遠超過這次控制分析的效益。

**但可以用多核心 CPU 平行化，這支程式已經內建這個功能。** 因為每組隨機
抽樣的參數組合彼此獨立，天生就是「embarrassingly parallel」的問題。
`run_morris_gsa.py` 用 Python 的 ProcessPoolExecutor，預設會抓你電腦的
全部邏輯核心數平行跑（可在檔案開頭手動指定 N_WORKERS）。

單次 ODE 模擬（300h）在雲端測試環境約 0.5-0.7 秒。預設 R_TRAJ=100，
代表 100×(58+1)=5900 次模擬：

| 核心數 | 預估時間 |
|--------|----------|
| 1 核心 | 約 55-65 分鐘 |
| 4 核心 | 約 15-18 分鐘 |
| 8 核心 | 約 8-10 分鐘 |
| 16 核心 | 約 4-5 分鐘 |

（你電腦的實際核心數可以在工作管理員「效能」頁看到，或執行
`python -c "import os; print(os.cpu_count())"`）

若想要更精細的信賴區間，可以把 R_TRAJ 調大（例如 200），時間會等比例
拉長；若只是想先看初步結果，調小到 30-50 也可以（信賴區間會較寬，但
排名的大方向通常已經穩定）。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
global_sensitivity_analysis/
├── gsa_common.py                   <- 共用函式（獨立命名，不會跟其他資料夾撞名）
├── edges_table_s2.py               <- Table S2 全部 106 條交互作用（已驗證節點名稱）
├── run_morris_gsa.py               <- 步驟 1：Morris 分析主程式
├── step2_topology_comparison.py    <- 步驟 2：拓樸比較
├── run_all.py                      <- 一鍵執行（依序跑步驟 1+2）
├── requirements.txt
└── models/
    └── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建（只需 NASH，
        與 manuscript Fig.6 local sensitivity 同一個 baseline，方便直接比較）

------------------------------------------------------------
環境設定
------------------------------------------------------------
    conda activate deseq2_env
    pip install -r requirements.txt

（deseq2_env 應該已經有 numpy/pandas/scipy；SALib 和 networkx 大機率需要
另外安裝。）

------------------------------------------------------------
執行方式（最簡單：一鍵跑完）
------------------------------------------------------------
打開 run_all.py，按 F5。會依序執行：
  1. Morris 抽樣 + 平行跑模擬（可續跑，中途中斷重新按 F5 會跳過已完成的）
  2. 分析 + 存 Excel（outputs/global_sensitivity_summary.xlsx）
  3. 拓樸比較 + 存 Excel（outputs/topology_comparison_summary.xlsx）

若想分開跑（例如先確認模擬跑完、隔天再跑拓樸比較），也可以：
    python run_morris_gsa.py
    python step2_topology_comparison.py

------------------------------------------------------------
如何調整參數
------------------------------------------------------------
在 run_morris_gsa.py 檔案開頭：
    R_TRAJ = 100          # 改這個數字調整精細度/時間
    N_WORKERS = None      # None = 用全部核心；也可以填數字強制指定
    RNG_SEED = 20260908   # 固定 seed，重跑會得到完全相同結果（可重現性）

------------------------------------------------------------
已知結果（雲端用小樣本 R_TRAJ=5 測試過，僅供驗證程式正確性，不是正式
結果——正式結果請以你本機用 R_TRAJ=100 跑出來的為準）
------------------------------------------------------------
- P_Hepatocyte_injury：top mu* 節點為 CASP3、Cytc、CASP7、Bax，與
  manuscript Fig.6 local sensitivity 排序高度一致（CASP7、CASP3 為
  兩大主要驅動者）
- P_Cell_death：top mu* 節點為 FasL、TNFa，與 manuscript 一致
- P_Inflammation：top mu* 節點為 IL_8，與 manuscript 一致
- 拓樸比較（小樣本初步結果，需要正式 R_TRAJ=100 跑完確認）：
  - out-degree 與 mu* 相關性普遍偏弱（rho≈0.15-0.26），代表排序不是
    單純「連結數多寡」能解釋
  - 但 P_Hepatocyte_injury 的 mu* 與「到 output 的距離」相關性極強
    （rho≈0.95，p<0.0001）——這點誠實地印證了 reviewer 的疑慮：對這個
    特定 output 而言，排序確實高度可以用「離 output 多近」預測，因為
    CASP3/CASP7 本來就是距離=1 的直接上游
  - P_Inflammation 則相關性較弱、不顯著——顯示不是每個 output 都有
    這個問題

------------------------------------------------------------
建議在 Response to Reviewers 中的寫法
------------------------------------------------------------
- 誠實報告：對 P_Hepatocyte_injury，global sensitivity 排名確實與拓樸
  距離高度相關（rho≈0.95），承認這個 output 的 sensitivity 結果很大程度
  上可以用「這些基因是直接上游」這個拓樸事實預測，這與 reviewer 的觀察
  一致
- 但同時報告：這個現象並非全面性的——P_Cell_death、P_Inflammation 的
  拓樸相關性明顯較弱，代表並非所有 output 的 sensitivity 排序都只是
  拓樸的平凡結果
- 補充：global sensitivity（Morris，大範圍、多變量同時擾動）與
  manuscript 原本的 local sensitivity（±1%、單變量）排序高度一致，
  說明原始結論在更嚴謹的方法下依然穩健，這是額外的正面佐證
- 附上所有 58 個節點的完整排名表（Excel）與信賴區間，回應 reviewer
  對「p-value/confidence」的具體要求

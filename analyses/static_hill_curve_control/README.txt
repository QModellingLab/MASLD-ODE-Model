Static Hill-curve Control (R3-1, second half)
================================================
用途：回應 npj SBA Reviewer 3, Major Comment 1 的第二個具體要求——
「a static Hill-curve calculation that predicts the effect without
integrating the ODEs」。edge_scramble_control 資料夾已經做了同一則意見
要求的第一個控制（degree-preserving edge-scramble）；這裡做的是「不依賴
AUC 時間積分窗口」的靜態/穩態版本控制。

------------------------------------------------------------
這支分析做兩件事（一次執行、幾十秒到一兩分鐘內完成，不需要跑大量
replicate，因為這是確定性計算，不是統計抽樣）
------------------------------------------------------------

PART A - 網路穩態（root-finding，不依賴積分窗口）
  對 NASH / NAFL 各自的「無藥 / silymarin ki=0.3」四種組合，用
  scipy.optimize.fsolve 直接對 dy/dt=0 求根，找出整個 70 節點網路的
  真正平衡點（穩態），而不是像 manuscript 那樣在 t=0-300h 這個特定窗口
  內做 AUC 積分。用一段較長時間的 ODE 模擬只是拿來當 fsolve 的初始猜測
  （純數值技巧，方便收斂），最終報告的是 root-finder 收斂後的解，並印出
  殘差 norm（||dy/dt||）以確認真的收斂到穩態。

  比較穩態算出的 %reduction / NAFL-NASH ratio，跟 manuscript AUC-based
  的 Table 4 數字是否方向、量級一致——如果穩態（跟積分窗口完全無關）也
  出現同樣的 NAFL>NASH 優勢，代表這個效應不是 AUC 窗口選擇的人工產物。

PART B - 單一反應的代數 Hill / receptor-reserve 檢驗
  完全不碰 ODE。針對三個核心 output 各自的直接上游 driver（Table S3：
  CASP3/CASP7 → P_Hepatocyte_injury；FasL/TNFa → P_Cell_death；
  IL_8/TGF_b1 → P_Inflammation），直接用代數 Hill 公式
  h(x)=x^n/(ksp^n+x^n)（n=ksp=2.0，跟 manuscript 一致）在該 driver 的
  NAFL vs NASH 初始 ratio（GENE_METADATA 內建值）上算「reserve」=1-h(x)，
  也就是這個 driver「自己的」活化反應在該疾病期還有多少未飽和的空間。
  這是 manuscript 已經引用的 receptor-reserve 文獻（Buchwald 2020;
  Wilder 1957）最單純、最字面意義上的代數版本，用來檢驗「NAFL 比較有效」
  這件事是否只是任兩個數字（NAFL ratio < NASH ratio）代入 Hill 公式後
  必然出現的結果。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
static_hill_curve_control/
├── common.py                  <- 共用函式（與 random_target_control 相同）
├── run_static_analysis.py     <- 主程式（單一檔案，一次執行完 Part A+B）
├── requirements.txt
└── models/
    ├── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
    └── ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py   <- 已內建

------------------------------------------------------------
執行方式
------------------------------------------------------------
    conda activate deseq2_env
    pip install -r requirements.txt     （通常 numpy/scipy 都已經有了）
    python run_static_analysis.py

直接在 Spyder 開啟 run_static_analysis.py 按 F5 也可以，不需要任何參數、
不需要批次執行，整個分析在雲端測試約 1-2 分鐘內完成（4 次 fsolve 求解 +
6 個代數計算）。

------------------------------------------------------------
已知結果（雲端測試過，供比對；本機重跑數字應該完全一致，因為整個
分析是確定性的，沒有隨機抽樣成分）
------------------------------------------------------------

PART A - 穩態 vs AUC-based 比較：
  Output                  SS ratio    AUC ratio(manuscript Table 4)
  P_Cell_death              2.036            2.44
  P_Hepatocyte_injury       1.782            2.10
  P_Inflammation            1.991            3.03

  → 三個 output 在穩態下的 NAFL/NASH 比值都 >1，方向與 AUC 結果一致，
    但量級略小於 AUC 值。代表：
    (a) NAFL>NASH 優勢不是 AUC 積分窗口（0-300h）的人工產物——即使看
        完全跟時間窗口無關的穩態，優勢依然存在，這是好消息，可以直接
        回應「Please add two controls」的字面要求。
    (b) 但穩態的量級普遍小於 AUC，代表 AUC 報告的優勢有一部分（並非
        全部）來自暫態動力學過程本身，而非單純穩態終點的差異——這點
        在 Response to Reviewers 應誠實說明，manuscript 用 AUC 而非
        穩態終值作為主要指標的理由（本來就是為了捕捉飽和速率的差異，
        而非終值本身）在這裡得到進一步支持。

PART B - 代數 reserve 檢驗（節錄）：
  Driver   reserve(NAFL)   reserve(NASH)   NAFL>NASH?
  CASP3       0.6947          0.7358          no
  CASP7       0.4773          0.4591          YES
  FasL        0.8684          0.7479          YES
  TNFa        0.6100          0.6273          no
  IL_8        0.3201          0.0609          YES（差異很大）
  TGF_b1      0.9578          0.9162          YES

  → 結果分歧、不是全面一致：6 個 driver 中 4 個支持 NAFL 有較大 reserve，
    2 個（CASP3、TNFa）方向相反。這代表「NAFL 比較有效」這件事**不是**
    任兩個 stage-specific 數字代入 Hill 公式後必然、trivial 就會出現的
    結果（如果是 trivial 的，6 個 driver 應該全部同方向）。
    IL_8 的 reserve 差異特別懸殊（0.32 vs 0.061），提示 P_Inflammation
    這個 output 的優勢有相當大一部分可能直接來自 CXCL8 這個節點自己的
    飽和程度差異，而不是網路動態放大的結果——這點值得在 Discussion
    誠實補充說明，呼應 manuscript 原本就有的 ascertainment-bias 討論。

------------------------------------------------------------
建議在 Response to Reviewers 中的寫法
------------------------------------------------------------
- 誠實呈現 Part A：window-free 的穩態計算證實方向性優勢不依賴 AUC 積分
  窗口選擇，但量級部分歸因於暫態動力學。
- 誠實呈現 Part B：driver 層級的代數 reserve 檢驗顯示結果並非全面一致
  地支持「trivial algebra」的解釋（6 個 driver 中 2 個方向相反），但也
  承認 IL_8/CXCL8 這個節點的效應主要可能來自其自身的飽和程度差異。
- 這樣的呈現方式跟 edge_scramble_control 的模式一致：坦承部分機制可以
  用簡單的 Hill 飽和 / receptor-reserve 解釋，但完整的三軸協調模式、
  以及某些 driver（如 CASP7、FasL）的行為，仍然需要完整的網路動態才能
  重現，不是單一代數步驟就能保證的。

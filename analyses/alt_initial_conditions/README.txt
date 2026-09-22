Alternative Initial-Condition Assumptions (R1-4)
===================================================
用途：回應 npj SBA Reviewer 1, Major Comment 4——
「為什麼 mRNA 表現量可以代表 inactive 蛋白濃度？為什麼在已經確立的
NAFL/NASH 組織中，所有訊號節點的初始活性都假設為 0？pathway output 的
inactive 值設為 100 也沒有生物學解釋。這些假設可能強烈影響模擬軌跡，
應該用替代的初始設定測試。」

------------------------------------------------------------
測試的兩個替代假設
------------------------------------------------------------
1. active_fraction：manuscript 預設所有分子節點的 active form 都從 0
   開始。這裡改成從「該節點 inactive form 初始值的某個比例」開始
   （0%、10%、30% 三種），模擬「疾病組織中訊號節點本來就有一定基礎活性」
   的替代假設。

2. p_output_inactive_level：manuscript 預設 12 個 pathway output 節點的
   inactive form 固定從 100 開始。這裡測試 50、100、200 三種替代值。

3×3=9 種組合，每組跑 NASH/NAFL × 有無 silymarin 共 4 次 ODE 模擬，用
manuscript 同一套 Table 4 指標（%reduction、NAFL/NASH advantage ratio）
評估結論穩健性。快速、確定性計算，單次執行約 30-40 秒。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
alt_initial_conditions/
├── altic_common.py       <- 共用函式
├── run_alt_ic.py         <- 主程式（單次執行，9 組合）
├── requirements.txt
└── models/
    ├── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
    └── ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py   <- 已內建

------------------------------------------------------------
執行方式
------------------------------------------------------------
    conda activate deseq2_env
    python run_alt_ic.py

（不需要參數，直接 F5 即可，約 30-40 秒跑完全部 9 組合。）

------------------------------------------------------------
已知結果（雲端測試，供比對）—— 對你有利的穩健性結果
------------------------------------------------------------
全部 9 個組合、全部 3 個核心 output，NAFL/NASH advantage ratio 全部
維持 >1（方向完全穩健，沒有任何翻轉案例）：

  P_Cell_death           全部 9 組皆 >1，range=[1.714, 2.680]
  P_Hepatocyte_injury    全部 9 組皆 >1，range=[1.502, 2.105]
  P_Inflammation         全部 9 組皆 >1，range=[2.121, 3.819]

值得注意的模式：
- active_fraction 越大（初始活性越高），效應量級（絕對 %reduction 和
  diff_pp）通常越小——這符合直覺：如果系統已經有一定基礎活性，藥物
  抑制的邊際效果自然較小，但「NAFL 比 NASH 更有效」這個方向性結論本身
  不受影響。
- p_output_inactive_level 越大，效應量級越大，但方向同樣不受影響。

------------------------------------------------------------
建議在 Response to Reviewers 中的寫法
------------------------------------------------------------
- 直接、正面回應：測試了 3×3=9 種替代初始條件組合（涵蓋兩個 reviewer
  質疑的假設：全零起始活性、pathway output 固定基準值 100），核心的
  stage-dependent early-intervention advantage 方向在所有測試組合下
  完全穩健維持（100% 的 27 筆 output×組合觀測皆為 ratio>1）
- 誠實補充：效應的絕對量級確實會隨這些假設變動而改變（尤其
  active_fraction 越高、量級越小），但這是預期中的、符合機制的模式，
  不影響核心的方向性 claim
- 可以引用這個結果直接反駁 reviewer 對「假設可能強烈影響模擬軌跡」的
  疑慮：軌跡的絕對數值確實會變，但支撐主要結論的相對比較（NAFL vs NASH）
  在這個測試範圍內是穩健的

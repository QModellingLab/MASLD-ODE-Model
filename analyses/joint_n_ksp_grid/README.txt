Joint n x ksp Parameter Grid (R3, minor point 3)
===================================================
用途：回應 npj SBA Reviewer 3, Minor Point 3——
「Vary n and ksp jointly rather than one at a time, and report how the
effect magnitude, not only its sign, responds.」

manuscript Supplementary Table S4 只做了「一次改一個參數」的掃描（n=1,2,4
固定 ksp=2；ksp=1,4 固定 n=2，共 5 組，含 primary）。這支程式做完整的
3×3 聯合網格（n∈{1,2,4} × ksp∈{1,2,4} = 9 組，涵蓋原本 5 組 + 4 組全新的
聯合組合），並報告效應「量級」（NASH/NAFL 各自的 %reduction、絕對差值
diff_pp），不只是方向（ratio>1 or <1）。

------------------------------------------------------------
技術實作：如何改動模型的 n 與 ksp
------------------------------------------------------------
manuscript 的 ODE 模型原始碼把所有 Hill 動力學常數（Vmax=ksp=n=kcat=2.0）
硬編碼成字面上的 "2.0"。每個 Hill 反應項的形式是：
    kcat * Regulator * Substrate**2.0 / (2.0**2.0 + Substrate**2.0)
也就是同時出現「2.0**2.0」（ksp 的 n 次方）跟裸露的「Substrate**2.0」
（受質本身的 n 次方項），而「2.0 * Regulator」這種乘法則是 kcat/Vmax，
維持不變（跟 Table S4 原本的做法一致，只動 n 和 ksp）。

程式用兩步驟字串取代：
  1. 先把 "2.0**2.0" 換成 "{ksp}**{n}"（處理分母的 ksp^n 項）
  2. 再把剩下所有 "**2.0" 換成 "**{n}"（處理裸露的受質/產物 n 次方項）

**這個替換方法已經過驗證**：重新跑 (n=2,ksp=2) 精確重現 manuscript
Table 4 的三個 output ratio（2.44/2.10/3.03）；重新跑 (n=1,ksp=2) 和
(n=2,ksp=1) 精確重現 Supplementary Table S4 已發表的數字
（1.006/2.316/4.366 與 2.646/2.584/2.397）——三組全部逐位吻合，證明
參數替換技術正確無誤。

------------------------------------------------------------
資料夾結構
------------------------------------------------------------
joint_n_ksp_grid/
├── nksp_common.py       <- 共用函式（含已驗證的 parameterize_nksp）
├── run_grid.py          <- 主程式（單次執行，9 組合 × 4 次 ODE 求解，約 20-30 秒）
├── requirements.txt
└── models/
    ├── ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py   <- 已內建
    └── ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py   <- 已內建

------------------------------------------------------------
執行方式
------------------------------------------------------------
    conda activate deseq2_env
    python run_grid.py

（不需要參數，直接 F5 即可，約 20-30 秒跑完全部 9 組合。）

------------------------------------------------------------
已知結果（雲端測試，供比對）
------------------------------------------------------------
**⚠️ 重要發現**：在 (n=1, ksp=1) 這個「只有聯合網格才會測到」的組合下，
P_Hepatocyte_injury 的 advantage ratio = **0.952（<1，方向翻轉！）**——
NASH 的介入效果反而略高於 NAFL，這跟 manuscript 主張的「NAFL 早期介入
優勢」方向相反。這是原本 Table S4 的「一次一個參數」掃描完全沒有測到的
組合（因為它只測了 n=1 配 ksp=2、以及 n=2 配 ksp=1，沒測過 n=1 配
ksp=1 這個聯合組合）。

其餘 8 個組合、以及 P_Cell_death、P_Inflammation 兩個 output 在全部
9 個組合中，方向都維持 ratio>1（NAFL>NASH）：
  P_Cell_death           全部 9 組皆 >1，range=[2.168, 2.917]
  P_Hepatocyte_injury    8/9 組 >1，但 (n=1,ksp=1) 例外 =0.952；range=[0.952, 2.646]
  P_Inflammation         全部 9 組皆 >1，range=[1.209, 4.956]

這與 manuscript Methods 已經寫的「P_Cell_death 最穩定、P_Hepatocyte_injury
對參數最敏感」高度一致，只是這次用聯合網格把這個「最敏感」的性質，
具體量化成一個真實會翻轉方向的組合。

------------------------------------------------------------
建議在 Response to Reviewers 中的寫法
------------------------------------------------------------
- 誠實揭露：完整 3×3 聯合網格顯示，P_Hepatocyte_injury 在一個特定的聯合
  參數組合（n=1, ksp=1，即最低 Hill 協同性 + 最低半飽和常數）下，方向
  確實會翻轉（ratio=0.952<1），而這是原本「一次一個參數」的 Table S4
  設計本來就無法偵測到的組合
- 同時報告：P_Cell_death 和 P_Inflammation 在全部 9 個聯合組合下依然
  穩健地維持方向，且 P_Hepatocyte_injury 其餘 8/9 組合也維持方向，只有
  這一個角落案例例外
- 建議在 manuscript 的 Limitations 或 Supplementary 中新增這個聯合網格
  的完整結果（本分析輸出的 Excel 含所有 9×3=27 筆數字），並在文字中
  明確承認這個特定參數組合下的方向不穩定性，同時說明 P_Cell_death 作為
  主要跨 cohort 驗證指標（manuscript 已有此定位）在所有測試過的參數
  設定下都是最穩健的 output，這個定位選擇因此得到進一步支持

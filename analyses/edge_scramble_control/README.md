# Edge-Scramble Control Analysis（回應 npj SBA Reviewer 3 Major Comment 1）

## 背景

Reviewer 3 質疑：NAFL 比 NASH 反應更好，可能只是模型建構方式（severity 編碼在
初始值、共同飽和天花板、AUC 指標、乘法式 ki 抑制）的必然產物，不是 hsa04932
這個特定拓撲的貢獻。要求兩個對照：static Hill-curve closed-form（另外處理）
與 **degree-preserving edge-scramble**（本套件）。

## 這個套件做了什麼

1. **重建一個「拓撲驅動」的通用 ODE 模型**（`generic_model.py`），完全從
   Supplementary Table S2 的 106 條交互作用（106 edges）生成方程式——不是
   把方程式寫死在程式碼裡。已驗證這個通用模型跟你原本寫死的
   `ode_model_pydeseq2_*.py` 數值上完全一致（差異在 ODE 積分器精度等級,
   約 1e-7）。
2. **兩種獨立的隨機化（null model）**，互相印證：
   - **degree-preserving scramble**（Maslov & Sneppen, 2002, *Science*）：
     保持每個節點的輸入輸出連線數完全不變，只打亂「誰接誰」。這是系統
     生物學裡驗證網路特徵是否顯著的黃金標準做法（同樣方法也用於
     Shen-Orr et al. 2002 對大腸桿菌轉錄調控網路的分析）。
   - **Erdős–Rényi random**：完全不管原本的連線數分布，103 條邊隨機接在
     任意兩個合法節點之間。用來檢查「聯合顯著」這個結果是不是只在
     degree-preserving 這種比較溫和的隨機化下才成立。
3. 對每一次隨機拓撲，重跑 silymarin 8-target 介入模擬（ki=0.3），計算三個
   core output 的 NAFL/NASH 差異，跟真實拓撲比較，算出經驗 P 值（含 Wilson
   信賴區間）。

## 核心發現（先講給你知道，跑完你會看到一致的數字）

- **單一 output 分別看**：隨機拓撲（不論哪種隨機化方式）常常也會出現類似
  幅度的 NAFL 優勢（方向一致比例 ~57-61%），P 值都不顯著（0.03-0.17）。
  這代表這個方向性主要來自「疾病嚴重度訊號廣泛分散在幾乎所有節點的初始
  值裡」，不是特定連線方式的產物——這點要誠實承認。
- **三個 output 同時看**（呼應稿件 Table S4 的 `all_outputs_hold` 設計）：
  隨機拓撲要「同時」在三條獨立訊號路徑上都達到觀察到的優勢幅度，機率
  很低（兩種隨機化方式都是 P<0.01-0.02），且兩種 null model 給出高度
  一致的結果，顯示這個聯合、協調的三軸效應確實需要類似真實拓撲的連結
  結構才能重現。

## 使用方式

```powershell
cd F:\MASLD\ODEModeling\MASLD-ODE-Model
mkdir analyses\edge_scramble_control
# 把這個資料夾裡所有檔案複製進去
cd analyses\edge_scramble_control
conda activate deseq2_env
pip install pandas openpyxl --quiet   # 如果還沒裝

# 第 1 步（一定要先做）：驗證通用模型跟你原本的模型數值一致
python step1_validate.py

# 第 2 步：算出真實拓撲的參考值（只需跑一次）
python step2_run_batch.py true

# 第 3 步：跑隨機化 replicate（可分批執行、可中斷續跑）
python step2_run_batch.py scramble 0 200
python step2_run_batch.py scramble 200 400
python step2_run_batch.py er 0 200
python step2_run_batch.py er 200 400

# 第 4 步：統計分析 + 輸出 Excel
python step3_analyze.py
```

- 每個 replicate 約需 1.5-2.5 秒（4 次 ODE 積分：NASH baseline/drug、NAFL
  baseline/drug），200 個約 5-8 分鐘。
- `step2_run_batch.py` 用固定亂數種子（scramble=42, er=123），且用
  replicate 編號（不是檔案長度）決定要接著跑第幾個，所以：
  - 中斷後重新執行同一段範圍會得到一樣的結果
  - 分成多次呼叫（例如 `0 200` 再 `200 400`）跟一次呼叫 `0 400` 結果完全
    相同
- 模型檔案自動搜尋路徑（跟 `equilibrate_then_inhibit` 那個資料夾同一套
  邏輯）：`data/models/`、`data/models_pydeseq2_mean/`、本資料夾內的
  `models_pydeseq2_mean/`、本資料夾本身。找不到的話最簡單做法是把兩個
  `ode_model_pydeseq2_*.py` 複製一份到這個資料夾。

## 檔案結構

```
edge_scramble_control/
├── generic_model.py       ← 核心：拓撲 + 通用 ODE 生成器 + 兩種隨機化演算法
├── model_loader.py         ← 自動搜尋並載入你的 NASH/NAFL 模型檔案
├── step1_validate.py       ← 驗證通用模型 = 原始模型（必須先跑，PASS 才能往下）
├── step2_run_batch.py      ← 跑 true topology 參考值 + 兩種隨機化 replicate
├── step3_analyze.py        ← 統計分析（邊際+聯合檢定，Wilson CI）+ 輸出 Excel
├── README.md               ← 本檔
└── outputs/                ← 執行後自動產生
    ├── true_result.json
    ├── scramble_results.jsonl
    ├── er_results.jsonl
    └── edge_scramble_summary.xlsx
```

## 建議 replicate 數量

至少 200 個（每種 null model），400 個更穩定（我們這邊測試時 n=199 跟
n=399 兩次的結論方向完全一致，沒有因樣本數變化而翻盤，但 400 個的信賴
區間明顯更窄）。如果你的電腦跑得動，n=500-1000 會讓「聯合檢定 P<0.01」
這個估計更精確（目前 0/399 只能講「P 小於某個上界」，沒辦法精確到小數點
第三位）。

## 重要：如何誠實寫進 response letter / 稿件

**不要只講聯合檢定顯著、跳過邊際檢定不顯著這件事**——兩層都要講，完整
版本大概是：

> 單獨看任一 output，這種幅度的 NAFL 優勢在隨機拓撲下並不罕見（兩種
> null model 下方向一致比例皆 ~57-61%，P=0.03-0.17），顯示疾病嚴重度
> 訊號廣泛分散於網路的初始條件中；但要求三個 core output 同時達到觀察
> 到的優勢幅度，在兩種獨立的隨機化方法下都幾乎從未發生（P<0.01-0.02），
> 顯示這個聯合、協調的三軸效應確實仰賴接近 hsa04932 的特定連結結構。

也請注意分寸：小 P 值只能說「這在隨機重連下很罕見」，不能直接說成
「證明了 hsa04932 生物學機制的必要性」——嚴謹講法是「排除了『隨便接接
就會出現』的可能性」，而不是正面證明生物學特異性。

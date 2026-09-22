#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step0_setup_environment.py  (GSE48452 pipeline)
================================================================
一鍵環境準備腳本。在跑 step1-4 之前先執行這支，會自動：

  1. 檢查必要 Python 套件是否安裝（缺的會列出 pip 指令）
  2. 檢查 4 個 pipeline script 是否都在資料夾內
  3. 嘗試下載 GSE48452_series_matrix.txt（若尚未下載）
  4. 檢查 DEG_pydeseq2_3comparisons.xlsx（GSE126848 DEG）是否就位
  5. 印出就緒狀態總表，告訴你下一步該做什麼

使用方式（Anaconda Spyder, deseq2_env）：
  1. 把這支與 step1-4 放在同一資料夾
  2. F5 執行
  3. 依照畫面提示補齊缺少的東西
================================================================
"""
import os
import sys
import gzip
import shutil
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# 設定
# ----------------------------------------------------------------------
REQUIRED_PACKAGES = {
    'numpy':       'numpy',
    'pandas':      'pandas',
    'scipy':       'scipy',
    'matplotlib':  'matplotlib',
    'openpyxl':    'openpyxl',
    'GEOparse':    'GEOparse',
}

REQUIRED_SCRIPTS = [
    'step1_parse_series_matrix.py',
    'step2_annotate_probes.py',
    'step3_DEG.py',
    'step4_concordance.py',
]

SERIES_MATRIX_FILE = 'GSE48452_series_matrix.txt'
SERIES_MATRIX_GZ    = 'GSE48452_series_matrix.txt.gz'
SERIES_MATRIX_URL   = (
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/'
    'GSE48452/matrix/GSE48452_series_matrix.txt.gz'
)

GSE126848_DEG_FILE = 'DEG_pydeseq2_3comparisons.xlsx'
# 若 DEG 檔在上層資料夾，這裡列出可能位置供自動搜尋
GSE126848_SEARCH_PATHS = [
    os.path.join(SCRIPT_DIR, GSE126848_DEG_FILE),
    os.path.join(SCRIPT_DIR, '..', GSE126848_DEG_FILE),
    os.path.join(SCRIPT_DIR, '..', '..', GSE126848_DEG_FILE),
]


# ----------------------------------------------------------------------
# 檢查函式
# ----------------------------------------------------------------------
def check_packages():
    print('\n[1/4] 檢查 Python 套件...')
    missing = []
    for mod, pip_name in REQUIRED_PACKAGES.items():
        try:
            __import__(mod)
            print(f'  OK   {mod}')
        except ImportError:
            print(f'  缺   {mod}')
            missing.append(pip_name)
    if missing:
        print('\n  >>> 請在 Anaconda Prompt 執行（先 conda activate deseq2_env）：')
        print(f'      pip install {" ".join(missing)}')
        return False
    print('  所有套件就緒。')
    return True


def check_scripts():
    print('\n[2/4] 檢查 pipeline scripts...')
    missing = []
    for s in REQUIRED_SCRIPTS:
        path = os.path.join(SCRIPT_DIR, s)
        if os.path.exists(path):
            print(f'  OK   {s}')
        else:
            print(f'  缺   {s}')
            missing.append(s)
    if missing:
        print('\n  >>> 以下 script 缺少，請從 Claude 提供的檔案複製進此資料夾：')
        for s in missing:
            print(f'      - {s}')
        return False
    print('  所有 script 就緒。')
    return True


def check_series_matrix():
    print('\n[3/4] 檢查 GSE48452 series matrix...')
    txt_path = os.path.join(SCRIPT_DIR, SERIES_MATRIX_FILE)
    gz_path  = os.path.join(SCRIPT_DIR, SERIES_MATRIX_GZ)

    if os.path.exists(txt_path):
        size_mb = os.path.getsize(txt_path) / 1e6
        print(f'  OK   {SERIES_MATRIX_FILE} 已存在 ({size_mb:.1f} MB)')
        return True

    # 若有 .gz 但沒解壓縮
    if os.path.exists(gz_path):
        print(f'  發現壓縮檔 {SERIES_MATRIX_GZ}，正在解壓縮...')
        try:
            with gzip.open(gz_path, 'rb') as f_in, open(txt_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f'  OK   已解壓縮為 {SERIES_MATRIX_FILE}')
            return True
        except Exception as e:
            print(f'  錯誤：解壓縮失敗 - {e}')
            return False

    # 都沒有，嘗試下載
    print(f'  未找到 series matrix，嘗試從 NCBI 下載...')
    print(f'  URL: {SERIES_MATRIX_URL}')
    try:
        urllib.request.urlretrieve(SERIES_MATRIX_URL, gz_path)
        size_mb = os.path.getsize(gz_path) / 1e6
        print(f'  下載完成 ({size_mb:.1f} MB)，正在解壓縮...')
        with gzip.open(gz_path, 'rb') as f_in, open(txt_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f'  OK   已解壓縮為 {SERIES_MATRIX_FILE}')
        return True
    except Exception as e:
        print(f'  下載失敗 - {e}')
        print('\n  >>> 請手動下載：')
        print(f'      1. 開啟瀏覽器到：{SERIES_MATRIX_URL}')
        print(f'      2. 下載 {SERIES_MATRIX_GZ}')
        print(f'      3. 解壓縮後將 {SERIES_MATRIX_FILE} 放入此資料夾')
        print(f'      4. 或直接放 .gz 檔，重跑本腳本會自動解壓縮')
        return False


def check_gse126848_deg():
    print('\n[4/4] 檢查 GSE126848 DEG 檔案...')
    found_path = None
    for path in GSE126848_SEARCH_PATHS:
        if os.path.exists(path):
            found_path = os.path.abspath(path)
            break

    local_path = os.path.join(SCRIPT_DIR, GSE126848_DEG_FILE)
    if found_path and os.path.abspath(found_path) == os.path.abspath(local_path):
        size_mb = os.path.getsize(found_path) / 1e6
        print(f'  OK   {GSE126848_DEG_FILE} 已在本資料夾 ({size_mb:.1f} MB)')
        return True
    elif found_path:
        print(f'  發現 {GSE126848_DEG_FILE} 於：')
        print(f'       {found_path}')
        print(f'  正在複製到本資料夾以保持自包含...')
        try:
            shutil.copy(found_path, local_path)
            print(f'  OK   已複製到本資料夾')
            return True
        except Exception as e:
            print(f'  複製失敗 - {e}')
            print(f'  >>> 請手動複製 {GSE126848_DEG_FILE} 到本資料夾')
            return False
    else:
        print(f'  缺   {GSE126848_DEG_FILE}')
        print('\n  >>> 請從以下位置複製 DEG 檔到本資料夾：')
        print('      C:\\Users\\USER\\Dropbox\\1_Research\\研討會\\ICEIB\\'
              'Rewrite 20260415 v1 轉投\\DEG_pydeseq2_3comparisons.xlsx')
        return False


# ----------------------------------------------------------------------
# 主程式
# ----------------------------------------------------------------------
def main():
    print('=' * 64)
    print('GSE48452 Validation Pipeline — 環境準備檢查')
    print(f'工作資料夾：{SCRIPT_DIR}')
    print('=' * 64)

    results = {
        '套件': check_packages(),
        'Scripts': check_scripts(),
        'Series matrix': check_series_matrix(),
        'GSE126848 DEG': check_gse126848_deg(),
    }

    print('\n' + '=' * 64)
    print('就緒狀態總表')
    print('=' * 64)
    all_ok = True
    for item, ok in results.items():
        status = 'OK  ' if ok else '待處理'
        print(f'  [{status}]  {item}')
        if not ok:
            all_ok = False

    print('\n' + '=' * 64)
    if all_ok:
        print('全部就緒！可以開始執行 pipeline：')
        print('  step1_parse_series_matrix.py  (F5)')
        print('  step2_annotate_probes.py      (F5)')
        print('  step3_DEG.py                  (F5)')
        print('  step4_concordance.py          (F5)')
        print('\n建議先跑 step1，把 console 輸出貼給 Claude 確認 group 偵測正確。')
    else:
        print('尚有項目待處理（見上方提示），補齊後重跑本腳本確認。')
    print('=' * 64)


if __name__ == '__main__':
    main()

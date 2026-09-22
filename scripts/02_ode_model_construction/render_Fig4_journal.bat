@echo off
REM ============================================================
REM render_Fig4_journal.bat
REM 用 Graphviz neato 引擎，把 Fig4_ODE_flat_v2.dot 渲染成
REM npj SBA 期刊雙欄寬度規格 (183mm @ 300dpi)
REM
REM 前提: 已安裝 Graphviz 並把 dot.exe / neato.exe 加入 PATH
REM       (https://graphviz.org/download/，Windows 安裝後勾選
REM        "Add Graphviz to the system PATH")
REM
REM 用法: 跟 Fig4_ODE_flat_v2.dot 放同一資料夾，直接雙擊執行
REM ============================================================

set DOTFILE=Fig4_ODE_flat_v2.dot
set OUTNAME=Fig4_journal_double_col

REM 雙欄寬度 183mm = 7.205 inch；高度上限給 100 inch（讓 Graphviz
REM 自動依內容決定實際高度，寬度固定 7.205 inch 等比例縮放）
neato -Tpng -Gdpi=300 -Gsize="7.205,100" %DOTFILE% -o %OUTNAME%.png
neato -Tpdf          -Gsize="7.205,100" %DOTFILE% -o %OUTNAME%.pdf
neato -Tsvg          -Gsize="7.205,100" %DOTFILE% -o %OUTNAME%.svg

echo.
echo 完成！輸出檔案：
echo   %OUTNAME%.png  (點陣，300dpi)
echo   %OUTNAME%.pdf  (向量，期刊投稿建議用這個)
echo   %OUTNAME%.svg  (向量，方便後續編輯)
echo.
echo 如果要改成單欄寬度(89mm)，把上面三行的 7.205 改成 3.504
pause

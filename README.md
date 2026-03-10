# Hamming Code C 語言實作

以 **C 語言**實作 **Hamming Code / Extended Hamming Code** 的編碼與解碼流程，包含從 `Sim.txt` 讀取 parity-check matrix、建立 generator matrix、產生測試資料、加入錯誤樣式、進行 syndrome decoding，並輸出估測結果 `u.txt`。

## 專案簡介

本專案根據作業規格，完成一套 Hamming Code 的端到端模擬流程：

1. 讀取 `Sim.txt` 中的 `n`、`r`、parity-check matrix `H` 與 error patterns  
2. 將 `H` 轉換為可用於編碼的 systematic form  
3. 根據關係式建立 generator matrix `G = [I | P]`  
4. 使用 LFSR 產生 `m × k` 個 information bits  
5. 執行編碼 `x = uG`  
6. 加入錯誤樣式 `y = x ⊕ z`  
7. 以 syndrome decoding 進行錯誤更正 / 偵測  
8. 將估測結果輸出為 `u.txt`

## 技術重點

- **C 語言實作**
- **pointer 與二維陣列操作**
- **動態記憶體配置與釋放**
- **bitwise / XOR / parity 邏輯**
- **matrix-based encoding / decoding**
- **syndrome 計算與錯誤定位**
- **輸入格式解析與結果輸出**
- **邊界條件檢查與除錯**

## LFSR 資料產生規則

專案使用以下 recursion 產生測試用 information bits：

`u_{l+6} = u_{l+1} ⊕ u_l`

初始條件為：

`u0 = 1, u1 = u2 = u3 = u4 = u5 = 0`

可產生週期為 63 的序列。

## 專案結構

```text
.
├─ README.md
├─ src/
│  └─ hamming_code.c
├─ data/
│  └─ Sim.txt
├─ report/
   └─ Hamming_Code_Report

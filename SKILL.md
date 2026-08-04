---
name: geoskill-climate-trend-analysis
description: '对温度/降水时序执行 Mann-Kendall 趋势检验与 Sen 斜率（Sen slope）估计，输出趋势斜率栅格、显著性（p 值）栅格与时序统计 JSON。Mann-Kendall trend test and Sen slope estimator for temperature/precipitation time series, outputting slope raster, significance (p-value) raster, and time-series JSON.'
---

# 气候趋势分析 | Climate Trend Analysis

Performs pixel-wise trend analysis on multi-temporal temperature / precipitation
time series (multi-band GeoTIFF cubes or composited series) to identify regions
with significant warming / cooling or wetting / drying trends. Suitable for
climate change monitoring, regional warming-rate assessment, and long-term
precipitation trend screening.

Three core algorithms:

- **Mann-Kendall trend test**: a non-parametric rank test with statistic S,
  standardized Z, and two-tailed p-value, including a tie correction. It does
  not require normality and is robust to outliers, making it the standard
  method for climate trend analysis.
- **Sen slope**: the median of all pairwise slopes (x_j − x_i)/(j − i),
  providing a robust estimate of trend magnitude (units / time step) that
  resists outliers.
- **OLS linear regression slope**: provided as a comparative reference; the
  correlation between the two slope estimates is reported to assess robustness
  differences.

The built-in `--synthetic` mode generates simulated series with spatially
varying warming trends, mild seasonality, and red noise, allowing the full
workflow to be validated without network access or real data.

## Dependencies / 依赖

```bash
pip install numpy rasterio scipy
```

## Usage / 使用方法

### Basic usage (synthetic data, offline)

```bash
python geoskill-climate-trend-analysis.py --bbox 116.0 39.0 117.0 40.0 --output-dir ./output
```

### Example 1: temperature warming trend (synthetic)

```bash
python geoskill-climate-trend-analysis.py --bbox 116 39 117 40 --variable temperature --n-dates 24 --output-dir ./temp_trend
```

### Example 2: precipitation trend (synthetic)

```bash
python geoskill-climate-trend-analysis.py --bbox 121 31 122 32 --variable precipitation --n-dates 30 --output-dir ./precip_trend
```

### Example 3: real multi-temporal rasters

```bash
python geoskill-climate-trend-analysis.py --input annual_temp_stack.tif --alpha 0.01 --output-dir ./real_trend
```

### Example 4: stricter significance level

```bash
python geoskill-climate-trend-analysis.py --bbox 116 39 117 40 --variable temperature --alpha 0.01 --output-dir ./strict --quiet
```

### Example 5: bbox-only auto-synthesis + quiet mode

```bash
python geoskill-climate-trend-analysis.py --bbox 110 30 111 31 --n-dates 20 --output-dir ./auto --quiet
```

## Output / 输出

| File | Format | Description |
|---|---|---|
| `trend_slope.tif` | GeoTIFF (float32, 2 band) | band1=Sen slope, band2=OLS slope, EPSG:4326 |
| `significance.tif` | GeoTIFF (float32, 1 band) | MK two-tailed p-value raster (p<alpha means significant) |
| `timeseries.json` | JSON | Per-date spatial mean series + Sen/OLS summary + significant fraction |
| `output-manifest.json` | JSON | Run manifest (inputs/outputs/QA/software versions) |

## Data Source / 数据源 / Source

- **Input mode**: local multi-temporal GeoTIFF (each band = one time step).
- **Synthetic mode**: generated locally, no external data source.

## Privacy / 隐私声明 / Privacy

- Runs fully offline by default; no network requests are made.
- `--synthetic` mode reads no external data.
- All computation is performed locally; user data is never uploaded.

## License / License

MIT

---

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-climate-trend-analysis
description: '对温度/降水时序执行 Mann-Kendall 趋势检验与 Sen 斜率（Sen slope）估计，输出趋势斜率栅格、显著性（p 值）栅格与时序统计 JSON。Mann-Kendall trend test and Sen slope estimator for temperature/precipitation time series, outputting slope raster, significance (p-value) raster, and time-series JSON.'
---

# 气候趋势分析 | Climate Trend Analysis

对温度 / 降水的多期时间序列（多波段 GeoTIFF 立方体或合成序列）做逐像元
趋势分析，识别显著变暖 / 变冷、变湿 / 变干区域。适用于气候变化监测、
区域增温速率评估、降水长期趋势筛查等场景。

核心算法三件套：

- **Mann-Kendall 趋势检验**：非参数秩检验，统计量 S、标准化 Z、双尾 p 值，
  含结（tie）修正。不要求正态分布、对异常值稳健，是气候趋势分析的标准方法。
- **Sen 斜率（Sen slope）**：所有点对斜率 (x_j − x_i)/(j − i) 的中位数，
  稳健估计趋势幅度（单位 / 时间步），抵抗离群点。
- **OLS 线性回归斜率**：作为对比参考，输出二者相关系数评估稳健性差异。

内置 `--synthetic` 模式生成含空间变化增温趋势 + 温和季节性 + 红噪声的模拟
序列，无需网络与真实数据即可验证流程。

## 依赖

```bash
pip install numpy rasterio scipy
```

## 使用方法

### 基本用法（合成数据，离线）

```bash
python geoskill-climate-trend-analysis.py --bbox 116.0 39.0 117.0 40.0 --output-dir ./output
```

### 示例 1：温度增温趋势（合成）

```bash
python geoskill-climate-trend-analysis.py --bbox 116 39 117 40 --variable temperature --n-dates 24 --output-dir ./temp_trend
```

### 示例 2：降水趋势（合成）

```bash
python geoskill-climate-trend-analysis.py --bbox 121 31 122 32 --variable precipitation --n-dates 30 --output-dir ./precip_trend
```

### 示例 3：真实多期栅格

```bash
python geoskill-climate-trend-analysis.py --input annual_temp_stack.tif --alpha 0.01 --output-dir ./real_trend
```

### 示例 4：更严格的显著性水平

```bash
python geoskill-climate-trend-analysis.py --bbox 116 39 117 40 --variable temperature --alpha 0.01 --output-dir ./strict --quiet
```

### 示例 5：仅 bbox 自动合成 + 静默

```bash
python geoskill-climate-trend-analysis.py --bbox 110 30 111 31 --n-dates 20 --output-dir ./auto --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `trend_slope.tif` | GeoTIFF (float32, 2 band) | band1=Sen 斜率，band2=OLS 斜率，EPSG:4326 |
| `significance.tif` | GeoTIFF (float32, 1 band) | MK 双尾 p 值栅格（p<alpha 即显著） |
| `timeseries.json` | JSON | 逐期空间均值序列 + Sen/OLS 汇总 + 显著比例 |
| `output-manifest.json` | JSON | 运行清单（输入/输出/QA/软件版本） |

## 数据源 / Source

- **输入模式**：本地多期 GeoTIFF（每波段 = 一个时间步）。
- **合成模式**：本地生成，无外部数据源。

## 隐私声明 / Privacy

- 默认完全离线运行，不发起任何网络请求。
- `--synthetic` 模式不读取任何外部数据。
- 所有计算在本地完成，不上传用户数据。

## License

MIT

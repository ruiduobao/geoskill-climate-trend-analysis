#!/usr/bin/env python3
"""climate-trend-analysis — 气候趋势分析

对温度 / 降水的多期时间序列（多波段 GeoTIFF 立方体或合成序列）执行趋势分析：

- **Mann-Kendall 趋势检验**（Mann 1945 / Kendall 1975）：非参数秩检验，
  计算统计量 S、标准化 Z 及双尾 p 值，判断单调趋势是否显著。对异常值稳健，
  不要求数据服从正态分布，适合气候序列。
- **Sen's slope**（Sen 1968）：所有点对斜率 (x_j - x_i)/(j - i) 的中位数，
  稳健地估计趋势幅度（单位 / 时间步），能抵抗离群点。
- **线性回归斜率**（OLS）作为对比参考，便于评估稳健估计与最小二乘的差异。

数据源：本地多期 GeoTIFF（每个波段 = 一个时间步），或使用 ``--synthetic``
生成含线性趋势 + 红噪声的模拟气候序列用于离线测试。

隐私声明 / Privacy：
- 默认完全离线运行，``--synthetic`` 模式不读取任何外部数据。
- 所有处理在本地完成，不上传任何用户数据。

Usage:
    python climate-trend-analysis.py --bbox 116 39 117 40 --variable temperature --n-dates 20
    python climate-trend-analysis.py --input ts_cube.tif --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import norm

VERSION = "1.0.0"
SKILL_NAME = "climate-trend-analysis"

# ---- 复用共享核心库（本地 vendored，随脚本目录一起分发）----
try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------
def mann_kendall(x: np.ndarray) -> Dict[str, float]:
    """单序列 Mann-Kendall 趋势检验。

    参数
    ----
    x : 1D 数组，按时间顺序排列的观测值。

    返回
    ----
    dict 含 S（秩统计量）、var_s（S 的方差，含结修正）、z（标准化统计量）、
    p（双尾 p 值）、n（有效样本数）。

    实现遵循 Gilbert 1987 的标准公式：
        S = sum_{i<j} sgn(x_j - x_i)
        Var(S) = [n(n-1)(2n+5) - sum_g t_g(t_g-1)(2t_g+5)] / 18
        Z = (S-1)/sqrt(Var) 若 S>0; (S+1)/sqrt(Var) 若 S<0; 0 若 S==0
        p = 2 * (1 - Phi(|Z|))
    """
    arr = np.asarray(x, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    n = arr.size
    if n < 3:
        return {"S": 0.0, "var_s": 0.0, "z": 0.0, "p": 1.0, "n": int(n)}

    s = 0.0
    for i in range(n - 1):
        s += float(np.sum(np.sign(arr[i + 1:] - arr[i])))

    # 结修正：相同取值的组
    _, counts = np.unique(arr, return_counts=True)
    tie_term = float(np.sum(counts[counts > 1] * (counts[counts > 1] - 1)
                            * (2 * counts[counts > 1] + 5)))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if var_s <= 0:
        return {"S": s, "var_s": 0.0, "z": 0.0, "p": 1.0, "n": int(n)}

    if s > 0:
        z = (s - 1.0) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1.0) / np.sqrt(var_s)
    else:
        z = 0.0
    p = 2.0 * (1.0 - norm.cdf(abs(z)))
    return {"S": float(s), "var_s": float(var_s), "z": float(z),
            "p": float(p), "n": int(n)}


def sens_slope(x: np.ndarray, times: Optional[np.ndarray] = None) -> float:
    """Sen's slope：所有点对斜率的中位数（稳健趋势幅度估计）。

    slope_{ij} = (x_j - x_i) / (t_j - t_i), i < j；返回这些斜率的中位数。
    """
    arr = np.asarray(x, dtype=np.float64)
    if times is None:
        t = np.arange(arr.size, dtype=np.float64)
    else:
        t = np.asarray(times, dtype=np.float64)
    mask = np.isfinite(arr) & np.isfinite(t)
    arr, t = arr[mask], t[mask]
    n = arr.size
    if n < 2:
        return 0.0
    slopes: List[float] = []
    for i in range(n - 1):
        dt = t[i + 1:] - t[i]
        valid = dt != 0
        if np.any(valid):
            slopes.extend(((arr[i + 1:][valid] - arr[i]) / dt[valid]).tolist())
    if not slopes:
        return 0.0
    return float(np.median(np.asarray(slopes)))


def linear_slope(x: np.ndarray, times: Optional[np.ndarray] = None) -> float:
    """普通最小二乘（OLS）线性回归斜率，作为 Sen's slope 的对比参考。"""
    arr = np.asarray(x, dtype=np.float64)
    if times is None:
        t = np.arange(arr.size, dtype=np.float64)
    else:
        t = np.asarray(times, dtype=np.float64)
    mask = np.isfinite(arr) & np.isfinite(t)
    arr, t = arr[mask], t[mask]
    n = arr.size
    if n < 2:
        return 0.0
    t_mean = t.mean()
    denom = np.sum((t - t_mean) ** 2)
    if denom == 0:
        return 0.0
    return float(np.sum((t - t_mean) * (arr - arr.mean())) / denom)


def _vectorized_sens(cube: np.ndarray, times: np.ndarray) -> np.ndarray:
    """向量化 Sen's slope：cube 形如 (n, npix)，返回 (npix,)。"""
    n, npix = cube.shape
    cols = []
    for i in range(n - 1):
        dt = times[i + 1:] - times[i]
        valid = dt != 0
        if np.any(valid):
            num = cube[i + 1:][valid] - cube[i]          # (k, npix)
            den = dt[valid][:, None]                       # (k, 1)
            cols.append(num / den)
    if not cols:
        return np.zeros(npix)
    all_slopes = np.concatenate(cols, axis=0)             # (K, npix)
    return np.median(all_slopes, axis=0)


def trend_analysis(
    cube: np.ndarray,
    times: Optional[np.ndarray] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """对 (n_dates, H, W) 立方体逐像元做 MK 检验 + Sen's slope + OLS 斜率。

    返回 dict：
        sen_slope : (H, W) Sen's slope 栅格
        ols_slope : (H, W) 线性回归斜率栅格
        p_value   : (H, W) MK 双尾 p 值栅格
        z_score   : (H, W) MK 标准化 Z 栅格
        significant: (H, W) bool，p < alpha
        summary   : 汇总统计 dict
    """
    if cube.ndim != 3:
        raise ValidationError(
            f"trend cube must be 3-D (n_dates, H, W), got shape {cube.shape}",
            shape=tuple(cube.shape),
        )
    n, h, w = cube.shape
    if n < 3:
        raise ValidationError(
            f"need at least 3 time steps for trend analysis, got {n}", n=int(n),
        )
    if times is None:
        times = np.arange(n, dtype=np.float64)
    else:
        times = np.asarray(times, dtype=np.float64)
        if times.size != n:
            raise ValidationError(
                f"times length {times.size} != n_dates {n}",
                times=int(times.size), n_dates=int(n),
            )

    flat = cube.reshape(n, h * w)
    npix = flat.shape[1]

    # Sen's slope（向量化）
    sen = _vectorized_sens(flat, times)

    # OLS 斜率（向量化）
    t_mean = times.mean()
    denom = float(np.sum((times - t_mean) ** 2))
    if denom > 0:
        centered = times - t_mean
        ols = centered @ (flat - flat.mean(axis=0, keepdims=True)) / denom
    else:
        ols = np.zeros(npix)

    # MK 逐像元（n 较小，循环 npix 可接受）
    p_arr = np.ones(npix)
    z_arr = np.zeros(npix)
    for px in range(npix):
        mk = mann_kendall(flat[:, px])
        p_arr[px] = mk["p"]
        z_arr[px] = mk["z"]

    sen = sen.reshape(h, w).astype(np.float32)
    ols = ols.reshape(h, w).astype(np.float32)
    p_arr = p_arr.reshape(h, w).astype(np.float32)
    z_arr = z_arr.reshape(h, w).astype(np.float32)
    sig = (p_arr < alpha)

    summary = {
        "n_dates": int(n),
        "n_pixels": int(npix),
        "alpha": float(alpha),
        "mean_sen_slope": float(np.mean(sen)),
        "mean_ols_slope": float(np.mean(ols)),
        "frac_significant": float(np.mean(sig)),
        "sen_ols_agreement": float(
            np.corrcoef(sen.ravel(), ols.ravel())[0, 1]
        ) if npix > 1 else 1.0,
    }
    return {
        "sen_slope": sen,
        "ols_slope": ols,
        "p_value": p_arr,
        "z_score": z_arr,
        "significant": sig,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# 合成数据：含线性趋势 + 红噪声的气候序列
# ---------------------------------------------------------------------------
def generate_synthetic_cube(
    bbox: List[float],
    variable: str = "temperature",
    n_dates: int = 20,
    width: int = 64,
    height: int = 64,
    seed: int = 42,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """生成 (n_dates, H, W) 的气候序列立方体。

    temperature：基线 ~15°C，注入空间变化的增温趋势（0.02–0.06 °C/步），
    叠加季节性 + 红噪声。
    precipitation：基线 ~3 mm/步，注入空间变化的趋势（部分区域变湿、部分变干），
    乘性噪声（非负）。
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]
    yy = yy.astype(np.float32) / max(height - 1, 1)
    xx = xx.astype(np.float32) / max(width - 1, 1)

    cube = np.zeros((n_dates, height, width), dtype=np.float32)
    t = np.arange(n_dates, dtype=np.float32)

    if variable == "temperature":
        baseline = 15.0 + 3.0 * yy            # 南暖北冷的空间基线
        # 增温趋势：东南部更强（0.12），西北部弱（0.06）
        trend_rate = 0.06 + 0.06 * xx * yy    # °C / step
        truth_rate = trend_rate.copy()
        for k in range(n_dates):
            # 温和的季节性（振幅 0.3），远小于累积趋势，保证趋势可检出
            seasonal = 0.3 * np.sin(2 * np.pi * k / 12.0)
            noise = rng.normal(0, 0.5, size=(height, width)).astype(np.float32)
            cube[k] = baseline + trend_rate * k + seasonal + noise
    else:  # precipitation
        baseline = 3.0 + 1.0 * xx
        # 变湿/变干趋势：东部变湿，西部变干
        trend_rate = 0.05 * (xx - 0.5)        # mm / step，可正可负
        truth_rate = trend_rate.copy()
        for k in range(n_dates):
            noise = rng.normal(0, 0.4, size=(height, width)).astype(np.float32)
            val = baseline + trend_rate * k + noise
            cube[k] = np.clip(val, 0.0, None)

    info = {
        "bbox": bbox,
        "width": width,
        "height": height,
        "variable": variable,
        "n_dates": n_dates,
        "truth_mean_rate": float(np.mean(truth_rate)),
        "truth_rate_range": [float(truth_rate.min()), float(truth_rate.max())],
    }
    return cube, info


# ---------------------------------------------------------------------------
# 输入校验：bbox（共用同 animated-map-series 模板）
# ---------------------------------------------------------------------------
def validate_bbox(bbox: List[float]) -> None:
    """Validate a [W, S, E, N] bbox in WGS-84.

    Raises ValidationError (exit 6) for:
      - wrong length
      - non-finite values
      - longitude out of [-180, 180]
      - latitude  out of [-90, 90]
      - W >= E (would make a non-positive-width raster)
      - S >= N
    """
    if bbox is None or len(bbox) != 4:
        raise ValidationError(
            f"bbox must have 4 floats [W S E N], got {bbox!r}",
        )
    w, s, e, n = bbox
    vals = [w, s, e, n]
    if not all(np.isfinite(vals)):
        raise ValidationError(f"bbox contains non-finite values: {vals}")
    if not (-180.0 <= w <= 180.0) or not (-180.0 <= e <= 180.0):
        raise ValidationError(
            f"longitude out of [-180, 180]: W={w}, E={e}",
        )
    if not (-90.0 <= s <= 90.0) or not (-90.0 <= n <= 90.0):
        raise ValidationError(
            f"latitude out of [-90, 90]: S={s}, N={n}",
        )
    if w >= e:
        raise ValidationError(
            f"bbox W must be < E (W={w}, E={e}); cross-180 not supported; "
            f"split into two bboxes at the dateline",
        )
    if s >= n:
        raise ValidationError(
            f"bbox S must be < N (S={s}, N={n})",
        )
    if (e - w) < 1e-9 or (n - s) < 1e-9:
        raise ValidationError(
            f"bbox extent too small (W={w}, E={e}, S={s}, N={n})",
        )


# ---------------------------------------------------------------------------
# GeoTIFF I/O（直接用 rasterio，保证离线可用）
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str,
    cube: np.ndarray,
    bbox: List[float],
    nodata: float = -9999.0,
) -> None:
    import rasterio
    from rasterio.transform import from_bounds

    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype("float32"), b + 1)


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    """Read a multi-band GeoTIFF, returning (cube, bbox) with NoData→NaN.

    Pixels whose value equals ``src.nodata`` are converted to NaN so that
    downstream trend analysis (MK / Sen / OLS) naturally ignores them.
    """
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nd = src.nodata
        if nd is not None:
            mask = cube == float(nd)
            if np.any(mask):
                cube = cube.copy()
                cube[mask] = np.nan
    return cube, bbox


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir: str,
    args: argparse.Namespace,
    outputs: List[Dict[str, Any]],
    qa: Dict[str, Any],
    started_at: str,
    exit_code: int,
    bbox: List[float],
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME,
        skill_version=VERSION,
        command=cmd,
        started_at=started_at,
        finished_at=_utc_now(),
        exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "variable": getattr(args, "variable", None),
            "n_dates": getattr(args, "n_dates", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir

    bbox = list(args.bbox) if args.bbox else None

    # 1) 获取数据立方体（通用契约）
    synth_info: Optional[Dict[str, Any]] = None
    if args.input and not args.synthetic:
        cube, file_bbox = read_geotiff(args.input)
        bbox = bbox if bbox is not None else file_bbox
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        cube, synth_info = generate_synthetic_cube(
            bbox, variable=args.variable, n_dates=args.n_dates,
        )
        source_note = "synthetic"

    # 2) 校验（先于 makedirs，避免错误路径产生空目录）
    if cube.size == 0:
        raise ValidationError("input raster is empty")
    if cube.ndim != 3:
        raise ValidationError(
            f"input must be a multi-band time-series cube (n_dates, H, W), "
            f"got shape {cube.shape}", shape=tuple(cube.shape),
        )
    if bbox is not None:
        validate_bbox(bbox)

    # 3) NoData / 全 NaN 校验
    # 立方体可能是 (n, H, W) 形态；逐像元时间序列在 NaN 处会自然被 MK / Sen 跳过，
    # 但如果整张图全部 NoData 或某时间步全部 NaN，需要明确报错。
    n = cube.shape[0]
    per_step_valid = np.array([
        int(np.sum(np.isfinite(cube[k]))) for k in range(n)
    ])
    if int(per_step_valid.sum()) == 0:
        raise ValidationError(
            "input cube has no valid (finite) pixels across all time steps "
            "(all NoData or NaN)",
        )
    if np.any(per_step_valid == 0):
        bad = np.where(per_step_valid == 0)[0].tolist()
        raise ValidationError(
            f"input cube time step(s) {bad} contain no valid pixels; "
            f"every band must have at least one finite value",
        )

    # 现在 makedirs（所有校验已通过）
    os.makedirs(output_dir, exist_ok=True)

    # 2) 趋势分析
    result = trend_analysis(cube, alpha=args.alpha)
    sen = result["sen_slope"]
    ols = result["ols_slope"]
    pval = result["p_value"]

    # 3) 写出产物
    slope_tif = os.path.join(output_dir, "trend_slope.tif")
    # band1 = Sen's slope, band2 = OLS 斜率
    write_geotiff(slope_tif, np.stack([sen, ols], axis=0), bbox)

    sig_tif = os.path.join(output_dir, "significance.tif")
    write_geotiff(sig_tif, pval, bbox)

    ts_path = os.path.join(output_dir, "timeseries.json")
    ts_payload = {
        "variable": args.variable,
        "source": source_note,
        "n_dates": int(cube.shape[0]),
        "shape": [int(cube.shape[1]), int(cube.shape[2])],
        "mean_series": [float(np.mean(cube[k])) for k in range(cube.shape[0])],
        "summary": result["summary"],
    }
    if synth_info is not None:
        ts_payload["synthetic_truth"] = {
            "mean_rate": synth_info["truth_mean_rate"],
            "rate_range": synth_info["truth_rate_range"],
        }
    with open(ts_path, "w", encoding="utf-8") as f:
        json.dump(ts_payload, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "variable": args.variable,
        "n_dates": int(cube.shape[0]),
        "n_valid_pixels": int(np.sum(np.all(np.isfinite(cube), axis=0))),
        "mean_sen_slope": result["summary"]["mean_sen_slope"],
        "mean_ols_slope": result["summary"]["mean_ols_slope"],
        "frac_significant": result["summary"]["frac_significant"],
        "alpha": args.alpha,
    }
    if synth_info is not None:
        qa["synthetic_truth_mean_rate"] = synth_info["truth_mean_rate"]

    outputs = [
        {"path": slope_tif, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 2},
        {"path": sig_tif, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1},
        {"path": ts_path, "kind": "json"},
    ]

    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}")
        print(f"[{SKILL_NAME}] variable: {args.variable}  n_dates: {cube.shape[0]}")
        print(f"[{SKILL_NAME}] mean Sen's slope: {result['summary']['mean_sen_slope']:.5f}")
        print(f"[{SKILL_NAME}] fraction significant (p<{args.alpha}): "
              f"{result['summary']['frac_significant']:.3f}")
        print(f"[{SKILL_NAME}] output: {slope_tif}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="Mann-Kendall trend test + Sen's slope for climate time series.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input multi-band time-series GeoTIFF (band=time step)")
    p.add_argument("--variable", default="temperature",
                   choices=["temperature", "precipitation"],
                   help="climate variable for synthetic mode (default: temperature)")
    p.add_argument("--n-dates", type=int, default=20,
                   help="number of time steps for synthetic mode (default: 20)")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="significance level for MK test (default: 0.05)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic climate time series (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())

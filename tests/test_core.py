"""Core algorithm tests for climate-trend-analysis."""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod


class TestMannKendall:
    def test_monotonic_increasing_significant(self):
        """严格单调递增序列 → 显著正趋势。"""
        x = np.arange(30, dtype=float)
        mk = mod.mann_kendall(x)
        assert mk["S"] > 0
        assert mk["z"] > 0
        assert mk["p"] < 0.05

    def test_monotonic_decreasing_significant(self):
        x = np.arange(30, 0, -1, dtype=float)
        mk = mod.mann_kendall(x)
        assert mk["S"] < 0
        assert mk["z"] < 0
        assert mk["p"] < 0.05

    def test_constant_series_not_significant(self):
        """常数序列 → S=0, Z=0, p=1。"""
        x = np.full(20, 5.0)
        mk = mod.mann_kendall(x)
        assert mk["S"] == 0.0
        assert mk["z"] == 0.0
        assert mk["p"] == pytest.approx(1.0)

    def test_short_series_safe(self):
        mk = mod.mann_kendall(np.array([1.0, 2.0]))
        assert mk["p"] == 1.0
        assert mk["n"] == 2

    def test_p_value_range(self):
        rng = np.random.default_rng(0)
        x = rng.normal(0, 1, 50)
        mk = mod.mann_kendall(x)
        assert 0.0 <= mk["p"] <= 1.0


class TestSensSlope:
    def test_exact_linear_recovery(self):
        """Sen's slope 精确恢复纯线性斜率 2.0。"""
        t = np.arange(25, dtype=float)
        x = 2.0 * t + 7.0
        assert mod.sens_slope(x, t) == pytest.approx(2.0, abs=1e-9)

    def test_negative_slope(self):
        t = np.arange(20, dtype=float)
        x = -0.5 * t + 100.0
        assert mod.sens_slope(x, t) == pytest.approx(-0.5, abs=1e-9)

    def test_robust_to_outlier(self):
        """加入大离群点，Sen's slope 仍接近真实斜率（稳健性）。"""
        t = np.arange(30, dtype=float)
        x = 1.0 * t + 5.0
        x[10] += 1000.0  # 极端离群点
        slope = mod.sens_slope(x, t)
        assert abs(slope - 1.0) < 0.2

    def test_default_times(self):
        x = 3.0 * np.arange(15, dtype=float)
        assert mod.sens_slope(x) == pytest.approx(3.0, abs=1e-9)


class TestLinearSlope:
    def test_exact_linear(self):
        t = np.arange(20, dtype=float)
        x = 2.0 * t + 7.0
        assert mod.linear_slope(x, t) == pytest.approx(2.0, abs=1e-9)

    def test_flat(self):
        assert mod.linear_slope(np.full(10, 3.0)) == pytest.approx(0.0, abs=1e-9)


class TestTrendAnalysis:
    def test_shape_preserved(self):
        rng = np.random.default_rng(1)
        cube = rng.normal(0, 1, (10, 16, 16)).astype(np.float32)
        res = mod.trend_analysis(cube)
        assert res["sen_slope"].shape == (16, 16)
        assert res["p_value"].shape == (16, 16)
        assert res["z_score"].shape == (16, 16)

    def test_detects_injected_trend(self):
        """注入强增温趋势 → 大部分像元显著且 Sen's 斜率接近真值。"""
        n, h, w = 20, 8, 8
        rng = np.random.default_rng(2)
        t = np.arange(n, dtype=np.float32)
        cube = np.zeros((n, h, w), dtype=np.float32)
        for k in range(n):
            cube[k] = 15.0 + 0.3 * k + rng.normal(0, 0.2, (h, w)).astype(np.float32)
        res = mod.trend_analysis(cube)
        # 几乎全部像元显著
        assert res["summary"]["frac_significant"] > 0.9
        # Sen's 斜率接近 0.3
        assert res["summary"]["mean_sen_slope"] == pytest.approx(0.3, abs=0.05)

    def test_too_few_dates_raises(self):
        cube = np.zeros((2, 4, 4), dtype=np.float32)
        with pytest.raises(mod.ValidationError):
            mod.trend_analysis(cube)

    def test_wrong_ndim_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.trend_analysis(np.zeros((4, 4), dtype=np.float32))


class TestSynthetic:
    def test_temperature_cube_shape(self):
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], variable="temperature", n_dates=15)
        assert cube.shape == (15, 64, 64)
        assert info["variable"] == "temperature"

    def test_temperature_has_positive_mean_trend(self):
        """合成的增温序列应被检出为正趋势。"""
        cube, info = mod.generate_synthetic_cube(
            [116, 39, 117, 40], variable="temperature", n_dates=24)
        res = mod.trend_analysis(cube)
        assert res["summary"]["mean_sen_slope"] > 0
        assert info["truth_mean_rate"] > 0

    def test_precipitation_nonnegative(self):
        cube, _ = mod.generate_synthetic_cube(
            [116, 39, 117, 40], variable="precipitation", n_dates=12)
        assert cube.min() >= 0.0


class TestGeoTiffIO:
    def test_roundtrip(self, tmp_path):
        cube = np.random.uniform(-1, 1, (2, 16, 16)).astype(np.float32)
        bbox = [116.0, 39.0, 117.0, 40.0]
        path = str(tmp_path / "t.tif")
        mod.write_geotiff(path, cube, bbox)
        assert os.path.exists(path)
        back, rb = mod.read_geotiff(path)
        assert back.shape == cube.shape
        np.testing.assert_allclose(rb, bbox, atol=1e-6)
        np.testing.assert_allclose(back, cube, atol=1e-5)

    def test_read_missing_raises(self):
        with pytest.raises(mod.UsageError):
            mod.read_geotiff("/nonexistent/x.tif")

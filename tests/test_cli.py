"""CLI argument parsing tests for climate-trend-analysis."""
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod, SCRIPT


def run_cli(args, timeout=120):
    return subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True, timeout=timeout,
    )


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--synthetic" in r.stdout
        assert "--variable" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        """无参数（无 --input 也无 --bbox）→ UsageError → exit 2。"""
        r = run_cli([])
        assert r.returncode == 2

    def test_bad_variable_exit_2(self):
        """非法 variable choice → argparse error → exit 2。"""
        r = run_cli(["--bbox", "116", "39", "117", "40", "--variable", "humidity"])
        assert r.returncode == 2

    def test_input_file_not_found_exit_2(self):
        """不存在的输入文件 → UsageError → exit 2。"""
        r = run_cli(["--input", "nonexistent_cube.tif"])
        assert r.returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        """--synthetic 无 --bbox → UsageError → exit 2。"""
        r = run_cli(["--synthetic"])
        assert r.returncode == 2


class TestCLISynthetic:
    def test_synthetic_runs(self):
        r = run_cli([
            "--bbox", "116", "39", "117", "40", "--synthetic",
            "--variable", "temperature", "--n-dates", "12", "--quiet",
            "--output-dir", "./_test_cli_synth",
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(
            os.path.dirname(SCRIPT), "_test_cli_synth", "output-manifest.json"))

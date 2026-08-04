"""Pytest fixtures for climate-trend-analysis tests."""
import importlib.util
import os
import sys
import pytest

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(SKILL_DIR, "geoskill-climate-trend-analysis.py")


def load_module():
    spec = importlib.util.spec_from_file_location("climate_trend_analysis", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["climate_trend_analysis"] = module
    spec.loader.exec_module(module)
    return module


mod = load_module()


@pytest.fixture
def skill_dir():
    return SKILL_DIR


@pytest.fixture
def script_path():
    return SCRIPT

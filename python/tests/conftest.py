"""Shared helpers for the tests."""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@lru_cache(maxsize=None)
def load_fixture(name: str) -> dict:
    """Return python/tests/fixtures/<name>.json, written by tools/make_fixtures.R."""
    with open(FIXTURE_DIR / f"{name}.json", encoding="utf-8") as f:
        return json.load(f)


def arr(x) -> np.ndarray:
    """Convert a JSON array to a float array; null becomes nan."""
    return np.asarray(x, dtype=float)


@pytest.fixture
def fixture():
    """Load a fixture by name, for example fixture("algebra")."""
    return load_fixture

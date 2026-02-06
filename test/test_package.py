from __future__ import annotations

import importlib.metadata

import openalea.archicrop as m


def test_version():
    assert importlib.metadata.version("openalea.archicrop") == m.__version__

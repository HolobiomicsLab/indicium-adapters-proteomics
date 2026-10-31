"""The package version is stated in three places, and they must agree.

``pyproject.toml`` is the source of truth. ``__version__`` is what callers read
at run time, and ``CITATION.cff`` is what a citation of the software records;
either one drifting means a user cites or reports a release that is not the one
they installed.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import indicium_adapters_proteomics

_ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version() -> str:
    with (_ROOT / "pyproject.toml").open("rb") as fh:
        return tomllib.load(fh)["project"]["version"]


def _citation_version() -> str:
    text = (_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*[\"']?([^\"'\s]+)[\"']?\s*$", text, re.MULTILINE)
    assert match, "CITATION.cff has no top-level version field"
    return match.group(1)


def test_dunder_version_matches_pyproject():
    assert indicium_adapters_proteomics.__version__ == _pyproject_version()


def test_citation_version_matches_pyproject():
    assert _citation_version() == _pyproject_version()

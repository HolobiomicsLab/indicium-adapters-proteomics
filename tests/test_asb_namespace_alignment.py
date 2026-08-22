"""Guard: the SHACL shapes must bind ``asb:`` to the base indicium targets.

The shape file hardcodes the ``asb:`` base as a literal ``@prefix``. If a future
base move updates indicium but misses the shapes file, the shape stops matching
the graph, pyshacl reports ``conforms=True`` against zero constraints, and every
invalid claim passes — a silent fail-open. ``test_shapes`` pins the synthetic
*data* graph to ``ASB_BASE`` but not the shapes, so this is the tripwire. The
metabolomics sibling carries the same guard.
"""

from __future__ import annotations

import re
from pathlib import Path

from indicium import ASB_BASE

_PKG = Path(__file__).resolve().parents[1] / "src" / "indicium_adapters_proteomics"
SHAPES_DIR = _PKG / "shapes"
SHAPE_FILE = "proteomics.shapes.ttl"
PREFIX_RE = re.compile(r"(?:@prefix|PREFIX)\s+asb:\s+<([^>]+)>")


def test_shape_binds_asb_to_indicium_base() -> None:
    text = (SHAPES_DIR / SHAPE_FILE).read_text()
    bases = PREFIX_RE.findall(text)
    assert bases, f"no asb: prefix declared in {SHAPE_FILE}"
    for base in bases:
        assert base == ASB_BASE, f"{SHAPE_FILE}: asb: bound to {base}, not {ASB_BASE}"


def test_no_shape_references_a_dead_base() -> None:
    """The base moved off two hosts nobody registered — neither may reappear."""
    # split so this guard scans itself without matching
    dead = ("asb.holobiomics" + ".org", "holobiomicslab" + ".eu/ns/asb")
    offenders = [
        f.name for f in SHAPES_DIR.glob("*.ttl") if any(host in f.read_text() for host in dead)
    ]
    assert offenders == []

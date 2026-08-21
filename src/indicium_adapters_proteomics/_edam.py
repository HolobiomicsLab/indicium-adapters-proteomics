"""Vendored EDAM label index — the concept spine for the proteomics adapter.

Generated offline from ``edamontology.org/EDAM.csv`` by ``scripts/build_index.py``
and shipped as a package resource, so enrichment is a pure function with no
network call at runtime. ``data/edam_index.json`` records the ontology roots the
index was cut from, alongside the label map itself.
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files
from typing import Final

EDAM_PREFIX: Final[str] = "http://edamontology.org/"


@lru_cache(maxsize=1)
def _payload() -> dict:
    raw = files("indicium_adapters_proteomics.data").joinpath("edam_index.json").read_text()
    return json.loads(raw)


def _normalise(text: str) -> str:
    """Case-fold and collapse whitespace so lookup is punctuation-stable."""
    return " ".join(str(text).lower().split())


def lookup(text: str) -> str | None:
    """Return the EDAM CURIE for ``text``, or None when it names no EDAM concept.

    Matching is exact after normalisation. Substring matching over an ontology
    manufactures confident-looking wrong CURIEs ("Peptide" inside "Peptide
    identification"), and a wrong identifier is worse than no identifier.
    """
    return _payload()["index"].get(_normalise(text))


def roots() -> tuple[str, ...]:
    """The EDAM classes this index was cut from — provenance for the vocabulary."""
    return tuple(_payload()["roots"])


def size() -> int:
    """Number of label/synonym entries in the index."""
    return len(_payload()["index"])

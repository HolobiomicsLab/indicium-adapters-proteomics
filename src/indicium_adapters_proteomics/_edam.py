"""Vendored EDAM label index — the concept spine for the proteomics adapter.

Generated offline from ``edamontology.org/EDAM.csv`` and shipped as a package
resource, so enrichment is a pure function with no network call at runtime.
``data/edam_index.json`` records the ontology roots the index was cut from
alongside the label map itself; unlike ``data/build_psimod_index.py`` no builder
for it ships here, so the JSON's own ``roots`` and ``source`` fields are the
provenance of record.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from importlib.resources import files
from typing import Final

EDAM_PREFIX: Final[str] = "http://edamontology.org/"

# A prefixed identifier, e.g. "GO:0006464" or "MOD:00000".
_CURIE: Final[re.Pattern[str]] = re.compile(r"^[A-Za-z][A-Za-z0-9_.]*:\S+$")


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

    Text that is already a CURIE is refused. EDAM lists the identifiers it
    cross-references among its own synonyms — ``GO:0006464`` and ``MOD:00000``
    both sit under Protein modifications — so an exact match would answer a slot
    that already names a term in another ontology with an EDAM *topic*, turning
    an identity into a subject heading.
    """
    if _CURIE.match(str(text).strip()):
        return None
    return _payload()["index"].get(_normalise(text))


def roots() -> tuple[str, ...]:
    """The EDAM classes this index was cut from — provenance for the vocabulary."""
    return tuple(_payload()["roots"])


def size() -> int:
    """Number of label/synonym entries in the index."""
    return len(_payload()["index"])

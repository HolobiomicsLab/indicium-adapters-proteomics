"""Maintainer-only builder for `_psimod_index.json` from the PSI-MOD OBO release.

Run on a PSI-MOD release, not at import time and not in CI:

    python -m indicium_adapters_proteomics.data.build_psimod_index PSI-MOD.obo

Source of truth is the HUPO-PSI controlled vocabulary
(https://github.com/HUPO-PSI/psi-mod-CV, file `PSI-MOD.obo`). Nothing in this
file invents a modification name: every index entry is a term's own `name:` or
one of its curated `synonym:` lines.

Two curation rules, both there to stop the index manufacturing wrong CURIEs:

* **RESID systematic names are excluded.** They are IUPAC chemistry strings
  ("(2R,3S)-2-amino-3-methylpentanoic acid") that no proteomics text or search
  engine emits, so they add collision surface without adding recall.
* **Ambiguous aliases are dropped, not arbitrated.** A synonym claimed by two
  MOD ids is removed from the index entirely; a term's own `name:` always wins
  over any synonym. A wrong identifier is worse than no identifier.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Final

# Curated synonym classes worth indexing — the labels that proteomics prose,
# search engines and the Unimod/UniProt exchange formats actually emit.
INDEXED_SYNONYM_TYPES: Final[frozenset[str]] = frozenset(
    {
        "PSI-MOD-label",
        "PSI-MOD-alternate",
        "PSI-MS-label",
        "Unimod-interim",
        "Unimod-description",
        "Unimod-alternate",
        "OMSSA-label",
        "DeltaMass-label",
        "RESID-name",
        "RESID-alternate",
    }
)

SOURCE_URL: Final[str] = (
    "https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo"
)


def _split_stanzas(obo_text: str) -> list[str]:
    """Split an OBO document into its `[Term]` stanzas, header discarded."""
    return obo_text.split("\n[Term]\n")[1:]


def _tag_values(stanza: str, tag: str) -> list[str]:
    """Return every value carried by `tag:` lines in one OBO stanza."""
    prefix = f"{tag}: "
    return [ln[len(prefix) :].strip() for ln in stanza.splitlines() if ln.startswith(prefix)]


def _synonym_label(line: str) -> tuple[str, str] | None:
    """Split an OBO `synonym:` value into (label, synonym type), or None.

    OBO shape is `"<label>" <scope> <type> []`; entries with no declared
    synonym type carry no curation provenance and are skipped.
    """
    if not line.startswith('"'):
        return None
    closing = line.find('"', 1)
    if closing < 0:
        return None
    trailing = line[closing + 1 :].split()
    if len(trailing) < 2:
        return None
    return line[1:closing], trailing[1]


def _normalise(text: str) -> str:
    """Case-fold and collapse whitespace, matching `_enrichment` lookup."""
    return " ".join(text.lower().split())


def _term_labels(stanza: str) -> tuple[str, list[str]]:
    """Return one term's own name plus its indexable synonym labels."""
    names = _tag_values(stanza, "name")
    synonyms = []
    for line in _tag_values(stanza, "synonym"):
        parsed = _synonym_label(line)
        if parsed and parsed[1] in INDEXED_SYNONYM_TYPES:
            synonyms.append(parsed[0])
    return (names[0] if names else ""), synonyms


class _IndexLayers:
    """The three accumulators `build_index` folds stanzas into."""

    def __init__(self) -> None:
        self.primary: dict[str, str] = {}
        self.secondary: dict[str, str] = {}
        self.ambiguous: set[str] = set()

    def add_synonym(self, key: str, mod_id: str) -> None:
        """Record a synonym, retiring the key entirely on a rival claim."""
        if key in self.ambiguous:
            return
        previous = self.secondary.get(key)
        if previous is not None and previous != mod_id:
            del self.secondary[key]
            self.ambiguous.add(key)
            return
        self.secondary[key] = mod_id

    def merged(self) -> dict[str, str]:
        """Primary names layered over surviving synonyms."""
        return {**self.secondary, **self.primary}


def _index_stanza(stanza: str, layers: _IndexLayers) -> None:
    """Fold one non-obsolete `[Term]` stanza into the primary/synonym layers."""
    if "\nis_obsolete: true" in f"\n{stanza}":
        return
    ids = _tag_values(stanza, "id")
    if not ids:
        return
    name, synonyms = _term_labels(stanza)
    if name:
        layers.primary.setdefault(_normalise(name), ids[0])
    for label in synonyms:
        layers.add_synonym(_normalise(label), ids[0])


def build_index(obo_text: str) -> dict[str, str]:
    """Build the normalised label to MOD CURIE index for a PSI-MOD OBO document.

    Primary `name:` labels win over synonyms; synonyms claimed by more than one
    term are dropped.
    """
    layers = _IndexLayers()
    for stanza in _split_stanzas(obo_text):
        _index_stanza(stanza, layers)
    return layers.merged()


def _data_version(obo_text: str) -> str:
    """Read the `data-version:` header of a PSI-MOD OBO document."""
    for line in obo_text.splitlines():
        if line.startswith("data-version: "):
            return line.split(": ", 1)[1].strip()
    return "unknown"


def main(obo_path: str, out_path: str) -> None:
    """Write `_psimod_index.json` from a local PSI-MOD OBO file."""
    obo_text = Path(obo_path).read_text(encoding="utf-8")
    payload = {
        "source": SOURCE_URL,
        "data_version": _data_version(obo_text),
        "synonym_types": sorted(INDEXED_SYNONYM_TYPES),
        "index": dict(sorted(build_index(obo_text).items())),
    }
    Path(out_path).write_text(json.dumps(payload, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

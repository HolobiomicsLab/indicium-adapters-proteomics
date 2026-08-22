"""The two identity lookups the proteomics adapter performs: UniProt and PSI-MOD.

Both are pure functions over vendored or published patterns — no network, no I/O
beyond a once-per-process read of the bundled PSI-MOD index.

**Protein accessions → UniProt.** The accession pattern is UniProt's own,
published in the UniProtKB manual (https://www.uniprot.org/help/accession_numbers)
and reproduced verbatim below. Recognising the *pattern* rather than a list of
accessions is what keeps this general: it accepts any UniProtKB accession ever
issued and invents none. Isoform suffixes (``P12345-2``) and the UniProtKB FASTA
header form (``sp|P12345|ALBU_HUMAN``) are unwrapped first, both being documented
UniProt conventions.

**PTMs → PSI-MOD.** Resolution is exact-after-normalisation against
``data/_psimod_index.json``, generated from the HUPO-PSI ``PSI-MOD.obo`` release
by ``data/build_psimod_index.py``. Labels claimed by more than one MOD term are
absent from that index by construction — "Phospho" is the PSI-MS label of ten
residue-specific terms in release 1.038.0, so it resolves to none of them rather
than to an arbitrary one.
``data/_ptm_aliases.json`` adds the process nouns papers actually write
("phosphorylation") on top, each pointing at a PSI-MOD *name* that is resolved
through the same index, so an alias whose target is renamed upstream disappears
instead of going stale.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from importlib.resources import files
from typing import Final

# UniProtKB accession grammar, from the UniProtKB manual. Wrapped in an optional
# isoform suffix, which UniProt appends as "-<n>".
UNIPROT_ACCESSION: Final[re.Pattern[str]] = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$"
)

# UniProtKB FASTA header: db|accession|entryName, db being sp (Swiss-Prot) or
# tr (TrEMBL).
_FASTA_HEADER: Final[re.Pattern[str]] = re.compile(r"^(?:sp|tr)\|([^|]+)\|", re.IGNORECASE)

# CURIE prefixes that already carry a UniProt accession, case-folded.
_UNIPROT_PREFIXES: Final[frozenset[str]] = frozenset({"uniprot", "uniprotkb", "up"})

_MOD_CURIE: Final[re.Pattern[str]] = re.compile(r"^MOD:\d{5}$", re.IGNORECASE)


def _normalise(text: str) -> str:
    """Case-fold and collapse whitespace, matching how the index was keyed."""
    return " ".join(str(text).lower().split())


@lru_cache(maxsize=1)
def _psimod_payload() -> dict:
    """Load the bundled PSI-MOD index once per process."""
    raw = files("indicium_adapters_proteomics.data").joinpath("_psimod_index.json").read_text()
    return json.loads(raw)


@lru_cache(maxsize=1)
def _ptm_index() -> dict[str, str]:
    """Return the PSI-MOD label index with the process-noun aliases folded in.

    An alias whose target name is not in the current PSI-MOD release is dropped
    silently — the alias file names targets, never MOD ids, so a renamed term
    costs recall rather than correctness.
    """
    index = dict(_psimod_payload()["index"])
    raw = files("indicium_adapters_proteomics.data").joinpath("_ptm_aliases.json").read_text()
    for alias, target in json.loads(raw).items():
        if alias.startswith("_"):
            continue
        resolved = index.get(_normalise(target))
        if resolved:
            index.setdefault(_normalise(alias), resolved)
    return index


def _strip_accession_wrapper(text: str) -> str:
    """Unwrap a FASTA header or CURIE prefix down to a bare accession candidate."""
    candidate = text.strip()
    header = _FASTA_HEADER.match(candidate)
    if header:
        return header.group(1).strip()
    prefix, _, rest = candidate.partition(":")
    if rest and prefix.lower() in _UNIPROT_PREFIXES:
        return rest.strip()
    return candidate


def uniprot_curie(text: str) -> str | None:
    """Return the ``uniprot:`` CURIE for a protein accession, or None.

    Accepts a bare accession (``P12345``), an isoform (``P12345-2``), an
    already-prefixed CURIE (``UniProtKB:P12345``) and the UniProtKB FASTA header
    form (``sp|P12345|ALBU_HUMAN``). Anything that is not a syntactically valid
    UniProtKB accession returns None rather than a guess.
    """
    if not text or not text.strip():
        return None
    candidate = _strip_accession_wrapper(text).upper()
    if not UNIPROT_ACCESSION.match(candidate):
        return None
    return f"uniprot:{candidate}"


def psimod_curie(text: str) -> str | None:
    """Return the ``MOD:`` CURIE for a modification name, or None.

    Accepts a PSI-MOD label, one of its curated synonyms, a process noun from
    the alias layer, or an already-formed ``MOD:#####`` CURIE. Matching is exact
    after case-folding and whitespace collapse: substring matching over a
    modification vocabulary manufactures confident-looking wrong CURIEs, and a
    wrong identifier is worse than no identifier.
    """
    if not text or not text.strip():
        return None
    candidate = text.strip()
    if _MOD_CURIE.match(candidate):
        return candidate.upper()
    return _ptm_index().get(_normalise(candidate))


def psimod_version() -> str:
    """The PSI-MOD release the bundled index was cut from."""
    return _psimod_payload()["data_version"]


def psimod_size() -> int:
    """Number of label entries in the bundled PSI-MOD index, aliases excluded."""
    return len(_psimod_payload()["index"])

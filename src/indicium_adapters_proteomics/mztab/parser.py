"""mzTab (proteomics) to asb:Claim projection — MTD, PRT, PEP and PSM sections.

Implementation note, inherited from the metabolomics sibling: we parse the
tab-separated sections directly rather than going through a general mzTab
loader. The sections are line-oriented and the semantics we need are per-row, so
a loader buys nothing and costs a dependency.

Claim shapes, one per row unless a row carries modification or abundance
columns, in which case those add claims of their own:

===============================  ================================================
mzTab source                     claim
===============================  ================================================
``MTD false_discovery_rate``     (file, ``fdr_controlled_at``, threshold)
``PRT`` row                      (file, ``identifies_protein``, accession)
``PRT`` search engine            (accession, ``identified_by``, engine)
``PRT`` database column          (accession, ``searched_against``, database)
``PEP`` row                      (accession, ``identifies_peptide``, sequence)
``PSM`` row                      (sequence, ``matched_to_spectrum``, spectra_ref)
``modifications`` column         (row subject, ``carries_modification``, CV id)
``*_abundance_assay[n]``         (row subject, ``has_abundance``, assay label)
===============================  ================================================

Protein accessions are resolved to ``uniprot:`` CURIEs when they are
syntactically valid UniProtKB accessions, and modification CV identifiers to
``MOD:`` when the file already names a PSI-MOD term. UNIMOD and CHEMMOD
identifiers are carried verbatim as the claim object and left unresolved: no
offline one-to-one mapping onto PSI-MOD exists, and inventing one would put a
wrong CURIE where none belongs.

**Commas are not separators on their own.** mzTab writes a CV parameter as
``[label, accession, name, value]`` (§5.10) and embeds parameters inside the
modifications column (§5.8) — as a per-position probability, and as a neutral
loss. Splitting either on a bare comma cuts one object into several, and the
fragments then read as further modifications; ``_split_top_level`` is what stops
a ``fragment neutral loss`` parameter being claimed as a PTM.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Final

from indicium_adapters_proteomics._claims import build_claim
from indicium_adapters_proteomics._enrichment import psimod_curie, uniprot_curie

# mzTab §5.9: a missing value MUST be reported as "null". The empty string is
# admitted alongside it only for the trailing columns a sparse row omits.
_NULL_VALUES: Final[frozenset[str]] = frozenset({"", "null"})

# Section header tag -> row tag, for the three identification tables.
_SECTIONS: Final[tuple[tuple[str, str], ...]] = (
    ("PRH", "PRT"),
    ("PEH", "PEP"),
    ("PSH", "PSM"),
)

# Quote characters mzTab uses to protect a param field containing a comma.
_QUOTES: Final[str] = "\"“”"
_BRACKET_DEPTH: Final[dict[str, int]] = {"[": 1, "]": -1}
_BRACKETED: Final[re.Pattern[str]] = re.compile(r"\[[^\]]*\]")
_PARAMETER_FIELDS: Final[int] = 4

# mzTab §5.8 closes the modification identifier namespaces to these three.
# Anything else in the column is a neutral-loss or scoring cvParam, which
# describes a fragment or a probability rather than naming a modification.
_MODIFICATION_PREFIXES: Final[str] = "MOD|UNIMOD|CHEMMOD"
# One site: a residue number or "null" when the site is unknown, optionally
# followed by its own cvParam.
_SITE: Final[str] = r"(?:\d+|null)(?:\[[^\]]*\])?"

# "3-MOD:00412", "3|4-UNIMOD:35" and the per-site parameter form
# "10[MS,MS:1001876, modification probability, 0.8]-MOD:00412".
_MODIFICATION: Final[re.Pattern[str]] = re.compile(
    rf"^(?P<positions>{_SITE}(?:\|{_SITE})*)"
    rf"-(?P<accession>(?:{_MODIFICATION_PREFIXES}):[^|,\s]+)",
    re.IGNORECASE,
)
# Positionless form, e.g. a CHEMMOD delta reported without a site.
_BARE_MODIFICATION: Final[re.Pattern[str]] = re.compile(
    rf"^(?P<accession>(?:{_MODIFICATION_PREFIXES}):[^|,\s]+)$", re.IGNORECASE
)

# Abundance value columns. Dispersion columns (_stdev_, _std_error_) describe a
# measurement rather than being one, so they are deliberately not matched.
_ABUNDANCE_COLUMN: Final[re.Pattern[str]] = re.compile(
    r"^(?:protein|peptide)_abundance_(?P<label>(?:assay|study_variable)\[\d+\])$"
)

# PRT columns naming how an identification was made, and the qualifier each maps to.
_PROVENANCE_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("search_engine", "identified_by"),
    ("database", "searched_against"),
)


def _split_top_level(text: str) -> list[str]:
    """Split `text` on the commas that sit outside every bracket and quote.

    mzTab embeds params in the modifications column (§5.8) and quotes any param
    field that itself contains a comma (§5.10), so a plain split would cut one
    object or one field into several.
    """
    parts: list[str] = []
    depth, start, in_quotes = 0, 0, False
    for index, char in enumerate(text):
        if char in _QUOTES:
            in_quotes = not in_quotes
        elif not in_quotes:
            depth = max(depth + _BRACKET_DEPTH.get(char, 0), 0)
        if char == "," and depth == 0 and not in_quotes:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return parts


def _parameter_fields(text: str) -> list[str] | None:
    """Return the four fields of an mzTab param, or None when `text` is not one.

    A param is written ``[CV label, accession, name, value]`` (§5.10) with any
    unavailable field left empty. Returns None for a plain column value, so a
    column that may hold either form passes non-params through untouched.
    """
    stripped = text.strip()
    if not (stripped.startswith("[") and stripped.endswith("]")):
        return None
    fields = [f.strip().strip(_QUOTES) for f in _split_top_level(stripped[1:-1])]
    return fields if len(fields) == _PARAMETER_FIELDS else None


def _parameter_value(text: str) -> tuple[str, str | None]:
    """Return an mzTab column value as (name, CV accession).

    A param contributes its name and accession; a plain value, such as a
    database name, is returned as written with no accession.
    """
    fields = _parameter_fields(text)
    if not fields or not fields[2]:
        return text, None
    return fields[2], fields[1] or None


def _value(row: dict[str, str], column: str) -> str | None:
    """Return a row's value for `column`, or None when absent or mzTab-null."""
    raw = (row.get(column) or "").strip()
    return None if raw.lower() in _NULL_VALUES else raw


def _parse_modification(segment: str) -> tuple[str, str] | None:
    """Split one mzTab modification object into (CV identifier, positions).

    Returns None for anything that does not name a modification: a neutral-loss
    cvParam, the "0" that means "this protein carries none", and "null". A
    per-site cvParam is stripped out of the positions, so a probability score
    never lands in the site list.
    """
    text = segment.strip()
    if text.lower() in _NULL_VALUES:
        return None
    placed = _MODIFICATION.match(text)
    if placed:
        positions = _BRACKETED.sub("", placed.group("positions"))
        return placed.group("accession"), "" if positions.lower() == "null" else positions
    bare = _BARE_MODIFICATION.match(text)
    return (bare.group("accession"), "") if bare else None


def _modification_claims(subject: str, row: dict[str, str], evidence_id: str) -> list[dict]:
    """Project a row's `modifications` column into `carries_modification` claims.

    One claim per modification listed. The site positions travel in evidence
    rather than in the object, so PTM localisation stays queryable without
    fragmenting the object into one value per site.
    """
    raw = _value(row, "modifications")
    if not raw:
        return []
    claims = []
    for segment in _split_top_level(raw):
        parsed = _parse_modification(segment)
        if not parsed:
            continue
        accession, positions = parsed
        claims.append(build_claim(
            subject, "carries_modification", accession,
            ontology_term=psimod_curie(accession),
            evidence={"positions": positions} if positions else None,
            evidence_id=evidence_id,
        ))
    return claims


def _abundance_claims(subject: str, row: dict[str, str], evidence_id: str) -> list[dict]:
    """Project a row's `*_abundance_assay[n]` columns into `has_abundance` claims."""
    claims = []
    for column, raw in row.items():
        match = _ABUNDANCE_COLUMN.match(column)
        if not match or not _value(row, column):
            continue
        claims.append(build_claim(
            subject, "has_abundance", match.group("label"),
            evidence={"abundance": raw.strip()},
            evidence_id=evidence_id,
        ))
    return claims


def _provenance_claims(accession: str, row: dict[str, str], evidence_id: str) -> list[dict]:
    """Project a PRT row's search_engine and database columns into claims.

    A search engine is reported as an mzTab param, so the claim object is the
    param's name and its PSI-MS accession rides along as the object's ontology
    term; a plain value such as a database name is claimed as written.
    """
    claims = []
    for column, qualifier in _PROVENANCE_COLUMNS:
        raw = _value(row, column)
        if not raw:
            continue
        name, term = _parameter_value(raw)
        claims.append(build_claim(
            accession, qualifier, name, ontology_term=term, evidence_id=evidence_id
        ))
    return claims


def _protein_claims(row: dict[str, str], file_id: str) -> list[dict]:
    """Project one PRT row into identification, provenance and quantity claims."""
    accession = _value(row, "accession")
    if not accession:
        return []
    evidence_id = f"mztab:{file_id}:PRT:{accession}"
    claims = [build_claim(
        f"mztab:{file_id}", "identifies_protein", accession,
        ontology_term=uniprot_curie(accession), evidence_id=evidence_id,
    )]
    claims += _provenance_claims(accession, row, evidence_id)
    claims += _modification_claims(accession, row, evidence_id)
    claims += _abundance_claims(accession, row, evidence_id)
    return claims


def _peptide_claims(row: dict[str, str], file_id: str) -> list[dict]:
    """Project one PEP row into a peptide-to-protein claim plus its quantities."""
    sequence = _value(row, "sequence")
    if not sequence:
        return []
    accession = _value(row, "accession") or f"mztab:{file_id}"
    evidence_id = f"mztab:{file_id}:PEP:{sequence}"
    claims = [build_claim(accession, "identifies_peptide", sequence, evidence_id=evidence_id)]
    claims += _modification_claims(sequence, row, evidence_id)
    claims += _abundance_claims(sequence, row, evidence_id)
    return claims


def _psm_claims(row: dict[str, str], file_id: str) -> list[dict]:
    """Project one PSM row into a spectrum match plus its modification claims."""
    sequence = _value(row, "sequence")
    if not sequence:
        return []
    psm_id = _value(row, "PSM_ID") or sequence
    evidence_id = f"mztab:{file_id}:PSM:{psm_id}"
    spectra_ref = _value(row, "spectra_ref") or f"mztab:{file_id}:PSM:{psm_id}"
    claims = [build_claim(sequence, "matched_to_spectrum", spectra_ref, evidence_id=evidence_id)]
    claims += _modification_claims(sequence, row, evidence_id)
    return claims


def _fdr_parameter_claims(value: str, file_id: str) -> list[dict]:
    """Project one `MTD false_discovery_rate` value into one claim per param."""
    claims = []
    for parameter in _BRACKETED.findall(value):
        fields = _parameter_fields(parameter)
        if not fields or not fields[3]:
            continue
        claims.append(build_claim(
            f"mztab:{file_id}", "fdr_controlled_at", fields[3],
            ontology_term=fields[1] or None,
            evidence={"fdr_scope": fields[2]} if fields[2] else None,
            evidence_id=f"mztab:{file_id}:MTD:false_discovery_rate",
        ))
    return claims


def _fdr_claims(metadata: list[list[str]], file_id: str) -> list[dict]:
    """Project `MTD false_discovery_rate` lines into `fdr_controlled_at` claims.

    The mzTab value is a param list such as
    ``[MS, MS:1001364, pep:global FDR, 0.01]``; each param becomes one claim
    whose object is the reported threshold and whose ontology term is the PSI-MS
    accession naming which FDR was controlled.
    """
    claims = []
    for fields in metadata:
        if len(fields) < 3 or fields[1] != "false_discovery_rate":
            continue
        claims += _fdr_parameter_claims(fields[2], file_id)
    return claims


def _read_file(path: Path) -> tuple[list[list[str]], dict[str, list[dict[str, str]]]]:
    """Read an mzTab file into its MTD lines and its three keyed row tables."""
    metadata: list[list[str]] = []
    headers: dict[str, list[str]] = {}
    tables: dict[str, list[dict[str, str]]] = {row: [] for _, row in _SECTIONS}
    header_of = {head: row for head, row in _SECTIONS}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            tag = fields[0]
            if tag == "MTD":
                metadata.append(fields)
            elif tag in header_of:
                headers[header_of[tag]] = fields
            elif tag in tables and tag in headers:
                # strict=False: sparse rows may omit trailing columns.
                tables[tag].append(dict(zip(headers[tag], fields, strict=False)))
    return metadata, tables


def parse_mztab_to_claims(mztab_path: Path | str) -> list[dict]:
    """Parse a proteomics mzTab file's MTD/PRT/PEP/PSM sections into claim dicts.

    Returns a list of claims in file order — FDR first, then proteins, peptides
    and PSMs. Raises FileNotFoundError if the file does not exist.
    """
    path = Path(mztab_path)
    if not path.exists():
        raise FileNotFoundError(f"mzTab file not found: {path}")
    metadata, tables = _read_file(path)
    file_id = path.stem
    claims = _fdr_claims(metadata, file_id)
    for row in tables["PRT"]:
        claims += _protein_claims(row, file_id)
    for row in tables["PEP"]:
        claims += _peptide_claims(row, file_id)
    for row in tables["PSM"]:
        claims += _psm_claims(row, file_id)
    return claims

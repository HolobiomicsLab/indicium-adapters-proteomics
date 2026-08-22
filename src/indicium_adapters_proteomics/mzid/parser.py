"""mzIdentML to asb:Claim projection — SequenceCollection and AnalysisData.

Parsed with the standard library's ElementTree using `{*}` wildcard tags, so the
same code reads mzIdentML 1.1 and 1.2 without pinning either namespace URI.

Claim shapes:

===================================  ============================================
mzIdentML source                     claim
===================================  ============================================
``Threshold`` FDR cvParam            (file, ``fdr_controlled_at``, value)
``PeptideEvidence`` (target)         (accession, ``identifies_peptide``, sequence)
distinct target-hit accession        (file, ``identifies_protein``, accession)
``Modification`` cvParam             (sequence, ``carries_modification``, CV id)
``SpectrumIdentificationItem``       (sequence, ``matched_to_spectrum``, spectrum)
===================================  ============================================

**What counts as an identified protein.** A ``DBSequence`` is a database entry,
not a result — a search database contributes thousands of them and identifies
none. This parser therefore emits ``identifies_protein`` only for accessions
reached through a ``PeptideEvidence`` whose ``isDecoy`` is not true, which is the
one rule that holds whether or not the file carries a ``ProteinDetectionList``.

**Quantification is absent by design.** mzIdentML describes identification;
abundances live in mzQuantML and in mzTab's abundance columns. Reading quantities
out of an identification file would mean inventing them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final
from xml.etree import ElementTree as ET

from indicium_adapters_proteomics._claims import build_claim
from indicium_adapters_proteomics._enrichment import psimod_curie, uniprot_curie

# PSI-MS names its false-discovery-rate terms with these words; the CV is the
# authority for the wording, so gating on it needs no vocabulary of our own.
_FDR_NAME_MARKERS: Final[tuple[str, ...]] = ("fdr", "false discovery rate")

_XSD_TRUE: Final[frozenset[str]] = frozenset({"true", "1"})


def _is_true(value: str | None) -> bool:
    """True iff an XSD boolean attribute is set and true."""
    return (value or "").strip().lower() in _XSD_TRUE


def _text_of(element: ET.Element, tag: str) -> str:
    """Return the stripped text of a child element, or an empty string."""
    child = element.find(f"{{*}}{tag}")
    return (child.text or "").strip() if child is not None and child.text else ""


def _peptide_sequences(root: ET.Element) -> dict[str, str]:
    """Map each Peptide element id to its amino-acid sequence."""
    return {
        peptide.get("id", ""): _text_of(peptide, "PeptideSequence")
        for peptide in root.iterfind(".//{*}Peptide")
    }


def _protein_accessions(root: ET.Element) -> dict[str, str]:
    """Map each DBSequence element id to its database accession."""
    return {
        entry.get("id", ""): (entry.get("accession") or "").strip()
        for entry in root.iterfind(".//{*}DBSequence")
    }


def _target_evidence(root: ET.Element) -> list[ET.Element]:
    """Return the PeptideEvidence elements that are not decoy hits."""
    return [
        evidence
        for evidence in root.iterfind(".//{*}PeptideEvidence")
        if not _is_true(evidence.get("isDecoy"))
    ]


def _evidence_span(evidence: ET.Element) -> dict[str, str]:
    """Return the start/end residue span of a PeptideEvidence, when reported."""
    return {key: evidence.get(key) for key in ("start", "end") if evidence.get(key)}


def _peptide_claims(evidence: ET.Element, names: tuple[dict, dict], file_id: str) -> list[dict]:
    """Project one non-decoy PeptideEvidence into an `identifies_peptide` claim.

    `names` is the (peptide id → sequence, DBSequence id → accession) pair built
    once per file.
    """
    sequences, accessions = names
    sequence = sequences.get(evidence.get("peptide_ref", ""), "")
    accession = accessions.get(evidence.get("dBSequence_ref", ""), "")
    if not sequence or not accession:
        return []
    return [build_claim(
        accession, "identifies_peptide", sequence,
        evidence=_evidence_span(evidence),
        evidence_id=f"mzid:{file_id}:PeptideEvidence:{evidence.get('id', '')}",
    )]


def _protein_claims(accessions: list[str], file_id: str) -> list[dict]:
    """Project the distinct target-hit accessions into `identifies_protein` claims."""
    seen: list[str] = []
    for accession in accessions:
        if accession and accession not in seen:
            seen.append(accession)
    return [
        build_claim(
            f"mzid:{file_id}", "identifies_protein", accession,
            ontology_term=uniprot_curie(accession),
            evidence_id=f"mzid:{file_id}:DBSequence:{accession}",
        )
        for accession in seen
    ]


def _modification_evidence(modification: ET.Element) -> dict[str, str]:
    """Return the site and mass shift a Modification element reports."""
    keys = ("location", "monoisotopicMassDelta")
    return {key: modification.get(key) for key in keys if modification.get(key)}


def _modification_claims(peptide: ET.Element, sequence: str, file_id: str) -> list[dict]:
    """Project a Peptide's Modification children into `carries_modification` claims.

    The modification site travels in evidence under `location`, keeping PTM
    localisation queryable without splitting the object across sites.
    """
    claims = []
    for modification in peptide.iterfind("{*}Modification"):
        for param in modification.iterfind("{*}cvParam"):
            accession = (param.get("accession") or "").strip()
            if not accession:
                continue
            claims.append(build_claim(
                sequence, "carries_modification", accession,
                ontology_term=psimod_curie(accession),
                evidence=_modification_evidence(modification),
                evidence_id=f"mzid:{file_id}:Peptide:{peptide.get('id', '')}",
            ))
    return claims


def _item_evidence(item: ET.Element) -> dict[str, str]:
    """Return the rank, charge and observed m/z a SpectrumIdentificationItem reports."""
    keys = ("rank", "chargeState", "experimentalMassToCharge")
    return {key: item.get(key) for key in keys if item.get(key)}


def _spectrum_claims(result: ET.Element, sequences: dict[str, str], file_id: str) -> list[dict]:
    """Project one SpectrumIdentificationResult into `matched_to_spectrum` claims.

    Only items that passed the search's own threshold are projected; a failing
    item is a rejected hypothesis, not a match.
    """
    spectrum_id = (result.get("spectrumID") or "").strip()
    claims = []
    for item in result.iterfind("{*}SpectrumIdentificationItem"):
        sequence = sequences.get(item.get("peptide_ref", ""), "")
        if not sequence or not _is_true(item.get("passThreshold")):
            continue
        claims.append(build_claim(
            sequence, "matched_to_spectrum", spectrum_id,
            evidence=_item_evidence(item),
            evidence_id=f"mzid:{file_id}:SpectrumIdentificationItem:{item.get('id', '')}",
        ))
    return claims


def _is_fdr_param(param: ET.Element) -> bool:
    """True iff a Threshold cvParam names a false-discovery-rate term with a value."""
    if not (param.get("value") or "").strip():
        return False
    name = (param.get("name") or "").lower()
    return any(marker in name for marker in _FDR_NAME_MARKERS)


def _fdr_claims(root: ET.Element, file_id: str) -> list[dict]:
    """Project FDR cvParams under any Threshold element into `fdr_controlled_at`."""
    return [
        build_claim(
            f"mzid:{file_id}", "fdr_controlled_at", (param.get("value") or "").strip(),
            ontology_term=(param.get("accession") or "").strip() or None,
            evidence={"fdr_scope": param.get("name")} if param.get("name") else None,
            evidence_id=f"mzid:{file_id}:Threshold:{(param.get('accession') or '').strip()}",
        )
        for param in root.iterfind(".//{*}Threshold/{*}cvParam")
        if _is_fdr_param(param)
    ]


def _identification_claims(root: ET.Element, file_id: str) -> list[dict]:
    """Project the SequenceCollection: proteins, peptides and their modifications."""
    names = (_peptide_sequences(root), _protein_accessions(root))
    evidence_elements = _target_evidence(root)
    hit_accessions = [names[1].get(e.get("dBSequence_ref", ""), "") for e in evidence_elements]
    claims = _protein_claims(hit_accessions, file_id)
    for evidence in evidence_elements:
        claims += _peptide_claims(evidence, names, file_id)
    for peptide in root.iterfind(".//{*}Peptide"):
        sequence = names[0].get(peptide.get("id", ""), "")
        if sequence:
            claims += _modification_claims(peptide, sequence, file_id)
    return claims


def parse_mzid_to_claims(mzid_path: Path | str) -> list[dict]:
    """Parse an mzIdentML file into ASB claim dicts.

    Returns claims in a stable order — FDR thresholds, then proteins, peptides,
    modifications and spectrum matches. Raises FileNotFoundError if the file does
    not exist, and `xml.etree.ElementTree.ParseError` if it is not well-formed
    XML.
    """
    path = Path(mzid_path)
    if not path.exists():
        raise FileNotFoundError(f"mzIdentML file not found: {path}")
    root = ET.parse(path).getroot()
    file_id = path.stem
    claims = _fdr_claims(root, file_id)
    claims += _identification_claims(root, file_id)
    sequences = _peptide_sequences(root)
    for result in root.iterfind(".//{*}SpectrumIdentificationResult"):
        claims += _spectrum_claims(result, sequences, file_id)
    return claims

"""Tests for the mzIdentML to asb:Claim projection."""

from __future__ import annotations

from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from indicium_adapters_proteomics.mzid import parse_mzid_to_claims


def _by_qualifier(claims: list[dict], qualifier: str) -> list[dict]:
    return [c for c in claims if c.get("qualifier") == qualifier]


def test_parses_fixture_without_error(mzid_path):
    claims = parse_mzid_to_claims(mzid_path)
    assert isinstance(claims, list)
    assert len(claims) >= 1


def test_missing_file_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        parse_mzid_to_claims(Path("/nonexistent/sample.mzid"))


def test_malformed_xml_raises_parse_error(tmp_path):
    broken = tmp_path / "broken.mzid"
    broken.write_text("<MzIdentML><SequenceCollection>", encoding="utf-8")
    with pytest.raises(ET.ParseError):
        parse_mzid_to_claims(broken)


def test_projection_is_deterministic(mzid_path):
    assert parse_mzid_to_claims(mzid_path) == parse_mzid_to_claims(mzid_path)


# --- protein identification -------------------------------------------------


def test_only_accessions_with_a_target_hit_are_identified(mzid_path):
    """A DBSequence is a database entry, not a result. P69905 is in the search
    database with no peptide evidence, so it must not be reported identified."""
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "identifies_protein")
    assert [c["object"] for c in claims] == ["P02768"]
    assert claims[0]["ontology_terms"]["object"] == "uniprot:P02768"


def test_decoy_evidence_produces_no_claims(mzid_path):
    claims = parse_mzid_to_claims(mzid_path)
    assert all("REV_Q9Y6K9" != c["object"] for c in claims)


def test_repeated_hits_yield_one_protein_claim(mzid_path):
    """Two peptides map to P02768; the protein must be claimed once."""
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "identifies_protein")
    assert len(claims) == 1


# --- peptide identification -------------------------------------------------


def test_peptide_evidence_links_peptide_to_protein_with_its_span(mzid_path):
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "identifies_peptide")
    assert {c["object"] for c in claims} == {"LVNEVTEFAK", "YICENQDSISSK"}
    first = next(c for c in claims if c["object"] == "LVNEVTEFAK")
    assert first["subject"] == "P02768"
    assert first["evidence"] == {"start": "66", "end": "75"}


# --- spectrum matches -------------------------------------------------------


def test_only_threshold_passing_items_become_spectrum_matches(mzid_path):
    """A SpectrumIdentificationItem with passThreshold=false is a rejected
    hypothesis, not a match."""
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "matched_to_spectrum")
    assert len(claims) == 1
    assert claims[0]["subject"] == "LVNEVTEFAK"
    assert claims[0]["object"] == "index=42"
    assert claims[0]["evidence"]["chargeState"] == "2"


# --- modifications ----------------------------------------------------------


def test_psi_mod_modification_resolves_and_keeps_its_location(mzid_path):
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "carries_modification")
    resolved = next(c for c in claims if c["object"] == "MOD:01060")
    assert resolved["subject"] == "YICENQDSISSK"
    assert resolved["ontology_terms"]["object"] == "MOD:01060"
    assert resolved["evidence"]["location"] == "3"


def test_unimod_modification_is_carried_verbatim_and_left_unresolved(mzid_path):
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "carries_modification")
    unresolved = next(c for c in claims if c["object"] == "UNIMOD:21")
    assert unresolved["ontology_terms"] == {}


# --- FDR --------------------------------------------------------------------


def test_threshold_fdr_param_becomes_an_fdr_claim(mzid_path):
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "fdr_controlled_at")
    assert len(claims) == 1
    assert claims[0]["object"] == "0.01"
    assert claims[0]["ontology_terms"]["object"] == "MS:1001364"


def test_valueless_threshold_params_are_not_fdr_claims(mzid_path):
    """'no threshold' (MS:1001494) carries no value and names no FDR, so it
    must not be read as a controlled error rate."""
    claims = _by_qualifier(parse_mzid_to_claims(mzid_path), "fdr_controlled_at")
    assert all(c["ontology_terms"]["object"] != "MS:1001494" for c in claims)


# --- namespace tolerance ----------------------------------------------------


def test_mzidentml_1_1_namespace_parses_the_same_way(tmp_path, mzid_path):
    """Wildcard tag matching means the 1.1 and 1.2 namespaces read alike."""
    text = mzid_path.read_text(encoding="utf-8").replace("mzIdentML/1.2", "mzIdentML/1.1")
    older = tmp_path / "v11.mzid"
    older.write_text(text, encoding="utf-8")
    assert len(parse_mzid_to_claims(older)) == len(parse_mzid_to_claims(mzid_path))


# --- dangling and incomplete references -------------------------------------


def test_evidence_pointing_at_an_unknown_protein_yields_nothing(tmp_path):
    """A dangling dBSequence_ref names no protein, so no link can be asserted."""
    mzid = tmp_path / "dangling.mzid"
    mzid.write_text(
        '<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.2">'
        "<SequenceCollection>"
        "<Peptide id=\"Pep_1\"><PeptideSequence>PEPTIDER</PeptideSequence></Peptide>"
        '<PeptideEvidence id="PE_1" dBSequence_ref="missing" peptide_ref="Pep_1"/>'
        "</SequenceCollection></MzIdentML>",
        encoding="utf-8",
    )
    assert parse_mzid_to_claims(mzid) == []


def test_modification_cvparam_without_an_accession_yields_nothing(tmp_path):
    mzid = tmp_path / "unnamed_mod.mzid"
    mzid.write_text(
        '<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.2">'
        "<SequenceCollection>"
        '<Peptide id="Pep_1"><PeptideSequence>PEPTIDER</PeptideSequence>'
        '<Modification location="2"><cvParam name="unknown modification"/></Modification>'
        "</Peptide></SequenceCollection></MzIdentML>",
        encoding="utf-8",
    )
    assert parse_mzid_to_claims(mzid) == []

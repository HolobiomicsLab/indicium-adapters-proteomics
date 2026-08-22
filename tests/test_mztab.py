"""Tests for the proteomics mzTab to asb:Claim projection."""

from __future__ import annotations

from pathlib import Path

import pytest

from indicium_adapters_proteomics.mztab import parse_mztab_to_claims


def _by_qualifier(claims: list[dict], qualifier: str) -> list[dict]:
    return [c for c in claims if c.get("qualifier") == qualifier]


def test_parses_fixture_without_error(mztab_path):
    claims = parse_mztab_to_claims(mztab_path)
    assert isinstance(claims, list)
    assert len(claims) >= 1


def test_missing_file_raises_filenotfound():
    with pytest.raises(FileNotFoundError):
        parse_mztab_to_claims(Path("/nonexistent/sample.mztab"))


def test_projection_is_deterministic(mztab_path):
    assert parse_mztab_to_claims(mztab_path) == parse_mztab_to_claims(mztab_path)


# --- identification ---------------------------------------------------------


def test_prt_row_yields_identifies_protein_with_uniprot_curie(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "identifies_protein")
    hits = {c["object"]: c for c in claims}
    assert set(hits) == {"P02768", "NOTANACC"}
    assert hits["P02768"]["ontology_terms"]["object"] == "uniprot:P02768"
    assert hits["P02768"]["evidence_id"].startswith("mztab:")


def test_non_uniprot_accession_carries_no_ontology_term(mztab_path):
    """An accession that fails the UniProtKB grammar still yields a claim — it
    just yields one with no CURIE attached."""
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "identifies_protein")
    other = next(c for c in claims if c["object"] == "NOTANACC")
    assert other["ontology_terms"] == {}


def test_prt_search_engine_and_database_become_provenance_claims(mztab_path):
    claims = parse_mztab_to_claims(mztab_path)
    engines = _by_qualifier(claims, "identified_by")
    databases = _by_qualifier(claims, "searched_against")
    assert [c["subject"] for c in engines] == ["P02768"]
    # A search engine is an mzTab param, so the claim carries the engine's name
    # and its PSI-MS accession, not the serialised param.
    assert engines[0]["object"] == "Mascot"
    assert engines[0]["ontology_terms"]["object"] == "MS:1001207"
    # The second PRT row reports a null search engine but a real database, so
    # it contributes a searched_against claim and no identified_by claim.
    assert [c["subject"] for c in databases] == ["P02768", "NOTANACC"]
    assert {c["object"] for c in databases} == {"UniProtKB"}


def test_pep_row_links_peptide_to_its_protein(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "identifies_peptide")
    assert len(claims) == 1
    assert claims[0]["subject"] == "P02768"
    assert claims[0]["object"] == "LVNEVTEFAK"


def test_psm_row_links_peptide_to_its_spectrum(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "matched_to_spectrum")
    assert len(claims) == 1
    assert claims[0]["subject"] == "LVNEVTEFAK"
    assert claims[0]["object"] == "ms_run[1]:index=42"


# --- FDR --------------------------------------------------------------------


def test_mtd_false_discovery_rate_yields_one_claim_per_parameter(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "fdr_controlled_at")
    assert [c["object"] for c in claims] == ["0.01", "0.05"]
    assert claims[0]["ontology_terms"]["object"] == "MS:1001364"
    assert claims[0]["evidence"]["fdr_scope"] == "pep:global FDR"
    # MS:1001214 is PSI-MS "protein-level global FDR"; the accession and the
    # scope it travels with must name the same thing.
    assert claims[1]["ontology_terms"]["object"] == "MS:1001214"
    assert claims[1]["evidence"]["fdr_scope"] == "protein:global FDR"


# --- modifications ----------------------------------------------------------


def test_psi_mod_modification_resolves_and_keeps_its_site(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "carries_modification")
    protein_mod = next(c for c in claims if c["subject"] == "P02768")
    assert protein_mod["object"] == "MOD:00696"
    assert protein_mod["ontology_terms"]["object"] == "MOD:00696"
    assert protein_mod["evidence"]["positions"] == "3"


def test_unimod_modification_is_carried_verbatim_and_left_unresolved(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "carries_modification")
    peptide_mod = next(c for c in claims if c["object"] == "UNIMOD:35")
    assert peptide_mod["subject"] == "LVNEVTEFAK"
    assert peptide_mod["ontology_terms"] == {}
    assert peptide_mod["evidence"]["positions"] == "3|5"


def test_multiple_modifications_on_one_row_yield_one_claim_each(tmp_path):
    mztab = tmp_path / "multi.mztab"
    mztab.write_text(
        "PSH\tsequence\tPSM_ID\tmodifications\tspectra_ref\n"
        "PSM\tPEPTIDER\t1\t3-MOD:00696,7-MOD:00696\tms_run[1]:index=1\n",
        encoding="utf-8",
    )
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "carries_modification")
    assert [c["evidence"]["positions"] for c in claims] == ["3", "7"]


def test_null_modification_column_yields_no_modification_claim(tmp_path):
    mztab = tmp_path / "null_mods.mztab"
    mztab.write_text(
        "PSH\tsequence\tPSM_ID\tmodifications\tspectra_ref\n"
        "PSM\tPEPTIDER\t1\tnull\tms_run[1]:index=1\n",
        encoding="utf-8",
    )
    assert not _by_qualifier(parse_mztab_to_claims(mztab), "carries_modification")


# --- quantification ---------------------------------------------------------


def test_abundance_columns_become_has_abundance_claims(mztab_path):
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "has_abundance")
    quantified = {c["subject"]: c for c in claims}
    assert set(quantified) == {"P02768", "LVNEVTEFAK"}
    assert quantified["P02768"]["object"] == "assay[1]"
    assert quantified["P02768"]["evidence"]["abundance"] == "1.2E6"


def test_dispersion_columns_are_not_read_as_abundances(mztab_path):
    """protein_abundance_stdev_assay[1] describes a measurement rather than
    being one, so it must not produce a has_abundance claim."""
    claims = _by_qualifier(parse_mztab_to_claims(mztab_path), "has_abundance")
    assert all(c["evidence"]["abundance"] != "0.11" for c in claims)


# --- section isolation ------------------------------------------------------


def test_small_molecule_sections_are_ignored(tmp_path):
    """mzTab-M SML rows belong to the metabolomics sibling; reading them here
    would mis-type a metabolite as a peptide."""
    mztab = tmp_path / "small_molecules.mztab"
    mztab.write_text(
        "SMH\tSML_ID\tchemical_name\tdatabase_identifier\n"
        "SML\t1\tD-Glucose\tCHEBI:17234\n",
        encoding="utf-8",
    )
    assert parse_mztab_to_claims(mztab) == []


def test_rows_without_a_header_are_skipped(tmp_path):
    mztab = tmp_path / "headerless.mztab"
    mztab.write_text("PRT\tP02768\tSerum albumin\n", encoding="utf-8")
    assert parse_mztab_to_claims(mztab) == []


# --- sparse and malformed rows ----------------------------------------------


def test_positionless_modification_is_claimed_without_a_site(tmp_path):
    """A CHEMMOD delta reported with no site still names a modification; only
    its localisation is unknown."""
    mztab = tmp_path / "chemmod.mztab"
    mztab.write_text(
        "PSH\tsequence\tPSM_ID\tmodifications\tspectra_ref\n"
        "PSM\tPEPTIDER\t1\tCHEMMOD:+15.9949\tms_run[1]:index=1\n",
        encoding="utf-8",
    )
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "carries_modification")
    assert [c["object"] for c in claims] == ["CHEMMOD:+15.9949"]
    assert "evidence" not in claims[0]


def test_prt_row_without_an_accession_yields_nothing(tmp_path):
    mztab = tmp_path / "no_accession.mztab"
    mztab.write_text("PRH\taccession\tdescription\nPRT\tnull\tunnamed\n", encoding="utf-8")
    assert parse_mztab_to_claims(mztab) == []


def test_psm_row_without_a_sequence_yields_nothing(tmp_path):
    mztab = tmp_path / "no_sequence.mztab"
    mztab.write_text("PSH\tsequence\tPSM_ID\nPSM\tnull\t1\n", encoding="utf-8")
    assert parse_mztab_to_claims(mztab) == []


def test_peptide_without_a_protein_is_attributed_to_the_file(tmp_path):
    """An unattributed peptide is still an identification; hanging it off the
    file keeps it in the graph instead of dropping it."""
    mztab = tmp_path / "orphan_peptide.mztab"
    mztab.write_text("PEH\tsequence\taccession\nPEP\tPEPTIDER\tnull\n", encoding="utf-8")
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "identifies_peptide")
    assert claims[0]["subject"] == "mztab:orphan_peptide"


def test_psm_without_a_spectra_ref_falls_back_to_its_own_evidence_id(tmp_path):
    mztab = tmp_path / "no_spectra_ref.mztab"
    mztab.write_text("PSH\tsequence\tPSM_ID\nPSM\tPEPTIDER\t7\n", encoding="utf-8")
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "matched_to_spectrum")
    assert claims[0]["object"] == claims[0]["evidence_id"] == "mztab:no_spectra_ref:PSM:7"


def test_fdr_parameter_without_a_value_yields_no_claim(tmp_path):
    """A parameter list that names an FDR term but reports no threshold states
    nothing about the error rate."""
    mztab = tmp_path / "empty_fdr.mztab"
    mztab.write_text(
        "MTD\tfalse_discovery_rate\t[MS, MS:1001364, pep:global FDR, ]\n", encoding="utf-8"
    )
    assert parse_mztab_to_claims(mztab) == []


# --- embedded cv parameters (mzTab §5.8, §5.10) ------------------------------


def _psm_modifications(tmp_path: Path, column: str) -> list[dict]:
    """Parse a one-row PSM file whose modifications column is `column`."""
    mztab = tmp_path / "mods.mztab"
    mztab.write_text(
        "PSH\tsequence\tPSM_ID\tmodifications\tspectra_ref\n"
        f"PSM\tPEPTIDER\t1\t{column}\tms_run[1]:index=1\n",
        encoding="utf-8",
    )
    return _by_qualifier(parse_mztab_to_claims(mztab), "carries_modification")


def test_a_neutral_loss_parameter_is_not_claimed_as_a_modification(tmp_path):
    """The spec reports a neutral loss as its own comma-separated cvParam. Its
    commas must not split it into fragments, one of which then reads as a
    further modification: MS:1001524 names a fragment, not a PTM."""
    claims = _psm_modifications(
        tmp_path, "3-UNIMOD:21,3-[MS, MS:1001524, fragment neutral loss, 63.998285]"
    )
    assert [c["object"] for c in claims] == ["UNIMOD:21"]


def test_a_neutral_loss_without_a_modification_yields_no_claim(tmp_path):
    claims = _psm_modifications(
        tmp_path, "[MS, MS:1001524, fragment neutral loss, 63.998285],7-UNIMOD:4"
    )
    assert [c["object"] for c in claims] == ["UNIMOD:4"]


def test_a_per_site_probability_parameter_keeps_the_modification_and_its_sites(tmp_path):
    """The spec's own PSM example. The probability param must neither destroy
    the modification nor leak into the site list."""
    claims = _psm_modifications(
        tmp_path,
        "3[MS,MS:1001876, modification probability, 0.8]"
        "|4[MS,MS:1001876, modification probability, 0.2]-MOD:00412,8-MOD:00412",
    )
    assert [c["object"] for c in claims] == ["MOD:00412", "MOD:00412"]
    assert [c["evidence"]["positions"] for c in claims] == ["3|4", "8"]


def test_zero_means_no_modifications_rather_than_one_named_zero(tmp_path):
    """mzTab writes "0" when a protein was searched and carries none."""
    assert _psm_modifications(tmp_path, "0") == []


def test_a_modification_with_an_unknown_site_is_still_claimed(tmp_path):
    """The site is mandatory but may be "null" for MS1 peptide mass
    fingerprinting; the modification is known even when its position is not."""
    claims = _psm_modifications(tmp_path, "null-MOD:00412")
    assert [c["object"] for c in claims] == ["MOD:00412"]
    assert "evidence" not in claims[0]


def test_a_quoted_comma_in_a_parameter_name_does_not_split_the_parameter(tmp_path):
    """mzTab quotes a param name containing a comma. Splitting on that comma
    would shift the fields and report the name fragment as the threshold."""
    mztab = tmp_path / "quoted_fdr.mztab"
    mztab.write_text(
        'MTD\tfalse_discovery_rate\t[MS, MS:1001364, "global FDR, peptide level", 0.01]\n',
        encoding="utf-8",
    )
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "fdr_controlled_at")
    assert [c["object"] for c in claims] == ["0.01"]
    assert claims[0]["evidence"]["fdr_scope"] == "global FDR, peptide level"


def test_a_user_parameter_search_engine_is_claimed_by_name(tmp_path):
    """A param with no CV label or accession still names the engine."""
    mztab = tmp_path / "user_param.mztab"
    mztab.write_text(
        "PRH\taccession\tsearch_engine\nPRT\tP02768\t[,,in-house search, ]\n",
        encoding="utf-8",
    )
    claims = _by_qualifier(parse_mztab_to_claims(mztab), "identified_by")
    assert claims[0]["object"] == "in-house search"
    assert claims[0]["ontology_terms"] == {}


def test_a_null_entry_among_several_modifications_is_skipped(tmp_path):
    """One unreported site does not invalidate the sites reported beside it."""
    claims = _psm_modifications(tmp_path, "3-MOD:00696,null,7-UNIMOD:35")
    assert [c["object"] for c in claims] == ["MOD:00696", "UNIMOD:35"]


def test_pep_row_without_a_sequence_yields_nothing(tmp_path):
    mztab = tmp_path / "no_pep_sequence.mztab"
    mztab.write_text("PEH\tsequence\taccession\nPEP\tnull\tP02768\n", encoding="utf-8")
    assert parse_mztab_to_claims(mztab) == []

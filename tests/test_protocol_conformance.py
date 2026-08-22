"""Protocol-conformance contract for the proteomics adapter.

The invariants every indicium DomainAdapter must hold, plus the two this
adapter exists to satisfy for ASB: it declares ``primary_domain`` (without it a
proteomics paper classifies as "other"), and it resolves concepts only by exact
match or by a published accession grammar (a wrong CURIE is worse than no CURIE).
"""

from __future__ import annotations

import copy
import re
from urllib.parse import urlparse

import pytest
from indicium_adapters import CURRENT_API_VERSION, DomainAdapter, SHACLProvider

SNAKE_CASE = re.compile(r"^[a-z][a-z0-9_]*$")
SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
# "=" is included because repository landing pages are query-parameter URLs
# (ProteomeXchange, MassIVE), where the accession appends after "ID=".
PREFIX_TERMINATORS = (":", "/", "_", "#", "=")


def test_satisfies_the_domain_adapter_protocol(adapter):
    assert isinstance(adapter, DomainAdapter)
    assert adapter.api_version == CURRENT_API_VERSION


def test_ships_shacl_shapes(adapter):
    assert isinstance(adapter, SHACLProvider)


def test_declares_primary_domain_so_asb_does_not_classify_as_other(adapter):
    assert adapter.primary_domain == "proteomics"
    assert adapter.domain_id == "proteomics"


def test_qualifiers_are_snake_case(adapter):
    """Perspicacité's asb:domainQualifier SHACL shape enforces this pattern."""
    bad = [q for q in adapter.qualifiers if not SNAKE_CASE.match(q)]
    assert not bad, f"non-snake_case qualifiers: {bad}"


@pytest.mark.parametrize(
    "qualifier",
    [
        "identifies_protein",
        "identifies_peptide",
        "matched_to_spectrum",
        "fdr_controlled_at",
        "carries_modification",
        "modified_at",
        "has_abundance",
        "quantified_by",
    ],
)
def test_qualifier_vocabulary_covers_the_projected_claim_shapes(adapter, qualifier):
    """Every qualifier the mzTab and mzIdentML parsers emit must be declared,
    or a projected claim would carry a qualifier the adapter does not admit."""
    assert qualifier in adapter.qualifiers


def test_vocabulary_terms_are_slugs(adapter):
    for term in (*adapter.subdomains, *adapter.techniques):
        assert SLUG.match(term), f"not a kebab-case slug: {term}"


def test_vocabulary_has_no_duplicates(adapter):
    assert len(set(adapter.subdomains)) == len(adapter.subdomains)
    assert len(set(adapter.techniques)) == len(adapter.techniques)


def test_ontology_prefixes_expand_to_curie_safe_iris(adapter):
    for prefix, base in adapter.ontology_prefixes.items():
        parsed = urlparse(base)
        assert parsed.scheme in ("http", "https"), f"{prefix} is not http(s)"
        assert base.endswith(PREFIX_TERMINATORS), f"{prefix} does not end CURIE-safely: {base}"


def test_extraction_context_names_the_domain_and_its_qualifiers(adapter):
    context = adapter.extraction_context()
    assert "proteomics" in context.lower()
    for qualifier in ("identified_by", "quantified_by", "digested_with"):
        assert qualifier in context


# --- enrichment ---------------------------------------------------------------


def test_enrich_claim_is_pure(adapter, base_claim):
    before = copy.deepcopy(base_claim)
    adapter.enrich_claim(base_claim)
    assert base_claim == before, "enrich_claim mutated its input"


def test_uniprot_accession_resolves_with_a_composite_justification(adapter, base_claim):
    """An accession is matched against a published grammar, not against a label,
    so the justification must not claim a lexical match."""
    out = adapter.enrich_claim(base_claim)
    assert out["ontology_terms"]["subject"] == "uniprot:P02768"
    assert out["ontology_term_justifications"]["subject"] == "semapv:CompositeMatching"


def test_ptm_name_resolves_to_psi_mod(adapter, base_claim):
    out = adapter.enrich_claim(base_claim)
    assert out["ontology_terms"]["object"] == "MOD:00696"
    assert out["ontology_term_justifications"]["object"] == "semapv:LexicalMatching"


def test_known_label_resolves_to_an_edam_curie(adapter):
    out = adapter.enrich_claim({"subject": "Proteomics", "object": "MaxQuant"})
    assert out["ontology_terms"]["subject"].startswith("edam:")
    assert out["ontology_term_justifications"]["subject"] == "semapv:LexicalMatching"


def test_lookup_is_case_and_whitespace_insensitive(adapter):
    spaced = adapter.enrich_claim({"subject": "  pRoTeOmIcS  "})
    exact = adapter.enrich_claim({"subject": "Proteomics"})
    assert spaced["ontology_terms"]["subject"] == exact["ontology_terms"]["subject"]


def test_unresolved_text_is_recorded_not_silently_dropped(adapter):
    out = adapter.enrich_claim({"subject": "a phrase that is not an EDAM concept"})
    assert out["ontology_terms_unmatched"]["subject"] == "a phrase that is not an EDAM concept"
    assert "ontology_terms" not in out


def test_substring_of_a_concept_does_not_resolve(adapter):
    """'Peptide ident' must not pick up 'Peptide identification' — a wrong CURIE
    is worse than no CURIE."""
    out = adapter.enrich_claim({"subject": "Peptide ident"})
    assert "ontology_terms" not in out


def test_an_already_resolved_slot_is_left_to_its_owner(adapter):
    claim = {"subject": "Proteomics", "ontology_terms": {"subject": "MS:1000031"}}
    out = adapter.enrich_claim(claim)
    assert out["ontology_terms"]["subject"] == "MS:1000031"


def test_empty_and_non_string_slots_are_skipped(adapter):
    out = adapter.enrich_claim({"subject": "", "object": None, "qualifier": 42})
    assert "ontology_terms" not in out
    assert "ontology_terms_unmatched" not in out


def test_a_resolved_slot_is_cleared_from_a_stale_unmatched_entry(adapter):
    """An upstream pass may have recorded a slot as unmatched; resolving it here
    must retract that record rather than leave the claim contradicting itself."""
    claim = {"subject": "P02768", "ontology_terms_unmatched": {"subject": "P02768"}}
    out = adapter.enrich_claim(claim)
    assert out["ontology_terms"]["subject"] == "uniprot:P02768"
    assert "ontology_terms_unmatched" not in out

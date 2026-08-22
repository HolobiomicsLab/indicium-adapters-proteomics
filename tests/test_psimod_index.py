"""Tests for the PSI-MOD index builder.

The builder runs offline against a PSI-MOD release, so the rules that matter —
what gets indexed, what gets dropped — are exercised here against small OBO
documents written inline rather than against the 1.8 MB release file.
"""

from __future__ import annotations

from indicium_adapters_proteomics.data.build_psimod_index import (
    INDEXED_SYNONYM_TYPES,
    build_index,
)

HEADER = "format-version: 1.2\nontology: mod\ndata-version: 1.038.0\n"


def _obo(*stanzas: str) -> str:
    """Assemble an OBO document from `[Term]` stanza bodies."""
    return HEADER + "".join(f"\n[Term]\n{body}\n" for body in stanzas)


def test_term_name_is_indexed_case_folded():
    index = build_index(_obo("id: MOD:00001\nname: Phosphorylated Residue"))
    assert index["phosphorylated residue"] == "MOD:00001"


def test_curated_synonyms_are_indexed():
    stanza = 'id: MOD:00002\nname: acetylated residue\nsynonym: "Acetyl" EXACT PSI-MOD-label []'
    index = build_index(_obo(stanza))
    assert index["acetyl"] == "MOD:00002"


def test_resid_systematic_synonyms_are_excluded():
    """IUPAC chemistry strings add collision surface without adding recall."""
    assert "RESID-systematic" not in INDEXED_SYNONYM_TYPES
    stanza = (
        "id: MOD:00003\nname: D-methionine\n"
        'synonym: "(2R)-2-amino-4-(methylsulfanyl)butanoic acid" EXACT RESID-systematic []'
    )
    index = build_index(_obo(stanza))
    assert "(2r)-2-amino-4-(methylsulfanyl)butanoic acid" not in index


def test_a_synonym_claimed_by_two_terms_is_dropped_entirely():
    """Arbitrating an ambiguous label would emit a confident wrong CURIE."""
    shared = 'synonym: "Phospho" RELATED PSI-MS-label []'
    index = build_index(
        _obo(
            f"id: MOD:00010\nname: O-phospho-L-serine\n{shared}",
            f"id: MOD:00011\nname: O-phospho-L-threonine\n{shared}",
        )
    )
    assert "phospho" not in index
    assert index["o-phospho-l-serine"] == "MOD:00010"


def test_a_synonym_repeated_by_one_term_survives():
    stanza = (
        'id: MOD:00012\nname: sulfated residue\nsynonym: "Sulfo" RELATED PSI-MS-label []\n'
        'synonym: "Sulfo" RELATED Unimod-interim []'
    )
    index = build_index(_obo(stanza))
    assert index["sulfo"] == "MOD:00012"


def test_a_term_name_wins_over_a_rival_synonym():
    index = build_index(
        _obo(
            "id: MOD:00020\nname: oxidized residue",
            'id: MOD:00021\nname: something else\nsynonym: "oxidized residue" EXACT RESID-name []',
        )
    )
    assert index["oxidized residue"] == "MOD:00020"


def test_obsolete_terms_are_not_indexed():
    stanza = "id: MOD:00030\nname: retired modification\nis_obsolete: true"
    assert "retired modification" not in build_index(_obo(stanza))


def test_synonyms_without_a_declared_type_are_skipped():
    """A synonym with no curation provenance carries no authority."""
    stanza = 'id: MOD:00040\nname: kept\nsynonym: "loose" []'
    index = build_index(_obo(stanza))
    assert "loose" not in index
    assert index["kept"] == "MOD:00040"

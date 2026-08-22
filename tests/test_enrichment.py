"""Tests for the two identity lookups: UniProtKB accessions and PSI-MOD labels.

The negative cases carry the weight here. An accession grammar that accepts
ordinary words, or a modification index that arbitrates an ambiguous label,
would emit confident wrong CURIEs — the failure this package is built to avoid.
"""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from indicium_adapters_proteomics import _edam
from indicium_adapters_proteomics._enrichment import (
    psimod_curie,
    psimod_size,
    psimod_version,
    uniprot_curie,
)

# --- UniProtKB accessions ---------------------------------------------------


@pytest.mark.parametrize(
    "text,expected",
    [
        ("P02768", "uniprot:P02768"),  # Swiss-Prot, [OPQ] grammar
        ("Q9Y6K9", "uniprot:Q9Y6K9"),
        ("A0A024R161", "uniprot:A0A024R161"),  # 10-character TrEMBL form
        ("P02768-2", "uniprot:P02768-2"),  # isoform suffix
        ("p02768", "uniprot:P02768"),  # case-normalised
        ("  P02768  ", "uniprot:P02768"),
        ("sp|P02768|ALBU_HUMAN", "uniprot:P02768"),  # UniProtKB FASTA header
        ("tr|A0A024R161|A0A024R161_HUMAN", "uniprot:A0A024R161"),
        ("UniProtKB:P02768", "uniprot:P02768"),  # already a CURIE
        ("uniprot:P02768", "uniprot:P02768"),
    ],
)
def test_uniprot_curie_recognises_published_accession_forms(text, expected):
    assert uniprot_curie(text) == expected


@pytest.mark.parametrize(
    "text",
    ["", "   ", "MaxQuant", "albumin", "REV_Q9Y6K9", "PXD000001", "GAPDH", "P0276"],
)
def test_uniprot_curie_returns_none_on_non_accessions(text):
    """Decoy prefixes, gene symbols, repository accessions and truncated
    accessions must all fail the UniProtKB grammar rather than be guessed at."""
    assert uniprot_curie(text) is None


# --- PSI-MOD modifications --------------------------------------------------


@pytest.mark.parametrize(
    "text,expected",
    [
        ("phosphorylated residue", "MOD:00696"),  # PSI-MOD preferred label
        ("Phosphorylated Residue", "MOD:00696"),  # case-insensitive
        ("  phosphorylated   residue ", "MOD:00696"),  # whitespace-collapsed
        ("phosphorylation", "MOD:00696"),  # process-noun alias
        ("MOD:00696", "MOD:00696"),  # already a CURIE
        ("mod:00696", "MOD:00696"),
        ("S-carboxamidomethyl-L-cysteine", "MOD:01060"),
        ("O-phospho-L-serine", "MOD:00046"),
    ],
)
def test_psimod_curie_resolves_labels_synonyms_and_aliases(text, expected):
    assert psimod_curie(text) == expected


@pytest.mark.parametrize("text", ["", "   ", "not a modification", "Peptide ident"])
def test_psimod_curie_returns_none_on_unknown_text(text):
    assert psimod_curie(text) is None


@pytest.mark.parametrize("text", ["Phospho", "Carbamidomethyl", "iTRAQ4plex"])
def test_ambiguous_short_labels_resolve_to_nothing(text):
    """These PSI-MS labels name several residue-specific MOD terms each. The
    index drops them rather than arbitrating, because a wrong identifier is
    worse than no identifier."""
    assert psimod_curie(text) is None


def test_unimod_identifiers_are_not_translated():
    """No offline one-to-one UNIMOD-to-PSI-MOD mapping exists, so a UNIMOD
    identifier must stay unresolved rather than be mapped by guesswork."""
    assert psimod_curie("UNIMOD:21") is None


# --- bundled index provenance -----------------------------------------------


def test_index_records_the_psimod_release_it_was_cut_from():
    assert psimod_version() != "unknown"
    assert psimod_size() > 1000


def test_every_ptm_alias_target_exists_in_the_index():
    """An alias names a PSI-MOD term, never a MOD id. If a release renames a
    target the alias silently stops resolving, so this test is the alarm."""
    raw = files("indicium_adapters_proteomics.data").joinpath("_ptm_aliases.json").read_text()
    aliases = {k: v for k, v in json.loads(raw).items() if not k.startswith("_")}
    unresolved = [alias for alias in aliases if psimod_curie(alias) is None]
    assert not unresolved, f"aliases with no PSI-MOD target: {unresolved}"


def test_aliases_and_labels_agree():
    """The alias layer must point at the same term the PSI-MOD label resolves to."""
    assert psimod_curie("acetylation") == psimod_curie("acetylated residue")
    assert psimod_curie("oxidation") == psimod_curie("oxidized residue")


# --- EDAM index provenance --------------------------------------------------


def _edam_index() -> dict[str, str]:
    """Read the bundled EDAM index straight from its package resource."""
    raw = files("indicium_adapters_proteomics.data").joinpath("edam_index.json").read_text()
    return json.loads(raw)["index"]


def test_edam_index_records_the_roots_it_was_cut_from():
    """The vocabulary's provenance has to be readable, or the EDAM branch this
    index was cut from becomes unanswerable after the fact. Non-emptiness is too
    weak a check: each declared root must be an EDAM topic, and must actually be
    represented in the index it supposedly generated."""
    assert _edam.size() > 0
    assert all(root.startswith(f"{_edam.EDAM_PREFIX}topic_") for root in _edam.roots())
    indexed = set(_edam_index().values())
    declared = [f"edam:{root.removeprefix(_edam.EDAM_PREFIX)}" for root in _edam.roots()]
    absent = [curie for curie in declared if curie not in indexed]
    assert not absent, f"declared roots absent from the index: {absent}"


def test_every_edam_entry_points_at_an_edam_topic():
    """The index is a topic spine. A value outside the edam: topic namespace
    would mean a concept from some other vocabulary had been folded in."""
    stray = sorted({v for v in _edam_index().values() if not v.startswith("edam:topic_")})
    assert not stray, f"index values outside the EDAM topic namespace: {stray}"


def test_a_foreign_curie_does_not_resolve_to_an_edam_topic():
    """EDAM lists the identifiers it cross-references among its own synonyms:
    GO:0006464 and MOD:00000 both sit under Protein modifications. Answering a
    slot that already names a GO or PSI-MOD term with an EDAM topic would
    rewrite an identity as a subject heading."""
    assert _edam.lookup("GO:0006464") is None
    assert _edam.lookup("MOD:00000") is None
    # The label those synonyms hang off still resolves, so the guard has not
    # simply switched EDAM lookup off.
    assert _edam.lookup("protein modifications") == "edam:topic_0601"

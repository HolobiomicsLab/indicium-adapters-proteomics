"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from indicium_adapters_proteomics import ProteomicsAdapter

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def adapter() -> ProteomicsAdapter:
    return ProteomicsAdapter()


@pytest.fixture
def mztab_path() -> Path:
    """A synthetic proteomics mzTab 1.0 file covering MTD, PRT, PEP and PSM."""
    return FIXTURES / "sample.mztab"


@pytest.fixture
def mzid_path() -> Path:
    """A synthetic mzIdentML 1.2 file with a decoy hit and a failing PSM."""
    return FIXTURES / "sample.mzid"


@pytest.fixture
def base_claim() -> dict:
    """A minimal claim whose subject is a UniProtKB accession."""
    return {
        "subject": "P02768",
        "qualifier": "carries_modification",
        "object": "phosphorylation",
        "ontology_terms": {},
        "evidence_id": "fixture-001",
    }

"""Tests for the hardcoded namespace constants.

The exhaustive-equality test is the point: a prefix cannot be added or a base
IRI edited without the change showing up here, where the authority for it has to
be named in review.
"""

from __future__ import annotations

import pytest

from indicium_adapters_proteomics import ProteomicsAdapter
from indicium_adapters_proteomics import _namespaces as ns


def test_proteomics_prefixes_complete():
    assert ns.PROTEOMICS_PREFIXES == {
        "MS": "http://purl.obolibrary.org/obo/MS_",
        "MOD": "http://purl.obolibrary.org/obo/MOD_",
        "PR": "http://purl.obolibrary.org/obo/PR_",
        "GO": "http://purl.obolibrary.org/obo/GO_",
        "OBI": "http://purl.obolibrary.org/obo/OBI_",
        "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
        "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
        "uniprot": "http://purl.uniprot.org/uniprot/",
        "PXD": "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=",
        "PRIDE": "https://www.ebi.ac.uk/pride/archive/projects/",
        "MSV": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=",
        "PASS": "http://www.peptideatlas.org/PASS/",
        "edam": "http://edamontology.org/",
    }


@pytest.mark.parametrize("prefix", ["PXD", "PRIDE", "MSV", "PASS"])
def test_repository_prefixes_present(prefix):
    """ProteomeXchange, PRIDE, MassIVE and PeptideAtlas deposits are the
    accessions a proteomics paper cites; the adapter must be able to expand
    the CURIEs the accession tier produces for them."""
    assert prefix in ns.PROTEOMICS_PREFIXES


@pytest.mark.parametrize("prefix", ["MS", "MOD"])
def test_psi_vocabularies_resolve_through_obo_purls(prefix):
    assert ns.PROTEOMICS_PREFIXES[prefix].startswith("http://purl.obolibrary.org/obo/")


def test_adapter_exposes_the_shared_prefix_map():
    """One prefix map, not two — the adapter must not carry a private copy."""
    assert ProteomicsAdapter().ontology_prefixes is ns.PROTEOMICS_PREFIXES

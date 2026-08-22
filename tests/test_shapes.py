"""End-to-end SHACL tests: validate synthetic asb:Claim graphs against the
bundled proteomics shapes.

The base IRI comes from `indicium.ASB_BASE` rather than a literal. A hardcoded
base goes silently dead after a base move: the synthetic Claim graphs would bind
a different base than the shapes, SHACL would match nothing, and every malformed
claim would wrongly "conform".
"""

from __future__ import annotations

import pyshacl
import rdflib
from indicium import ASB_BASE

from indicium_adapters_proteomics import ProteomicsAdapter
from indicium_adapters_proteomics._namespaces import PROTEOMICS_PREFIXES

ASB = ASB_BASE


def _claim_graph(predicate: str, value: str) -> rdflib.Graph:
    """Build a one-triple asb:Claim graph carrying `predicate` = `value`."""
    graph = rdflib.Graph()
    graph.bind("asb", ASB)
    subject = rdflib.URIRef(f"{ASB}claim-sample")
    graph.add((subject, rdflib.RDF.type, rdflib.URIRef(f"{ASB}Claim")))
    graph.add((subject, rdflib.URIRef(f"{ASB}subject"), rdflib.Literal("P02768")))
    graph.add((subject, rdflib.URIRef(f"{ASB}{predicate}"), rdflib.Literal(value)))
    return graph


def _conforms(predicate: str, value: str) -> bool:
    """Validate a single-term claim graph against the adapter's shapes."""
    conforms, _, _ = pyshacl.validate(
        _claim_graph(predicate, value),
        shacl_graph=ProteomicsAdapter().shacl_shapes(),
        advanced=True,
    )
    return conforms


def test_shacl_shapes_returns_non_empty_graph():
    graph = ProteomicsAdapter().shacl_shapes()
    assert isinstance(graph, rdflib.Graph)
    assert len(graph) > 0


def test_recognised_prefixes_conform():
    assert _conforms("subject_ontology_term", "uniprot:P02768")
    assert _conforms("object_ontology_term", "MOD:00696")
    assert _conforms("object_ontology_term", "edam:topic_0121")


def test_unrecognised_prefix_is_a_violation():
    assert not _conforms("subject_ontology_term", "CHEBI:17234")


def test_shape_covers_every_declared_prefix():
    """The shape's prefix alternation and PROTEOMICS_PREFIXES must stay in step,
    or a newly declared prefix would be rejected by the adapter's own shapes."""
    rejected = [
        prefix
        for prefix in PROTEOMICS_PREFIXES
        if not _conforms("subject_ontology_term", f"{prefix}:X1")
    ]
    assert not rejected, f"prefixes missing from proteomics.shapes.ttl: {rejected}"


def test_claims_without_ontology_terms_conform():
    """The ontology-term slots are optional; a claim that resolves nothing is
    not thereby invalid."""
    graph = rdflib.Graph()
    subject = rdflib.URIRef(f"{ASB}claim-bare")
    graph.add((subject, rdflib.RDF.type, rdflib.URIRef(f"{ASB}Claim")))
    graph.add((subject, rdflib.URIRef(f"{ASB}subject"), rdflib.Literal("P02768")))
    conforms, _, _ = pyshacl.validate(
        graph, shacl_graph=ProteomicsAdapter().shacl_shapes(), advanced=True
    )
    assert conforms

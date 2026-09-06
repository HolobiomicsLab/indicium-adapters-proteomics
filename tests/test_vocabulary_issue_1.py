"""Issue #1: the EDAM topics the harvest carries must be declared.

866 of 2,179 harvested proteomics rows carried only EDAM topics this adapter
never uttered, so ``domain_admission`` could name no field for them. These
tests pin the terms that closed the gap.
"""

from __future__ import annotations

import re

SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")

ADDED_SUBDOMAINS = (
    "protein-interactions",
    "protein-modifications",
    "proteins",
    "protein-folds-and-structural-domains",
    "protein-folding-stability-and-design",
    "protein-structural-motifs-and-surfaces",
    "protein-secondary-structure",
    "protein-binding-sites",
    "protein-sites-features-and-motifs",
    "protein-targeting-and-localisation",
    "immunoproteins-and-antigens",
    "membrane-and-lipoproteins",
)

ADDED_TECHNIQUES = (
    "proteomics-experiment",
    "protein-structure-analysis",
)


def test_the_measured_subdomains_are_declared(adapter):
    missing = [term for term in ADDED_SUBDOMAINS if term not in adapter.subdomains]
    assert not missing, f"issue #1 subdomains still undeclared: {missing}"


def test_the_measured_techniques_are_declared(adapter):
    missing = [term for term in ADDED_TECHNIQUES if term not in adapter.techniques]
    assert not missing, f"issue #1 techniques still undeclared: {missing}"


def test_every_declared_term_is_a_lowercase_hyphen_slug(adapter):
    bad = [t for t in (*adapter.subdomains, *adapter.techniques) if not SLUG.match(t)]
    assert not bad, f"not lowercase hyphen-separated slugs: {bad}"


def test_no_term_is_declared_twice_anywhere_in_the_adapter(adapter):
    terms = [*adapter.subdomains, *adapter.techniques]
    duplicates = sorted({t for t in terms if terms.count(t) > 1})
    assert not duplicates, f"declared more than once: {duplicates}"

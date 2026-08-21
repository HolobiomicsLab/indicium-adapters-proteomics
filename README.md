# indicium-adapters-proteomics

A basic `proteomics` `DomainAdapter` for the [indicium](https://github.com/HolobiomicsLab/indicium) ecosystem.

It owns **tier-3 concerns only** (ADR 0001): the classification vocabulary and
claim/ontology enrichment. It does **not** recognise repository accessions —
those belong to the canonical accession table in `asb-domain-core`, and a
second copy of those patterns here would be the duplicate-then-drift failure
ADR 0001 exists to prevent.

## What it contributes

| Attribute | Source |
|---|---|
| `primary_domain` | `proteomics` — without this, a proteomics paper classifies as `other` in AgenticScienceBuilder |
| `subdomains` | see the class docstring: derived from EDAM where EDAM has a hierarchy to read, curated where it does not |
| `techniques` | curated — EDAM models assays as operations and data, not as a technique hierarchy |
| `qualifiers` | curated domain claim qualifiers, all snake_case |
| `ontology_prefixes` | prefix → base IRI, for expanding the CURIEs the accession tier produces |
| `enrich_claim` | exact-match lookup against a vendored EDAM label index |

## Concept resolution is deliberately exact-match

Substring and fuzzy matching over an ontology manufacture confident-looking
wrong CURIEs, and a wrong identifier is worse than no identifier. A slot that
carries text but resolves to nothing is recorded under
`ontology_terms_unmatched` rather than passed over, so "this text is not an
EDAM concept" stays distinguishable from "no attempt was made".

## Install

```bash
uv pip install -e .
```

The `indicium.adapters` entry point makes it discoverable via
`indicium_adapters.discover_adapters()` on install — no configuration.

## Regenerating the EDAM index

`data/edam_index.json` is generated offline from `edamontology.org/EDAM.csv`
and records the ontology roots it was cut from. Regenerate on an EDAM release.

## Licence

MIT.

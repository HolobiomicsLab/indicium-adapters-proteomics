# indicium-adapters-proteomics

Proteomics-domain `DomainAdapter` for the [indicium](https://github.com/HolobiomicsLab/indicium) ecosystem.

- **`ProteomicsAdapter`** — mass-spectrometry-based proteomics claims: protein and peptide identification, quantification, PTM localisation, FDR control, interaction and localisation. Enrichment resolves protein accessions to UniProtKB, modifications to PSI-MOD, and remaining free text to EDAM.

It registers under the `indicium.adapters` entry-point group and is discoverable via `indicium_adapters.discover_adapters()`.

Also bundles two PSI-format projections into ASB claim dicts, mirroring how the
metabolomics sibling projects mzTab-M SML/SMF rows:

- `indicium_adapters_proteomics.mztab.parse_mztab_to_claims` — mzTab 1.0 MTD, PRT, PEP and PSM sections.
- `indicium_adapters_proteomics.mzid.parse_mzid_to_claims` — mzIdentML 1.1/1.2 search results.

## What it contributes

| Attribute | Source |
|---|---|
| `primary_domain` | `proteomics` — without this, a proteomics paper classifies as `other` in AgenticScienceBuilder |
| `subdomains` | curated: EDAM's proteomics topic has no sub-topic hierarchy to derive one from |
| `techniques` | curated: EDAM models assays as operations and data, not as a technique hierarchy |
| `qualifiers` | curated domain claim qualifiers, all snake_case |
| `ontology_prefixes` | prefix → base IRI, for expanding the CURIEs the accession tier produces |
| `enrich_claim` | UniProtKB accession grammar, then PSI-MOD labels, then EDAM labels |
| `shacl_shapes` | CURIE-format validation for `asb:Claim` graphs |

## Namespaces

| Prefix | Authority |
|---|---|
| `MS`, `MOD` | HUPO Proteomics Standards Initiative (PSI-MS, PSI-MOD), via OBO Foundry PURLs |
| `PR` | Protein Ontology consortium |
| `GO` | Gene Ontology Consortium |
| `OBI` | Ontology for Biomedical Investigations |
| `NCBITaxon`, `UBERON` | NCBI Taxonomy, Uberon |
| `uniprot` | UniProt Consortium |
| `PXD` | ProteomeXchange Consortium |
| `PRIDE` | EMBL-EBI PRIDE Archive |
| `MSV` | MassIVE (UCSD) |
| `PASS` | PeptideAtlas submission system (Institute for Systems Biology) |
| `edam` | EDAM (ELIXIR) |

## Concept resolution is deliberately exact-match

Substring and fuzzy matching over an ontology manufacture confident-looking
wrong CURIEs, and a wrong identifier is worse than no identifier. Three
consequences, all intentional:

- A label claimed by more than one PSI-MOD term — `Phospho` is the PSI-MS label
  of ten residue-specific terms in release 1.038.0 — resolves to **none** of
  them.
- UNIMOD and CHEMMOD identifiers in an mzTab or mzIdentML file are carried
  verbatim as the claim object and left unresolved: no offline one-to-one
  mapping onto PSI-MOD exists.
- A slot that carries text but resolves to nothing is recorded under
  `ontology_terms_unmatched` rather than passed over, so "this text is not a
  concept we know" stays distinguishable from "no attempt was made".

Protein accessions are the one non-lexical case: they are matched against
UniProtKB's own published accession grammar, so the hit is rule-based and is
justified as `semapv:CompositeMatching` rather than `semapv:LexicalMatching`.

## Accessions are not this package's job

Per ADR 0001 this adapter owns **tier-3 concerns only**: classification
vocabulary and claim/ontology enrichment. It does not recognise repository
accession patterns (PXD, PRIDE, MassIVE, PASS) — those belong to the canonical
accession table in `asb-domain-core`, and a second copy of those patterns here
would be the duplicate-then-drift failure ADR 0001 exists to prevent. What it
contributes for identifiers is `ontology_prefixes`, so a consumer can expand
the CURIEs that tier produces.

## Install

```bash
uv pip install -e .
```

The `indicium.adapters` entry point makes it discoverable on install — no
configuration.

## Usage

```python
from indicium_adapters import discover_adapters

adapter = discover_adapters()["proteomics"]
enriched = adapter.enrich_claim({
    "subject": "P02768", "qualifier": "carries_modification",
    "object": "phosphorylation", "ontology_terms": {}, "evidence_id": "abc",
})
# enriched["ontology_terms"] == {"subject": "uniprot:P02768", "object": "MOD:00696"}
```

```python
from indicium_adapters_proteomics import parse_mzid_to_claims, parse_mztab_to_claims

claims = parse_mztab_to_claims("results.mztab") + parse_mzid_to_claims("results.mzid")
```

## Regenerating the bundled indexes

`data/edam_index.json` was generated offline from `edamontology.org/EDAM.csv`
and records the ontology roots and source it was cut from. No builder for it
ships here — unlike the PSI-MOD index below, it cannot be regenerated from this
repository, so its JSON `roots` and `source` fields are the provenance of
record. `tests/test_enrichment.py` checks that every entry points at an EDAM
topic and that each declared root is represented.

`data/_psimod_index.json` is generated from the HUPO-PSI `PSI-MOD.obo` release
and records the `data-version` it was cut from:

```bash
curl -LO https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo
python -m indicium_adapters_proteomics.data.build_psimod_index \
    PSI-MOD.obo src/indicium_adapters_proteomics/data/_psimod_index.json
```

`data/_ptm_aliases.json` maps process nouns ("phosphorylation") onto PSI-MOD
term *names*, never onto MOD ids, so a release that renames a term costs recall
rather than correctness. `tests/test_enrichment.py` fails if an alias stops
resolving.

## Development

```bash
git clone https://github.com/HolobiomicsLab/indicium-adapters-proteomics
cd indicium-adapters-proteomics
pip install -e .[dev]
pytest
```

## Licence

MIT — see `LICENSE`. Bundled PSI-MOD data: CC-BY 4.0 (HUPO-PSI). Bundled EDAM
data: CC-BY-SA 4.0 (ELIXIR).

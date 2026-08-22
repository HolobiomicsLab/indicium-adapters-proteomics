# Contributing to indicium-adapters-proteomics

Thanks for your interest!

## Setup

```bash
git clone https://github.com/HolobiomicsLab/indicium-adapters-proteomics
cd indicium-adapters-proteomics
pip install -e .[dev]
```

## Development workflow

1. Create a feature branch from `main`.
2. Write tests first (TDD). The suite is `pytest`.
3. Run `ruff check .` and `pytest` before pushing.
4. Open a PR describing the change.

## Adding a term to a vocabulary

Every term in `subdomains`, `techniques`, `ontology_prefixes` and
`data/_ptm_aliases.json` must be an independently verifiable domain entity, and
the PR must name the authority it comes from — a standards body, a database, a
vendor, a published grammar. Terms scraped from whichever paper happened to be
failing are the failure mode this rule exists to prevent: they make a detector
look healthy while it silently matches nothing outside that paper. If you
cannot name the authority, leave the term out.

`tests/test_namespaces.py` asserts the prefix map exhaustively, so a new prefix
cannot land without the change being visible in review.

## Refreshing bundled data

See [README § Regenerating the bundled indexes](README.md#regenerating-the-bundled-indexes).
`data/build_psimod_index.py` runs against a PSI-MOD release downloaded
separately and is maintainer-only.

## Code of conduct

By participating you agree to abide by the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).

"""indicium-adapters-proteomics — the proteomics-domain adapter for indicium.

Exposes :class:`ProteomicsAdapter` via the ``indicium.adapters`` entry point, so
``indicium_adapters.discover_adapters()`` finds it on install and ASB's
classification vocabulary picks up ``primary_domain="proteomics"``.

Also bundles two PSI-format projections — ``parse_mztab_to_claims`` (mzTab
MTD/PRT/PEP/PSM sections) and ``parse_mzid_to_claims`` (mzIdentML search
results) — and the two identity lookups the adapter enriches with,
``uniprot_curie`` and ``psimod_curie``, re-exported here for bridge layers that
need to normalise accessions or modification names before emitting claims.
"""

from __future__ import annotations

from indicium_adapters_proteomics._enrichment import psimod_curie, uniprot_curie
from indicium_adapters_proteomics.mzid import parse_mzid_to_claims
from indicium_adapters_proteomics.mztab import parse_mztab_to_claims
from indicium_adapters_proteomics.proteomics import ProteomicsAdapter

__version__ = "0.2.1"
__all__ = [
    "ProteomicsAdapter",
    "parse_mzid_to_claims",
    "parse_mztab_to_claims",
    "psimod_curie",
    "uniprot_curie",
]

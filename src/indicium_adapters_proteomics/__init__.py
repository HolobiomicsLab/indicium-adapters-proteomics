"""Proteomics DomainAdapter for the indicium ecosystem.

Exposes :class:`ProteomicsAdapter` via the ``indicium.adapters`` entry point, so
``indicium_adapters.discover_adapters()`` finds it on install and ASB's
classification vocabulary picks up ``primary_domain="proteomics"``.
"""

from .proteomics import ProteomicsAdapter

__all__ = ["ProteomicsAdapter", "__version__"]

__version__ = "0.1.0"

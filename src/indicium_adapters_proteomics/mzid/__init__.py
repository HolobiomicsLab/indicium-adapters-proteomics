"""mzIdentML parser submodule — project search results into asb:Claim dicts.

PSI mzIdentML standard: https://hupo-psi.github.io/mzIdentML/

Emits `fdr_controlled_at`, `identifies_protein`, `identifies_peptide`,
`carries_modification` and `matched_to_spectrum` claims. Quantification is not
covered: mzIdentML does not carry abundances.
"""

from __future__ import annotations

from indicium_adapters_proteomics.mzid.parser import parse_mzid_to_claims

__all__ = ["parse_mzid_to_claims"]

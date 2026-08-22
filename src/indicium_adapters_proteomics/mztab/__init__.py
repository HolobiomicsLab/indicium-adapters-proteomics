"""mzTab parser submodule — project MTD/PRT/PEP/PSM sections into asb:Claim dicts.

PSI mzTab 1.0 standard: https://hupo-psi.github.io/mzTab/

Six claim types are produced: `fdr_controlled_at`, `identifies_protein`,
`identified_by`, `searched_against`, `identifies_peptide`,
`matched_to_spectrum`, plus `carries_modification` and `has_abundance` for rows
that carry modification or abundance columns.

mzTab-M (metabolomics, SML/SMF/SME sections) is out of scope here — that is the
`indicium-adapters-metabolomics` sibling's file, and the small-molecule sections
are silently ignored by this parser rather than mis-typed as peptides.
"""

from __future__ import annotations

from indicium_adapters_proteomics.mztab.parser import parse_mztab_to_claims

__all__ = ["parse_mztab_to_claims"]

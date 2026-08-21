"""ProteomicsAdapter — the proteomics DomainAdapter (api_version=1)."""

from __future__ import annotations

from typing import Final

from . import _edam

# Claim slots worth resolving. ``qualifier`` is included because a qualifier
# literal often carries an assay name ("quantified by label-free quantification")
# which is itself an EDAM concept.
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", "qualifier")

# Every CURIE this adapter emits comes from exact-after-normalisation matching
# against vendored EDAM labels, so this is the only justification it can honestly
# claim.
_LEXICAL: Final[str] = "semapv:LexicalMatching"


class ProteomicsAdapter:
    """Mass-spectrometry-based proteomics.

    Owns the tier-3 concerns only (ADR 0001): classification vocabulary and
    claim/ontology enrichment. It deliberately does **not** recognise repository
    accessions (PXD, PRIDE, MassIVE) — those belong to the canonical accession
    table in ``asb-domain-core``, and a second copy of those patterns here would
    be the duplicate-then-drift failure ADR 0001 exists to prevent. What this
    adapter contributes for identifiers is ``ontology_prefixes``, so a consumer
    can expand the CURIEs that tier produces.

    **Where the vocabulary comes from, and where it does not.** Concept
    resolution runs against EDAM, vendored and cut from the Proteomics topic plus
    the Proteomics-experiment, Protein-modifications and Protein-expression
    branches (``_edam.roots()``). ``subdomains`` and ``techniques``, by contrast,
    are **curated here rather than derived**: EDAM's proteomics topic has no
    sub-topics to read, so unlike the astronomy adapter — which reads its
    vocabulary out of the UAT — this one cannot claim the field's own hierarchy.
    They are an editor's list, and should be replaced by a derived vocabulary if
    PSI ever publishes a topic hierarchy suitable for the purpose.

    Every qualifier matches ``^[a-z][a-z0-9_]*$`` (snake_case), required because
    Perspicacité routes adapter qualifiers onto the open ``asb:domainQualifier``
    slot whose SHACL shape enforces that pattern.
    """

    domain_id = "proteomics"
    api_version = 1

    # ASB's classification layer reads ``primary_domain`` in preference to
    # ``domain_id`` (agentic_science_builder.classification._adapter_vocab) and
    # unions it into the controlled vocabulary.
    primary_domain = "proteomics"

    qualifiers: frozenset[str] = frozenset(
        {
            # --- Identification and search ------------------------------------
            "identified_by",  # object is a search engine or scoring method
            "searched_against",  # object is a sequence database
            "fdr_controlled_at",  # object is the false-discovery-rate threshold
            "digested_with",  # object is the protease
            "cleaved_by",  # object is an enzyme, in a biological sense
            "fragmented_by",  # object is the dissociation method (HCD, CID, ETD)
            # --- Quantification -----------------------------------------------
            "quantified_by",  # object is the quantification method or tool
            "labelled_with",  # object is TMT, iTRAQ, SILAC, dimethyl
            "normalized_by",  # object is the normalisation strategy
            "differentially_abundant_in",  # object is the condition or group
            # --- Modification and enrichment ----------------------------------
            "modified_at",  # object is the residue or site
            "carries_modification",  # object is the PTM
            "enriched_by",  # object is the enrichment chemistry or affinity step
            # --- Interaction and localisation ----------------------------------
            "co_purifies_with",  # object is an interaction partner
            "interacts_with",
            "localised_to",  # object is a compartment or tissue
            # --- Validation ----------------------------------------------------
            "validated_by",  # object is an orthogonal assay
        }
    )

    # PSI owns the mass-spectrometry vocabularies; UniProt and PR own protein
    # identity; PRIDE and MassIVE own the deposits the accession tier resolves.
    ontology_prefixes: dict[str, str] = {
        "MS": "http://purl.obolibrary.org/obo/MS_",
        "MOD": "http://purl.obolibrary.org/obo/MOD_",
        "PR": "http://purl.obolibrary.org/obo/PR_",
        "GO": "http://purl.obolibrary.org/obo/GO_",
        "OBI": "http://purl.obolibrary.org/obo/OBI_",
        "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
        "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
        "uniprot": "http://purl.uniprot.org/uniprot/",
        "PRIDE": "https://www.ebi.ac.uk/pride/archive/projects/",
        "edam": "http://edamontology.org/",
    }

    # Curated, not derived — see the class docstring.
    subdomains: tuple[str, ...] = (
        "bottom-up-proteomics",
        "top-down-proteomics",
        "targeted-proteomics",
        "phosphoproteomics",
        "glycoproteomics",
        "immunopeptidomics",
        "chemoproteomics",
        "structural-proteomics",
        "interactomics",
        "spatial-proteomics",
        "single-cell-proteomics",
        "clinical-proteomics",
        "plasma-proteomics",
        "metaproteomics",
        "proteogenomics",
        "degradomics",
        "secretomics",
        "protein-turnover",
        "thermal-proteome-profiling",
        "ptm-proteomics",
        "computational-proteomics",
    )

    techniques: tuple[str, ...] = (
        "data-dependent-acquisition",
        "data-independent-acquisition",
        "parallel-reaction-monitoring",
        "selected-reaction-monitoring",
        "label-free-quantification",
        "tandem-mass-tag",
        "itraq",
        "silac",
        "lc-msms",
        "maldi-tof",
        "ion-mobility-ms",
        "orbitrap-ms",
        "cross-linking-ms",
        "hydrogen-deuterium-exchange-ms",
        "native-ms",
        "affinity-purification-ms",
        "de-novo-peptide-sequencing",
        "spectral-library-search",
    )

    # ---- protocol methods --------------------------------------------------

    def extraction_context(self) -> str:
        return (
            "Domain: mass-spectrometry-based proteomics. Name subjects and objects "
            "with the terminology of EDAM and the PSI vocabularies where one "
            "applies (for example 'peptide identification', 'protein "
            "quantification'), using the singular preferred label rather than a "
            "paraphrase. Prefer a concrete tool, instrument or database name over "
            "a generic word: 'MaxQuant' rather than 'the search engine', "
            "'UniProtKB/Swiss-Prot' rather than 'the database'. Keep a protein and "
            "the peptide it was inferred from distinct. Domain qualifiers "
            "available in addition to the base set: "
            + ", ".join(sorted(self.qualifiers))
            + "."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Attach EDAM CURIEs to the slots whose text names an EDAM concept.

        Pure: returns a new dict and never mutates the input.

        A slot that carries text but resolves to nothing is recorded under
        ``ontology_terms_unmatched`` rather than passed over. Silence would make
        "this text is not an EDAM concept" indistinguishable from "no attempt was
        made", which is the distinction a downstream mapping set needs.
        """
        out = dict(claim)
        terms: dict[str, str] = dict(out.get("ontology_terms") or {})
        justifications: dict[str, str] = dict(out.get("ontology_term_justifications") or {})
        unmatched: dict[str, str] = dict(out.get("ontology_terms_unmatched") or {})

        for slot in _ENRICHABLE_SLOTS:
            if slot in terms:
                continue  # another adapter already resolved this slot
            value = out.get(slot)
            if not isinstance(value, str) or not value.strip():
                continue
            curie = _edam.lookup(value)
            if curie:
                terms[slot] = curie
                justifications[slot] = _LEXICAL
                unmatched.pop(slot, None)
            else:
                unmatched[slot] = value

        if terms:
            out["ontology_terms"] = terms
            out["ontology_term_justifications"] = justifications
        if unmatched:
            out["ontology_terms_unmatched"] = unmatched
        return out

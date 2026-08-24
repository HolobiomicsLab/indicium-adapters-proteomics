"""ProteomicsAdapter — the proteomics DomainAdapter (api_version=1)."""

from __future__ import annotations

from collections.abc import Callable
from importlib.resources import files
from typing import Final

import rdflib

from indicium_adapters_proteomics import _edam
from indicium_adapters_proteomics._enrichment import psimod_curie, uniprot_curie
from indicium_adapters_proteomics._namespaces import PROTEOMICS_PREFIXES

# Claim slots worth resolving. ``qualifier`` is included because a qualifier
# literal often carries an assay name ("quantified by label-free quantification")
# which is itself an EDAM concept.
_ENRICHABLE_SLOTS: Final[tuple[str, ...]] = ("subject", "object", "qualifier")

# Resolution ladder, most specific first, each paired with the SEMAPV
# justification that honestly describes it. UniProt is a published accession
# grammar, so its hits are rule-based rather than lexical; PSI-MOD and EDAM are
# exact-after-normalisation label matches.
_RESOLVERS: Final[tuple[tuple[Callable[[str], str | None], str], ...]] = (
    (uniprot_curie, "semapv:CompositeMatching"),
    (psimod_curie, "semapv:LexicalMatching"),
    (_edam.lookup, "semapv:LexicalMatching"),
)


def _slot_text(claim: dict, slot: str) -> str | None:
    """Return a claim slot's text when it is a non-empty string, else None."""
    value = claim.get(slot)
    if not isinstance(value, str) or not value.strip():
        return None
    return value


def _store(claim: dict, key: str, value: dict) -> None:
    """Set `key` on a claim when `value` is non-empty, else drop any stale entry.

    Dropping matters: a slot resolved on this pass must not leave behind an
    `ontology_terms_unmatched` record from an earlier one, or the claim
    contradicts itself.
    """
    if value:
        claim[key] = value
        return
    claim.pop(key, None)


def _resolve(text: str) -> tuple[str, str] | None:
    """Resolve free text to a (CURIE, SEMAPV justification) pair, or None.

    Walks the resolution ladder in order and returns the first hit, so a
    UniProt accession is never re-read as a modification label.
    """
    for resolver, justification in _RESOLVERS:
        curie = resolver(text)
        if curie:
            return curie, justification
    return None


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
    resolution runs against three vendored, independently published sources:
    UniProtKB's own accession grammar, the HUPO-PSI PSI-MOD release, and EDAM
    cut from the Proteomics topic plus the Proteomics-experiment,
    Protein-modifications and Protein-expression branches (``_edam.roots()``).
    ``subdomains`` and ``techniques``, by contrast, are **curated here rather
    than derived**: EDAM's proteomics topic has no sub-topics to read, so unlike
    the astronomy adapter — which reads its vocabulary out of the UAT — this one
    cannot claim the field's own hierarchy. They are an editor's list, and
    should be replaced by a derived vocabulary if PSI ever publishes a topic
    hierarchy suitable for the purpose.

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
            "identifies_protein",  # subject is the file or study, object an accession
            "identifies_peptide",  # subject is the protein, object the peptide
            "matched_to_spectrum",  # subject is the peptide, object the spectrum ref
            "identified_by",  # object is a search engine or scoring method
            "searched_against",  # object is a sequence database
            "fdr_controlled_at",  # object is the false-discovery-rate threshold
            "digested_with",  # object is the protease
            "cleaved_by",  # object is an enzyme, in a biological sense
            "fragmented_by",  # object is the dissociation method (HCD, CID, ETD)
            # --- Quantification -----------------------------------------------
            "quantified_by",  # object is the quantification method or tool
            "has_abundance",  # object is the assay or study variable measured
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

    ontology_prefixes: dict[str, str] = PROTEOMICS_PREFIXES

    # Set by the 2026-08-25 cross-domain vocabulary audit; a term two
    # sciences both declare cannot tell them apart, so it was dropped.
    subdomains: tuple[str, ...] = (
        "quantitative-proteomics",
        "interactomics",
        "phosphoproteomics",
        "glycoproteomics",
        "ptm-proteomics",
        "proteogenomics",
        "metaproteomics",
        "clinical-proteomics",
        "single-cell-proteomics",
        "spatial-proteomics",
        "immunopeptidomics",
        "chemoproteomics",
        "peptidomics",
    )

    # Set by the 2026-08-25 cross-domain vocabulary audit; a term two
    # sciences both declare cannot tell them apart, so it was dropped.
    techniques: tuple[str, ...] = (
        "shotgun-proteomics",
        "top-down-proteomics",
        "de-novo-peptide-sequencing",
        "peptide-mass-fingerprinting",
        "two-dimensional-gel-electrophoresis",
        "affinity-purification-ms",
        "cross-linking-ms",
        "native-ms",
        "hydrogen-deuterium-exchange-ms",
        "ion-mobility-ms",
        "limited-proteolysis-ms",
        "thermal-proteome-profiling",
        "proximity-labeling",
        "phosphopeptide-enrichment",
        "proximity-extension-assay",
        "protein-interaction-experiment",
        # Restored: the first audit pass deleted this while resolving the
        # gene-expression collision, but the two are different terms and no
        # other science declares it. Without it proteomics had no word for
        # the thing it measures.
        "protein-expression",
    )

    # ---- protocol methods --------------------------------------------------

    def extraction_context(self) -> str:
        """Return the prompt fragment appended to Perspicacité's extraction call."""
        return (
            "Domain: mass-spectrometry-based proteomics. Name subjects and objects "
            "with the terminology of EDAM and the PSI vocabularies where one "
            "applies (for example 'peptide identification', 'protein "
            "quantification'), using the singular preferred label rather than a "
            "paraphrase. Prefer a concrete tool, instrument or database name over "
            "a generic word: 'MaxQuant' rather than 'the search engine', "
            "'UniProtKB/Swiss-Prot' rather than 'the database'. Give a protein as "
            "its UniProtKB accession where the paper reports one, and a "
            "modification by its PSI-MOD name ('phosphorylated residue') or the "
            "process noun ('phosphorylation'). Keep a protein and the peptide it "
            "was inferred from distinct. Domain qualifiers available in addition "
            "to the base set: " + ", ".join(sorted(self.qualifiers)) + "."
        )

    def enrich_claim(self, claim: dict) -> dict:
        """Attach UniProt, PSI-MOD or EDAM CURIEs to the slots that name one.

        Pure: returns a new dict and never mutates the input. Slots another
        adapter already resolved are left alone.

        A slot that carries text but resolves to nothing is recorded under
        ``ontology_terms_unmatched`` rather than passed over. Silence would make
        "this text names no concept we know" indistinguishable from "no attempt
        was made", which is the distinction a downstream mapping set needs.
        """
        out = dict(claim)
        terms: dict[str, str] = dict(out.get("ontology_terms") or {})
        justifications: dict[str, str] = dict(out.get("ontology_term_justifications") or {})
        unmatched: dict[str, str] = dict(out.get("ontology_terms_unmatched") or {})
        for slot in _ENRICHABLE_SLOTS:
            text = _slot_text(out, slot)
            if text is None or slot in terms:
                continue
            resolved = _resolve(text)
            if resolved is None:
                unmatched[slot] = text
                continue
            terms[slot], justifications[slot] = resolved
            unmatched.pop(slot, None)
        _store(out, "ontology_terms", terms)
        _store(out, "ontology_term_justifications", justifications)
        _store(out, "ontology_terms_unmatched", unmatched)
        return out

    def shacl_shapes(self) -> rdflib.Graph:
        """Return ASB-claim SHACL shapes (CURIE format validation only)."""
        graph = rdflib.Graph()
        turtle = files("indicium_adapters_proteomics.shapes").joinpath("proteomics.shapes.ttl")
        graph.parse(turtle.open("rb"), format="turtle")
        return graph

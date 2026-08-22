"""Hardcoded namespace URI constants for the proteomics adapter.

Every prefix here is owned by a named authority, listed alongside it. Nothing is
inferred from a paper or a corpus: if an upstream namespace ever changes URI,
update it in one place — this file.
"""

from __future__ import annotations

PROTEOMICS_PREFIXES: dict[str, str] = {
    # HUPO Proteomics Standards Initiative controlled vocabularies, served as
    # OBO Foundry PURLs. MS carries instruments, acquisition methods and
    # search-engine scores; MOD carries protein modifications.
    "MS": "http://purl.obolibrary.org/obo/MS_",
    "MOD": "http://purl.obolibrary.org/obo/MOD_",
    # Protein Ontology (PRO consortium) — protein forms and complexes.
    "PR": "http://purl.obolibrary.org/obo/PR_",
    # Gene Ontology Consortium — function, process, cellular component; the
    # target of every enrichment-analysis claim a proteomics paper makes.
    "GO": "http://purl.obolibrary.org/obo/GO_",
    # Ontology for Biomedical Investigations — assay and protocol classes.
    "OBI": "http://purl.obolibrary.org/obo/OBI_",
    # NCBI Taxonomy and Uberon, for the species and tissue a sample came from.
    "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
    "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
    # UniProt Consortium — the identity namespace for a protein accession.
    "uniprot": "http://purl.uniprot.org/uniprot/",
    # ProteomeXchange Consortium — the PXD identifier a paper cites; PRIDE
    # (EMBL-EBI) and MassIVE (UCSD) are two of the deposition repositories
    # behind it. PASS is PeptideAtlas's own submission namespace (Institute for
    # Systems Biology); PeptideAtlas builds have no stable per-build PURL, so
    # only the submission namespace is represented here.
    "PXD": "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=",
    "PRIDE": "https://www.ebi.ac.uk/pride/archive/projects/",
    "MSV": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=",
    "PASS": "http://www.peptideatlas.org/PASS/",
    # EDAM (ELIXIR) — the topic/operation/data spine this adapter enriches
    # free-text claim slots against.
    "edam": "http://edamontology.org/",
}

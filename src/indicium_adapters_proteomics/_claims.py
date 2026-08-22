"""The claim-dict constructor shared by the mzTab and mzIdentML projections.

Both parsers emit the same 5-slot ASB claim shape, so the shape lives in one
place: a change to how ontology terms or evidence are attached cannot drift
between the two formats.
"""

from __future__ import annotations

from typing import Any

# Every CURIE these parsers attach comes from matching a name or an identifier
# a file already carries, so this is the only justification they can claim.
LEXICAL_MATCHING = "semapv:LexicalMatching"


def build_claim(
    subject: str,
    qualifier: str,
    obj: str,
    *,
    ontology_term: str | None = None,
    evidence: dict[str, Any] | None = None,
    evidence_id: str = "",
) -> dict[str, Any]:
    """Build one 5-slot ASB claim dict.

    `ontology_term` is attached to the object slot together with a lexical
    matching justification; `evidence` is included only when non-empty, so a
    claim never carries an empty evidence bag.
    """
    claim: dict[str, Any] = {
        "subject": subject,
        "qualifier": qualifier,
        "object": obj,
        "ontology_terms": {},
        "evidence_id": evidence_id,
    }
    if ontology_term:
        claim["ontology_terms"]["object"] = ontology_term
        claim["ontology_term_justifications"] = {"object": LEXICAL_MATCHING}
    if evidence:
        claim["evidence"] = evidence
    return claim

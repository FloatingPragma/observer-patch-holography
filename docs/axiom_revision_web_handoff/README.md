# Axiom Revision Web Handoff

Bundle for external agents and authors migrating public web surfaces to the
three-axiom basis. The bundle is self-contained: an external agent works
from these files without repository access beyond the two normative
sources named below.

## Scope

In scope: rewriting, repairing, or swapping copy on public learning and
website routes so that every surface states exactly three axioms, uses the
canonical renderings in `axiom_copy.json`, respects the claim boundaries
in `claim_contract.md`, and passes the acceptance tests in
`route_matrix.json`. The route matrix enumerates the known high-risk
routes; routes discovered during the work are added to the route change
manifest with a proposed rewrite class.

Out of scope: any change to theory content, papers, ledgers, registries,
Lean sources, simulators, or claim statuses. Copy questions that the
bundle does not answer are returned as questions; they are never resolved
by invention.

## Authority order

1. `claims/axiom_registry.yaml` (machine-readable normative source)
2. `docs/AXIOM_REFERENCE.md` (reader-facing reference)
3. This bundle (`docs/axiom_revision_web_handoff/` and
   `docs/AXIOM_REVISION_LEARNING_HANDOFF.md`)

Where two levels disagree, the higher level wins and the disagreement is
reported in the return package. The bundle files carry sha256 hashes in
`handoff_manifest.json`; an agent verifies hashes before use.

## Permission boundary

The external agent has no license to invent, extend, weaken, or
reclassify theory content. Formal axiom text is consumed verbatim from
`axiom_copy.json` and never paraphrased. Status labels move only by dated
artifact in the source repository, never by web edit. Staging previews
are permitted; production deployment requires owner approval and is
outside the agent's authority. User progress data on learning routes is
preserved across the migration.

## Required return artifacts

1. Route change manifest: one row per touched route with rewrite class,
   copy assets consumed, diff summary, and any routes added beyond the
   route matrix.
2. Preview URL for the complete staged migration.
3. Test reports: stale-token scan output, forbidden-claim scan output,
   and formal-copy-verbatim check output, each with pass/fail per route
   and the exact scan configuration used.

## Bundle contents

| File | Role |
|---|---|
| `README.md` | Scope, authority order, permission boundary, return artifacts. |
| `claim_contract.md` | Permitted statements, required premises, forbidden promotions, status-class wording rules. |
| `axiom_copy.md` | Human-readable canonical renderings of the three axioms. |
| `axiom_copy.json` | Machine-readable mirror of `axiom_copy.md`; the file sites consume. |
| `glossary_migration.csv` | Term migration and glossary table with search aliases. |
| `route_matrix.json` | High-risk routes with rewrite classes, copy assets, claim dependencies, and acceptance tests. |
| `EXTERNAL_AGENT_PROMPT.md` | Complete brief for the implementing coding agent. |
| `handoff_manifest.json` | Bundle manifest with sha256 hashes and build notes. |
| `../AXIOM_REVISION_LEARNING_HANDOFF.md` | Source packet for learning-material authors: narrative, cards, toy models, countermodels, dependency map, lessons, FAQ, web copy blocks, migration table. |

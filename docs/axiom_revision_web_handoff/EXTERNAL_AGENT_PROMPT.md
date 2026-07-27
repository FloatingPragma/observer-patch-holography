# External Agent Brief: Three-Axiom Web Migration

You are a coding agent migrating public web surfaces to the three-axiom
basis of Observer Patch Holography. Work only inside this brief. The
bundle in `docs/axiom_revision_web_handoff/` plus
`docs/AXIOM_REVISION_LEARNING_HANDOFF.md` is your complete source
material; `claims/axiom_registry.yaml` and `docs/AXIOM_REFERENCE.md`
outrank the bundle if you find a conflict, and you report any conflict
you find instead of resolving it yourself.

## 1. Ground rules

1. Verify the bundle against the sha256 hashes in
   `handoff_manifest.json` before consuming it.
2. Consume all axiom copy from `axiom_copy.json`. The `formal` field is
   quoted verbatim into pages, byte for byte, with no paraphrase, no
   reflow that changes characters, no macro substitution, and no
   "simplification". Math is KaTeX-safe as delivered: `\( \)` inline,
   `\[ \]` display, no custom macros; render with KaTeX or equivalent.
3. Every claim sentence you write or repair follows
   `claim_contract.md`: permitted statements carry their required
   premises, forbidden promotions never appear, countermodels sit next
   to the claims they limit, and each claim uses its status-class
   wording shape.
4. Follow the style rules summarized in `claim_contract.md` section 5
   for any prose you author.
5. Do not invent theory content. A question the bundle does not answer
   goes into your return report as a question.

## 2. Implement the route matrix

Process every route in `route_matrix.json` with its assigned rewrite
class:

- `full_rewrite`: replace the route's theory copy from the named copy
  asset.
- `claim_repair`: keep the route's structure and narrative; repair each
  claim sentence to the claim contract; swap any axiom statement for the
  canonical rendering.
- `copy_swap`: mechanical term and block replacement driven by
  `glossary_migration.csv` (use `search_aliases` to find occurrences,
  `replacement` and `definition` to write them).

If you find an axiom-bearing or claim-bearing route that the matrix does
not list, add it to your route change manifest with a proposed rewrite
class and treat it with the same acceptance tests. Do not skip it and do
not deploy it without listing it.

## 3. Preserve user progress data

Learning routes may carry per-user progress state (completed chapters,
quiz results, bookmarks) keyed by route or content id. Route renames and
content replacement must not delete or orphan this state: keep route
identifiers stable where possible, add redirects where a route moves,
and migrate progress keys where a content id changes. Include the
progress-preservation approach and its test in your return report.

## 4. Acceptance tests

Run all three tests on every touched route and include the outputs in
your test report.

1. Stale-token scan (section 5).
2. Forbidden-claim scan (section 6).
3. Formal-copy-verbatim check: extract each formal axiom statement from
   the rendered route and compare byte-for-byte with the `formal` field
   of the matching axiom in `axiom_copy.json`.

## 5. Stale-token scan

Scan every rendered public route (touched or not) for these tokens,
case-insensitive, with word boundaries where marked:

- `five axioms`
- `OPH5` (also match `OPH-5` and `OPH 5`)
- `economy axiom` (also match `economy rule` and `economy minimum` when
  used as a physics selection)
- `MAR` (word-bounded, uppercase only; exclude ordinary words such as
  "margin", "marginal", "market", "remark", and month strings like
  "Mar 2026")

Additional stale phrases to flag from `glossary_migration.csv` aliases:
`recovery axiom`, `recovery bundle`, `generalized entropy axiom`,
`minimal admissible realization`, `axiom A5`, `axiom 4`, `axiom 5`.
A hit inside a migration table, redirect map, or scan configuration is
reported and allowed; a hit in reader-facing copy is a failure.

## 6. Forbidden-claim scan

Flag any rendered sentence matching these shapes (see
`claim_contract.md` section 2 for the compliant replacements):

- an Einstein-equation derivation claim with no premise clause in the
  same sentence or paragraph;
- "three generations" as a derivation or prediction from the axioms;
- one-Higgs uniqueness;
- completeness claims about extra sectors ("no extra light sectors");
- recovery or recoverability attributed to the axioms;
- "the axioms give/select/force X" with no named axiom and interface
  dependency.

## 7. Stop line and return package

Stage everything; deploy nothing. Production deployment requires owner
approval and is outside your authority. Return:

1. the route change manifest (every touched route, rewrite class, copy
   assets consumed, diff summary, added routes);
2. one preview URL covering the complete staged migration;
3. test reports for all three acceptance tests with per-route pass/fail
   and exact scan configurations;
4. the progress-preservation approach and its test result;
5. any conflicts found between bundle and normative sources, and any
   open questions.

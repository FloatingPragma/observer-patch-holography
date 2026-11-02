# OPH Writing Style Guide

Binding rules for all OPH prose: papers, READMEs, docs, book, and blog source
material. Code comments follow the same two constraints: state the exact
research status and avoid machine-generated prose mannerisms.

## State-only language

Material states the exact current state with no reference to past or future
states of the research.

- Banned: "now", "already", "previously", "no longer", "recently", "used to",
  "currently", "presently", "latest", "still", "not yet", "going forward",
  "in the future", "will be added", "an earlier version", "has been updated",
  and "new" when it describes a research artifact or result.
- In papers, also banned: "work in progress", "stays open", "remains open",
  "remains partial", and equivalent live-project phrasing. State the
  scientific boundary directly: "the construction supplies no physical clock"
  or "no continuum theorem is proved here."
- Scientific provenance may identify dataset editions, publication dates,
  frozen-protocol timestamps, artifact versions and hashes. Maintained guides
  and research prose omit editorial dates, update stamps and audit chronology.
  Papers give the theorem, premises, evidence and nonclaims directly.
- Keep `/docs` small: update the owning page, consolidate overlapping material,
  and place implementation-specific mathematics beside its code. Do not add
  individual run reports, duplicate status summaries or migration-only pages.
  The root holds only high-level programs, references and ledgers. No dated
  file, progress report, audit result, integration record or PR handoff is
  added anywhere under `docs/`; that history belongs in the GitHub pull
  request or issue. A `docs/research/` note states a current result its code
  verifies, without history. Release papers do not cite `docs/` files, and
  code, receipts and registers do not bind them as evidence. The
  [docs index](../README.md#rules-for-this-directory) states the full rules.
  Preserve substantive proofs, counterexamples, evidence and source commitments
  when consolidating documents.
- Scientific classifications such as conditional, unproved, refuted, or
  outside the theorem's hypotheses are allowed when they classify a claim
  rather than a work item.

## AI giveaways

Remove or refactor on sight:

- Em-dashes. Use commas, colons, parentheses, or separate sentences.
- "not X, but Y" and "not only X, but also Y" sentence shapes.
- Short punchy intro sentences that tee up a paragraph ("The result?",
  "Here is the catch.", "This matters.").
- Stock intensifiers and connectives: "crucially", "importantly", "notably",
  "moreover", "furthermore", "in essence", "essentially", "arguably",
  "delve", "robust", "comprehensive", "seamless", "landscape" (figurative),
  "tapestry", "journey", "unpack", "It's worth noting", "In conclusion".
- Rule-of-three flourishes ("fast, simple, and powerful").
- Anthropomorphized documents ("this paper aims to", "the section seeks to").
  The paper states, proves, reports.
- Bullet lists where prose carries the argument better. Bullets are for
  genuinely enumerable items.

## Banned wording

- "honest", "honestly", "honesty". The word never appears in any surface.
  Honesty is the default assumption and needs no label.
- Other moral-character adjectives used to advertise a surface or an accounting
  convention. Good faith is the default assumption and needs no label.

## Voice

Prose matches Bernhard's own register: varied sentence lengths, dry, direct,
first person used sparingly. Reference sample:
https://muellerberndt.medium.com/building-a-secure-nft-gaming-experience-a-herdsmans-diary-1-91aab11139dc

## Standalone papers

Each paper is a publication that an outside reader finishes without access to
the repository, the ledgers, or the other papers.

- A paper never reads as a list of project-internal labels. Internal tier and
  rule names (QFT-Q1 through QFT-Q4, MGNS-1, Q0, FJ, DAG as a tier graph)
  carry no meaning outside the project, so the prose states the thing itself:
  "exact finite quantization", "formal perturbative quantization", "the modular
  algebra-state reconstruction data".
- Premise-register identifiers, observation-ledger identifiers, GitHub issue
  numbers, issue-lane codes, project phases, and live status summaries never
  appear in a paper. Examples include PR-65, OL-C6, issue B19, and ``lane
  \#743.'' Put that information in the owning ledger, plan, or issue.
- This prohibition covers every rendered part of a paper, including titles,
  abstracts, body text, footnotes, captions, tables, appendices, and source or
  provenance notes. A stable scientific citation or artifact path may support
  reproducibility. It never carries task ownership, open/closed state,
  completion history, or a project update.
- A paper records the scientific result and its boundary. State the theorem,
  hypotheses, evidence, counterexample, conditionality, or nonclaim directly.
  Work queues, roadmaps, milestone narration, and descriptions of what a
  revision closed belong in the workspace plan or the owning issue.
- Where a paper depends on a companion result, it states the result in words
  and cites the companion paper. It does not hand the reader a label.
- After any paper-prose edit, run `python3 tools/check_reader_style.py` from
  the repository root. When changing this policy or its checker, also run
  `python3 -m pytest -q tools/test_gates_actually_fail.py`. Do not add an
  exception merely to retain project-control prose in a paper.

## Abstracts and informal surfaces

- Abstracts stay short and informal, and are usually left untouched. An
  abstract that has grown past roughly 200 words is a candidate for tightening,
  not for another clause.
- Abstracts, informal descriptions, and the book NEVER carry code references
  or internal identifiers (D12, sigma_ref, CL-3, GAP-A5, DK-01, QFT-Q2, MAR,
  MGNS-1, Q0 and similar). Those live in technical sections, docs, and ledgers.
- Abstracts carry the main results of the paper without calling them
  achievements.
- An abstract reads on its own. A reader who has seen none of the other papers
  follows every sentence in it.

## Comparing research programs

Make OPH's case through its concrete constructions: bounded self-reading
patches, protected overlap records, and the derived structures connected to
them. Compare named formulations using primary sources and the same criteria:
explicit inputs, strength of derivation, controlled approximation, physical
identification and prospective tests. Credit established results and shared
ideas. Numerical claims of superiority require a defined metric and evidence;
theorem counts, author reputation and publicity do not measure physical truth.
Keep detailed comparisons in overview discussions. Specialist papers retain
their scientific focus, and READMEs give the short positive case.

## Repository READMEs

- A README is an informal front door for interested non-specialists. It gives
  the physical idea, the supported result, and the main remaining step in
  plain language. Proof inventories and complete premise lists belong in the
  linked papers and ledgers.
- A highlighted result uses a short heading and one compact paragraph. As a
  working ceiling, keep it below 100 words and link at most three primary
  destinations. Several code files or theorem modules are grouped behind one
  paper, proof, receipt, or ledger link.
- Headings name familiar physics directly: "the four laws of thermodynamics",
  "the Standard Model gauge group", and "three-dimensional space". They do
  not use internal status classes such as "conditional finite theorem
  package".
- Technical boundaries are summarized in one sentence about the main
  remaining work. The README does not reproduce a referee report, issue
  checklist, theorem dependency chain, mutation inventory, or every caveat.
- Accuracy remains binding. A conditional result is called conditional, a
  postdiction is not called a prediction, and an abstract mathematical object
  is not called physical before its attachment is established. The wording
  presents the open work as the research route rather than as a pile of
  negative status labels.

## Acronyms

- ALL-CAPS acronyms are used as sparingly as the sentence allows. Where the
  spelled-out words fit, they are used instead.
- An acronym that survives is spelled out at its first use in every paper
  independently, including acronyms that are standard in the field.
- Project-internal acronyms and named rules (CFQ, KMS in the OPH-specific
  sense) never appear in an abstract, in the book, or in any informal surface.
  In technical sections they are defined before first use.
- Proper names of outside tools and results keep their usual form
  (Isabelle/HOL, Lean-QIT, Mathlib, CHSH after one spelled-out use).

## Book

- The book is pop-science and standalone. It explains the physics to a reader
  who has never opened a paper.
- The book NEVER carries internal identifiers, tier labels, rule names, or
  repository paths.
- Quantities the book itself defines and explains (a named capacity, a named
  constant) stay, because the book gives the reader their meaning on the page.
- A prose paragraph in the book stays at or below 150 words. A longer argument
  is divided where the idea changes, so the reader gets somewhere to breathe.
- A technical stretch begins with the physical question or a concrete image.
  Equations are introduced in words, and the prose after them says what the
  reader should carry forward.
- The basis stays visible through the narrative: bounded observers make
  records and recover the specification of their world, agreement fixes what
  counts as shared, and since nothing unrecovered may be assumed, maximum
  randomness fills what the recovered constraints leave free. This does not
  imply that preference or consciousness chooses physical outcomes.

## Basis and premises

- A technical surface that states the basis distinguishes the closure
  hypothesis, operational agreement, the executable specification and the
  adopted state-selection and response-endogeneity principles, and cites
  `docs/BASIS.md`. Internal realization motivates the latter principles;
  closure alone does not entail an entropy minimizer or inner action on the
  connected response algebra. A surface that
  states a theorem premise A1, A2 or A3 pairs it with both forms in the same
  local context: one plain-language sentence, then the concise formal
  statement, then what the premise constrains, then what it does not imply,
  then a citation to `docs/AXIOM_REFERENCE.md`. A surface that only consumes a
  premise names the exact dependency and cites the reference. No surface
  counts A1-A3 as axioms of OPH or calls them its basis.
- The symbol \(A_5\) names the alternating group on five letters, never an
  axiom. Where confusion is possible, the group is typeset with a subscript
  and introduced as a group.
- Retired principles are never mentioned by name in active prose. A former
  consequence appears only through its surviving premises: a conditional
  window, a declared completion, a classified interface, or a countermodel.

## Consistency

- Values, status labels, and claim boundaries agree across every paper and
  public surface. The ledger is the source of truth; papers cite it.
- One name per object. A renamed object is renamed everywhere in the same
  release.

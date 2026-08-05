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
- Allowed: "is work in progress" for unclosed derivations and similar open
  lanes. "Is open" for open obligations. Dated artifact names (a certificate
  carries its date in its filename; that is provenance, not narrative).
- Ledgers and audits identify artifact versions, dates, hashes, and exact
  classifications without narrating the evolution of the research. Papers
  state the theorem, premises, and status label that apply.
- Claim-status idiom such as "stays on record as a display packet only" is a
  classification of an artifact, and stays legal.

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
- Where a paper depends on a companion result, it states the result in words
  and cites the companion paper. It does not hand the reader a label.

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
- The three-axiom intuition remains visible through the narrative: bounded
  observers make records, agreement constrains shared meaning, and maximum
  randomness fills what those constraints leave free. This does not imply that
  preference or consciousness chooses physical outcomes.

## Axiom basis

- A surface that states or counts the axiom basis gives exactly three axioms
  and pairs every axiom with both forms in the same local context: one
  plain-language sentence, then the concise formal statement, then what the
  axiom constrains, then what it does not imply, then a citation to
  `docs/AXIOM_REFERENCE.md`. A surface that only consumes an axiom names the
  exact dependency and cites the reference.
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

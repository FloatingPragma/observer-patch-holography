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
  rule names (QFT-Q1 through QFT-Q4, MAR, MGNS-1, Q0, FJ, DAG as a tier graph)
  carry no meaning outside the project, so the prose states the thing itself:
  "exact finite quantization", "formal perturbative quantization", "the economy
  minimum", "the modular algebra-state reconstruction data".
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

## Acronyms

- ALL-CAPS acronyms are used as sparingly as the sentence allows. Where the
  spelled-out words fit, they are used instead.
- An acronym that survives is spelled out at its first use in every paper
  independently, including acronyms that are standard in the field.
- Project-internal acronyms and named rules (MAR, CFQ, KMS in the OPH-specific
  sense) never appear in an abstract, in the book, or in any informal surface.
  In technical sections they are defined before first use.
- Proper names of outside tools and results keep their usual form
  (Isabelle/HOL, Lean-QIT, Mathlib, CHSH after one spelled-out use).

## Book

- The book is pop-science and standalone. It explains the physics to a reader
  who has never opened a paper.
- The book NEVER carries internal identifiers, tier labels, rule names, or
  repository paths. Concepts appear under plain descriptions: "the economy
  rule", not an acronym.
- Quantities the book itself defines and explains (a named capacity, a named
  constant) stay, because the book gives the reader their meaning on the page.

## Consistency

- Values, status labels, and claim boundaries agree across every paper and
  public surface. The ledger is the source of truth; papers cite it.
- One name per object. A renamed object is renamed everywhere in the same
  release.

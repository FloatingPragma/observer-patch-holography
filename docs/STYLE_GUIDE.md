# OPH Writing Style Guide

Binding rules for all OPH prose: papers, READMEs, docs, book, blog source
material. Code comments follow the spirit. The voice is Bernhard's: plain,
direct, concrete, technically confident, dry humor where it fits, zero fluff.

## State-only language

Material states the exact current state with no reference to past or future
states of the research.

- Banned: "now", "already", "previously", "no longer", "recently", "used to",
  "going forward", "in the future", "will be added", "an earlier version",
  "has been updated", "new" (as in "the new certificate").
- Allowed: "is work in progress" for unclosed derivations and similar open
  lanes. "Is open" for open obligations. Dated artifact names (a certificate
  carries its date in its filename; that is provenance, not narrative).
- History and withdrawals belong in ledger and audit surfaces (closure ledger,
  proof spine, falsification program, issue threads), never in paper prose.
  A paper states the current theorem, the current premises, the current
  status label. Nothing else.
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

- moral-character adjectives used to advertise a surface or an accounting
  convention. Good faith is the default assumption and needs no label.

## Abstracts and informal surfaces

- Abstracts stay short and informal, and are usually left untouched.
- Abstracts, informal descriptions, and the book NEVER carry code references
  or internal identifiers (D12, sigma_ref, CL-3, GAP-A5, DK-01 and similar).
  Those live in technical sections, docs, and ledgers.
- Abstracts carry the main results of the paper without calling them
  achievements.
- Acronyms and named rules (MAR, CFQ, KMS in the OPH-specific sense) are
  defined before first use, in every paper independently.

## Consistency

- Values, status labels, and claim boundaries agree across every paper and
  public surface. The ledger is the source of truth; papers cite it.
- One name per object. A renamed object is renamed everywhere in the same
  release.

## Bernhard's voice, in brief

Reference: the public posts at muellerberndt.medium.com. Sentences are short
and declarative. Claims come with the number and the artifact. Jokes are dry
and land in passing, without a setup. The reader is treated as smart and
busy. When something is unproven, the text says exactly that in one clause
and moves on.

# Solved Physics Problems

This folder holds standalone markdown writeups for OPH applications to named
physics problems. They act as supplemental papers for public reading and OPH
Sage ingestion. They do not enter the TeX/PDF rendering or publication
pipeline.

Each writeup follows the same rule: an OPH explanation must name the
observer-like self-reading structure. The relevant object is a bounded physical
or software patch with local state, visible boundaries or ports, readback,
records, repair moves, and public evidence receipts. Normal forms classify what
survives repair; selection, yield, power, or material performance requires a
branch-specific source law, action, repair ledger, or receipt bundle.

Use GitHub-safe math delimiters in this folder: inline math as `$...$` and
display math as fenced `math` blocks. Avoid `\(...\)` and `\[...\]`, because
GitHub Markdown can treat equation lines such as `=` as headings before math
rendering.

## Writeups

- [Fractional Quantum Hall States](fractional_quantum_hall.md): fractional Hall
  phases as edge/holonomy normal forms, Abelian $K$ matrix recovery, and the
  $5/2$ selector no-go theorem.
- [High-Temperature Superconductivity](high_temperature_superconductivity.md):
  high $T_c$ as charge $2e$ repair plus $U(1)_Q$ phase confluence, with an
  inverse-design protocol and no recipe claim.
- [Plasma Fusion and Confinement](plasma_fusion.md): fusion confinement as a
  quotient-visible boundary-repair ledger, with H-mode, ELM, Lawson, DD, and
  net-power claims kept in separate receipt tiers.

## Scope Rule

A writeup can close a physics problem at the OPH theorem-package level. Device,
material, and experiment claims need their own public evidence bundles. A failed
applied branch falsifies that branch without automatically falsifying the compact
recovered-core paper.

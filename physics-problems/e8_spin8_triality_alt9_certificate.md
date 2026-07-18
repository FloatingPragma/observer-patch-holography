# Finite E8/Spin8 Triality Claim Statement For An Alt(9) Double Cover

## Motivating Result

This certificate note starts from a finite mathematical construction. Griess
and Lam obtain a $2\!\cdot\!\mathrm{Alt}(9)$ subgroup in exceptional-lattice
and moonshine gluing through triality for groups of type $D_4$
([Griess and Lam, "A moonshine path from E8 to the
Monster"](https://arxiv.org/abs/0910.2057)).
Their construction uses finite $\mathrm{Spin}^{+}(8,3)$, $E_8$-lattice, and
Leech-lattice data. The narrower question here is whether a related real
$\mathrm{Spin}(8)$ vector/half-spin construction can be recorded as an exact,
public matrix certificate without being promoted into physical $E_8$
speculation.

**Status:** claim statement; exact verification pending. The
repository contains only the claim statement under
`code/e8_triality_claim_statement/`, not the promised Sage programs, matrices, lattice
bases, check output, or stable hash bundle. The statements below are therefore
verification targets, not a proved OPH theorem.

Date: 2026-07-08
Audit revision: 2026-07-11

## Problem, Standard Mathematics, and OPH Contribution

**Standard mathematics.** The Weyl group of an $A_8$ root system is
$\mathrm{Sym}(9)$, with even subgroup $\mathrm{Alt}(9)$. The deleted
permutation representation gives an orthogonal eight-dimensional model. Its
orientation-preserving subgroup can be pulled back through
$\mathrm{Spin}(8)\to\mathrm{SO}(8)$, and triality permutes the vector and two
half-spin representations. A positive-definite even unimodular lattice of rank
eight is isometric to the $E_8$ lattice.

**Unresolved certificate target.** The missing work is not the naming of those
standard facts. It is an exact common construction that verifies the chosen
$A_8\subset E_8$ embedding, the nonsplit spin pullback, the invariant
half-spin lattice, both mod-2 orbit partitions, and a compatible triality map in
one basis-controlled bundle.

**OPH contribution.** OPH contributes the evidence discipline: finite visible
state, exact records, quotient invariants, repairable checks, and public
receipts. It does not make this group-theoretic calculation easier or unique;
the calculation belongs to standard Clifford, lattice, and computational group
theory. Its OPH role is as a local representation-closure sidecar only.

## Certificate Specification

**Target 1** (root subsystem and vector presentation). Supply an explicit
$A_8$ simple-root chain inside a fixed $E_8$ Gram basis and verify its Cartan
matrix. Its Weyl group must act as $\mathrm{Sym}(9)$, with
$\mathrm{Alt}(9)$ acting through the deleted permutation representation

```math
V=\left\{x\in\mathbb R^9:\sum_{i=1}^{9}x_i=0\right\}.
```

**Target 2** (nonsplit double cover). In this declared representation,
$(12)(34)$ has a two-dimensional $-1$-eigenspace. The standard spin-lift
criterion then gives square $-1$ for either lift. Together with the Schur
multiplier of $\mathrm{Alt}(9)$, this identifies the pullback as the nonsplit
Schur double cover $2\!\cdot\!\mathrm{Alt}(9)$, rather than
$2\times\mathrm{Alt}(9)$. The certificate must check this using exact
Clifford generators rather than treating the group name as an input label.

**Target 3** (half-spin $E_8$ lattice). Under the positive half-spin
representation $\Delta^+$, exhibit a rank-eight lattice basis $B_+$ and
verify, for every declared generator $g$,

```math
[g]_{B_+}\in\mathrm{GL}(8,\mathbb Z),
\qquad
[g]_{B_+}^{\mathsf T}G_{E_8}[g]_{B_+}=G_{E_8}.
```

The Gram matrix must be positive definite, integral, even, and unimodular. Only
after those checks may the lattice be identified with $E_8$ and the image
placed in $\mathrm{Aut}(E_8)=W(E_8)$.

**Target 4** (mod-2 fingerprints). Reduce both integral presentations in the
same declared quadratic space $E_8/2E_8\setminus\{0\}$ and verify the proposed
orbit partitions

```math
\mathrm{Alt}(9)_{\rm vec}:\{9,36,84,126\},
\qquad
(2\!\cdot\!\mathrm{Alt}(9))_{\Delta^+}:\{120,135\}.
```

Both lists sum to $255$, as required. Since central $-I$ becomes the identity
modulo two, the effective mod-2 action of the spin copy is again an
$\mathrm{Alt}(9)$ action. If exact orbit enumeration produces these different
partitions, the two images cannot be conjugate in $O_8^+(2)$.

**Target 5** (triality compatibility). Exhibit an outer triality automorphism
of $\mathrm{Spin}(8)$ that carries the declared vector representation to the
positive-half-spin representation, together with compatible lattice bases and
mod-2 reductions. Triality relates or fuses the presentations under an outer
automorphism; it does not make the nonconjugate $O_8^+(2)$ images conjugate.

## Evidence Labels and Verification Obligations

The current evidence labels are:

- **standard background:** $W(A_8)\cong\mathrm{Sym}(9)$, the deleted
  permutation representation, the spin-lift order criterion, rank-eight even
  unimodular uniqueness, and $D_4$ triality; the final receipt must cite or
  prove each imported fact it uses;
- **model-derived target:** the chosen exact matrices, invariant half-spin
  lattice, orbit representatives and stabilizers, and integral triality
  intertwiner;
- **public evidence:** pending until the raw bundle is committed and its checks
  can be rerun from a clean environment;
- **physical evidence:** none, and none is implied by this algebraic sidecar.

For public verifier status, a future certificate bundle must contain the Sage source,
exact matrix data, lattice bases, mod-2 orbit computation, machine-readable
check receipts, and stable hashes named in the claim statement. The checks must derive
the group orders, relations, Gram invariance, orbit exhaustion, and triality
intertwining from the matrices; a manifest that merely repeats the expected
answers is not a certificate.

## Claim Boundary

Passing this bundle would establish a finite algebraic representation
certificate. It would support exceptional representation-closure and triality
bookkeeping. It would not prove OPH, derive the Standard Model quotient, prove a
physical $E_8$ realization, close the heterotic critical-edge CFT, or count as
a hardware receipt.

## Remaining Extension

The proposed certificate concerns the nonsplit
$2\!\cdot\!\mathrm{Alt}(9)$ subgroup and its triality-related vector and
positive-half-spin presentations. Extending it to a double cover of
$\mathrm{Sym}(9)$ requires orientation-reversing permutation generators. Such
generators lift naturally through a Pin group, exchange the two half-spin
modules, and require an explicit normalization convention. The claimed
$\sqrt2$ normalization and its relation to the Griess-Lam construction remain
separate verification obligations rather than consequences of the present
specification.

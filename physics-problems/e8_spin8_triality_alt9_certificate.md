# Finite E8/Spin8 Triality Certificate For An Alt(9) Double Cover

## Motivating Result

This certificate note starts from a finite mathematical surprise. The group
`2.Alt(9)` appears naturally in exceptional lattice and moonshine bookkeeping
through triality and gluing data
([Hohn, "A moonshine path from E8 to the Monster"](https://arxiv.org/abs/0910.2057)).
The OPH question is whether the same kind of finite, public, representation
certificate can be kept separate from physical `E8` speculation.

**Status:** solved as a finite certificate statement and OPH claim-boundary
note. Public verifier promotion remains gated on the raw Sage scripts, exact
matrix data, check outputs, and stable hashes under `code/e8_triality/`.

Date: 2026-07-08

## Role In The OPH Stack

OPH uses finite observer-visible data, records, repair maps, quotient invariants, and public receipts as primitive evidence carriers. The present certificate belongs to that finite exact discipline on the exceptional-symmetry side:
``` math
\mathrm{Spin}(8)\ \text{triality} + E_8\ \text{lattice preservation}
\quad\Rightarrow\quad
\text{finite exceptional representation certificate}.
```
It supports the $E_8$-type representation-closure lane used by compact-gauge and heterotic-local bookkeeping. It does not replace the recovered-core theorem stack, the Minimal Admissible Realization (MAR) selection of the Standard Model quotient, the critical-edge CFT gate, or any physical carrier evidence rule.

## Certificate Statement

**Theorem 1** (Finite $E_8/\mathrm{Spin}(8)$ triality certificate). *The certificate target is an exact finite matrix construction with the following data.*

1.  *An $A_8$ root subsystem inside the $E_8$ root lattice. Its Weyl group is $W(A_8)\cong\mathrm{Sym}(9)$, and the even subgroup is $\mathrm{Alt}(9)$.*

2.  *The permutation involution $(12)(34)\in\mathrm{Alt}(9)$ lifts to $\mathrm{Spin}(8)$ with square $-1$. Hence the preimage of $\mathrm{Alt}(9)$ is the nonsplit Schur double cover $2\!\cdot\!\mathrm{Alt}(9)$, not $2\times\mathrm{Alt}(9)$.*

3.  *Under the positive half-spin representation $\Delta^+$, the lifted group preserves an even unimodular determinant-one lattice, hence an $E_8$ lattice. The spin-side copy therefore lands in $\mathrm{Aut}(E_8)=W(E_8)$.*

4.  *The vector and positive-half-spin mod-2 orbit fingerprints on $E_8/2E_8\setminus\{0\}$ are different:
    ``` math
    \mathrm{Alt}(9)_{\mathrm{vec}}:\{9,36,84,126\},
    \qquad
    (2\!\cdot\!\mathrm{Alt}(9))_{\Delta^+}:\{120,135\}.
    ```
    Thus the two copies are not conjugate in $O_8^+(2)$.*

5.  *The outer triality automorphism of $\mathrm{Spin}(8)$ permutes the vector and two half-spin eight-dimensional representations and identifies these otherwise nonconjugate presentations.*

*Remark 1* (Notation). Here $A_8$ denotes the root subsystem and $\mathrm{Alt}(9)$ denotes the alternating group. The shorthand $A_9$ is avoided for the group because $A_9$ also denotes a root system.

*Remark 2* (Claim boundary). This is a finite algebraic certificate. It supports exceptional representation-closure and triality bookkeeping. It does not prove OPH, derive the Standard Model quotient, prove physical $E_8$ realization, close the heterotic edge CFT, or count as a hardware receipt.

*Remark 3* (Public receipt gate). For public verifier status, the repository bundle must include the Sage source, exact matrix data, lattice bases, mod-2 orbit computation, stdout or machine-readable check receipts, and stable hashes under `code/e8_triality/`. Until that bundle is populated, this note records the certificate statement and its OPH claim placement rather than a standalone public reproduction bundle.

## Remaining Extension

This certificate concerns the nonsplit $2\!\cdot\!\mathrm{Alt}(9)$ subgroup and its triality-fused vector and positive-half-spin presentations. The full Griess-Lam $2.\mathrm{Sym}(9)$ double-cover construction requires adjoining odd Clifford lifts. Those odd lifts exchange the two half-spin modules and carry the associated $\sqrt 2$-normalization. That is a separate certificate gate.

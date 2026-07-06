---
title: "Arithmetic Chain from mod-6 Splitting to Eisenstein, Coxeter-Todd, Leech, and Monster Structures"
status: "background-note"
claim_level: "expository"
oph_dependency: false
oph_relevance:
  - finite residue classes
  - Eisenstein/A2 geometry
  - coding-theoretic glue
  - exceptional lattice background
  - sporadic group structure
not_claimed:
  - OPH theorem
  - RH proof
  - derivation of Leech lattice from mod-6 alone
  - derivation of Standard Model Z6 quotient
  - derivation of Monster from arithmetic alone
---

# An Expository Note on the Arithmetic Chain

$$(\mathbb{Z}/6\mathbb{Z})^\times \to \mathbb{Z}[\omega] \to K_{12} \to \Lambda_{24} \to \mathbb{M}$$

*Background note / expository map, May 2026*

---

### Status

This is a **background contribution / expository map**, not an OPH theorem and not a derivation of the Leech lattice or the Monster from mod-6 arithmetic. The note records a chain of standard mathematical structures that may be useful as background for OPH work on finite symmetry, edge codes, exceptional lattices, and implementation-facing regulator charts.

### Abstract

We present an expository account of a five-step chain of standard mathematical structures, starting from the multiplicative group $`(\mathbb{Z}/6\mathbb{Z})^\times = \{1, 5\}`$ and reaching the Monster sporadic simple group $`\mathbb{M}`$:

$$(\mathbb{Z}/6\mathbb{Z})^\times \to \mathbb{Z}[\omega] \to K_{12} \to \Lambda_{24} \to \mathbb{M}.$$

Each bridge requires substantial independent mathematical data: algebraic number theory for Bridge 1, the ternary Golay code for Bridge 2, the binary Golay code for Bridge 3, and vertex operator algebra theory for Bridge 4. All results cited are standard, due to Coxeter, Todd, Conway, Sloane, Griess, Frenkel-Lepowsky-Meurman, Borcherds, and others. **No new theorems are claimed**, and **the chain is ancestral, not constructive**: $`(\mathbb{Z}/6\mathbb{Z})^\times`$ does not generate the later structures; it indexes a single arithmetic dichotomy that initiates the chain. For OPH, the relevance is as background: the note records how finite residue classes, Eisenstein/$`A_2`$ geometry, code glue, exceptional lattices, and sporadic symmetry form a standard mathematical ladder that OPH can reference when studying finite boundary codes and exceptional symmetry patterns.

---

## 1. Introduction and Scope

The group of units modulo 6 is $`(\mathbb{Z}/6\mathbb{Z})^\times = \{1, 5\}`$, a cyclic group of order 2. This elementary object turns out to index the first link in a chain of genuine mathematical structures of increasing depth.

**What this note does:** describe each bridge precisely, with the correct mathematical content for each construction.

**What this note does not do:** prove the Riemann Hypothesis, derive the Standard Model $`\mathbb{Z}_6`$ quotient, claim that mod-6 arithmetic generates the Leech lattice or the Monster, or assert any novel theorem.

---

## 2. Bridge 1: From $`(\mathbb{Z}/6\mathbb{Z})^\times`$ to the Eisenstein Integers

**Definition 2.1.** Let $`\omega = e^{2\pi i/3} = -\frac{1}{2} + \frac{\sqrt{3}}{2}i`$. The Eisenstein integers are

$$\mathbb{Z}[\omega] = \{a + b\omega : a, b \in \mathbb{Z}\} \subset \mathbb{C}.$$

As a real lattice, $`\mathbb{Z}[\omega]`$ is the $`A_2`$ root lattice, a hexagonal lattice. The norm of $`a + b\omega`$ is

$$N(a + b\omega) = a^2 - ab + b^2.$$

**Theorem 2.2** (Splitting of rational primes in $`\mathbb{Z}[\omega]`$). Let $`p`$ be a rational prime with $`p \geq 5`$. Then:

1. If $`p \equiv 1 \pmod{6}`$, then $`p = \pi \bar{\pi}`$ in $`\mathbb{Z}[\omega]`$ for a non-unit Eisenstein prime $`\pi`$. We say $`p`$ **splits**.
2. If $`p \equiv 5 \pmod{6}`$, then $`p`$ remains prime in $`\mathbb{Z}[\omega]`$. We say $`p`$ is **inert**.

*Proof.* Standard, via Dedekind–Kummer applied to $`\mathbb{Q}(\omega) = \mathbb{Q}(\sqrt{-3})`$, whose discriminant is $`-3`$. A prime $`p \neq 3`$ splits iff $`-3`$ is a square mod $`p`$, which by quadratic reciprocity occurs exactly when $`p \equiv 1 \pmod{3}`$. For $`p \geq 5`$, this is equivalent to $`p \equiv 1 \pmod{6}`$. $`\blacksquare`$

**Remark 2.3** (Verification).

For example, $`7 \equiv 1 \pmod{6}`$ splits as $`7 = (3 + \omega)(3 + \bar{\omega})`$, while $`11 \equiv 5 \pmod{6}`$ remains inert. Likewise, $`13`$ and $`19`$ split, while $`17`$ and $`23`$ remain inert.

**Connection to $`(\mathbb{Z}/6\mathbb{Z})^\times`$.** The splitting/inertness dichotomy is governed by the residue of $`p`$ modulo 3, since the quadratic character of $`\mathbb{Q}(\sqrt{-3})`$ is the Kronecker symbol $`\chi_{-3}`$, which has conductor 3. The lift to modulus 6 is via the natural isomorphism

$$(\mathbb{Z}/6\mathbb{Z})^\times \cong (\mathbb{Z}/2\mathbb{Z})^\times \times (\mathbb{Z}/3\mathbb{Z})^\times \cong (\mathbb{Z}/3\mathbb{Z})^\times \cong \mathbb{Z}/2\mathbb{Z},$$

where the mod-2 component is trivial. So $`(\mathbb{Z}/6\mathbb{Z})^\times \cong C_2`$ **indexes** the split/inert dichotomy for rational primes $`p \notin \{2, 3\}`$; it does not generate or determine the structure of $`\mathbb{Z}[\omega]`$, which exists independently.

This is the precise sense in which $`(\mathbb{Z}/6\mathbb{Z})^\times`$ enters the chain: as an indexing label on an arithmetic dichotomy, not as a generating algebraic object.

---

## 3. Bridge 2: From $`\mathbb{Z}[\omega]`$ to the Coxeter–Todd Lattice $`K_{12}`$

This bridge requires more than a single congruence condition. The Coxeter–Todd lattice $`K_{12}`$ is constructed from $`\mathbb{Z}[\omega]^6`$ using the **ternary Golay code** as glue.

**Definition 3.1** (Ternary Golay code). The ternary Golay code $`\mathcal{C}_6 = [6, 3, 4]_3`$ is a 3-dimensional subspace of $`\mathbb{F}_3^6`$ with minimum Hamming distance 4. It is the unique (up to equivalence) self-dual $`[6, 3]`$ code over $`\mathbb{F}_3`$ with these parameters. Its 27 codewords each lie at Hamming distance at least 4 from every other.

**Definition 3.2** (Coxeter–Todd lattice via the ternary Golay code). Let $`\theta = 1 - \omega`$, so that $`N(\theta) = 3`$ and $`\mathbb{Z}[\omega]/(\theta) \cong \mathbb{F}_3`$. The Coxeter–Todd lattice is the set of $`(z_1, \ldots, z_6) \in \mathbb{Z}[\omega]^6`$ such that the residue vector

$$(z_1 \bmod \theta, \ldots, z_6 \bmod \theta) \in \mathcal{C}_6,$$

where the reduction $`\mathbb{Z}[\omega]^6 \to \mathbb{F}_3^6`$ is applied componentwise. The constraint requires the residue vector mod $`\theta`$ to be a codeword of the ternary Golay code. The zero-sum condition is insufficient by itself.

The index of $`K_{12}`$ in $`\mathbb{Z}[\omega]^6`$ is $`|\mathbb{F}_3^6 / \mathcal{C}_6| = 3^{6-3} = 27`$.

**Proposition 3.3** (Properties of $`K_{12}`$; Coxeter–Todd 1953, Conway–Sloane 1983). The lattice $`K_{12}`$ has:

- **Dimension**: 12 (real).
- **Minimal vector norm**: 4 (with standard scaling).
- **Kissing number**: 756.
- **Discriminant**: $`3^6 = 729`$.
- **Even integral**, with no norm-2 vectors.

It is the densest known 12-dimensional lattice packing.

**Remark 3.4.** A naive sublattice such as $`\{(z_1, \ldots, z_6) : \sum z_i \equiv 0 \pmod{\theta}\}`$ has discriminant $`3^7 = 2187`$, not $`3^6 = 729`$, and so is *not* the Coxeter–Todd lattice. The full ternary Golay constraint is essential.

**Remark 3.5.** The norm $`N(\theta) = 3`$ is the same "3" that appears in the conductor of $`\chi_{-3}`$ and in the ramification of 3 in $`\mathbb{Z}[\omega]`$: the prime 3 ramifies as $`3 = -\omega^2 \theta^2`$.

---

## 4. Bridge 3: From $`K_{12}`$ to the Leech Lattice $`\Lambda_{24}`$

The third bridge also requires Golay code data, this time the **binary** Golay code.

**Definition 4.1** (Extended binary Golay code). The extended binary Golay code $`G_{24}`$ is the unique (up to isomorphism) binary linear code with parameters $`[24, 12, 8]`$: length 24, dimension 12, minimum Hamming distance 8.

**Theorem 4.2** (Leech lattice; standard construction). The Leech lattice $`\Lambda_{24}`$ consists of the vectors $`\frac{1}{\sqrt{8}}(a_1, \ldots, a_{24})`$ with $`a_i \in \mathbb{Z}`$ satisfying:

1. All $`a_i`$ are congruent modulo 2 (i.e., they share a common parity, denoted $`m \in \{0, 1\}`$).
2. The reduction $`(a_i \bmod 4)/2 \in \mathbb{F}_2`$ for each $`i`$ produces a vector in $`\mathbb{F}_2^{24}`$ that is a codeword of $`G_{24}`$.
3. The sum satisfies $`\sum_i a_i \equiv 4m \pmod{8}`$.

This is the standard construction from Conway–Sloane [1], Chapter 4. Several equivalent formulations exist; see Remark 4.3.

**Remark 4.3** (Why naive Construction A is insufficient). The naive Construction A applied to $`G_{24}`$, namely $`\{x \in \mathbb{Z}^{24} : x \bmod 2 \in G_{24}\}`$ scaled by $`\frac{1}{\sqrt{8}}`$, gives a related but distinct 24-dimensional even unimodular lattice that *does* contain norm-2 vectors. The Leech lattice $`\Lambda_{24}`$ is *rootless* (no norm-2 vectors), which requires the additional mod-4 and parity conditions above. Equivalent constructions exist via:

- Construction $`A_4`$ applied to a lift of $`G_{24}`$ to $`\mathbb{Z}/4\mathbb{Z}`$, with a parity check appended (Conway–Sloane).
- Construction B applied to $`G_{24}`$ followed by point-doubling (Forney's cubing construction).
- Gluing from any Niemeier lattice via the appropriate modifying procedure.
- The Weyl vector construction in the Lorentzian lattice $`\mathrm{II}_{25,1}`$.

All of these produce the same lattice.

**Properties of $`\Lambda_{24}`$:**

- **Dimension**: 24.
- **Minimal vector norm**: 4 (rootless: no vectors of norm 2).
- **Kissing number**: $`196{,}560`$.
- **Determinant**: 1 (unimodular, even).

It is the unique even unimodular lattice in $`\mathbb{R}^{24}`$ with no roots (Niemeier 1973).

**Connection back to $`K_{12}`$.** Conway–Sloane (1983) show that $`K_{12}`$ is the sublattice of $`\Lambda_{24}`$ fixed by a specific order-3 automorphism, realizing $`\Lambda_{24}`$ as a 12-dimensional $`\mathbb{Z}[\omega]`$-lattice that doubles $`K_{12}`$ in a controlled way. This embedding is the precise sense in which $`K_{12}`$ "leads to" $`\Lambda_{24}`$; the full Leech lattice requires the binary Golay code data on top.

---

## 5. Bridge 4: From $`\Lambda_{24}`$ to the Monster Group $`\mathbb{M}`$

This bridge differs from the previous three. The Monster arises as the symmetry group of a vertex operator algebra constructed from the Leech lattice, rather than as a lattice or code.

### 5.1 The Conway groups

**Theorem 5.1** (Conway 1968). The automorphism group of $`\Lambda_{24}`$ is the Conway group $`\mathrm{Co}_0 = \mathrm{Aut}(\Lambda_{24})`$, with order

$$|\mathrm{Co}_0| = 2^{22} \cdot 3^9 \cdot 5^4 \cdot 7^2 \cdot 11 \cdot 13 \cdot 23.$$

Its quotient $`\mathrm{Co}_1 = \mathrm{Co}_0 / \{\pm I\}`$ is one of the 26 sporadic simple groups, of order approximately $`4.16 \times 10^{18}`$. Two further sporadic simples $`\mathrm{Co}_2`$, $`\mathrm{Co}_3`$ arise as stabilizers in $`\mathrm{Co}_1`$ of vectors of specific minimal norms.

### 5.2 The Moonshine module and the Monster

The Monster does *not* act directly on $`\Lambda_{24}`$. It arises as follows.

**Construction 5.2** (Frenkel–Lepowsky–Meurman 1988). Starting from the Leech lattice vertex operator algebra $`V_\Lambda`$ (a bosonic conformal field theory of central charge 24), the **Moonshine module** is defined as a $`\mathbb{Z}/2`$-orbifold:

$$V^\natural = V_\Lambda^+ \oplus V_\Lambda^{T,+},$$

where $`V_\Lambda^+`$ is the involution-invariant subspace of $`V_\Lambda`$, and $`V_\Lambda^{T,+}`$ is the involution-invariant subspace of the unique twisted module $`V_\Lambda^T`$ for the lattice involution $`x \mapsto -x`$.

**Theorem 5.3** (Frenkel–Lepowsky–Meurman 1988; Borcherds 1992). The automorphism group of $`V^\natural`$ is the Monster sporadic simple group:

$$\mathrm{Aut}(V^\natural) = \mathbb{M}, \qquad |\mathbb{M}| \approx 8.08 \times 10^{53}.$$

The graded character of $`V^\natural`$ is $`j(\tau) - 744`$, and for each $`g \in \mathbb{M}`$, the McKay–Thompson series $`T_g(\tau) = \mathrm{Tr}_{V^\natural}(g \cdot q^{L_0 - 1})`$ is a Hauptmodul for a genus-zero subgroup of $`\mathrm{SL}_2(\mathbb{R})`$. This is **monstrous moonshine** (Conway–Norton 1979, proved by Borcherds 1992).

**Remark 5.4** ($`\mathrm{Co}_1`$ inside $`\mathbb{M}`$). The Conway group $`\mathrm{Co}_1`$ from Theorem 5.1 appears inside the Monster as a quotient of a maximal subgroup: $`2^{1+24}.\mathrm{Co}_1 \leq \mathbb{M}`$, the centralizer of a $`2A`$-involution in $`\mathbb{M}`$. This makes precise the sense in which the Leech-lattice symmetries are a "piece of" the Monster's structure, without being its full source. The Monster is fundamentally larger and richer than $`\mathrm{Co}_0`$.

**Remark 5.5** (Other orbifolds, multiple constructions). Carnahan (2017) proved that for each fixed-point-free automorphism of $`\Lambda_{24}`$ satisfying a no-massless-states condition, the corresponding cyclic orbifold of $`V_\Lambda`$ is isomorphic to $`V^\natural`$. This yields 51 distinct orbifold constructions of the Moonshine module, in bijection with the 51 non-Fricke conjugacy classes of $`\mathbb{M}`$.

---

## 6. The Complete Chain

| Step | Object | Bridge | Tool required |
|---:|---|---|---|
| 0 | $`(\mathbb{Z}/6\mathbb{Z})^\times`$ | (starting point) | none |
| 1 | $`\mathbb{Z}[\omega]`$ | $`B_1`$ | Algebraic number theory |
| 2 | $`K_{12}`$ | $`B_2`$ | Ternary Golay code over $`\mathbb{F}_3 \cong \mathbb{Z}[\omega]/(\theta)`$ |
| 3 | $`\Lambda_{24}`$ | $`B_3`$ | Binary Golay code with mod-4 / parity conditions |
| 4 | $`\mathbb{M}`$ | $`B_4`$ | Vertex operator algebra orbifold ($`V^\natural`$ from $`V_\Lambda`$) |

Each bridge requires a genuinely new ingredient that is not contained in the previous step:

1. **Bridge 1** (Number theory): The Kronecker symbol $`\chi_{-3}`$ and the Dedekind–Kummer splitting theorem.
2. **Bridge 2** (Coding theory + lattice theory): The ternary Golay code $`[6,3,4]_3`$ as algebraic glue inside $`\mathbb{Z}[\omega]^6`$.
3. **Bridge 3** (Coding theory + lattice theory): The binary Golay code $`[24,12,8]`$ together with mod-4 and parity conditions to ensure rootlessness.
4. **Bridge 4** (Vertex operator algebras + group theory): The orbifold construction of $`V^\natural`$, requiring the FLM construction (1988), Borcherds' proof of monstrous moonshine (1992), and Griess's original construction of $`\mathbb{M}`$ (1982).

---

## 7. What This Note Does and Does Not Establish

1. **The chain is real, standard mathematics.** Every step is a known theorem, properly attributed. The structural ancestry from mod-6 arithmetic to the Monster sporadic group is a genuine feature of the mathematical landscape.

2. **The Riemann Hypothesis is not proved.** The connection between the splitting behavior of primes modulo 6 and the locations of zeros of $`\zeta(s)`$ is *not* what this chain captures. The Weil explicit formula relates prime distributions to zero locations, but splitting behavior and zero locations are distinct phenomena. Nothing in this chain bears on RH.

3. **The chain is ancestral, not generative.** $`(\mathbb{Z}/6\mathbb{Z})^\times`$ does not build, generate, or determine the later structures. It indexes a single arithmetic dichotomy. Each subsequent structure requires substantial independent mathematical data: the Eisenstein integers exist independently of mod-6 arithmetic, the Coxeter–Todd lattice requires the ternary Golay code as additional input, the Leech lattice requires the binary Golay code with mod-4 conditions, and the Monster requires the entirety of vertex operator algebra theory.

4. **The Monster's relationship to $`\Lambda_{24}`$ is genuinely indirect.** Unlike Bridges 1-3, where each output is constructed directly from the input plus auxiliary code data, Bridge 4 produces a *VOA from the lattice*. The Monster is the symmetry group of that VOA, not of $`\Lambda_{24}`$ itself. The Monster does not act on the Leech lattice; $`\mathrm{Aut}(\Lambda_{24}) = \mathrm{Co}_0`$, and $`\mathrm{Co}_1`$ appears inside $`\mathbb{M}`$ only as a maximal subgroup of an involution centralizer.

5. **No connection to physics is claimed.** Some physics models use a quotient of the form

$$G_{\mathrm{phys}} = \mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)/\mathbb{Z}_6,$$

   which is a quotient of a Lie group by a central cyclic subgroup. It is *not* the same object as the multiplicative group $`(\mathbb{Z}/6\mathbb{Z})^\times`$ studied here, which is a finite abelian group on its own. Any apparent rhyme between these two appearances of the symbol "$`\mathbb{Z}_6`$" is not made precise in this note and should not be inferred to be such.

6. **No new theorems are claimed.** Every result cited is standard and decades old.

---

## 8. Conclusion

The arithmetic of $`(\mathbb{Z}/6\mathbb{Z})^\times = \{1, 5\}`$ enters the foundational ancestry of five major mathematical structures: the Eisenstein integers (via the splitting dichotomy), the Coxeter–Todd lattice (via the ternary Golay code), the Leech lattice (via the binary Golay code), the Conway groups (as lattice automorphisms), and the Monster sporadic simple group (via the Moonshine module).

The chain crosses four fields: number theory, lattice theory, coding theory, and vertex operator algebra theory. Each crossing is mathematically substantial. The connections are real and well-known to specialists, appearing in standard references including Conway-Sloane [1] and Frenkel-Lepowsky-Meurman [6].

The value of understanding this chain is the precise mathematical map it provides: a tiny arithmetic object, the units modulo 6, indexes the first link of a sequence that culminates, five non-trivial steps later, in the largest sporadic simple group. The role of $`(\mathbb{Z}/6\mathbb{Z})^\times`$ throughout is that of an indexing label, not a generator.

---

## References

[1] Conway, J. H. and Sloane, N. J. A. *Sphere Packings, Lattices and Groups*, 3rd ed. Springer-Verlag, 1999.

[2] Coxeter, H. S. M. and Todd, J. A. "An extreme duodenary form." *Canadian Journal of Mathematics*, **5** (1953), 384–392.

[3] Conway, J. H. "A characterisation of Leech's lattice." *Inventiones Mathematicae*, **7** (1969), 137–142.

[4] Conway, J. H. and Sloane, N. J. A. "The Coxeter–Todd lattice, the Mitchell group and related sphere packings." *Math. Proc. Cambridge Philos. Soc.*, **93** (1983), 421–440.

[5] Davenport, H. *Multiplicative Number Theory*, 3rd ed. Springer-Verlag, 2000.

[6] Frenkel, I., Lepowsky, J., and Meurman, A. *Vertex Operator Algebras and the Monster*. Academic Press, 1988.

[7] Borcherds, R. E. "Monstrous moonshine and monstrous Lie superalgebras." *Inventiones Mathematicae*, **109** (1992), 405–444.

[8] Conway, J. H. "A perfect group of order 8,315,553,613,086,720,000 and the sporadic simple groups." *PNAS*, **61** (1968), 398–400.

[9] Conway, J. H. and Norton, S. P. "Monstrous moonshine." *Bull. London Math. Soc.*, **11** (1979), 308–339.

[10] Griess, R. L. "The Friendly Giant." *Inventiones Mathematicae*, **69** (1982), 1–102.

[11] Neukirch, J. *Algebraic Number Theory*. Springer-Verlag, 1999.

[12] Carnahan, S. "51 constructions of the Moonshine module." arXiv:1707.02954 (2017).

[13] Duncan, J. F. R. and Mack-Crane, S. "The Moonshine module for Conway's group." *Forum of Mathematics, Sigma*, **3** (2015), e10.

[14] Niemeier, H.-V. "Definite quadratische Formen der Dimension 24 und Diskriminante 1." *J. Number Theory*, **5** (1973), 142–178.

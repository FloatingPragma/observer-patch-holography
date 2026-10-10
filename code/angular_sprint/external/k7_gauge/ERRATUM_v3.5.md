# Erratum — GIFT / K₇ framework v3.5

**Issued** 2026-08-22. **Applies to** the v3.5 release
(Zenodo concept `10.5281/zenodo.16891489`, version record `21296168`,
deposited 2026-07-10) and to the copies of those sources distributed in this
repository under `publications/papers/`.

**This erratum does not trigger a new deposit.** The v3.5 record stands as
published; the corrections below are recorded here and will be integrated
into **v3.6**. Nothing in `publications/papers/` has been regenerated — the
sources are left as deposited, so that this erratum can be read against them,
with two exceptions, both limited to links and identifiers: the live Lean links
repointed to `Arithmon/K7-Lean` on 2026-07-17 (one line each in the main text
and supplements S1, S2, markdown only), and the three [E] citations corrected
on 2026-10-01 (section 4).

---

## 1. The closed-form Calabi–Yau residual bound does not bound what it says

**Where.** `publications/papers/markdown/k7_framework_3_5_main.md` lines 131
and 760; `…/k7_framework_3_5_S1_foundations.md` line 528; and the
corresponding passages in `publications/papers/tex/k7_framework_3_5_main.tex`
and `…_S1.tex`.

**What is written.** That the Ricci-flat Kähler (Calabi–Yau) metric on the
Z₂³-equivariant K3 fibre admits an explicit closed-form approximation whose
order-3 residual (667 parameters) satisfies an *interval-rigorous,
assumption-free, machine-checked* bound

    Var(log R) ≤ ε₃ = 1309/10⁷ ≈ 1.309 × 10⁻⁴ < 10⁻³

on a frozen 4000-point Krawczyk-certified witness, formalised in Lean 4
(`K3ClosedFormWitness`).

**What is wrong.** The 667-parameter witness was fitted **and evaluated**
with the contraction `V† H V`, which is not the holomorphic pullback of the
ambient (1,1) form to the K3; the correct pullback is `Vᵀ H V̄`. That metric
convention was retracted on **2026-07-13**, and the defect was measured
**directly on this witness** on **2026-07-19**, over the same 4000 points:

| quantity | value |
|---|---:|
| `Var(log R)` over all 4000 points, code convention — *what ε₃ bounds* | 1.309 × 10⁻⁴ |
| `Var(log R)` on the 3478 points positive-definite in **both** conventions, code | 1.445 × 10⁻⁴ |
| `Var(log R)` on those same points, **correct pullback** | **5.315 × 10⁻¹** |
| ratio true / code | **3677** |
| `log det_code − log det_true` | mean 0.534, std 0.730 |
| positive-definite census | 4000 (code) vs **3478** (correct pullback) |

So ε₃ bounds an **auxiliary quantity**, not the Calabi–Yau residual of a
Kähler metric on the K3. And in the correct convention the ansatz is **not**
close to Ricci-flat: 5.3 × 10⁻¹ is not a small residual, and 522 of the 4000
points do not even carry a positive-definite form.

**What still stands.** The arithmetic and the surface-level certification are
untouched:

* the Lean theorems are `native_decide` statements about integers and about
  the aggregate of the 4000 serialised rational endpoints — they remain valid;
* the Krawczyk–Rump certification that each box contains an **exact** common
  zero of the defining equations is a statement about the **surface**,
  independent of any metric convention;
* the *method* — the frame-invariant determinant identity, the rank-one
  caching, forward interval arithmetic, and the Lean re-computation of the
  variance envelope from rational endpoints — transfers unchanged.

The public Lean modules carry the matching correction
(`GIFT/Foundations/K3ClosedFormWitness.lean` and
`…/K3ClosedFormBoxEnclosures.lean` in the `K7-Lean` repository, 2026-08-22).

**What v3.6 must do.** Restate the passage as a bound on the auxiliary
quantity, or re-establish the certificate on the current metric datum (the
degree-3 witness v2, `holomorphic_pullback_VT`, gauge det M = 1, frozen
2026-07-16). The honest target for a re-run on the current datum is of order
10⁻³ – 10⁻², not 1.3 × 10⁻⁴.

---

## 2. The datum-level instantiation of [E] is no longer realisable as posed

**Where.** The passages of `k7_framework_3_5_main.md` (§ around line 131 and
the two-layer claim boundary of §9) that discharge the datum-level analytic
existence scheme at the normalised datum 𝒟₀ of the companion paper [E],
with `R₀(𝒟₀) ≤ 4.9 × 10³`.

**What is wrong.** On **2026-08-19** a Stokes-type obstruction (internally
`P3-R1`, ratified) established that the **global single-root** datum 𝒟₀ of
[E] is not realisable as posed.

**What still stands.** The **conditional theorems** of [E] are unaffected —
they are conditional on the two-slot external structure pack, and that
conditionality is exactly what the obstruction touches at the level of
*instantiation*, not of proof. Corollary B remains valid **as a computation**;
what fails is that its datum is no longer available in the single-root global
form. A multi-root salvage is in progress.

**What v3.6 must do.** Present the datum-level layer as conditional on a
datum that is *not currently instantiated*, rather than as discharged at 𝒟₀.

---

## 3. Citation hygiene noted for v3.6 (not an error of substance)

The v3.5 sources cite the **version** DOI `10.5281/zenodo.19893371` of the
companion spectral-geometry paper [B] 31 times, and never its **concept** DOI
`10.5281/zenodo.18920367`. The concept DOI is the one that should be cited:
it resolves to the latest version. To be corrected in v3.6.

---

## 4. Paper [E] is cited by the DOI of its dataset, not of the preprint (added 2026-10-02)

**Where.** Three places in the v3.5 main text: the "Established at the
datum" item of the claim boundary, the "Companion analytic paper [E]"
paragraph, and the bibliography entry **[E]**, in
`k7_framework_3_5_main.md` and in `k7_framework_3_5_main.tex`, and therefore
in the deposited PDF.

**What is written.** "concept DOI 10.5281/zenodo.21209413" for the paper [E].

**What is wrong.** The two Zenodo records of [E] were swapped in our
citations. According to the Zenodo API (`resource_type`, files):

| record | concept DOI | version DOI | type | content |
|---|---|---|---|---|
| **[E], the preprint** | **10.5281/zenodo.21215956** | 10.5281/zenodo.21215957 | Preprint | main paper (68 pp) + supplement (25 pp), PDF |
| [E], its dataset | 10.5281/zenodo.21209412 | 10.5281/zenodo.21209413 | Dataset | data + verification scripts, tar.gz |

So the number printed in v3.5 resolves to the **dataset**, and it is a
*version* DOI, not a concept DOI. The paper should be cited as
**10.5281/zenodo.21215956**. The swap was found on 2026-10-01 by
cross-checking our DOI registry against the author's ORCID record.

**Already corrected.** On 2026-10-01 (commit `366f69e`), the repository
surfaces (README, CITATION, STRUCTURE, wiki, `publications/papers/README.md`)
and the three passages of the distributed copies
`publications/papers/markdown/k7_framework_3_5_main.md` and
`publications/papers/tex/k7_framework_3_5_main.tex` were changed to cite
21215956. Apart from the Lean links of 2026-07-17 (see the header), these
are the only changes to the copies since the deposit. The deposited PDF and
the Zenodo record of v3.5 are immutable and keep 21209413.

**What v3.6 must do.** Cite [E] by its preprint concept DOI 21215956. Where
the dataset is meant, cite it explicitly as such (concept DOI 21209412).

---

## 5. The fine-structure relation is labelled at the wrong scale (added 2026-10-02)

**Where.** Supplement S2, §20b (relation #18), whose formula is written
α⁻¹(M_Z) = (248 + 8)/2 + 99/11 + (65/32)(1/61), in
`k7_framework_3_5_S2_derivations.md` and `k7_framework_3_5_S2.tex`, and the
row α⁻¹(M_Z) of the observable reference (wiki, corrected on 2026-10-02).

**What is wrong.** The value is compared with 137.036, which is the inverse
fine-structure constant **at low energy**, in the Thomson limit, α⁻¹(0). The
running coupling at M_Z is a different quantity (about 127.93 in the
five-flavour MS-bar scheme, PDG). The main text is already consistent with
the low-energy reading: §5.4 calls 137.033 "the topological low-energy value
α_em⁻¹" and keeps it distinct from the M_Z values of the RGE chain. The label
in S2 contradicts it.

**What still stands.** The relation itself, 267489/1952 = 137.033299…, and
its classification (STRUCTURAL: it uses the imposed normalization
det(g) = 65/32).

**How to read the agreement.** Against α⁻¹(0) = 137.035999177(21) (CODATA
2022), the deviation is −19.7 ppm (0.002%). That is a numerical proximity,
not an agreement within experimental uncertainty: the measurement is about
10⁵ times more precise than the gap. The internal observable record also
carried an unsourced uncertainty of ±0.0001. It now cites CODATA 2022.

**What v3.6 must do.** Write the relation as α⁻¹(0), and state the
comparison as above.

---

## Provenance

Sections 1–3 were surfaced by an internal audit of the corpus on
2026-08-22 and verified against the primary measurements before being written
here. The measurement underlying §1 is dated 2026-07-19 and predates this
erratum by five weeks; the delay is a process failure of ours, recorded as
such, and it is what prompted the audit.

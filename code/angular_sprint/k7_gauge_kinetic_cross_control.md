# K7 / OPH gauge kinetic cross-control

**Frozen baseline after rebase:** OPH `0660c94573c8f79c1b86955fbd1d7d927a0e4fa2`; K7 `210480b2ba8a13c98e6a6f146fbfd965e4904c00`; K7-Lean `2553ed170d04db8777f84fb2f9ef2b66b386f7fb`. All four frozen OPH scientific source blobs remain unchanged from preregistration.

## Result

`NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES`.

The frozen OPH matter-trace statistic is reconstructed from `k=(10/3,2,2)` and the declared `(nG,nH)=(3,1)` beta column `b=(41/6,-19/6,-7)`. Direct cofactor expansion gives `det(x,k,b)=(-23/3)x1+37x2-(218/9)x3`, with integer zero locus `69x1-333x2+218x3=0`. The vector role is `(alpha_Y^-1, alpha_2^-1, alpha_3^-1)` in the declared hypercharge normalization. OPH's physical kinetic-form selector remains open.

| K7 candidate | Provenance typing | OPH common-object gate |
|---|---|---|
| Type-I alpha relation plus weak angle and strong coupling | `267489/1952` is explicitly `alpha_em^-1(0)` in the pinned erratum; `sin²(theta_W)=3/13` and `alpha_s=sqrt(2)/12` have MZ comparison/interpretation in the paper. | Fails: mixed scale; scheme/threshold bridge unknown. The constructed vector is a negative control only. |
| B-test / holonomy ray | GUT-normalized ratio `14:7:2`; exact normalization gives GUT ray `(42,21,6)sqrt(2)` and coherent Y ray `(70,21,6)sqrt(2)`. The pinned paper conflicts internally: the theorem display labels `91 sqrt(2)` at `M_Z`, while its proof note and holonomy-sequence discussion call it the GUT/M_GUT exact scale. | Fails: scale conflict, MSSM branch, conditional boundary role, and unresolved scheme/threshold contract differ from frozen OPH SM one-loop plane. |
| Type-III RGE | The pinned paper publishes a common-`M_Z` numerical output row: `sin²θ_W=0.2377`, `alpha_em^-1=131.19`, and `alpha_s=0.1224` for split spectrum (`0.1038` for all-MSSM). It describes two-loop MSSM running from `alpha_GUT^-1=25.3`, `sin²θ_W=3/13` at `M_GUT`, with `tan(beta)=2`, MSSM above 3165 GeV and SM plus gauginos below. The exact vector and reproducer are not established by the four frozen source surfaces consumed here; experimental-input status is UNKNOWN. | Fails: rounded numerical output, two-loop MSSM/split-spectrum branch, and incomplete source-clean run provenance. No exact determinant verdict is computed. |

No candidate is established by the frozen source surfaces consumed here as closing the observable, common scale, compatible scheme, U(1) convention, field content, beta branch, running direction/boundary role, and source-clean derivation checks. The Type-III table gives a numerical common-scale output candidate, but it is conditional and rounded; it does not close the exact common-object gate. This verdict is scoped to the vendored K7 erratum, main paper and K7-Lean gauge-sector file; it is not a claim that no other K7 repository artifact could supply a bridge.

## Exact controls

The forbidden mixed-scale assembly uses `A0=267489/1952` and computes

```text
xY = (10 A0/13, 3 A0/13, 6 sqrt(2))
D  = -82654101/25376 + 1308 sqrt(2) != 0.
```

It is labeled `NONPHYSICAL_MIXED_SCALE_CONTROL` and cannot be promoted to an off-plane verdict. For the B-test ray, forced substitution into the OPH SM statistic gives `D=-855 sqrt(2) != 0`; it is labeled `RGE_BRANCH_MISMATCH_CONTROL` and also cannot be promoted.

The producer and independent verifier use exact rational arithmetic and a canonical pair `(a,b)` for `a+b sqrt(2)`. Zero is exactly `a=0 and b=0`; no floating tolerance is used. The Y-to-GUT check rescales all first-row entries `(x1,k1,b1)` by `3/5`; the determinant gets the same nonzero row factor, preserving both zero and nonzero status. Partial row conversions are hostile controls.

## Custody and firewalls

Four minimal pinned K7/K7-Lean source files are vendored under `external/k7_gauge/`. The receipt carries repository, commit, path, Git blob, SHA-256, and byte count. The independent verifier checks the pinned repository/commit/path/blob metadata, recomputes source hashes and byte counts, and rejects contradictory control vectors, determinant statuses and final comparison fields. Its JSON loader rejects duplicate keys and nonfinite constants; integer coefficients cannot be supplied through lossy float-to-integer coercion. CI is offline. The OPH repository field records the Arithmon fork carrying the same pinned OPH commit.

The finite source-typing record also has a separately, manually reviewed semantic digest, binding all annotations, rejection reasons, duplicate values and field names. The digest omits the OPH columns/cofactors/integer locus, the two exact control determinants, and source SHA-256/byte counts: those are independently recomputed by the verifier. It never imports or calls the producer. Changing a source interpretation requires a new review of this semantic binding; matching the digest alone proves neither arithmetic nor source interpretation and selects no physical model.

The vendored main paper contains experimental-comparison text, but no public measurement value is parsed or consumed as a cross-control input; the sealed OPH comparison column remains unopened. Historical construction provenance is a separate question: the pinned sources do not establish a target-blind origin for the target-exposed Type-I agreements. The receipt therefore records construction-input provenance as `UNKNOWN` for those relations and for the conditional B-test; it does not assert that a particular formula was fitted. The `false` measurement-input flags refer only to this cross-control. The K7 historical/conditional epistemic labels are retained. The OPH matter-trace branch is not physically selected, the port-response branch is not identified with K7, and even a future valid off-plane result would not falsify OPH globally.

Reopen only when a pinned K7 artifact supplies a reproducible common-scale inverse-coupling triple with explicit scheme, normalization, field content and boundary role compatible with the frozen OPH branch, or OPH derives a branch matching a frozen K7 vector.

## Observer scale and conditional transport contract

The frozen inputs are OPH `879ed37ad2c1933ecb6dd8f54271362ee3f4d41c`, K7 `26d18fdd2fce61cc0231ec3485f8f5fe7efd5fc6`, and K7-Lean `70dc589bd7c3f77db54a2acb6e67986a7ef69468`; the contract replays on integration baseline `7b0997b4a8a425679d9e2946f82c86cd92289059`. The historical receipt above retains its own pins. The current K7 erratum, §8.2, withdraws the reading of `3/13` as a GUT boundary and the dependent RGE chain. Its electromagnetic relation is at `q²=0`; the weak and strong relations are interpreted at `M_Z`, with no claim of a common-scale triplet. These are source-status statements, not measurement inputs.

The native axial counter has exact covariance `(4/99) I₃` per uniform attempt and diffusion generator `(2/99) Δ`. Its binary current mark has one-attempt law `(28/33,5/33)`. Assigning arbitrary positive length and time units `ℓ,τ` preserves these native records and assigns diffusion `2ℓ²/(99τ)`. This gives a countermodel to identification of physical units from this declared signature. It proves no universal nonselection theorem for every possible OPH construction. The native radius does not itself supply physical continuum position, an RG subtraction scale, or a laboratory energy.

The missing gauge datum is narrower than an absolute common scale. On the frozen one-loop matter branch, multilinearity gives exactly

```text
D(x-t b) = D(x),  D(x)=det(x,k,b),  k=(10/3,2,2), b=(41/6,-19/6,-7).
```

The vanishing-plane test therefore needs no chosen absolute `μ` **conditional on** one common subtraction prescription and compatible branch and threshold transport. Independent component shifts do not obey this identity; neither do arbitrary scheme changes, thresholds, or transport on a different beta branch. In particular, transporting only the first component of `k` by `-b₁` gives `D=943/18`; using the MSSM column `(11,1,-3)` for common transport gives drift `-76/3`. These controls cannot be interpreted as physical K7 comparison results.

The irreducible missing datum at these pins is source-clean transport of the heterogeneous K7 relations into one common subtraction prescription with compatible running, field content, normalization and threshold matching. The matter kinetic branch also lacks a physical selector. No K7 vector is admitted and `NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES` remains the interface verdict. The exact native-scale and conditional transport certificate is `k7_gauge_scale_contract.py`, with receipt `runtime/k7_gauge_scale_contract.json`, independent verifier `verify_k7_gauge_scale_contract_independent.py`, and hostile controls `test_k7_gauge_scale_contract.py`.

## Commands

```bash
python3 code/angular_sprint/k7_gauge_kinetic_cross_control.py --check
python3 code/angular_sprint/verify_k7_gauge_kinetic_cross_control_independent.py \
  --receipt code/angular_sprint/runtime/k7_gauge_kinetic_cross_control.json
```

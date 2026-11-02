# Native energy readout, neutral spectrum and repaired spin dynamics

This records the preceding construction stage. The [next integrated results](GAPLESS_NATIVE_SPATIAL_RESPONSE.md)
add an ordered gapless spin sector, native local observers, lattice stress
refinement and actual spatial particle labels. Earlier scope statements below
refer to this stage's proofs.


This campaign certifies the centered cube's **full finite spectral gap**, executes preparation and a four-bin charged-energy measurement of its actual ground, constructs an interacting spin limit with complete field-repair records, and supplies nine independent matter responses for the existing Einstein theorem stack. The measured-particle dictionary, an interacting relativistic continuum with massive and massless sectors, and the complete gravity requirements are not established.

The source records are the [campaign](native_neutral_spin_response_campaign.json) and [independent review](native_neutral_spin_response_review.json). The four attachments are FA88–FA91. Their finite and continuum statements use the reference states, clocks and apparatus specified below; they cannot be combined by silently changing those inputs.

## Five issues and the resulting evidence

| Issue | Concrete result | Requirement that remains |
| --- | --- | --- |
| [#1138 / CP-02](https://github.com/FloatingPragma/observer-patch-holography/issues/1138), continuing regime and reference | FA89 executes preparation of the centered ground; FA90 constructs source-derived ground states and continuing repaired spin dynamics along a specified volume/coupling trajectory. | Join reference preparation, observer access and continuing dynamics on one admitted process; finite prepared states do not establish every faithful-reference or modular-clock premise. |
| [#1139 / CP-03](https://github.com/FloatingPragma/observer-patch-holography/issues/1139), protecting invariant | FA90 proves that the full two-flavour C/E field instrument preserves every internal spin state when each site has one fermion, including all coin/raw records. | The stated field extension and large-volume observable dictionary need their native observer realization; this is not a proof for an arbitrarily reinterpreted one-bit carrier. |
| [#1142 / CP-06](https://github.com/FloatingPragma/observer-patch-holography/issues/1142), interacting chiral field description | FA90 controls a growing-volume limit with noncommuting Heisenberg interactions, complete repairs, and positive-energy ground-state dynamics. | Its spatial spacing is fixed. A jointly controlled relativistic continuum with the required charged, massive and massless sectors is not established. |
| [#1143 / CP-07](https://github.com/FloatingPragma/observer-patch-holography/issues/1143), particle dictionary and spectrum | FA88 resolves neutral spin/density measures and the full finite gap; FA89 executes the centered charged-energy POVM. | Source energies and hidden-spin labels require relativistic dispersion, physical species identification and a common laboratory-unit calibration. |
| [#1145 / CP-09](https://github.com/FloatingPragma/observer-patch-holography/issues/1145), joint gravity premises | FA91 constructs a common positive-state variation space, stress tangents and finite modular/entropy responses in nine independent directions. | Geometric stress balance, area/Newton normalization, entropy stationarity, joint geometric refinement and metric dynamics remain distinct obligations. |

These results update the evidence for the five issues; they do not by themselves satisfy all five acceptance contracts.

## FA88: complete finite neutral gap and observable spectra

The reference is the unique centered eight-site ground of

```
K = H - N/2,     c=1/8, U=1, hbar=1.
```

It has N=8, hidden spin zero and eta zero. This is different from the original uncentered N=4 ground. All following finite energies are dimensionless source-clock excitation energies.

Exact spin/eta and signed-translation projectors cover all 4900 balanced states and all 12870 N=8 states. Complete integer eigenbasis proposals are admitted by rational Gram and residual bounds. Together with the existing all-other-charge separation greater than 1/8, the new certificate proves the gap of the entire 65536-state source:

```
0.015021254176 <= Delta_finite <= 0.015021254177.
```

The first excited level has multiplicity three, hidden spin one, eta zero and spatial character 111. This closes the missing quantitative **finite** neutral-gap premise. It also bounds finite Euclidean connected correlations by `exp(-0.015021254176*t)` times the corresponding vacuum-vector norms.

| Probe | Certified lowest visible excitation | Positive lines | Exactly dark candidates |
| --- | --- | ---: | ---: |
| Local hidden Sz | [0.01502125417612, 0.01502125417614] | 135 | 1 |
| Local density n−1 | [0.81180294910720, 0.81180294910723] | 135 | 1 |
| Local charged addition/removal, inherited FA85 | [0.36168271180597, 0.36168271180599] | 228 per branch | No unresolved local weight |

Each neutral zero-weight candidate has source energy K=−2. An exact integer eigenvector and its exact energy-−2 probe pullback prove orthogonality to the ground; the zero is not inferred from a small floating amplitude. Every other listed neutral weight is strictly positive, including density weights of order `10^-20`.

Every spatial residue entry is enclosed. Each spin component and density reaches 357 dimensions; the three spin components together with density reach 1428 orthogonal neutral dimensions. The other 11441 excited N=8 states are outside this specified probe-cyclic range. The complete-sector gap certificate includes them.

Full-Fock continuity gives the current-divergence measures by multiplying each density or spin weight by `E^2`. Independent transverse edge-current channels are not resolved by this operation. The checks include the exact sum rules

```
<Sz^2> + <(n-1)^2>/4 = 1/4,
M1(density) = 4 M1(Sz).
```

The [full proof](../centered_neutral_spectrum/README.md), [receipt](../centered_neutral_spectrum/results.json) and [complete neutral line table](../centered_neutral_spectrum/lines.csv) retain these distinctions. A positive fixed-cube gap is not a relativistic rest mass or a uniform continuum gap. The existing charged finite-gap dictionary exclusion also remains: its energy spread is below eleven, which cannot reproduce the charged-lepton ratios using one positive scale.

## FA89: executed centered-ground preparation and four-bin readout

The source-derived ground witness is compiled into 262137 native amplitude-tree pulses, followed by 184 pulses implementing the full-Fock inverse frame. A 13-pulse charge instrument and a 607071-pulse centered clock give **869405 pulses** in a complete continuing history. The preparation retains all 65536 field rows; the clock retains 262144 quantum rows and two classical charge columns.

The clock explicitly emits the number correction `exp(+itN/2)` and the certified ground-energy phase. Combined with the original H evolution, these implement the requested excitation clock `K-E_ground`. Histories have durations 0, 1, 2 and 3, followed by a positive four-point Fourier readout.

The executed joint probabilities below include the charge tag; each column sums to approximately one half. They are rounded numerical outputs of the native simulation, not conditional probabilities obtained by renormalizing success.

| Clock bin | Addition and bin | Removal and bin |
| ---: | ---: | ---: |
| 0 | 0.358982913 | 0.358982913 |
| 1 | 0.089230748 | 0.089230748 |
| 2 | 0.023106883 | 0.023106883 |
| 3 | 0.028679456 | 0.028679456 |

The proved complete-state and total-variation error is approximately **0.0287933351168**, below 0.03. This bound combines the target, source controls, frame, number/ground phases and full-Fock product-formula errors. Numerical probability-accounting residuals and export-failure weights are separate diagnostics.

All 216 raw charge/clock pointer combinations are listed, including 208 proper-code zeros, and every earlier export failure remains an absorbing recorded outcome. The four broad bins measure the certified 228-line charged distribution at finite resolution; they do not resolve 228 separate peaks. No random laboratory shots are claimed.

The [native preparation and readout proof](../centered_native_preparation/README.md) specifies blank registers, coherent memory, control precision, schedule and feedforward. This experiment has no interleaved bare field repair during its preparation or clock evolution. It does not prove autonomous ground selection or continuing protection under that different law.

## FA90: interacting repaired-spin limit and finite spacing collapse

On the one-fermion-per-site, two-flavour carrier, the complete C/E repair instrument preserves the internal spin, with a fair coin and raw record zero. All other branches are retained for inputs outside that carrier. For m free intervals of total duration T, the full recorded-state error obeys

```
error <= min(2, 8m*tau/U + 6T*tau^3/U^2),    tau=c|E|.
```

This holds for arbitrary spin inputs, untouched references and record-adaptive finite histories. With fixed exchange J and independent edge-clock rate Gamma, the explicit joint trajectory

```
c=J|E|^4,    U=J|E|^8
```

gives the Poisson-averaged bound

```
error <= min(2, 8(1+Gamma|E|T)/|E|^3 + 6JT/|E|) -> 0.
```

A spin locality estimate completes the limit to interacting Heisenberg dynamics on the cubic lattice. Actual source grounds supply limiting spin ground states and positive-energy dynamics. The lattice spacing stays fixed, charged excitation energies diverge, and the admitted two-flavour repair extension is explicit. None of those facts supplies the missing relativistic charged continuum.

For a separately specified periodic source family with odd transverse width W, length `L=2W^4`, and the same coupling trajectory, a twist argument gives

```
E1-E0 <= 5J/W^2 + J/W^6 -> 0.
```

This is the first spacing **counting ground multiplicity**. Degeneracy can make that spacing zero without closing the gap above the whole ground space. The bound therefore does not establish an infinite-volume massless particle. See the [complete repaired-spin proof](../protected_spin_limit/README.md).

## FA91: nine independent stress and finite entropy responses

The input is the same centered free eight-Weyl continuum used in FA87. Explicit shell-supported CAR operators define positive normal states on a common local algebra. Their local stress tangents span all nine real symmetric tracefree tensor directions at a point. One Weyl summand is perturbed; the other seven remain in the common vacuum.

Nine actual boosted-ball modular pairings distinguish those directions. For

```
0 < R/L <= 2^-28,
```

the normalized response matrix is invertible with inverse infinity norm at most `272/15`. The same holds for finite vacuum-relative entropy differences with parameter step `(R/L)^6`. Both normalized error channels vanish as `O((R/L)^2)`. Positivity and a uniform quadratic relative-entropy remainder prove the first law across the whole nine-dimensional variation space.

The [stress-response construction](../weyl_gravity_response/README.md) supplies explicit matter data for the existing Einstein theorem stack. It does not supply the missing geometric null balance, edge/area normalization, positive Newton calibration, generalized-entropy stationarity, joint lattice/collar/chart limit or responsive metric. Those conditions must be derived on the same process rather than replaced with matter identities.

## Reproduction and custody

Run the strict package tests and each cold replay from the repository root:

```sh
python -W error -m pytest -q model/centered_neutral_spectrum model/centered_native_preparation model/protected_spin_limit model/weyl_gravity_response
python -W error -m model.centered_neutral_spectrum --check
python -W error -m model.centered_native_preparation --check
python -W error -m model.protected_spin_limit --check
python -W error -m model.weyl_gravity_response --check
```

The native replay executes the complete preparation and clock and is substantially more expensive than the symbolic checks. Each package binds its source, preserves all declared outcomes and reports historical dependency custody separately. New mathematical checks cannot silently repair an older receipt or promote a physical claim. The [campaign](native_neutral_spin_response_campaign.json) and [review](native_neutral_spin_response_review.json) record the integrated status.

The four packages pass **96 strict tests**. An independent integration audit
passes 809 checks while preserving the previous 87 attachments, all 84
canonical premise rows, 346 claim classifications and historical accepted
evidence. All eight primary paper roots and four new proof fragments compile;
the flagship remains **53 pages**. CI workflow files are unchanged. The full
completion-program and attachment-map validators still reject inherited
source-hash mismatches, which remain visible rather than being recertified
without their original executions.

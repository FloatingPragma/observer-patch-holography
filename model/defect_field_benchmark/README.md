# Static charged quasiparticles on the actual colored repair carrier

This package joins the existing calibrated Weyl–Higgs action to the
**same per-particle internal states** used by
[`defect_field_bridge`](../defect_field_bridge/REPRESENTATION.md). The
native hard-core repair branches then exchange two such particles by
`−SWAP`, and the companion source program executes their internal time
evolution. The chosen observer and memory design is a legitimate model
input. Its consequences and the measured calibration inputs are kept
separate.

The spatial carrier is an explicit extension: each of twelve native ports
receives the 90-dimensional internal quasiparticle space. Its full CAR
space has 1,080 modes. The colored two-particle space is a sector of this
extended model, not literally a sector of the original single-cell Fock90.

Run the finite certificate with:

```sh
python3 -m model.defect_field_benchmark --check
python3 -m pytest -q -W error model/defect_field_benchmark
```

## One action, all canonical directions

The input is the existing `calibrated_matter.pairing_kernel()`:

\[
\Delta=\epsilon_{\rm spin}\otimes\overline K,\qquad
B_0=\begin{pmatrix}0&\Delta\\\Delta^\dagger&0\end{pmatrix}.
\]

Here `K` is the actual 45-component symmetric Higgs kernel, including
complex CKM coefficients in the established all-left convention. We do
not replace it with a diagonal list of desired masses. Electric charge
commutes with `B_0`, so each charge block is diagonalized separately.
`quasiparticles()` returns the 90-column matrix `W`, all 180 canonical
columns `[W,CW]`, the energies, the full electric observable and invariant
energy/charge projectors. `C(u,v)=(conj(v),conj(u))`. It checks both CAR
completeness and reconstruction of the **entire** action and charge.

There are 84 positive charged directions and six neutral zero-energy
annihilators. Their six conjugates are also retained. Choosing the original
annihilator half in the zero block fixes a static vacuum convention; it is
not a positive-energy continuation to nonzero momentum. The minimal law
here has massless neutrinos and does not reproduce measured neutrino
oscillations. The existing optional calibrated Weinberg extension is a
different, available action choice whose canonical join is not asserted
by this certificate.

The many-body normal form is

\[
H=E_{\rm vac}+\sum_a E_a\gamma_a^\dagger\gamma_a,
\qquad E_{\rm vac}=-\tfrac12\sum_a E_a.
\]

The bare-particle number is not conserved by pairing. The static
quasiparticle number is conserved, so its two-particle sector legitimately
supplies the two internal colors of the native carrier. Every crossing
acts by `−SWAP`; a global even Fock parity is not substituted for exchange.
All 534,600 hard-core states are present. The 48,060 same-port collision
states are counted as a complement outside this particular repair law.

`validate_chart` accepts arbitrary unitary rotations inside degenerate
energy/charge spaces. Receipts serialize only invariant quantities and
compare floating values with explicit tolerances. No raw eigenvector
equality, visibility cutoff or favorable degenerate vector is used.

## Common clock and source action

The companion [`dynamics.py`](../defect_field_bridge/dynamics.py) uses the
specified native-record feedback

\[
U_j=\exp[-i\tau(1+j/4)(E+E_{\rm vac}/2)+i\kappa jQ],
\quad j\in\{-1,0,1\}.
\]

The scalar vacuum energy occurs once in `U_j⊗U_j`. There is one supplied
clock `τ` in inverse GeV, with `ℏ=1`, for all species. The positive factor
`1+j/4` rescales the whole static action and therefore preserves the common
canonical basis. The native exchange action itself is
`π P_(SWAP=+1)`, whose exponential is `−SWAP`. The full branch combines
this action with the internal phase and the native partial injection;
its original Kraus factor is `1/sqrt(2)`.

Replicating the internal action at twelve ports would give
`H_spatial=12 E_vac+Σ_(x,a) E_a n_(x,a)`. The compiled convention explicitly
subtracts the common scalar `11 E_vac I`, leaving one global `E_vac`.
This reference choice is specified by the model; it does not follow from
the single-cell diagonalization and is not an absolute vacuum-energy
measurement. Restricting that chosen spatial action to the native
hard-core two-particle sector gives precisely the displayed phase law.

This is a constructive operator action on the repair/memory process. Its
mass coefficients are the specified calibrated inputs. It is not an
identity between energy and the number of native attempts, movements or
compiler pulses. The preceding [`defect_energy`](../defect_energy/README.md)
certificate constrains precisely such scalar count proposals.

The source Ramsey witness uses two different charged-lepton energies with
equal charge and an opposite-charge spectator. Its current-conditioned
probability is
`cos²[τ(1+j/4)(E_b−E_a)/2]`. It measures an energy **difference** in the
declared clock units. The common electric phase cancels on this neutral
pair, so the witness does not itself measure charge or an electric
response. A separate neutral QND charge pointer in the same companion
process gives `p0=(1−sin q)/2` and distinguishes the allowed charge values.
Absolute common energy offsets also require a further coherent
reference. The full source outputs, failure and padding channels are in
the companion dynamics receipt; this static certificate does not replace
that source execution.

Once the companion output has been generated, run
`python3 -m model.defect_field_benchmark.readout --check` to replay the
inverse observation. It validates every raw source branch and recovers a
gap interval from the nonzero branches, using the explicit first-lobe
assumption `0≤τ_eff ΔE≤π`. It also inverts the separate charge pointer.
The intervals describe numerical execution error, not finite-sample
experimental uncertainty.

## Calibration and current scope

The nine Yukawas, CKM matrix and Higgs background are the repository's
replayed SMDR 1.3 common-`M_Z` calibration inputs. We independently compare
the diagonalized charged blocks to `y v/√2`: spin multiplicity two for
each lepton sign, and spin-times-color multiplicity six for each quark
sign. These are running action parameters, not direct pole-mass
predictions. The input table has a pinned fit, scheme, scale and external
matching provenance. This package produces zero new independently
predicted mass values and performs no additional fit or scan.

Every one of the twelve original current generators is transformed as
the full Nambu matrix `diag(T,−T.T)`. The primitive source basis mixes
physical axes, so every individual basis current has a nonzero anomalous
block. Exact linear algebra on the actual fundamental generators acting
on the neutral Higgs direction proves a three-dimensional orbit and a
nine-dimensional stabilizer. The corresponding nine exact combinations
have vanishing pair blocks: color plus electric charge. The other three
directions retain their pair-creating blocks and cannot act solely inside
the static two-quasiparticle space. The tests also check the rest spin
`su(2)` algebra and spin Casimir `3/4`; a spacetime spin-statistics theorem
is not inferred from this finite check.

## The mass–inertia constraint

With the **same** calibrated `K`, the existing Weyl symbol uses the Nambu
pair `(p,−p)`, lower block `−h(−p).T`, and supplied speed `c=1`. The finite
Clifford and anticommutator checks yield

\[
E(p)^2=E(0)^2+c^2|p|^2,
\qquad \partial_{p_i}^2E(0)=c^2/E(0),\quad m=E(0)/c^2.
\]

Thus one cannot independently fit an inertial curvature once this action
and clock are fixed. This is the same-form obligation appearing in
[`InternalEnergyInertia.lean`](../../Lean/Geometry/InternalEnergyInertia.lean).
The native port motion has not been mapped to this momentum or speed;
the rest basis is not held fixed at finite momentum, and finite-momentum
source execution requires both momentum sectors. Source propagation,
particle-changing interactions, collision states and physical detector
observables remain the concrete joins needed for particle identification.

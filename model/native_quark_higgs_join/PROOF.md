# Quark and quantum Higgs interference on one native source

The source is the complete eight-cell Hamiltonian `h=H0+2^-34 T` with the
actual thirty-mode fermion, five-state scalar and PW119 link carriers.
Every Yukawa coefficient, scalar penalty and link penalty is one. This
construction uses the original physical two-source preparation, its retained
spin references and energy-agreement pointers. The observable source clock
is `H_total=h_active+h_partner+4N_ref`, where the full thirty-mode static
reference has its actual native currents and initial vacuum. No reference
energy, source sector or scalar register is removed.

## Exact operators and their target images

For each native quark annihilator `A_i=c_i`, define `C_i=[H_local,c_i]`.
The scalar penalty commutes with `c_i`. The literal quantum Yukawa
polynomial gives

`C_i=-adjoint(sum_channel screened_eom(channel,i))`.

`operators.py` checks this identity for all twenty-four quark components,
together with all twelve gauge, all three spin and both raw lepton/baryon
current actions. Each `C_i` is the sum of two scalar partial isometries
times unit-norm CAR operators, so `||C_i||<=2`. This is an operator
identity on all occupations and scalar states. It uses no scalar background
substitution or canonical CAR assertion for the composite.

The source annihilator and its composite transform in the conjugate of
the named native family. Their hypercharge changes are `2/3`, `-1/3` and
`-1/6` for `u^c`, `d^c` and `Q`, respectively. Their raw lepton change is
zero, and three times their baryon change is `1`, `1`, and `-1`.
Each source probe is odd. These finite internal labels do not assert
physical Lorentz chirality.

For a family of internal dimension `d_r`, the compensated operators on the
vacuum domain of the reference are

`A=sum_(spin,a in r) r_(spin,a)^dagger c_(spin,a)/sqrt(2d_r)`,
`C=sum_(spin,a in r) r_(spin,a)^dagger C_(spin,a)/sqrt(2d_r)`.

They are zero on the orthogonal reference-input complement. Orthogonality
of the one-particle reference vectors gives `||A||<=1`, `||C||<=2`.
The reference creation representation contracts the annihilator
representation. Thus both compensated probes commute with all joint
Gauss generators, spin and total raw charges, and have even total parity.
The reference has energy four on success and zero on failure. Full finite
unitary completions retain all other reference occupations. The scalar
contraction over both spin components requires no selected spin direction.
The apparatus addresses the declared origin; no translation-covariant
averaging over apparatus positions is asserted.

The actual particle-hole transformation maps each original polynomial to
the complete 224-dimensional quark target module. The literal composite
map is constructed using all four elementary maps, every scalar row and
the actual sixteen-dimensional source block. Full image inclusion in the
236-dimensional weight carrier and its complete 26/26/52-dimensional
Hamiltonian kernels is checked. In those actual coordinates,

`C_reduced=H_target A_reduced-A_reduced H_source`.

This check fixes the source orientation and the Higgs image directly;
matching target eigenvalues alone would not fix either.

Let `e` be the exact local ground root and `epsilon_j,w_j` the complete
parent quark measure. For every line, including every zero residue,

`elementary weight=w_j`,
`real elementary/composite cross weight=(epsilon_j-e) w_j`,
`composite weight=(epsilon_j-e)^2 w_j`.

The exact source moments `m_k=tr(P_e A^dagger H_target^k A)` give total
cross weight `m_1-e m_0` and composite norm squared
`m_2-2e m_1+e^2 m_0`. The parent proves `m_0=1/2`. All roots and residues
are retained. The local `u^c` and `d^c` measures coincide for these equal
couplings; no physical up/down mass equality or exact equality of their
full interacting measures is inferred.

## Full interacting bin and a positive shared spectral contribution

Write `V=2^-34 T`, `v=360*2^-34`, and let `g` be the actual unique neutral
physical ground. The reference vacuum commutes with `h`, so

`[h,A]=C+[V,A]`, `||[V,A]||<=2v`.

The reference commutator is separately `[4N_ref,A]=4A`. It is not part
of the source excitation energy.

The parent quark package fixes a radius-`128v` bin of the actual operator
`h_active+4N_ref` around the dominant atomic line. Let `P` be its exact
projector, `[L,U]` its source excitation enclosure after subtracting four
and the actual neutral ground, and `[w_lo,w_hi]` its exact-ground
unnormalized elementary probability enclosure. Then

`P C g=(h-E_g) P A g-P[V,A]g`,
`Re <P A g,P C g> >= L w_lo-2v = gamma > 0`.

Cauchy-Schwarz gives `||P(C/2)g||^2 >= gamma^2/(4w_hi)`.
Because the finite spectral sum of real cross weights is positive, at
least one actual eigenspace in this bin carries a positive common
contribution. A selector maximizing that positive cross weight is defined
on the finite set. It differs from a maximum-elementary-residue selector;
no equivalence between these selectors is assumed. The result neither
isolates a single pole nor identifies the full compensated charged ground.

## Complete interference records and preparation errors

For an independent fair classical setting `s=+1,-1`, use the contraction

`M_s=(A+s C/2)/2`, `||M_s||<=1`,
`K_s=sqrt(I-M_s^dagger M_s)`.

A Halmos unitary supplies the complete pump. The energy service is
controlled by pump success. On failure the energy pointer remains zero;
on success it records `P` or `I-P`. On an exact ground input,

`p_(+,hit)-p_(-,hit)=Re <P A g,P C g>/2 >= gamma/2`.

The fair-setting joint contrast is at least `gamma/4`. These are
unnormalized pump-and-bin probabilities; neither line normalization nor
an additional factor one half is inserted.

Let `q` be the inherited physical preparation success probability and
`rho` its conditional active-source state, with
`D(rho,|g><g|)<=delta`. The complete flagged detector instrument contracts
trace distance. Its signed hit statistic has values `+1,-1,0`, so the
joint contrast changes by at most `delta`: equivalently, the difference
of the two fair-setting effects is half the difference of two effects
in `[0,I]`, with spectral range at most one. The original-attempt contrast
is therefore at least `q_lo (gamma/4-delta)>0`. No exact-ground projection
is inserted into the preparation.

For each sign the full conditional law is

`Tr(K_s rho K_s^dagger)`,
`Tr(P M_s rho M_s^dagger)`,
`Tr((I-P) M_s rho M_s^dagger)`.

Their sum is one. Multiplication by `q/2`, the complete fine preparation
failure law of total weight `1-q`, and the two explicit zero
failed-pump/energy-hit labels give the full instrument. The setting,
pump flag and energy flag have eight labels after successful preparation,
with six potentially nonzero labels. Every native preparation reference,
failed filter branch and source complement is retained.

`bounds.py` also encloses every individual interference record. Its atomic
line weight is

`w_(j,s)=w_j (1+s(epsilon_j-e)/2)^2/4`.

The complete sum equals the number-field expression
`(m_0+(m_2-2e m_1+e^2 m_0)/4+s(m_1-e m_0))/4`.
The local moment operators of the spin-summed compensated probes are
scalars on the atomic ground doublet, so these are also the probabilities
in the actual spin-entangled cube trial `u`.

For the selected atomic vector, whose interval lies at distance `d` from
the bin complement, the projected norm is at least
`sqrt(w_(0,s)(1-(v/d)^2))`. Each other line contributes at most
`v sqrt(w_(j,s))/d_j` to its projection into this bin. All other lines,
including tiny and zero weights, enter that sum. Subtracting the leakage
sum and the certified `||g-u||<=eta` gives the clipped lower amplitude.
The upper amplitude adds these errors. This uses `||M_s||<=1`.
The total pump probability differs from its atomic value by at most
`eta`, since it is an effect and the pure-state trace distance is bounded
by `||g-u||`. Subtraction gives the outside-bin interval, and complement
gives pump failure. The actual prepared probability of each effect then
changes by at most `delta`. Multiplication by `q_lo/2` gives each
original-attempt joint setting lower bound. All arithmetic is outward.

## Clock, apparatus and scope

Pump and spectral readout each take one source unit. For each finite
unitary service `U=exp(-iL)`, the declared drive
`V(t)=exp(-it H_total)Lexp(it H_total)` implements full evolution
`exp(-iH_total)U` over that interval. Common autonomous source evolution
continues on every failure and no-operation branch. The energy projectors
commute with that evolution, so the record law includes these durations
without a long-time perturbative replacement of the actual clock.
The conditional prepared state need not be stationary.

The finite apparatus includes the exact full-source spectral-bin service
and programmable unitary completion. These are declared services; no
efficient synthesis, primitive repair implementation or enormous
full-matrix execution is claimed. The retained storage upper bound is
`D^2 2^30 2^32 Mprep^2 2^4`, where
`D=2^240 5^8 119^24`, `Mprep=1197721483502453176152084669`, and the last
four bits include the retained preparation flag, setting, pump and energy.
Controller synthesis workspace is additional. The full bare Hamiltonian
norm is bounded by `1080+720*2^-34`, so work per completed service is at
most twice that norm.

These are finite quantum-Higgs interference records with native quark
current labels on the original source. They supply no quark propagation,
asymptotic colored state, physical rest mass, measured weak chirality or
Standard Model matching theorem.

## Normalized admitted states and common detector histories

A separate setting chooses the contraction `X=A` or `X=C/2`, using the
same full reference, preparation, source Hamiltonian and absolute bin.
Let `a=PAg`, `b=PCg`, `e_c=(L+U)/2`, `d=(U-L)/2` and `s=sqrt(w_lo)`.
The full EOM gives

`||b-e_c a|| <= d||a||+2v <= r||a||`, `r=d+2v/s`.

The executed bounds have `0<r<e_c`, so `||b||>=(e_c-r)||a||>0`.
The distance from `e_c a` to the line spanned by `b` is at most
`r||a||`. Hence the normalized exact-ground states have trace distance
at most `r/e_c` and infidelity at most `(r/e_c)^2`. This geometric
inequality is the pure `filtered_bounds` lemma in
`model/native_higgs_filtered_state`; only that inequality is reused.
Its lepton source, selected channel, extended port reference and
preparation depth are not inputs here.

The actual original two-source preparation certifies a conditional pair
fidelity `F_pair`. Therefore its active-source ground fidelity is at least
`F_pair`. For a purification `Psi` of the complete conditional preparation,
there is a normalized common environment vector `chi` with

`||Psi-g tensor chi|| <= epsilon=sqrt(2(1-F_pair))`.

This environment includes the physical partner and retained preparation
references, with any mathematical purification registers added for the
comparison. Both settings use the same `chi`. Writing
`a_Psi=(PA tensor I)Psi` and `b_Psi=(PC tensor I)Psi` gives

`||a_Psi|| >= s-epsilon`,
`||b_Psi|| >= (e_c-r)s-2epsilon`,
`||b_Psi-e_c a_Psi|| <= r||a_Psi||+(r+2+e_c)epsilon`.

Thus the normalized accepted states, including their common retained
preparation environment and static reference, obey

`D(rho_A,rho_C) <= [r+(r+2+e_c)epsilon/(s-epsilon)]/e_c`.

The upper bound is below `5.193e-6` for every family. The comparison
conditions on successful preparation, pump and admission in the same
absolute bin. The classical setting label is excluded from the compared
state. A subsequent common detector policy may use the source, retained
reference, preparation environment and identical accepted flags; it may
not read the different setting labels. Every finite common instrument
history, including adaptive policies on its own subsequent records,
contracts trace distance and therefore has total-variation difference at
most this bound. Common actual-H waits of arbitrary duration do not add
a perturbative clock error. The instrument chosen using the initial pump
setting is outside this comparison.

There is an independent source/reference-only check using the actual
preparation trace bound `delta`. For one fixed complete flagged
pump/bin instrument, let `A_s,B_s` be its subnormalized accepted outputs
on the prepared and exact-ground inputs, with traces `p_s,q_s`. Output
contractivity and the complementary flag imply

`||A_s-B_s||_1+|p_s-q_s| <= 2delta`.

Since `A_s/p_s-B_s/q_s=(A_s-B_s)/q_s+(1/p_s-1/q_s)A_s`, the normalized
trace distance is at most `delta/q_s`. Applying this separately to `A`
and `C/2`, then adding the exact-ground line-distance bound, gives the
recorded independent marginal comparison. This check does not replace
its conditioning denominators with an unnormalized error.

The complete experiment retains, for each fair setting, preparation
failure, `K_X=sqrt(I-X^dagger X)`, successful admission and its complement.
The failed-pump energy-hit labels are exactly zero. Its eight post-
preparation setting/pump/energy labels have six potentially nonzero
outcomes. The code gives outward intervals for all six, original-attempt
probability lower bounds, and distinct total-pump and admitted-bin success
rates. No closeness of the unconditional setting laws is asserted. The
same storage and clock services stated above apply, with this binary
setting replacing the interference setting.

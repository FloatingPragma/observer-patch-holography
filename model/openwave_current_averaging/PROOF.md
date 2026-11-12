# Uniform averaging of the unfiltered native weak current

## Source and free-mode transformation

The source is the complete 270-mode fermionic Hamiltonian of
`model/openwave_weak_decay`: fifteen native internal coordinates, three
families, two spins and three momenta `{-q,0,q}`. Its 31464 canonical quartic
current words are unchanged, including all quark/lepton currents, spectators,
scattering and reference response. The free reference fills negative bands
and leaves zero bands vacant. Charged masses are `(m,1,M)=(e_1/e_2,1,e_0/e_2)`
and `q=(1-m)/2`. The source is `H=H0+G V`, with declared comparison coupling
`G=2^-72`, in the finite-box convention of the native contact. No physical
Fermi normalization or continuum rate is assigned to this coefficient.

No occupation penalty or selected transition Hamiltonian is added.
Free-clock averaging compares this unfiltered source with its complete
resonant current. A small invariant orbit is derived from that averaged
current after retaining every source mode.

For charged left mode `(p,s,f)` set `epsilon=-sigma_s p`,
`E=sqrt(p^2+m_f^2)`. Its actual Higgs singlet partner is at `(-p,1-s,f)`.
With the native Higgs spin sign, define

`u=sqrt((E+epsilon)/(2E))`, `v=sigma_s sqrt((E-epsilon)/(2E))`,
`c_L=u a_L-v a_R†`, `c_R=u a_R+v a_L†`.

These transformations preserve every canonical anticommutator and diagonalize
the same free pairing. Rest coefficients are signed `1/sqrt(2)`. For each
unpaired massless mode, `c=a` on a vacant nonnegative-energy mode and `c=a†`
on a filled negative-energy mode. All 270 modes are transformed.

Every quartic word becomes ordered quasiparticle words, with their actual
contractions retained when acting on states. There are 55575 ordered terms
and 299 distinct Bohr-frequency labels. A word's doubled frequency is an
integral linear combination of

`(m,1,M,sqrt(q^2+m^2),sqrt(q^2+1),sqrt(q^2+M^2))`.

The massless frequency `(1-m)/2` is included exactly. Native analyzer rational
intervals and integer square-root enclosures certify every nonzero label
away from zero. Their minimum modulus exceeds `1/50000`, numerically about
`2.35028e-5`. The zero label contains 14675 ordered terms. No near-resonant
word is silently made resonant.

## Exact complete resonant orbit

Let `Vbar` sum all zero-frequency words. Start with the middle charged rest
mode of spin minus one half and an idle opposite-charge reference. Apply
every resonant word to it, then to each additional state its exact column
generates. For nonzero momenta the Bogoliubov coefficients are represented
by `x_f,y_f`, with `x_f^2+y_f^2=1`. Polynomial reduction by these identities
uses exact rational arithmetic, retaining `sqrt(2)` separately. Every
coefficient outside the resulting orbit cancels identically, without a
floating-point zero threshold.

The orbit contains precisely the parent and a light charged rest state of
opposite spin, a middle-family neutral at `+q` and a light-family antineutral
at `-q`. Both states have free energy one. In canonical quasiparticle phases,

`Vbar|orbit = [[173/2,1],[1,85]]`.

A receiver phase change reverses both off-diagonal signs and changes no
record. The diagonal entries are actual current expectations and differ by
`3/2`; deleting them would remove a source-induced detuning. Thus

`p_receiver(t)=(16/25) sin^2(5Gt/4)`.

The first peak is `16/25` at `Gt=2pi/5`; parent probability is `9/25`.
Every other averaged source output has probability zero on this preparation
by the proved invariant orbit. The unfiltered source retains all those
outputs with the bound below.

## Uniform full-source comparison

Each canonical quartic word is a contraction. The native coefficient sum
`L=41952` bounds `||V||`. Cesaro averaging of the full unitary free-clock
conjugates gives `||Vbar||<=L`. Each transformed canonical factor has
coefficient absolute sum at most `sqrt(2)`, so the expanded word coefficient
sum is at most `4L`. With `delta=1/50000`, the oscillatory primitive

`F(t)=int_0^t [exp(iH0s)V exp(-iH0s)-Vbar] ds`

obeys `||F(t)||<=2*(4L)/delta=:A=16780800000` for every real time. This uses
all words and their certified frequencies, independently of Fock dimension.

Let `U_I(t)` be the exact full interaction-picture propagator and
`Ubar(t)=exp(-iG Vbar t)`. Duhamel's formula, integrated by parts with
`F'=V_I-Vbar`, gives

`||U_I(t)-Ubar(t)|| <= G A (1+2G L t)`.

The endpoint contributes `G A`; derivatives of the two unitary factors
contribute `G^2 A t(||V||+||Vbar||)`. Unitarity removes any exponential norm
factor. Tensoring any external quantum reference preserves this bound.

At `G=2^-72` and `0<=Gt<=2pi/5`, the rational bound `pi<22/7` gives a uniform
vector error `epsilon<1/2500000`. Parent and receiver probabilities differ
by at most `2epsilon<8e-7`. The full orthogonal remainder has amplitude at
most `epsilon`, hence probability below `1.6e-13`. At the first peak the
unfiltered receiver probability lies in `[0.6399992,0.6400008]`; the receipt
retains sharper rational endpoints. This is a full-source theorem at a long
declared clock time, not enumeration of `2^270` amplitudes. The exact derived
orbit is executed at five scaled times and compared with its analytic law.

At `G=1/64` this bound exceeds one and gives no useful approximation. It does
not license replacing that stronger source by its resonant orbit or by the
separate continuum decay member.

## Records, energy and controls

The complete instrument consists of parent, receiver and the full orthogonal
source remainder, without conditional renormalization. Initial full energy
is `1+(173/2)G`; the pure receiver expectation is `1+85G`. Autonomous total
energy is conserved. Complete readout dephases the interaction and can do
work. The three projectors commute with `H0`, so absolute mean readout work
is bounded by `2G L`. The executed averaged records retain their actual work
divided by `G`.

Zero coupling gives zero receiver probability. Deleting the computed
detuning changes the averaged peak from `16/25` to one and falsifies that
altered effective model as a comparison of this source. The resonant box,
weak coupling, reference preparation and clock are declared; no irreversible
decay, observed lifetime or physical neutrino map is inferred.

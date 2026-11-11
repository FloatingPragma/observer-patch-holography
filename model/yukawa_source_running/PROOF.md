# Finite Higgs thresholds and conditional running

## Source and domain

The source is the complete charged-lepton Higgs join in
[`interacting_extension/higgs_join`](../interacting_extension/higgs_join/PROOF.md).
It retains both full native gauge links, all Gauss states in sectors
`Q*=N_F-2N_H=0,1`, scalar occupancy `N_H=0,1`, and six compensated
charged-addition probes. Their dimensions are 438 and 1,596. Colored matter
occupation and the quark Higgs vertices are zero in this source domain.
The static compensating reference is part of each charged probe.

The family matrix `Y` is reconstructed from the complete 128-record native
law, its 12-state memory transition, its stationary probabilities and the
native ray projectors. The producer calls `exact_source_packet`; the retained
expression is checked against that reconstruction. Its singular values are
approximately `(0.000264834486,0.051673720427,0.948061445086)`.
The second hopping matrix is the source ordered-moment expression
`sigma=Y+(1-cos(log(5/4)))D+sin(log(5/4))J`.
The preparation and feedback services in that producer retain their declared
architecture. This source is distinct from the five-state memory in
`native_repair_mass`; its dressing bound is not transferred here.

The Hamiltonian has gauge-field term `H_G`, hopping term `hT`, Higgs term
`g(W+Wdagger)`, scalar energy `sN_H` and reference energy `rQ*`.
This package derives a finite consequence for `h=0,s>0`; it computes
complete numerical controls at nonzero `h`. The coefficients `g,s,h,r`,
the ground-sector preparation and a physical clock are not fixed by the
native family record law. `build_join` exposes these choices explicitly.
An addition frequency is an energy difference between the two prepared
reference sectors. It is not a dimensionless Yukawa coefficient.

## Exact blocking calculation

Write `S=Tr(Ydagger Y)` and order the right-handed creation modes by
family and spin. The literal Higgs vertex creates one scalar, one weak
doublet fermion and one charged singlet. Each scalar component has a unit
native invariant coefficient. Both antisymmetric spin terms are retained.
The independently enumerated canonical anticommutation rules give

```
Wdagger W |vacuum> = 4 S |vacuum>,
Wdagger W |one e> = [4 S I6 - 2 (Ydagger Y)^T tensor I2] |one e>.
```

The transpose follows from the conjugate-right all-left convention.
`exact_local_car()` verifies all 36 entries for nine arbitrary complex
symbols. It does not assume that `Y` is real or diagonal.

Every column of the exact rational retained carrier is also checked.
For all nine family entries, `H_G W=W H_G`; for each of the six additions,
`H_G Bdagger=Bdagger H_G`. The 72 identities separately include the
electric operator, loop and loop adjoint. Every `W` entry raises scalar
occupation by one. Consequently the scalar parity `(-1)^N_H` conjugates
`g` to `-g`, without changing the spectrum. All scalar and complementary
states in the finite join remain included.

Fix a normalized vacuum gauge eigenvector of energy `E_G`. Its cyclic
Higgs subspace has matrix

```
E_G I2 + [[0, 2 g sqrt(S)], [2 g sqrt(S), s]].
```

Let `y_i` be the three singular values. The right singular family direction
and either spin give a charged cyclic subspace with coupling squared
`c_i^2=4S-2y_i^2`. Define

```
d0 = sqrt(s^2+16 g^2 S),
di = sqrt(s^2+16 g^2 S-8 g^2 y_i^2).
```

Prepare the lower eigenstate of the vacuum cyclic subspace. The two
charged-addition frequencies, including an arbitrary reference offset, are

```
delta_i = (d0-di)/2+r = 4 g^2 y_i^2/(d0+di)+r,
omega_i = (d0+di)/2+r.
```

These are exact eigenstate differences on this reducing cyclic subspace.
The full-carrier calculation separately checks that the chosen state is the
unique `Q*=0` ground state for every listed nonzero-coupling control. No
claim that it is the absolute ground state for all possible scalar and
gauge parameters is required by the cyclic-subspace theorem.

At `g=s=1,r=0` the three lower frequencies are approximately
`(3.57176314e-8,0.00136026487,0.52897528448)`.
For a small `y_i`, `delta_i=2g^2 y_i^2/d0+O(y_i^4)`.
This is an interaction-induced addition cost from Pauli blocking. It
amplifies an inherited hierarchy without generating a Yukawa vertex from
zero. If `Y=0`, all lower addition costs are exactly zero.

For `0<y_i<y_j`, `g!=0` and `r=0`, rationalization proves

```
delta_j/delta_i
 = (y_j/y_i)^2 (d0+di)/(d0+dj)
 >= (y_j/y_i)^2.
```

Exact leading-principal-minor inertias bracket the source eigenvalues.
The resulting adjacent lower bounds are `106502400/2809>37914` and
`89870400/267289>336`. They hold for every positive `s` and nonzero `g`
in the no-hopping, zero-offset domain. The comparison values are not used
to select this formula, its parameters or these bounds.

## Complete spectral weights and coupling inverse

Take `g>0`, without loss of spectral generality. Put

```
a0=sqrt((1+s/d0)/2), b0=sqrt((1-s/d0)/2),
ai=sqrt((1+s/di)/2), bi=sqrt((1-s/di)/2),
ki=sqrt((4S-2y_i^2)/(4S)).
```

The two creation probabilities per family and spin and the blocked
probability are

```
w_i^-=(a0 ai+b0 bi ki)^2,
w_i^+=(-a0 bi+b0 ai ki)^2,
p_i^blocked=b0^2 y_i^2/(2S).
```

They sum to one. The zero-matrix limit is taken continuously. The fixed
coordinate probes have residue matrices `w_i^+/- Pi tensor I2`, where
`Pi` are spectral projectors of `(Ydagger Y)^T`. Thus no eigenvector is
inserted into a preparation or detector. The analytic record law is
compared against all 1,596 final states of the complete carrier at five
times, with the blocked effect checked separately. All dark rows are
retained in the exported spectral arrays. A visibility threshold controls
only the displayed summary.

For the primary source, the smallest upper-line weight per spin is about
`4.9e-17`; discarding it changes the exact inverse problem. No finite-shot
precision certificate is claimed.

If `r=0`, the exact identity
`delta_i omega_i=2g^2 y_i^2` recovers the linear Yukawa ratios from both
lines. A complete failure record also removes an unknown common `r`.
Let `B` be the trace of the six-probe blocked effect and assume all three
lower/upper line pairs have been resolved by their common family
projectors. Then

```
B=b0^2,  kappa=1-2B=s/d0,
di=omega_i-delta_i,
d0^2=2 sum_i di^2/(5+kappa^2),
r=(omega_i+delta_i-d0)/2,
g^2 y_i^2=(d0^2-di^2)/8.
```

The identity for `d0` follows by summing
`di^2=d0^2-8g^2y_i^2` and using
`16g^2S=d0^2-s^2`. `reconstruct_coupling` takes only the paired lines
and blocked effect. Controls recover `|g|y_i`, `s` and the reference
offset, including a nonzero offset. It does not separate `g` from `Y`,
select a source coupling, reconstruct unobserved left-family phases or
establish a physical Yukawa matching convention.

This gives a constructive route to the nonzero-hopping source: reconstruct
the complete matrix response and its Schur self-energy, retain the measured
wave-function residues, and extract the Higgs vertex with a specified
external-leg normalization. The exact no-hopping inverse is a control for
that calculation. A source-selected interacting preparation, reference and
propagating field identification are separate required premises.

## Hopping and reference controls

The computation retains the complete spectra and addition records at
`h=1/4,1/2,1`, using the actual `Y,sigma` hopping matrices. At `h=g=s=1`
the first three spin-paired addition frequencies are approximately
`(2.58004e-8,0.00100665351,0.452014700)`. With the same hopping and
`g=1/4`, negative-frequency charged-addition records occur. No branch is
discarded when it crosses zero. The finite calculation therefore supplies
no parameter-independent positive mass identification.

Adding `rQ*` shifts every addition frequency by `r`, retaining all residues
and blocking probabilities. The opposite-Higgs control verifies the parity
identity. A common identity family coefficient gives equal family lines;
zero coefficient gives zero lower lines at zero hopping. A genuinely
complex family control checks the conjugate-right convention against the
complete finite matrices. All controls are declared, with no empirical
ranking or fitted selection.

## Conditional one-loop running

The minimal Standard Model is supplied as a separate continuum comparison
theory, with no right-handed neutrino Yukawa. Use GUT-normalized `g1` and
`t=log(mu/mu0)`. The implemented matrix equations are the one-loop
equations (3) and (4) of Xing and Zhou,
[arXiv:1804.01925](https://arxiv.org/html/1804.01925v2), with `Ynu=0`.
Their beta functions have the form `dY_f/dt=A_f(t)Y_f`. Whenever the
solution is finite, the fundamental matrix solving `dF_f/dt=A_f F_f`,
`F_f(0)=I`, is invertible. Hence

```
Y_f(t)=F_f(t)Y_f(0), rank Y_f(t)=rank Y_f(0).
```

Zero input does not generate a vertex. The coupled polynomial flow is
locally unique and reversible on every nonsingular finite interval.
It cannot erase all boundary dependence in finite running time. Family
rotations act covariantly. With identity matrices as every family input,
the equal-family subspace is invariant, so no hierarchy emerges from a
flavour-symmetric boundary. This statement is about the stated equations,
not arbitrary nonperturbative ultraviolet completions.

For charged leptons, the left singular vectors are fixed and

```
16 pi^2 d log(y_j/y_i)/dt = (3/2)(y_j^2-y_i^2).
```

Thus every positive ordered charged-lepton ratio contracts toward the
infrared, independently of the running quark and gauge values. A larger
middle/small ratio cannot be obtained by minimal one-loop infrared
evolution from a smaller source ratio. Two-loop corrections, extra neutrino
couplings and additional flavour-breaking thresholds are outside this
theorem.

The numerical control uses `(Yu,Yd,Ye)=(Y,sigma,Y)`, gauge values
`(1/2,1/2,1/2)` and `t=-40`. These are explicit matching and scale choices,
without a GeV interpretation. The source lepton ratios change from
`(195.1170,18.3471)` to `(195.0282,16.0571)`. The calculation checks exact
gauge running, covariance, rank-one and zero support, and inverse flow.
An analytic scalar attractive-fixed-point control also retains strictly
positive sensitivity to its nonzero boundary value at finite time. No
fixed point or dimensional-transmutation scale is selected by this source.

## Empirical separation and scope

The producer reads no empirical mass or mixing packet. `compare.py`
separately reads the named empirical common-`M_Z` Yukawa packet after the
complete source calculation. Its ratios are comparison quantities; their
mass-to-Yukawa matching was calibrated. Neither those numbers nor their
normalization enter the threshold or running producer.

The exact no-hopping threshold ratios differ by orders of magnitude from
the charged-lepton comparison. The primary source `y2/y1` lies below its
comparison ratio and cannot be raised by the stated infrared flow.
The actual source coefficient, the finite addition spectrum and a
canonically normalized continuum Yukawa vertex are three different
quantities. The inverse proves how to distinguish the first two within
the declared finite architecture. It does not close their physical
continuum matching, select the scalar controls or derive laboratory masses.

# Observations define the reported response

## Interpretation

Within this OPH investigation, a particle identification is an interpretation
of committed detector records. A Hamiltonian eigenvalue is an intermediate
calculation, not by itself a particle. The implementation therefore supplies
a complete instrument, its back-action and outcome probabilities, and
simulated retained records. It does not claim that the finite model or its
ground state is selected by OPH.

## Interacting addition response

Let \(H_0,H_1\) act on the complete physical no-reference and compensated
charged sectors of the joined gauge–matter–Higgs control. The prepared source
state is the unique ground \(\Omega\) of the declared \(H_0\), with energy
\(E_0\). A charged creation probe \(b_j:\mathcal H_0\to\mathcal H_1\) creates
the fermion and the compensating detector reference together. All six
family/spin probes are retained.

The response inferred from measurements is

\[
C_{jk}(t)=\langle\Omega|
b_j^\dagger e^{-i(H_1-E_0)t}b_k|\Omega\rangle.
\]

For a complete \(H_1\) eigenbasis,

\[
C_{jk}(t)=\sum_n e^{-i(E_n-E_0)t}
\langle b_j\Omega|n\rangle\langle n|b_k\Omega\rangle.
\]

Every frequency, including any nonpositive frequency, must be retained.
Taking absolute values of addition energies could hide an inconsistent
reference calibration or an unstable prepared sector. A vanishing residue
cannot be promoted into an observed feature.

## Complete forward-only Ramsey instrument

The literal CAR maps in this control are partial isometries. On
\(\mathcal H_0\oplus\mathcal H_1\), define

\[
G_j=\begin{pmatrix}0&b_j^\dagger\\b_j&0\end{pmatrix},
\qquad
W_{s,j}=e^{-is\pi G_j/2}=I-G_j^2-isG_j,\quad s=\pm1.
\]

The unitary identity follows from \(G_j^3=G_j\). Define
\(x_j=b_j\Omega\) and \(a_j=(I-b_j^\dagger b_j)\Omega\); the latter is the
blocked branch and is **not discarded**.

Prepare an ancillary detector in \(|+\rangle\). In path zero apply forward
evolution \(U(t)\), then \(W_{+,j}\). In path one apply \(W_{s,k}\), then
the same forward \(U(t)\). The cross overlap is

\[
F_{s,jk}(t)=A_{jk}(t)+sC_{jk}(t),
\quad
A_{jk}(t)=\langle a_j|e^{-i(H_0-E_0)t}|a_k\rangle.
\]

Measuring the ancillary \(X\) or \(Y\) gives the complete laws

\[
p_{\pm}^{X}=\frac{1\pm\operatorname{Re}F_{s,jk}}2,\qquad
p_{\pm}^{Y}=\frac{1\pm\operatorname{Im}F_{s,jk}}2.
\]

Both outcomes of every setting are recorded. Taking the difference between
the two sign settings reconstructs

\[
C_{jk}=(F_{+,jk}-F_{-,jk})/2.
\]

There is no postselection and no reversal of the system's time evolution.
The complete postmeasurement maps themselves are

\[
M^X_\pm=(W_{+,j}U\pm U W_{s,k})/2,\qquad
M^Y_\pm=(W_{+,j}U\mp iU W_{s,k})/2.
\]

They obey \(\sum_\pm M_\pm^\dagger M_\pm=I\) on the full direct-sum
carrier. The public ramsey_kraus_actions API returns both unnormalized
postmeasurement states for arbitrary vectors or a purification. Sample
postmeasurement vectors for the prepared input and an independent arbitrary
input are retained. Numerical probability clipping, if any, is logged.

The quarter-turn probe itself, controlled detector ancilla and initial
ground preparation are declared laboratory services; their autonomous
native implementation has not been derived.

The compensating reference carries the appropriate gauge charge and is
assigned odd fermion parity, making the compound probe even. That parity
assignment is part of the measurement model. A static charge label alone
does not automatically supply a physical detector with those properties.

More explicitly, order the reference CAR mode after all matter modes and
define \(b_j=c_j^\dagger d_{\rm ref}^\dagger\),
\(b_j^\dagger=d_{\rm ref}c_j\). Reference creation acts first on an even
matter state; reference annihilation acts after the matter annihilation has
left an even state. Both reference signs are therefore positive. Reversing
the creator order introduces a minus sign. The separate exact carrier
audit checks every retained probe/basis case.

The implementation checks \(X/Y\) signs by direct joint ancilla evolution
on an eight-dimensional control and checks the branch overlaps independently
on the full 2,034-dimensional joined system, including cross-family probes.

## Statistical records

The simulation retains packed bits for every shot of every diagonal
detector setting. For \(M\) such Bernoulli settings and \(S\) shots each,
Hoeffding's inequality and a union bound give simultaneous error radius

\[
\epsilon=\sqrt{\frac{2\log(2M/\alpha)}{S}}
\]

for all measured \(\pm1\) means with failure probability at most \(\alpha\).
Half-differences inherit the same conservative bound. The files label these
as simulated records, not experimental observations. Sampling uncertainty
is not identified with an intrinsic particle linewidth.

The finite record window is also explicit. A very small calculated
frequency need not accumulate an identifiable phase over that window.
The existence of a numerical response line is therefore distinguished
from successful finite-data identification of that line.

## An exact reference-calibration ambiguity and its gauge cancellation

The joined carrier conserves \(Q_*=N_F-2N_H\). Thus

\[
\epsilon_F N_F+\epsilon_H N_H
=\epsilon_F Q_*+(2\epsilon_F+\epsilon_H)N_H.
\]

Changing \(\epsilon_F\) by \(\delta\) and \(\epsilon_H\) by \(-2\delta\)
adds the constant \(\delta q\) in charge sector \(q\). All within-sector
record laws are unchanged, while every \(q=0\to1\) addition frequency
shifts by \(\delta\). The reference detector's own energy also shifts
these addition frequencies. This comparison holds the reference energy
fixed.

There is an essential qualification. In this restricted sector
\(N_L=N_H\), \(N_{e^c}=Q_*+N_H\), so the total native matter-plus-scalar
hypercharge is

\[
q_{6Y}=6N_{e^c}-3N_L-3N_H=6Q_*.
\]

If the offset is a common gauge potential, the compensating reference
shifts by \(-\delta\). Its contribution cancels the matter shift exactly.
The implementation checks that the resulting **complete** Hamiltonians,
and therefore all record laws, are unchanged in both sectors.

This exposes an incompletely specified reference calibration. It is
**not** a theorem that intrinsic particle masses are arbitrary, and it
does not license fitting masses with a gauge artifact. A common source
model of matter and measuring reference must fix the physical addition
protocol and implement consistent gauge cancellations.

Likewise \(H\mapsto\kappa H\) together with \(t\mapsto t/\kappa\)
leaves this complete dimensionless record law unchanged. A common physical
clock is required to infer an energy readout. A numerical record-capacity
bound alone does not fix that clock or the conserved-charge offset.

## Information available in the finite shot schedule

A short observation interval alone is not a statistical impossibility theorem.
For the actual independently prepared Ramsey shots, a quantitative bound is
available. Compare the stated Hamiltonian with
`H1' = H1 - delta P`, where `P` is an orthogonal spectral projector of `H1`.
The initial state, both probes, the neutral Hamiltonian and the detector
calibration are identical in this comparison. Write
`w = ||P b_k Omega||²`; partial isometry gives `0 <= w <= 1`.
Only the charged part of the second ancilla path acquires the additional
phase. The two normalized pre-readout pure states have overlap

```
z(t) = 1 + (exp(i delta t) - 1) w/2,
1 - |z(t)|² = w(1-w/2)(1-cos(delta t)) <= delta² t²/4.
```

For independent trials indexed by `i`, pure-state fidelity multiplies.
The inequality `1-product(1-x_i) <= sum(x_i)` and contraction under the
complete detector measurement imply

```
TV(all retained records) <= min(1, |delta| sqrt(sum_i t_i²)/2).
```

This bounds every inference procedure applied to those records. It does not
assume a Fourier estimator or discard unsuccessful detector outcomes.
For the supplied schedule there are `24*4096` shots at each of the rational
times `i/10`, `i=0,...,128`. At `|delta| <= 3/10^8`, the upper bound is
strictly below `1/2500`. Thus a commuting charged-sector energy shift of
this size changes any event probability in the entire recorded experiment
by less than 0.0004. A prior restriction that excludes the comparison
Hamiltonian is additional information, outside this bound.
An observer allowed to read the simulator specification may obtain that
information and compute the frequency. This is a bound on the finite detector
record alone, not on every inference available to an OPH observer.

The theorem is exact, and the squared bound is evaluated with rational
arithmetic. Its application to an exactly vanishing smallest line requires
an enclosure of that line; the separate 50/80-digit numerical audit is not
an interval certificate. The theorem does not promote that numerical audit
into an enclosure or identify the line with a physical particle mass.

## Executable domain checks

The record API requires finite Hermitian Hamiltonians and genuine partial
isometries, since a general contraction does not make `I-G²-iG` unitary.
Supplied eigenpairs are checked for full dimension, ordering, orthonormality
and their residual against the actual Hamiltonian before use. Both Kraus
outputs act on arbitrary full-carrier vectors or purification columns.
The tests compare their complete matrices with independent exponentials,
including complex probes and noncommuting Hamiltonians.

The source contract checks consumed bytes and their static project import
closure against the recorded source commit. It does not infer provenance
from the working tree's commit label. Cached source data revalidate their
inputs before returning. Floating residuals are diagnostics checked against
explicit tolerances; the stored acceptance windows are not interval proofs.

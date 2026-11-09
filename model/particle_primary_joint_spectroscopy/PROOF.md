# Complete joint energy and spatial-phase spectroscopy

## One finite experiment on the primary source

The source is exactly the [primary member](../particle_common_contract/PROOF.md),
its stationary physical polynomial density rho, unchanged H, Gauss
constraints and canonical conjugate reference with idle Hamiltonian zero.
The [spatial Fourier instrument](../particle_primary_momentum/PROOF.md)
supplies every K_(b,k), including complementary branches. Its labels are
finite-aperture wavevectors, without an assumed conserved momentum PVM.
The instrument can use either actual finite filtered probes or their
separate unsmeared versions. These experiments are distinguished below.

This construction reuses the complete Fourier clock and summed Fejer
tail of [spectral-tube readout](../spin_shell_readout/PROOF.md) and its
`bounds.py`. That theorem's scalar kernel bounds are source independent.
Its periodic spin Hamiltonian, ground, momentum and particle assignments
are not imported. The complete
[source clock error](../particle_filtered_probe/PROOF.md) and
[transfer-to-excitation bound](../particle_primary_mass_readout/PROOF.md)
apply to the same primary H. A finite, error-certified branch comparison
does not require an infinite continuum theorem.

For J clock values and step tau define

```
F_a(H)=J^-1 sum_(r=0)^(J-1) exp(2pi i ar/J) exp(-ir tau H).
```

Finite Fourier orthogonality proves sum_a F_a*F_a=I. First record this
energy outcome, then apply the complete spatial instrument, then record
another energy outcome. All fine Kraus maps are

```
L_(a,b,k,c)=F_c(H) K_(b,k) F_a(H).                          (1)
```

Summing in reverse instrument order gives sum L*L=I. There are exactly
2M J² records, including unused bins, failed spatial probes and the
initial energy record. Every source operator preserves Gauss support.
Their neutral apparatus does not act on the retained reference.
The two energy registers each have dimension J and the spatial apparatus
has dimension 2M. These registers are separate from the two-bit compact
Ramsey apparatus; that service's carrier bound does not cover this
complete spectroscopy experiment.

In a common H/rho eigenbasis, their original probability is

```
p(a,b,k,c)=sum_(i,j) p_i |K_(b,k),ji|²
                     |F_a(E_i)|² |F_c(E_j)|².              (2)
```

This is a positive coupling of the exact stationary joint response to
two clock kernels. No success probability is divided out. The first
measurement leaves the source marginal rho unchanged because rho
commutes with every F_a. Its complete joint source/reference output
does change, and the initial pointer is retained. Equation(1) is not
declared equal to a Ramsey instrument or to a protocol omitting that
record. Source stationarity does not imply stationarity or unchanged
coherence of the purification.

## Signed energy labels, alias control and complete errors

Let R be the certified norm bound for the complete source H. Choose
tau=1/(4R). Use signed centers

```
e_a=2pi s_a/(J tau),   s_a=a if 2a<J, otherwise a-J.
```

The whole source spectrum has phase |tau E|<=1/4. For an energy
resolution eta with tau eta<=1, every real error |e_a-E|>=eta has
circular phase distance at least tau eta. Indeed,
wrapping a signed-bin separation greater than pi leaves distance at
least pi-1/4>1. Thus there is no alias contribution. No part of the
source spectrum or probability has been removed.

For address spacing h the spatial labels are q_k=2pi k/(Mh), modulo
2pi/h. Their spacing is 2pi/(Mh), an aperture convention rather than a
conserved momentum resolution. The inherited open-aperture identity
is K'_k-exp(iq_k h)K_k=exp(iq_k h)(B_M-B_0)/M, whose norm is at most
2/M. Dividing by h gives 2/(Mh), which need not shrink at fixed physical
aperture. Boundary mismatch is therefore explicit; closing the addresses
into a periodic source is not part of this experiment.

The inherited complete tail is

```
ell(J,d)=min(1,1/[2(Jd)²]+1/[2Jd]).                         (3)
```

Use the rational lower turn margin d=tau eta/7, since 2pi<7. Then a
complete energy readout has error probability at most ell(J,d), uniformly
in the true source energy. A power-of-two J>=1/(d beta) suffices for
ell<=beta<1. For both energy records, the event that either error
exceeds eta has original probability at most 2ell. More sharply, a
specified spatial channel of original total weight s loses at most
2ell s. The implementation uses the normalized cap s<=1.

Execute each base clock tick by the inherited r-substep Cayley service.
With this tau its operator error is at most 1/(768r²). Each controlled
clock power has error at most (J-1)/(768r²), so the complete clock
instrument has that half-trace bound for any reference. Composing the
two clocks gives

```
epsilon_clock <=min(1,2(J-1)/(768r²)).                       (4)
```

These rational approximations are compared to the unchanged H clock;
their Cayley eigenphases are not declared a replacement energy spectrum.
The spatial instrument is exact for its declared finite contractions.
Comparison with its ideal triangular filters adds the existing
sqrt(delta²+2delta) all-reference allowance, without an address-count
factor. Compiler, preparation and sampling errors require explicit
additional allowances. Failure maps are compared only for the same
instrument construction. No full primary histogram or pulse tape is
evaluated by these bounds.

## From actual joint records to excitation intervals

Suppose an original complete-record event, in a specified set of spatial
labels, has probability at least w and reports
e_c-e_a in [l,u]. Outside the two clock tails, its true transfer lies in
[l-2eta,u+2eta]. Let epsilon be the selected mean ground excess and r0>0.
The positive source coupling in the mass-readout theorem bounds the
original weight with initial excitation above r0 by epsilon/r0, since
the sum of the admitted spatial effects is at most I. Therefore

```
final excitation in [max(0,l-2eta),u+2eta+r0]
has original weight >=max(0,w-2ell-epsilon_clock-e_extra-epsilon/r0).
                                                                    (5)
```

This conditional adapter does not invent w or execute the histogram.
The supplied event must retain all other trials, and e_extra includes
its actual sampling/compilation allowance. No numerical E0 is needed
for (5). Alternatively, a source-certified ground interval [g_l,g_u]
and a final-energy record [a,b] give final excitation in
[max(0,a-eta-g_u),b+eta-g_l], with one-clock tail and its complete
event errors. A variational Rayleigh quotient supplies a ground upper
bound, not an excited pole or a ground lower bound. Producer validation
must identify the same full physical H and member, not just matching
string labels.

`ground_energy_bracket` regenerates the primary member and the complete
physical packet variational producer. Its Rayleigh quotient gives g_u.
The printed upper endpoint is rounded upward on a 2^-32 source-energy
grid, adding less than one grid step to that bound.
It uses either the declared full source norm floor or an authenticated
physical endpoint lower-bound receipt for g_l. Both bounds refer to the
same coupling-one H and full physical carrier. The producer's trial
does not replace rho, and its residual interval does not identify a
ground excitation. `observed_final_energy_to_excitation` subtracts this
verified bracket from the supplied final-energy event, keeping the final
clock tail, both implemented-clock errors and supplied event error.
It needs no initial cooling tail. The adapter computes no probability
from the interval supplied to it.

A resolved energy/wavevector ridge is an operational structure-factor
result. Before the mass-readout interval for E²-v²p² applies, the same
finite experiment additionally needs an attachment of its wavevector
windows to physical final-momentum intervals, initial-momentum control
where transfers are used, a common branch and a justified dispersion
error. Those finite-domain premises can be stated with bounded errors;
they do not require an infinite limit. Defining the detector, displaying
broad support or fitting a ridge alone proves none of them. The selected
source has no derived rest-energy or GeV mass in this package.

## A positive broad joint band on the actual source

There is a source result beyond defining the detector. For the separate
unsmeared experiment, the k=0 blank map is the Hermitian aperture mean
Abar=M^-1 sum_j A_j. Its distinct plaquette commutator vectors on the
physical blank are orthogonal by the inherited Haar-loop argument.
Consequently

```
||[H,Abar]b||²=K/(M a²),      K=4949/20.
```

The blank b is a matrix-element witness only; the input remains the
selected rho>=2^-N P_G. The original second transfer moment of Abar
is at least 2^-N K/(M a²). The actual polynomial density is decreasing
in energy. Since Abar is Hermitian, its original positive-frequency
weight dominates the reflected negative weight. All source transfers
lie in [-W,W], with W the certified full width. Put

```
delta=a^-1 2^-ceil(N/2),
w=[K/M-1] 2^-N/(2a²W²)>0.                                (6)
```

For the selected M=4, delta<W and the second-moment/passivity estimate
proves mu_(blank,k=0)([delta,W])>=w. The label k=0 is not identified
with a physical rest frame. This is an actual broad finite response
bound, without a resolved line, ridge or surviving physical band.

Choose eta=delta/4 and each clock's tail and complete implementation
allowance at most w/16. The recorded event with blank k=0 and transfer
label in [delta/2,W+delta/2] then has probability at least 3w/4: the
two tails cost w/8 and the two implementations cost w/8. Every other
record remains in the instrument. The code supplies this sufficient
plan without allocating the enormous registers. If w=c_w 2^-N and
eta=(1/(4a))2^-h, h=ceil(N/2), it takes

```
J=2^j,  j>=N+h+ceil(log2[16/(c_w tau/(28a))]),
r=2^s,  2s>=j+N+ceil(log2[1/(48c_w)]).                    (7)
```

The clock counts and precision exponents are exact symbolic resource
bounds. An uncounted compiler or empirical sampling error cannot be
added while preserving the same 3w/4 claim. These potentially enormous
costs are not an efficient or executed experiment. Along refinement
delta and w can vanish and W grows. Neither (6) nor (7) is a branch
resolution theorem or a nonzero limiting mass.

## Exact independent control

Take H=diag(0,1), rho=diag(3/4,1/4), two spatial probes
(sigma_x/2,-sigma_x/2), and the exact source step
U=diag(1,-i)=exp(-i pi H/2). Four energy bins have signed centers
(0,1,-2,-1) and resolve both source levels exactly. All 64 joint rows
are retained. Their four nonzero weights are the original upward record
3/16, downward record1/16, and complementary no-change records9/16 and
3/16. Their sum is one. Direct density-matrix evaluation and the positive
kernel law (2) agree. Initial energy measurement preserves rho but changes
the joint purification, which is checked explicitly. These matrices are
a control of the theorem and instrument, not the primary spatial source.

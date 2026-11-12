# Coherent compilation of the native lepton rest pairing

The source is the complete analyzer and its full fermion pairing V in
[lepton_mass_compiler](../lepton_mass_compiler/PROOF.md). Its physical
vacuum embedding J obeys J*VJ=H_Y, with Y the actual three-family analyzer
effect. The rest-site kinetic term is zero. Each site has 90 physical
fermion modes, 324000 enlarged modes and ||V||=6. The supplied neutral Higgs
background, analyzer phase 2/9 and retrospective family assignment are
unchanged. This construction changes the controller: it uses coherent
ancillary phase rotations and one final complete measurement, without
intermediate resets. Those operations are declared apparatus services.

## Exact selection on every Fock sector

Let N_A count every ancillary fermion and R=exp(i pi N_A/2). The vacuum
embedding satisfies RJ=J. Every creation-pair monomial of V has zero, one
or two ancillary creators, and annihilation terms have the opposite count.
The four conjugates V_r=R^(-r) V R^r therefore satisfy

    V_bar = (V_0+V_1+V_2+V_3)/4,
    V_bar J = J H_Y.

The finite root-of-unity sum removes both the one-ancillary and the
two-ancillary terms, leaving only physical pairs. This is an operator
identity on the entire enlarged Fock carrier. In particular, it controls
all physical occupations and arbitrary untouched references. Two phases
remove the odd terms but leave ancillary pair creation; the executable
retains this failed control. R is a phase on internal-representation
copies and commutes with every native internal current. On the neutral
Higgs background the retained conserved current is the unbroken electric
one, with the same boundary as the original source.

The parity S of the complete left lepton-doublet leg conjugates V to -V,
since every pairing contains exactly one fermion from that leg. It
commutes with R. Thus a negative signed V segment can be implemented by S,
positive-duration V evolution and S. No reverse physical time is required.
There is no claim that this conjugation reverses a nonzero kinetic term.

## A fourth-order finite product

Write A_r=V_r/4. Define the symmetric word

    S2(h)=exp(-ih A0/2) exp(-ih A1/2) exp(-ih A2/2)
          exp(-ih A3) exp(-ih A2/2) exp(-ih A1/2) exp(-ih A0/2).

For p=1/(4-cuberoot(4)), set

    S4(h)=S2(ph)^2 S2((1-4p)h) S2(ph)^2.

The executable proves equality of every noncommuting word through degree
four with exp(-ih V_bar): 341 word coefficients, reduced exactly modulo
60p^3-48p^2+12p-1. No matrix example substitutes for this identity. Changing
p to 2/5 violates it.

Since 1/4<p<5/12, the sum of the absolute source durations is
(8p-1)|h|<7|h|/3. On s rest sites the sum of segment norm bounds is at most
14s|h|. For 14s|h|<=1, the difference between the two Taylor tails is at
most 2 exp(14s|h|)(14s|h|)^5/5! <= (14s|h|)^5/20.
Unitary telescoping over n cycles, and the norm bound for two isometries
with arbitrary references, give the unhalved complete-output trace bound

    epsilon <= (14s)^5 |t|^5 / (10 n^4),   n >= 14s|t|.

Final ancillary measurement and retention of its quantum complement are
channels, so the same bound applies to their complete records. No vacuum
postselection is used. Ancillary leakage probability is also bounded by
the squared isometry error. The exact bound is independent of numerical
matrix exponentiation and of the physical input state.

For one rest site, t=1 and epsilon<=10^-12, n=15229 suffices: 533015
Hamiltonian segments and one initial ancillary vacuum, compared with
108000000000000 reset rounds under the original source's conservative
bound. These counts compare two declared protocols with the same target
and trace budget. Each segment additionally requires analyzer and phase
services. The total native interaction duration is (8p-1)t; neither that
resource duration nor the signed target time is identified with a passive
laboratory clock.

## Executed complete records and source masses

The source's exact rank-one active-span reduction gives four fermion modes
per spin-oriented family: two physical and two ancillary modes. On all
sixteen Fock states V^3=V, so each source segment has the exact polynomial
form I-i sin(h)V+(cos(h)-1)V^2. High-precision evaluation and repeated
squaring execute the entire cycle power, not a loop with omitted records.
All four physical input columns and an arbitrary retained reference are
controlled by the proof; the numerical check uses a maximally entangled
four-dimensional reference and retains every final ancillary output.
The source coefficients are taken from the midpoint of its rational
intervals. These high-precision checks are diagnostics, not outward-rounded
floating certificates; the complete-output theorem uses the exact source
coefficients and its exact norm bound.

A compensated Ramsey preparation uses the target ground and a single
charged excitation, with an idle opposite-charge reference. The four
complete effects (I+X)/4, (I-X)/4, (I+Y)/4, (I-Y)/4 vanish in X,Y on the
orthogonal complement, which is retained. The source's phase-readout lemma
gives energy error at most 4 epsilon at t=1. Applying its certified electron
denominator floor 1/20000 gives guarded mass-ratio and electron-calibrated
GeV intervals. No empirical muon or tau value selects this compiler.

The three target coefficients and their approximate calibrated values are
unchanged: 0.00051099895069, 0.105659414496 and 1.776984970813 GeV. The first
is an empirical anchor. Finite-shot noise, anchor uncertainty, a selected
Higgs vacuum, passive clock attachment, interacting lepton identification
and the physical muon discrepancy are separate from the compiler bound.
This rest-site result does not assert a kinetic propagation theorem or an
exact finite-cycle Hamiltonian spectrum.

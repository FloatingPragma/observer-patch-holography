# A complete attenuator for the coherent family action

The configured coherent analyzer of `model/standard_model_source_map/KOIDE_SOURCE.md`
can attenuate its two oriented family amplitudes while retaining every loss
record. This supplies an explicit action map, including its full fermionic
extension. It does not select the attenuation from the primitive C/E law.
The supplied analyzer, loss rotation, phase reference, blank flags and clock
are part of the enlarged apparatus.

## Source weights and complete dilation

The actual 128-label C/E instrument, including its 104 zero labels, gives the
five disjoint probabilities `(10,5,5,2,50)/72`. The parent coherent service
prepares their square roots, rather than measuring those labels and restoring
interference. In the parent's verified cyclic frame let

    R^3=I, R*=R^2, Tr R=0,
    a=sqrt(10/72), b=sqrt(5/72), |z|=1,
    C=aI+b(zR+conjugate(z)R^2).

The old work flag has eight coordinates, including three padding coordinates.
Add a two-state loss flag. On work flags j=1,2 only, rotate the loss coordinate
by `[[r,-ell],[ell,r]]`, where `0<=r<=1` and `ell=sqrt(1-r^2)`. On all other
work flags use identity. This defines a unitary L_r on all sixteen flags,
including nonblank loss inputs. Use the same controlled family operators on
both values of the loss flag. Apply the original three-point Fourier readout
on both loss sectors. With the old preparation V, Fourier F and controlled
family unitary B, the complete analyzer is

    U_r=(I2 tensor F tensor I3)(L_r tensor I3)
          (I2 tensor B)(I2 tensor V tensor I3).             (1)

All factors are unitary. The code includes the parent's fourth family padding
coordinate and checks the resulting 64-column matrix. The physical family
restriction has dimension 48. The sixteen blank-input Kraus maps are, for
m=0,1,2 and omega=exp(2pi i/3),

    K_(0,m)=[aI+r b(omega^(-m)zR+omega^(-2m)conjugate(z)R^2)]/sqrt(3),
    K_(1,m)=ell b(omega^(-m)zR+omega^(-2m)conjugate(z)R^2)/sqrt(3),
    K_(0,3)=sqrt(2/72)I,  K_(0,4)=sqrt(50/72)I,              (2)

with every other map zero. The loss outputs are Fourier-combined outputs;
they are not two separately measured oriented branches. Fourier orthogonality
and r^2+ell^2=1 give sum K* K=I on every family input, including references.
Deleting the loss rows violates this identity when r<1.

For a tracial family input each of the three pass probabilities is
`5(1+r^2)/108`, each of the three loss probabilities is `5(1-r^2)/108`,
and the other nonzero probabilities are `1/36,25/36`. Their sum is one.
The original total loss weight is `5(1-r^2)/36`. It is never renormalized away.
These formulas are not asserted for every input; its full record law is (2).

## The action map and its selection boundary

The selected output *effect*, used coherently inside a complete controller, is

    Y_r=K_(0,0)*K_(0,0)=C_r^2/3,
    C_r=aI+r b(zR+conjugate(z)R^2).                       (3)

It is not a success-conditioned mass. If C is positive, then
`C_r=(1-r)aI+rC` is positive. On the parent's positive chamber this gives

    Tr Y_r=5(1+r^2)/36,
    Q_r=Tr Y_r/(Tr sqrt(Y_r))^2=(1+r^2)/3.              (4)

For any declared chi>=0 set r=exp(-chi); then (3) is exactly the requested
bare-to-effective amplitude map. In the MCPR convention
`D=[I+sqrt(2)r cos(delta+2pi k/3)]^2`, the three coefficients of Y_r are
`(5/108)D_k`. Reproducing the declared MCPR common coefficient A additionally
requires the supplied common gain `108A/5`. Phase, attenuation and gain are
separate controls. The original primitive weights alone determine none of
those additional choices.

This is demonstrated without an inverse fit: the same primitive C/E labels,
geometry and coherent-service architecture admit r=1,3/5,0. They give
Q=2/3,34/75,1/3 respectively, with complete normalized instruments and the
same event counts. Thus those fixed data do not entail a unique attenuation.
A fully specified controller, including r, has a definite output law.
This counterfamily does not claim ambiguity of a fixed program or replace the
actual correlated-memory dynamics by an arbitrary attenuator.

In particular, evaluating the stipulated scalar fixed point
`chi=(u^2/512)/(1-u/77)` for `0<=u<77` supplies an input to (1), not a proof
that geometry or a finite closure theorem selects that equation. The
existing MCPR realization realizes chosen spectral coefficients; it is not
used to derive (1), nor does it select chi here. The actual default unread
memory diffusion is a different channel, as checked by the owning atlas
`COMPATIBILITY.md` and `memory_diffusion.py`.

## Native controls, errors and clock

Write the work-flag bits little endian as j0,j1,j2 and the loss bit as l.
The projector onto j=1 or2 is

    Pi=(I+Z_j2-Z_j0 Z_j1-Z_j0 Z_j1 Z_j2)/4.

For theta=acos(r), the added stage is `L_r=exp(-i theta Y_l Pi)`.
Its four commuting Pauli terms are compiled with the existing native
`observer.source_scalar.compiler.pauli_rotation`. This uses the actual
six-level source controls and their proper-code load/export services.
Every primitive retains its complement failure labels; on the invariant
proper code their exact weights are zero. No global phase is fitted out.
The executable control uses r=3/5 and tests every input column of this loss
stage. All zero-angle requests are retained in the r=0,1 count controls.

The exact unitary theorem allows real theta; it does not assert a finite
exact digital word for every transcendental exp(-chi). If the implemented
angle differs by at most e, the loss-unitary norm error is at most e, as is
the complete analyzer error; its compressed Y error is at most 2e. Certified
native primitive/service errors must be added separately. Floating matrix
comparisons are execution diagnostics, not interval certificates for e.

The parent analyzer has a compiled Fourier suffix. Inserting the loss stage
immediately before that suffix constructs a native tape for (1); the full
64-column replay is a separate bounded execution check. The original C/E
statistics, coherent preparation, phase controls, complete measurement and
fresh-flag services remain explicit. A native pulse count is not a physical
time normalization. Repeating a fixed-time analyzer at arbitrarily small
weak-kick angle does not make its service overhead vanish.

## Complete weak feedback and full CAR

Let Q project onto output loss0/work0, including all three family coordinates,
and put P_r=U_r* Q U_r. It is a projection of rank three. Apply the parent
logical weak controller with P_r: prepare the blank flags, apply U_r, kick
`exp(-ih X tensor Q)`, uncompute, and record all sixteen return flags. Then

    A0=I+(cos h-1)I tensor Y_r-i sin h X tensor Y_r,
    sum_(j!=0) A_j* A_j=4sin^2(h/2)I tensor(Y_r-Y_r^2)<=h^2 I/4. (5)

The loss bit remains coherent inside U_r and is read only after uncomputation.
The standalone analyzer probabilities following (2) are therefore distinct
from the return probabilities in (5). Measuring the intermediate loss bit
would implement a different controller. This remains a complete instrument
on all target/reference inputs. The
parent contraction telescope, with N rounds and h=t/N, gives unhalved
complete recorded-output distance at most `9t^2/(4N)` from evolution by
`X tensor Y_r` with the all-zero archive, and nonzero-archive weight at most
`t^2/(4N)`. A per-round certified dilation error e adds `2Ne` and `Ne`
respectively. Every other outcome is retained; no conditional probability
or renormalized generator enters this limit. The model clock t is the
accumulated supplied kick angle, not the total native service duration.

The full-CAR theorem in `model/charged_lepton_response/CAR_KOIDE.md` depends
only on a full analyzer unitary and its compressed rank-three projection,
not on eight rather than sixteen flags. Extend the ninety physical modes
to `C2_spin tensor M15 tensor C16_flag tensor C3_family`: 1,440 modes,
including 1,350 fresh ancillary modes. The flag coordinate is gauge neutral;
its fermionic copies carry the ordinary M15 charges. Apply U_r on the entire
lepton doublet and conjugate(U_r) on its singlet conjugate, and identity
elsewhere. Both support projectors commute with all twelve native currents.
The exterior lift is therefore a full-current-preserving unitary.

Conjugating the existing rank-one charged Higgs creation pairing gives
oriented coefficient U_r* Q U_r and its transpose. Compression by the
ancillary vacuum gives Y_r and Y_r^T on every physical occupation sector.
This is the same CAR functor identity as the owning theorem; it never uses
Gamma(sum A)=sum Gamma(A). Its kick norm is `L=6|v_H|`, since Q
has rank three, and its compressed pairing norm is
`2|v_H|Tr Y_r=5(1+r^2)|v_H|/18`.

Retain the complete quantum nonvacuum outputs with invariant vacuum/nonvacuum
flags and supply fresh ancillary vacua each round. The owning full-CAR bound
is then `3L^2t^2/N=108|v_H|^2t^2/N` in unhalved complete trace distance,
with original nonvacuum probability at most `36|v_H|^2t^2/N`. A per-round
dilation error e adds `2Ne` and `Ne`. Before classical flag commitment the
coherent-archive bound is the weaker `sqrt(2)L|t|/sqrt(N)` in isometry norm.
A common H0 requires the invariant-vacuum and conjugation hypotheses stated
in CAR_KOIDE; no arbitrary interacting H0 is imported. The 1,440-mode native
execution, physical clock selection and measured pole-mass map are not
asserted. The finite controller is a constructive bare-to-effective action
map, with its required apparatus and unresolved source selection explicit.

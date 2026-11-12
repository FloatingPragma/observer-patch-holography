# Double-cover rotation and odd-composite exchange

Use four weak-doublet modes indexed `(spin, weak component)` and two
charge-six singlet modes indexed by spin. The corresponding native
central charges are `(-3,-3,-3,-3,6,6)`. Occupying two weak modes and one
singlet mode gives `binom(4,2)*2=12` basis states, all centrally neutral.
The colour generators vanish because both field representations are
colour singlets. The three weak generators are the exterior lifts of

```text
g_a = diag(I_spin tensor sigma_a/2, 0_2).
```

Their simultaneous kernel has dimension six. The script computes its
basis rather than installing six basis vectors. The spin coefficient
generators are the exterior lifts of

```text
j_a = diag(sigma_a/2 tensor I_weak, sigma_a/2).
```

If `B` is the unnormalized kernel basis and `G=B*B` its positive Gram
matrix, the reduced matrices are `J_a=G^-1 B* dGamma(j_a) B`. The script
checks invariance, all three commutators and Gram-adjointness. Their
Casimir has eigenvalues `3/4` with multiplicity two and `15/4` with
multiplicity four. The two exact orthogonal projectors are

```text
P_(1/2)=(15I/4-J^2)/3,       P_(3/2)=I-P_(1/2).
```

For each of the 48 binary-cube quaternions, the elementary spinor is
`U(q)=q0 I-i sum_a q_a sigma_a`. The actual three-creation exterior
matrix is computed from its `3x3` minors and restricted to the gauge
kernel. Unitarity uses the physical Gram matrix. The traces on the two
Casimir sectors are `2q0` and `(2q0)^3-4q0`, the complete `G1` and `H`
characters. All 48 spatial rotations permute the eight cube sites, so a
uniform orbital preparation retains these labels.

Four successive quarter turns give `-I_6` and eight give `I_6`, including
the rank-two `G1` component. The sign is not removed as an irrelevant
global phase: it distinguishes two coherent rotation paths in an
interference experiment, whose physical implementation is outside this
script.

## Complete disjoint exchange test

Each kernel column becomes a polynomial of degree three in creators,
with its exact exterior coefficients. Copies at distinct sites act on
twelve active modes. The program applies both orders of all 36 column
pairs to all `2^12=4096` occupation states, yielding 147,456 exact basis
comparisons, and verifies

```text
F_(x,v) F_(y,w) = -F_(y,w) F_(x,v),       x != y.
```

At least one product must be nonzero, excluding a vacuous check. The
algebraic extension is the `(-1)^(3*3)` sign from nine interchanges of
disjoint CAR factors; it holds with spectator modes as well. The script
does not assert canonical free CAR for composite operators at the same
site. The elementary CAR structure is an input to this result.

| Equation or check | Implementation |
| --- | --- |
| 12 occupation masks | [spin.py:29](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L29), `basis_masks` |
| Literal exterior signs and bilinear lift | [spin.py:35](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L35), `apply_word`; [spin.py:47](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L47), `one_body_lift` |
| Gauge kernel, spin and Casimir | [spin.py:68](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L68), `local_kernel` |
| All 48 quaternions and multiplication | [spin.py:98](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L98), `quaternion_elements`; [spin.py:116](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L116), `quaternion_product` |
| Spinor and three-particle exterior rotation | [spin.py:123](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L123), `spinor`; [spin.py:127](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L127), `exterior_rotation` |
| Full and double turn | [spin.py:143](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L143), `verify_turns` |
| Full two-site exchange | [spin.py:166](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L166), `verify_exchange` |
| Physical cube orbit and character restrictions | [spin.py:182](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L182), `_uniform_orbit_check`; [spin.py:193](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/spin.py#L193), `run` |

**Scope**: this calculation constructs the weak-epsilon, charge-six composite
kernel in an exact canonical chart, checks its complete binary-cube
rotation labels, and executes disjoint odd exchange on every state of its
two-site active-mode Fock carrier. The spinor action and canonical
anticommutation relations are declared inputs. This is no derivation of
fermionic statistics from a classical field or of the relativistic
spin-statistics theorem.

The source is public RER commit
`3bcacf1276ddb3ead35be2b9adcd63b6643e2d35`,
[`model/particle_primary_composite_spin/PROOF.md`](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/particle_primary_composite_spin/PROOF.md).
The portable calculation retains its positive native increment
`(0,2,0,0,1)`. It does not execute the full primary Hamiltonian.

The optional pinned-source audit maps the six portable modes to the
actual positive-beta source modes `[6,10,21,25,11,26]`, with the original
filled negative-beta sea. Its exact signed exterior map intertwines all
three spin generators and maps the whole kernel projector to the source
projector; all twelve native Gauss generators annihilate it. Reversing
one exterior sign makes the native Gauss check fail. This checks the
chart equivalence directly rather than inferring it from dimensions.


## Inputs, negative controls and omissions

The native representations, spinor action, CAR quantization, selected
occupation sector and uniform orbital preparation are supplied. There is
no measured mass, fitted angle or empirical normalization. The finite
cubic labels do not uniquely determine continuum spin. A stable electron,
an isolated mass pole, Lorentz symmetry and a laboratory spin-interference
readout are not established. These omissions keep `full_criterion_met`
false for OpenWave row 7.

The negative controls erase the central sign by multiplying the
quarter-turn matrix by `exp(i*pi/4)`, and replace anticommuting creators
by commuting hard-core creators. The same turn and exchange checks reject
them. Both failures are required for a successful result.

```sh
python -m oph_openwave.spin
python -m unittest discover -s tests -p test_spin.py
python audits/check_spin_source_chart.py --rer /path/to/pinned/RER
```

The exact computed matrices, all characters and the full comparison count
are in the emitted payload. The [independent adversarial review](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/evidence/charge_spin_review.json) reproduces the Casimir and all projected characters using a vector-tensor-spinor construction, then checks all 147,456 exchanges through permutation inversions. All three broken controls are detected; [the audit script](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/audits/charge_spin_review.py) is executable. The imported spinor and CAR premises retain their qualification.

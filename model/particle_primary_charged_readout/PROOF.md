# Charge transfer and energy on one complete source history

The selected member is the eight-cell open Wilson/Gauss source. Its
polynomial density, idle conjugate reference, interaction, source clock and
energy origin are unchanged. The
[charged endpoint construction](../particle_primary_charged_string/PROOF.md)
supplies 48 link-dressed successes and one failure, with a populated
original energy band. A coherent finite counter attaches an operational
charge-transfer record to that same experiment without an initial charge
measurement. This is a record of charge moved by the intervention; it is
not a single-particle charge or mass identification.

## 1. Actual local charges and a finite meter

Let the five native integer hypercharges and internal dimensions be

\[
 (q_j,d_j)=(-4,3),(-3,2),(1,6),(2,3),(6,1).
\]

Each block has four Dirac rows, so its cell occupation satisfies
\(0\le N_{j,z}\le4d_j\). The local matter charge
\[
                  Q_z=\sum_jq_jN_{j,z}
\]
has integer spectrum in \([-72,72]\). The native beta sea has charge zero.
These number operators commute with all local gauge actions and with the
product link projection. Their spectral projectors therefore act on the
same physical Gauss carrier. This is the native integer charge convention;
no electroweak vacuum or laboratory electromagnetic identification is used.

Each meter has cyclic radix \(M=5\). The admitted probe transfers at
one cell are exactly \(-6,0,+6\); their residues are \(4,0,1\).
These three values are distinguished even though the unknown initial
charge is not determined modulo any larger register. Radices one through
four collide on this alphabet, so five is the smallest ordinary cyclic
native-number-addition meter for these transfers. This is not a minimum
over arbitrary quantum encodings.

There is one meter per cell, all prepared in \(|0\rangle\), giving
meter dimension \(5^8=390625\). Their Hamiltonian is zero. Every meter
outcome is retained; values two and three have zero weight in the exact
experiment. Decoding depends on the proved transfer alphabet. A different
intervention with a shift such as +1 would alias with +6, so this meter is
not a readout of arbitrary charge differences in \([-144,144]\).

Write \(S\) for the cyclic shift. The physical controlled-number unitary
is
\[
 W=\sum_{\mathbf q}P_{\mathbf q}\otimes
                         \bigotimes_z S_z^{q_z}.                 \tag{1}
\]
The joint projectors are complete; hence \(W\) is unitary and commutes
with the source Gauss action. The finite cyclic-pointer construction is
the same one used in the
[physical boundary-state readout](../particle_physical_boundary_state/PROOF.md).
An explicit implementation factors (1) into one controlled addition of
\(q_j\) per occupied native CAR mode. There are 480 such controlled shifts
per coupling, and 960 for the two couplings below. These counts specify
the apparatus gates, not an executed large-source simulation or a cost
bound for the underlying charged probe and energy clocks.

## 2. Coherent cancellation of the initial charge

For a directed bond from \(y\) to \(x\), the success operator is
\[
 K_b=B_{xy\alpha}/\sqrt{48},\qquad
 B_{xy\alpha}=c^\dagger(p_{x\alpha})c(m_{y\alpha})
                                  \chi_6(U_{xy})P_{0,xy}.
\]
The reverse bond has the conjugate character. On the full physical source,
\[
 [Q_z,K_b]=d_{b,z}K_b,\qquad
 d_{b,z}=6(\delta_{zx}-\delta_{zy}).                       \tag{2}
\]
The complementary failure is
\(K_f=(I-\sum_bK_b^\dagger K_b)^{1/2}\). Every success effect commutes
with every \(Q_z\), so \([Q_z,K_f]=0\) by functional calculus. Set
\(\mathbf d_f=0\). All 49 operators consequently have definite transfer
vectors, and their effects sum to the identity.

Couple the blank meters with \(W^\dagger\), apply the unchanged complete
probe instrument, then couple with \(W\). No charge or meter is read
between the couplings, and no source evolution is silently inserted there.
For every branch, (2) gives the exact operator identity
\[
 W(K_b\otimes I)W^\dagger
       =K_b\otimes\bigotimes_zS_z^{d_{b,z}}.              \tag{3}
\]
For example, a source charge component \(\mathbf q\) first writes
\(-\mathbf q\) coherently. The branch changes its source charge to
\(\mathbf q+\mathbf d_b\). The second coupling removes the unknown
initial component and leaves the meter at \(\mathbf d_b\), independently
of that component. Linear superpositions and all source-reference
coherences obey the same identity. The meter is measured only after this
cancellation.

If the original probe is retained coherently as an isometry with label b,
the enlarged output is exactly its isometric embedding
\(|b\rangle|0\rangle\mapsto|b\rangle|\mathbf d_b\rangle\).
This copies a function of an orthogonal pointer label, not an arbitrary
quantum state. On classical complete records the source/reference maps
are exactly the original maps with the additional deterministic transfer
label. Thus the counter does not add independent information about the
unknown input charge; it realizes the operator charge-transfer identity
as a physical source coupling and a readable record.

Reading the first meter would instead dephase the unknown initial charge.
It is a different apparatus. No transfer of an energy bound to that
dephased preparation is made.

## 3. Same-clock energy records and an optional terminal charge readout

Let \(F_i(H)\) and \(F_j(H)\) be the two complete finite source-energy
instruments from the
[joint spectroscopy construction](../particle_primary_joint_spectroscopy/PROOF.md).
Place the two coherent couplings around the probe, after the first energy
clock and before the second. Equation (3) gives every fine-history Kraus
operator:
\[
 L_{i,b,j,\mathbf r}=
        1_{\mathbf r=\mathbf d_b}\,F_j(H)K_bF_i(H).        \tag{4}
\]
All meter labels, both clocks, probe failures and the idle reference are
retained. The formal meter-refined alphabet contains
\(49J^2\,5^8\) outcomes; (4) identifies its exact zero-probability
outcomes. The nonzero transfer labels sum to zero across the cube and
carry \(+6,-6\) at opposite endpoints. They obey Gauss constraints at
the intervention, without a nonneutral ambient-state replacement.

The original charged-response producer is regenerated, with its actual
source/state/reference/clock identifiers. The positive event includes one
of the 48 success labels as well as the stated energy band. The exact
success-summed energy measure has
weight above \(3/10\) in \([39/4,135/4]\). Finite clock tails and
implementation bounds give original reported transfer weight greater
than \(1/4\) in \([77/8,271/8]\). Equation (4) proves that every record
in that joint success/energy event also has the operational opposite-charge
transfer label,
with exactly the same original weight. No probability is divided by
successful charge or energy selection. Failure records are retained, and
their energy transfers are not included in this success-summed event.
The inherited
`Lean/QFT/GroundRecordTail.lean` consumer converts the same-source
positive transfer to final ground-relative excitation with its cooling
allowance; the coherent meter changes none of its premises.

A complete absolute-charge readout can be appended after the final energy
record. Its additional maps are \(P_{\mathbf q}L_{i,b,j,\mathbf r}\),
with every one of the \(145^8\) charge tuples retained. Since
\(\sum_{\mathbf q}P_{\mathbf q}=I\), summing the terminal charge record
preserves every earlier joint energy/probe/transfer probability, for an
arbitrary quantum reference. The postmeasurement source/reference state
has the stated charge-readout disturbance.

This terminal absolute charge need not equal the transfer label. The
initial absolute charge was not measured, and the intervening energy
instrument does not commute with local charge. The counter records charge
moved by the probe event, not an assertion that the final energy-resolved
state has permanently localized endpoint charges. Measuring two absolute
charges without retaining their disturbance is not an alternative proof
of (4).

## 4. Finite errors and independent quantum controls

Suppose each cell coupling approximates its exact unitary with operator
error at most e. The product coupling error is at most \(8e\); the two
couplings give total isometry error at most \(16e\). The complete output
half trace distance is therefore at most \(\min(1,16e)\), uniformly for
the full source and any retained reference. Unitary gates, energy
instruments and terminal measurements cannot increase this distance.
Subtract this event allowance from the inherited clock-corrected band
weight. It has no factor equal to the number of outcomes or the inverse
success probability. Exact controlled shifts have e=0 in the declared
finite architecture.

The executable small control uses three physical neutral states:
a filled-minus pair with constant link, a link-dressed pair with charges
\((+6,-6)\), and its reversed pair with charges \((-6,+6)\). Their Gauss operators vanish on the explicitly
tested CAR/link embedding. On these states take
\[
 H_c=\begin{pmatrix}0&2&2\\2&2&0\\2&0&2\end{pmatrix},
 \quad Q_c=\operatorname{diag}(0,6,-6),\quad
 \rho_c=(5I-H_c)^2/59.
\]
The clock spectrum is \((-2,2,4)\), the density is stationary and
faithful, and \([H_c,Q_c]\ne0\). This is a control of the readout
identities, not the selected primary source. Its two raising operators
divided by \(\sqrt2\) and the excited-space failure are a complete
instrument. A radix-five coherent meter executes (3), including every
zero and failure record. Exact two-energy and terminal-charge records
check (4) with the canonical entangled reference. Reading the meter
prematurely changes the energy histogram, and the final absolute charge
can differ from the probe's transfer label. Those countercontrols
prevent an accidental nondisturbing two-charge interpretation.

The primary full-source probability bound uses the actual authenticated
producer and operator identity (3); it does not use control frequencies
as primary data. No large-source histogram, isolated charged particle,
spacetime spin, scattering law or physical rest mass is calculated here.

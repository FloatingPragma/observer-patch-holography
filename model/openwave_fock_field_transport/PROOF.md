# Complete Fock transport and the dynamical field connection

## An exact map on the declared native source

Take the complete 192-mode CAR algebra and its Fock representation of
`openwave_common_mass_source`. The 144 colourful modes are unchanged.
Pair each left charged coordinate with its conjugate coordinate of opposite
spin and the same cell and family; do the same for each neutral coordinate
and its conjugate singlet. This partitions the 48 lepton modes into 24 pairs.
All family mixing remains in the Hamiltonian. It is not suppressed to make
the map diagonal.

For a pair (a,b) and its spin sign sigma, set

    c_a=(a-sigma b*)/sqrt(2),  c_b=(b+sigma a*)/sqrt(2).

On the ordered complete Fock basis (empty,a,b,ab), the even-sector rotation
with entries c=s=1/sqrt(2), and identity on the two odd states, is an explicit
unitary implementer. The executable checks its CAR identities entry by
entry. Changing sigma changes the off-diagonal signs. Implementers of
disjoint pairs are even and commute, so their product implements the map
on the entire 2^192-dimensional Fock space. No enumeration or occupation
restriction is needed for this finite-algebra theorem.

The exporter transports every actual linear or bilinear factor of all 114
native source terms. It checks the exact inverse on each supplied factor
using rational arithmetic and sqrt(2), with full CAR normal ordering.
Products, positive squares, scalar weights and their order are retained.
The receipt exports every canonical normal-ordered factor coefficient,
in addition to its hash and inverse check, so the complete finite polynomial
can be reconstructed without an implicit effective-Hamiltonian replacement.
The colourful factors are fixed. Each of the 24 actual semileptonic lepton
currents has nonzero pair terms in its image; their omission changes the
source. The supplied floating sparse snapshots are checked as exact dyadic
coefficients for this identity, without identifying their rounded values
with the analytic calibrated constants. The algebraic substitution theorem
also applies to the stated analytic coefficients themselves.

This map gives a canonical CAR polynomial for the same interacting finite
source. It does not replace that polynomial by a free diagonal Hamiltonian.
In particular, the neutral family mixing and all conjugate sectors remain.

## Charge, the longitudinal field, reference and work

Both modes in each charged pair have opposite electric charge; neutral
pairs have charge zero. Thus each local native charge satisfies
Q_x(c)=Q_x(a), with no reference constant. The producer checks this exact
identity over all 192 modes. Consequently the pair implementer commutes
with the full longitudinal Gauss dressing in `openwave_common_electromagnetism`.
For each native charge occupation, that dressing has the unique labels

    R=-3 sum_x Q_x,  E_x=3 V Q_x+R,  sum_x E_x=0.

Every matter occupation has integral labels. All integer rotor and quantum
reference spaces remain in the kinematic carrier. The physical constraint
kernel is unitarily equivalent to the entire native Fock space. Its
orthogonal constraint-failure projector remains an outcome with zero
probability on a physical preparation. H, each dressed current, the electric
square and its declared self-energy subtraction all transform by the same
unitary. Because the native pair transform is field independent, its
connection is exactly zero; there is no omitted derivative term in this
particular longitudinal construction.

The pure Dirac reference becomes its exact unitary image, including neutral
flavour mixing. It is not replaced by a convenient bare vacuum. If K_r is
any complete instrument and rho a state with any external quantum reference,
conjugate H, rho and every K_r by the same unitary tensored with the identity
on that reference. Then each full joint output, outcome probability and
mean work Tr(H(Phi(rho)-rho)) is identical. Completeness follows from
sum K_r* K_r=I. The executable checks a complete four-outcome instrument,
a nonzero exact work value and an entangled-reference joint output. This
is a direct use of the finite complete-instrument theorem, without a
postselection or discarded spectator sector.

## Full-Fock lift of the field-dependent rotor diagnostic

The signed-band diagnostic `openwave_dirac_field_connection` has a different
carrier: an independent full L2(U1) rotor. Its complete two-mode Fock lift is

    H=alpha E^2 + m(N_a+N_b-1)
      + t(exp(-iq theta)a* b*+exp(iq theta)b a) + epsilon I,
    epsilon=sqrt(m^2+t^2).

Here alpha=1/32 and t=1/16 are the inherited diagnostic settings. The two
odd states have energy epsilon before the electric term; they are retained.
The empty/pair block is exactly the signed diagnostic, with pair flux n
and empty flux n+q for every integer n. Every odd state at every integer
flux is also retained. Its energies are alpha n^2+epsilon. The complete
carrier is the orthogonal direct sum of these even two-dimensional blocks
and odd one-dimensional blocks.

Let c^2=(1+m/epsilon)/2 and s=t/(2 epsilon c). The full Fock implementer U
has U_00=U_33=c, U_03=s exp(iq theta), U_30=-s exp(-iq theta), and identity
on the odd subspace. Its exact connection is

    A=-i U* U'
      =q[s^2(N_a+N_b-1)
         +cs(exp(-iq theta)a* b*+exp(iq theta)b a)].

Thus U*HU=alpha(E+A)^2+epsilon(N_a+N_b). The electric square contains
-i A', -2i A derivative_theta and A^2. In particular
A^2=q^2 s^2(N_a+N_b-1)^2. Smooth periodic U preserves the Sobolev domains
of E and E^2, so these identities hold for the full self-adjoint rotor
operators, beyond a finite Fourier test. The executable independently
checks all four Fock columns and the separate differential terms at the
exact rational control c=4/5,s=3/5,q=3. The pair-creation component of A
on the empty state has squared norm 1296/625, so omitting it fails exactly.

For all three charged inherited mass parameters and both rotor charges,
the producer checks exact rational enclosures of every even block with
-16<=n<=16 and retains both odd modes. Every unlisted block has energy
at least (1/32)14^2. This exceeds the enclosed ground energy, proving that
the reported ground interval concerns the full integer-flux operator.
The complete original-input records retain preparation failure, both even
bands, both odd sectors and every other flux block. Initial preparation
work uses this full-source ground interval; terminal work is retained.
The displayed propagation probabilities are floating diagnostics, with
the analytic reducing-block and completeness identities stated separately.

## The exact transfer obstruction

The longitudinal source has a finite physical carrier of dimension 2^192.
The independent-rotor diagnostic has a countably infinite physical carrier.
There is no unitary onto both complete carriers. Its addition therefore
requires an explicitly extended native architecture; the pair map alone
cannot turn the existing longitudinal source into that rotor source.

On each diagnostic even block the connection-omitted field also differs by
the nonzero off-diagonal alpha cs((n+q)^2-n^2). At n=0 its square is

    alpha^2 q^4 t^2 / (4(m^2+t^2)) > 0.

The executable gives this exact number for every inherited charged family
and both signs. The positive intraband Pauli-root member does not become the
signed-Dirac field conjugate by dropping these terms. Conversely this control
does not prove that no enlarged native construction can implement a
relativistic particle theory.

All inherited longitudinal source energies and calibration ancestry are
unchanged by the exact native unitary. No additional interacting pole or
physical GeV mass is identified by the diagnostic. A common transverse
interacting source and its continuum/particle transfer remain separate
mathematical requirements.

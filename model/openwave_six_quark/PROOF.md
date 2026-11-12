# Complete native six-quark binding and intrinsic graph quadrupole

## One declared current source and its selection

Use the full sixty-CAR two-cell native matter carrier. Its active quark
module is the actual Q weak doublet with three colours and two spins at each
cell. All conjugate quarks, leptons and complementary occupations remain in
the source. Fixed species occupations are invariant. On the active module,

    H=mN+t sum_(edge,a)(c_0a-c_1a)*(c_0a-c_1a)
        +sum_x C_x²+(1/4)sum_x S_x²+(1/4)sum_x I_x²,
    m=1, t=1/8.

C is the actual native colour current with its orthonormal Cartan metric;
S is spin and I is the native weak current restricted to colourful Q. Leptons
do not acquire this extra I cost. The full action is positive and H>=mN.
Equal spin and weak coefficients impose a declared spin/weak exchange
symmetry;1/4 is inherited from the existing spin cost. This equality is an
architectural input, not a derived physical coupling or measured nuclear fit.
The member without the I cost is retained as a control. No composite-state
projector is inserted in either Hamiltonian.

The native colour blocks on both Q components are exactly equal. The code
checks all81 entries of their Fierz identity

    sum_a(T_a)_ij(T_a)_kl=(delta_il delta_jk)/2-(delta_ij delta_kl)/6.

With integer CAR columns X_ab=E_ab+E_ba, Y_ab=E_ab-E_ba,
D1=n_r-n_g and D8=n_r+n_g-2n_b, the full occupation matrix is exactly

    48 C_x²=12sum_(a<b)(X_ab²-Y_ab²)+12D1²+4D8².

For spin or weak Z=2Jz and its raising column,
48(1/4)J²=3(Z²+2Z)+12Jplus*Jplus. Mass and positive hopping also have
integer numerators over48. The complete11700-word N=6,Nu=3,Nspin-up=4
sector is therefore represented by an integer matrix divided by48. Every
fermionic sign is obtained from occupation-word parity.

## Exact complete reduction to twenty-seven states

A colour singlet has two quarks of each colour. Restriction to that exact
Cartan weight leaves1800 occupation words and is invariant. Let C3=3C²,
S2 be total spin Casimir and I2 total weak Casimir. On this sector form

    P_colour=product_(c=3,6,8,12,18)(1-C²/c),
    P_spin=(S2-6)(S2-12)/40,
    P_weak=product_(c=2,6,12)(1-I2/c).

Their product is computed in exact integer sparse arithmetic, reducing a
common rational denominator after each multiplication. The final denominator
is2880. Every multiplication checks a conservative integer-overflow bound.
Exact matrix identities verify P*=P, P²=P, C²P=0, S2P=2P, I2P=0 and HP=PH.
Conversely every colour-zero,S=1,I=0 vector is fixed by each factor, so this
is the entire target sector, not a selected subset. Its exact trace is27.

Twenty-seven independent integer columns of P have a positive exact LDL
Gram decomposition. They therefore span the complete target sector. This
reduces the11700-dimensional source to an exact27-dimensional Hermitian
matrix whose entries are rational numbers divided by square roots of rational
Gram factors. Both the balanced and no-weak-cost actions use the same exact
space. The omitted orthogonal sectors are retained in the source and in the
complete failure instrument.

## Full spectral certificate and thresholds

For each exact reduced matrix the producer encloses all algebraic entries,
rounds a proposed eigenframe and eigenvalues to dyadics, and computes its
Gram and residual products with exact integers. Let eta bound the eigenframe
Gram defect, rho its residual and K the matrix/eigenvalue norm. Polar
unitarization and the self-adjoint spectral comparison put every ordered
true eigenvalue within2rho+4K eta of its corresponding dyadic proposal.
Outward algebraic-entry errors are included in rho and K. Thus the complete
spectrum, rather than an eigensolver's selected low value, is certified.

The same balanced baryon source supplies its complete120-dimensional
certificate. Its minimum exceeds3.733595788226506. An exact normalized
integer six-quark trial has Rayleigh energy about7.220709704169545. Hence

    2E_baryon,min-E_sixquark,min >0.24648187228.

For a colour-nonsinglet three-quark fragment, local3+0 occupation has a
nontrivial colour cost at least3, while local2+1 occupation costs at least8/3.
Thus every nonsinglet three-quark fragment lies at least3+8/3. The complete
120-dimensional colour-singlet baryon minimum is therefore the actual lowest
three-quark fragment energy. A colour-zero pair of nonsinglet fragments has
conjugate nontrivial colours and also lies above the six-quark trial.

Every nontrivial partition other than3+3 includes compensating
nonzero-triality fragments. Their colour Casimir costs are at least4/3 each,
so the threshold is at least6+8/3. Spin, weak and kinetic costs are positive.
The six-quark energy lies below this threshold as well. It is therefore a
stationary bound state relative to every quark-fragment partition in this
fixed-source sector. Each fragment uses a separate noninteracting copy of
the same finite action; no asymptotic continuum scattering theorem is implied.
For perturbations preserving the admitted occupation sectors, uniform
generator perturbations of norm epsilon in the six-quark and each
baryon copy preserve the gap when3epsilon is smaller than the certified
binding margin.

The same complete spectral enclosure proves the unbalanced S=1,I=0 minimum
is above two unbalanced baryons by more0.00099626457. It is an unbound control
in this target channel. It is not a statement that every other six-quark
channel of that action is unbound. Bare unequal up/down number energies add
the same3u+3d scalar to the six-quark state and a proton/neutron pair; this
common shift cancels in their binding comparison.

## Intrinsic charge quadrupole

Take graph coordinates z=0,1 and the number centre X_N=(sum z_i N_i)/6.
The complete charge observable is

    Qzz=2sum_i q_i(z_i-X_N)².

On an I=0 state every isovector charge expectation multiplied by a number
scalar vanishes. Since q=1/6+I3 on Q, the exact expectation identity is

    <Qzz>=<N0 N1>/18,   0<=<Qzz><=1/2.

This is an intrinsic finite-graph charge moment. It does not identify a
physical length, a three-dimensional nuclear shape or spatial Noether spin.

Let P be the co-location subspace N0=0 or6 within the target sector. There,
PHP=29/4. Positivity gives QHQ>=6 for Q=1-P, and the hopping norm gives
||PHQ||<=6t=3/4. If E<=E_trial<29/4 and p=<Q>, then

    E>=29/4-(5/4)p-(3/2)sqrt(p)
      >=29/4-(11/4)sqrt(p).

Thus every ground vector satisfies

    <Qzz> >=(5/18)[(29/4-E_trial)/(11/4)]² >0.

An exact integer trial also gives its quadrupole and residual. The enclosed
next eigenvalue separates the ground line; the residual/gap bound transfers
the full state and all bounded observables to that exact stationary line.
The conditional graph quadrupole is approximately0.00263818496 and the
trial has approximately99.0634percent co-location. It is a tightly clustered
six-quark candidate, not an established extended physical deuteron.

## Complete records, additional partner and calibration

The admitted input is the maximally mixed complete11700-word sector, with an
idle compensating charge−1 reference. The exact rational rank-one trial
filter and its complement have success1/11700 and failure11699/11700.
Every success and failure occupation probability is retained, including
zero rows. The full-Fock complement has an identity failure outcome and zero
probability for this preparation. Since the initial state is maximally mixed,
the complete nonselective filter does zero mean work; its preparation energy
is explicitly recorded. The trial is labelled variational. Its trace-distance
bound to the stationary ground is uniform under the source's unitary evolution
for every time, and transfers to arbitrary complete observer records.

Spin/weak exchange maps S=1,I=0 to a degenerate S=0,I=1 partner. Its charge
triplet0,1,2 is retained. These extra finite-source states are not declared
new observed particles. One unchanged electron unit displays a conditional
six-quark stationary energy and binding energy; it supplies no physical
strong scale, deuteron mass or length. All absolute mass and species
identifications remain separate from the finite binding theorem.

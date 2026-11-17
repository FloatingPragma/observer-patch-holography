# CCF-2 signed lepton kinetics: exact construction and finite witness

The full generator in `dirac.py::dirac_terms` returns native CAR/link
monomials for all charged and neutral leptons, both spins and all three
families. It uses the common 96-mode-per-vertex convention and is consumed
by the common-source assembler. The signed Dirac/Wilson operator acts
on the full Nambu spinor, including pair terms. [RESULTS.json](RESULTS.json)
contains its executed algebra, finite-carrier controls and source attachment.

## Algebra and charge

The exact pair map is a CAR Bogoliubov transformation on each complete
four-state pair; disjoint even pair transformations commute. It preserves
`-3 N_e + 3 N_p` exactly. The neutral map has zero electromagnetic charge.
The lower spinor components include the antisymmetric spin-conjugation
sign, rather than interpreting a raw antiparticle index as particle spin.

Every component of a Nambu column changes local charge by `-q` on
annihilation, including the antiparticle creation components. Therefore
`Psi_v^dagger M Psi_u` changes charges by `(-q,+q)` at `(u,v)` and its
link shift `-q` has that same incidence boundary. Its adjoint reverses
both. Each derivative and Wilson monomial thus commutes with every local
Gauss operator. The native inverse map preserves these identities.
The exact per-monomial audit passes on the edge and triangle, with 144 and
360 combined native monomials respectively. These terms act on the lepton
coordinates; the common source retains all other native modes.

No Weyl inverse or cancellation of hard-wall shifts is used by the finite
matrix check. It applies the complete monomial and its link shift to every
Fock/flux column, setting a shift beyond the wall to zero.

## Dirac, Wilson and spin identities

The matrices obey

    {alpha_i,alpha_j}=2 delta_ij I,  {alpha_i,beta}=0, beta^2=I,
    alpha_i alpha_j = i epsilon_ijk Sigma_k  (i != j).

Consequently, for the uniform orthogonal specialization with a compatible
scalar Dirac mass,

    h(p)^2 = [c^2 sum_i sin(a p_i)^2/a^2
             +(m+r sum_i(1-cos(a p_i))/a)^2] I.

For positive m and r the positive energy has its unique minimum at p=0
modulo reciprocal periods: the Wilson mass is at least m, with equality
only when every cosine is one. Other momentum bands remain in the law.
This analytic observation is not a theorem about the full interacting
source's three physical charged minima, nor a claim of a chiral continuum
limit. Both charged and neutral derivative terms use the same principal
coefficient c; the supplied chart and observer-limit identification require their own proof.

A directional derivative need not commute with a transverse
spin component. The Clifford-square identity retains the commutators of
covariant momenta. If a physical momentum/field identification supplies
`[Pi_i,Pi_j]=i q epsilon_ijk B_k`, its spin term is `-q Sigma.B`.
With independently identified kinetic mass and spin `S=Sigma/2`, this is
the algebraic factor-two precursor. This component has not calculated an
intrinsic dipole, a physical g factor or a Larmor curve.

The retained neutral mass factors contain their supplied rounded flavour
mixing. `mass_check.py` checks the actual 49-monomial local lepton
restriction of the compiled snapshot. Both complete 12-by-12 charged and
neutral Nambu matrices close, are Hermitian, carry the correct charges,
and anticommute exactly with all three spatial alpha matrices. The check
also rejects a charge-preserving, spin-summed native-density mutation.
These identities establish local mass and principal-kinetic compatibility.
They do not determine the spectrum of the full interacting or geometry-
coupled source.

## Pair conversion on a complete finite carrier

In the canonical polynomial the charged derivative contains 24
pair-current monomials on one edge. Its total canonical-occupation
commutator has 48 nonzero terms when neutral terms are included. These
counts are edge diagnostics in both the edge and triangle receipt
sections; they do not count the triangle's full canonical polynomial.

The z-edge test keeps active family-zero up-particle and paired
down-antiparticle modes at both sites. All 16 active Fock occupations
remain. Inactive channels are put in specified odd pair states:
charged up particles/down antiparticles occupied, and neutral particles
occupied. The complete added derivative and Wilson terms preserve these
spectators by Pauli blocking. The code checks this by acting every term
on every retained column and rejecting any nonzero leakage, including
the Wilson scalar and onsite terms. It is a reduction of the full
generated lepton kinetic, not deletion of unwanted pair terms.

The remaining spectator charge is +3 at each site. A retained root
reference of -6 compensates it. Flux cap 6 gives all 208 active
Fock/flux states; Gauss-error states remain in the matrix and instrument.
The initial active electron at site0 and positron at site1 has E=-6;
the active vacuum has E=-3. The literal generated Hamiltonian element is

    <active vacuum,E=-3 | H | pair,E=-6> = -i/2.

The finite witness adds a declared unit active mass and electric energy
`E^2/16` to this generated kinetic. It does not assert that its unit mass
is the unchanged binary parent coefficient, or that the entire parent
source is proved to preserve these spectators. The master
source must check its own full-parent attachment. All other local mass
spectators would contribute an explicitly retained scalar if that
attachment is proved.

At source time 1/4:

| Quantity | Complete-carrier result |
| --- | ---: |
| Active-vacuum probability | 0.01363257479013916 |
| Same probability with pair-changing derivative removed | 0 |
| Gauss-mismatch probability | 0 |
| Closed-source energy residual | 9.93e-16 |
| Terminal occupation/flux meter work | -0.0948604429141433 |

All 208 occupation/flux outcomes are saved, including zeros and
complements. The control keeps Wilson hopping, mass, electric stiffness,
the source state and detector while removing the number-changing
derivative; it is the same-alphabet control for the particular pair
conversion claim. The exact nonzero matrix element establishes a
short-time channel independently of floating evolution. The numerical
run demonstrates that channel at the declared finite time; its residuals
are not interval-certified floating-error bounds.

This is active-pair conversion into the same longitudinal field/background,
not an on-shell annihilation rate, photon identification or free outgoing
radiation result. The common magnetic face source provides additional
degrees of freedom, whose scattering/energy-shell properties remain to be
tested. Charge is conserved while separate canonical particle and
antiparticle occupations intentionally are not.

## Constants, reference and clock boundary

Every Wilson normal-ordering constant remains. A constant added to an
isolated source changes no source-only probabilities, but coupling that
source energy to an autonomous lapse makes the reference convention
physical for the combined geometry response. No constant is discarded on
the ground that it was an unobservable phase before that coupling.

No stationary electron eigenstate is claimed to have a self-starting
observable mass clock. Its global phase remains unobservable. The
separate retained-reference interference protocol and its actual measured
energy gap must state what they compare; a pair/vacuum gap is not renamed
the electron mass frequency. Physical clocks, spectrum, continuum,
and observer geometry require separate physical comparisons.

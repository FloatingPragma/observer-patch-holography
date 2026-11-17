# Two-cell colour, spin and flavour source

This is a newly declared24-CAR-mode source. Its cells, spins, flavours, currents and Hamiltonian are explicit; evidence from the earlier12-mode spinful source is not transferred by renaming coordinates. The two spatial cells each carry two physical spins, two weak flavours and three colours. Mode `3*(4*cell+alpha)+colour` uses `alpha=(up u,up d,down u,down d)`.

The source is

`H=m(N_left-N_right) -(1/8) sum_alpha,colour(c_left† c_right+h.c.) + kappa C_left²`,

with zero exterior colour flux, complete globalSU3 singlet sectors, kappas1/2,1,2, and the same declared432-analyzer light coefficient as the earlier signed-vacuum experiments. The common energy conversion remains conditional: it uses the electron rest anchor but has not been joined to a common spatial Dirac clock. Reusing this small coefficient as a quark mass slot is a declared architecture, not an inherited quark-Higgs theorem or a target fit.

The exact character product `product_colour,alpha(1+x_colour*y_alpha)^2` reconstructs all16,777,216 full Fock states. Multiplication by the SU3 and SU4 Weyl denominators gives every colour-singlet/spin-flavour highest-weight multiplicity. Colour-singlet dimensions in N=0,3,...,24 are1,120,2520,14112,24696,14112,2520,120,1. The largest multiplicity block is45, while the N12 SU4 singlet has multiplicity35. These counts are independently reconstructed from the exact coefficient character, not assumed from the resulting matrices.

Each source block is generated inside its actual CAR Cartan-weight basis. Polynomial spectral projectors in actual quadratic colour and flavour Casimirs extract the requested highest weights. Cyclic action of the actual mass, hopping and electric operators supplies a basis whose modular rank equals the independent exact multiplicity. Actual CAR raising operators annihilate every column; every compressed Hamiltonian component stays in the same colour kernel and highest-weight eigenspace. Thus the reduced blocks are exhaustive. Integer congruences preserve the Gram matrix and every source component. Signed10-bit limb products compute the same arbitrary-precision matrix products without machine-integer overflow, with direct-bigint controls.

All competing irreps with N<=12 are included. The explicit CAR transformation `c_left -> c_right†`, `c_right -> -c_left†` preserves the mass and hopping terms. It maps `C_left²` to `C_right²`, which equals `C_left²` on total colour singlets. Thus it intertwines the same Hamiltonian at N and24-N, with SU4 highest weight mapped to `(6-lambda4,6-lambda3,6-lambda2,6-lambda1)`. No separate energy normalization is introduced for the conjugate sectors.

The actual physical-spin current is `S_i=sum c†(sigma_i/2 tensor I_flavour tensor I_colour)c`; weak isospin uses the analogous flavour generator. Relative to the half-filled reference, `B=(N-12)/3` and `Q=(2/3)N_u-(1/3)N_d-2=I3+B/2`. These are finite conserved operator identities, not a transferred electromagnetic field-coupling theorem.

The exact globalSU4 currents commute with every source component. The symmetric three-index SU4 representation20 contains both `(S,I)=(1/2,1/2)` and `(3/2,3/2)`, of dimensions4 and16, with equal energy within each multiplicity eigenstate. The complete calculation also includes all other SU4 representations, including the4-dimensional representations that contain only spin/isospin1/2. A low-spin state must therefore be identified from the actual full spectrum, not inferred solely from a chosen20-dimensional representation.

A spin-unitary bond on a tree cannot supply a new spectral splitting in this source class. For bond matricesU, the recursive local change of coordinates `W_child=U_parent_child† W_parent` removes every bond matrix while leaving scalar onsite terms and colour electric Casimirs invariant. Consistent transport of states and currents also leaves measured matrix elements unchanged. This statement is limited to this tree architecture and does not exclude a genuinely new Dirac or magnetic construction.

The source does not yet supply a four-component Dirac spinor, a physical relativistic axial current, a native full preparation/readout, a continuum transfer, or a completed physical nucleon mass result. The separate local auxiliary-current source has its own state and vacuum obligations; its corrections cannot be added as if they had already been derived from this Hamiltonian.

The complete calculation certifies every ground interval in all43 blocks at all three kappas,129 block/member cases. An independent exact integer inertia implementation is checked against rational congruence on64 matrices. Floating spectra propose endpoints only; exact inertia supplies the bounds. The computational mass is rounded to denominator2^36, and the uniform operator error24|m_rounded-m_actual| is added to every nonzero-sector interval. Both the computational inertia interval and the resulting true-source interval are recorded separately.

The source-wide ground state is the N12 SU4 singlet for all three members. The B=+1 ground state is the symmetric SU4 representation20; its B=-1 conjugate has highest weight(3,3,3,0). Thus spin/isospin1/2 states now exist at the lowest baryon energy, but they remain degenerate with spin/isospin3/2 states. Isospin symmetry also leaves neutron and proton energies equal.

| kappa | Certified conditional baryon mass, GeV | Certified conditional two-baryon binding, GeV |
| --- | --- | --- |
| 1/2 | [1.4624348232, 1.4624349507] | [0.1692370159, 0.1692372710] |
| 1 | [0.8033762163, 0.8033763438] | [0.1456337564, 0.1456340115] |
| 2 | [0.3709269864, 0.3709271139] | [0.0732302779, 0.0732305330] |

Binding uses2(E_B1-E_vac)-(E_B2-E_vac), so a positive interval establishes a finite-source two-baryon bound threshold. All three intervals are strictly positive. The B2 ground transforms in SU4[0,3,0], dimension50. No attraction coefficient was added: this is the same complete hopping/electric source and the same source-wide vacuum. It is not a physical deuteron, nuclear-spectrum or halo identification, and the baryon masses miss their physical targets.

Five actual CAR state frames independently replay every Gram, mass, hopping and electric matrix exactly. Six highest-weight triquark polynomial probes at each kappa retain all target spectral lines and their failure probabilities. The four symmetric20 probes have nonzero lowest-line residues; the two mixed20 probes retain separate target sectors. Probabilities are numerical diagnostics with vacuum-frame orthogonality residual below7.7e-12, not certified probability intervals or a complete native apparatus execution. Full SU4-multiplet apparatus records have not been executed.

The frozen local auxiliary source has separately received independent full-state current checks in `../species_leptons/auxiliary_current/`. Its actual dressing changes the axial current, but its positive-gap boundgA>25/18 excludes the target near1.2753. That result is not a correction already joined to the present source.

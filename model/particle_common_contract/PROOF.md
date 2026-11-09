# One declared source, state, reference and clock

## Construction fixed by this contract

`oph.open_wilson_gauss_observer.v1` selects one spatial construction.
Its choices are part of the program, with no claim of architectural
uniqueness or selection from measured masses. The mission permits such a
choice. This closes the specification requirement: for every integer
n>=0 there is one source operator, one prepared density, one retained
reference and one time convention. It does not prove that their records
match nature, or that the narrower C/E repair interpreter implements this
extended quantum-field program.

The source is an open cubic box with side `L_n=2^(2n+1)`, spacing
`a_n=2^(-n)`, native-word link degree `ell_n=n+1`, all sixty CAR modes
per cell, electric coefficient105 and magnetic coefficient1/7.
Every cell owns three outgoing links. Impose zero external charge at every
endpoint, including ghost endpoints. Retain all outgoing electric terms,
interior Wilson hopping bonds and interior plaquettes. There are3L³
owned links, L³+3L² gauge vertices,3L²(L−1) hopping bonds and
3L(L−1)² plaquettes. The Fock order is site/Dirac/internal row; the link
order is owner/direction. Native word spaces use length-then-lexicographic
order, exact Haar inner products and positive-norm Gram-Schmidt, skipping
dependent words. These prescriptions fix matrix coordinates without a
chosen eigenbasis or a ground-state oracle.

The Hamiltonian is the complete parent Wilson expression followed by the
product of the per-link word projections. Its onsite one-body coefficient
is3beta_D/a; forward hopping is
`(-beta_D+beta_D gamma_i)R(U_i)/(2a)` plus its adjoint. Electric energy is
`105 sum C_trace/(2a)` and magnetic energy is
`-sum ReTr_M15 R(U_p)/(105a)`. Keep that energy zero. The native representation
and implementation files are bound inside the source identifier. The
source identity includes the transitive project import closure and verified
frozen data ancestry. This conservatively binds comparison source bytes;
it does not import their states or conclusions. The
electric and magnetic multipliers satisfy `e0*kappa0=15`; this matches
their quadratic principal characteristic to the matter symbol, without
asserting an interacting Lorentz limit.

The source-only physical carrier is the joint kernel of the Gauss
generators. Its projector P is the product of commuting local Haar
projectors. Lie components at one site need not commute. The empty
Fock/constant-link vector is a nonzero physical blank. The
[Gauss construction](../particle_gauss_ground_state/PROOF.md) proves
that P preserves the Hamiltonian and supplies the finite invariant-space
algorithm. No charged source-only vector is substituted for a physical
state; a charged probe requires its explicitly compensated apparatus.

## One preparation

Let `[L,U]` and the log-dimension upper bound b be the exact source bounds
returned by `gauss.spatial_plan`, including coefficient105. Set
`epsilon_n=2^(-n-4)`, `beta=2b/epsilon_n`, choose the least dyadic integer
s>=max(1,beta(U−L)), and the least nonnegative Taylor order m whose
canonical half-trace bound is at most
`min(1/4,epsilon_n/(2(U−L)))`. The scalar-width case uses1/4.
Then

```
Q=P[T_m(beta(U I-H)/(2s))]^s,
rho=Q²/Tr Q²,    Psi=vec(Q)/sqrt(Tr Q²).
```

This is a definition by finite algebraic operations on the declared
carrier. The canonical Householder mapping the source/reference blank to
Psi is the exact preparation primitive of this program. It is identity
when the target equals the blank. Its normalized positive overlap fixes
the Householder phase. The construction retains every outcome and has no
success filter. The prior finite entropy theorem and positive-polynomial
bound certify the energy tolerance relative to the **restricted** ground.
`member(n)` consumes that actual planner rather than copying its bounds.

The full spatial matrix and preparation tape are not executed by the
contract exporter. The exact algebraic program and a finite-precision
compiled realization are different objects. A compiler error tolerance
is recorded separately; it does not select a unique approximate state.
Any claimed implementation must bind its deterministic compiler, tape and
complete output error. Exact source stationarity is asserted for rho;
it is not asserted for an arbitrary state within the tolerance.

## Reference and time

The reference is the full conjugate source carrier, initially the conjugate
blank, retained in every complete output. Its idle Hamiltonian is zero.
For this source-only Gauss preparation both source and reference are
separately gauge invariant. The joint purification generally evolves
under H tensor I even though the source marginal is stationary.

The sole source time convention is
`U(t)=exp(-itH)`, `tau_t(A)=exp(itH)Aexp(-itH)`, in source units with hbar=1.
Exact propagators specify the mathematical program. A finite approximation
has its separately certified clock error; its generator is not silently
identified with H. Counter-evolution on the reference is an explicit
apparatus gate when used, and does not change the idle reference Hamiltonian.

On the positive support, the modular Hamiltonian is exactly
`K=log(Tr Q²)I−2log Q`. This polynomial state is not exactly Gibbs, so
stationarity does not identify modular flow with tau at time −beta*s.
That affine relation applies to the exact Gibbs target, a different state.
No proper-time or GeV conversion is selected here. A fixed-coframe strain
experiment differentiates H at fixed bare a and fixed input rho; preparing
a different rho(V) defines a different experiment.

## Compatibility and claim domains

`member(n)` binds the source, state, reference and clock separately and
also binds their full combined member. Verification regenerates the
literal member and rejects substitutions even when a caller rewrites
the supplied hashes. Consequence receipts bind both that member and
their result content. A typed consumer regenerates the actual dimension,
preparation or characteristic calculation and compares its entire payload;
unregistered claims and substituted results are rejected even with correctly
recomputed hashes. This checks common-object compatibility, rather than
arbitrary mathematical truth. Each consequence also needs its owning proof.

The generated receipt attaches the actual source dimension/energy bounds,
Gauss preparation plan, quadratic clock calculation, regional support floor,
[complete boundary-charge construction](../particle_physical_boundary_state/PROOF.md)
and [complete local record comparison](../particle_primary_locality/PROOF.md).
The boundary consumer exhibits actual Wilson-loop vectors in eleven charge
sectors of this physical state and retains its full finite charge menu.
It resolves finite preparation including boundary charges; it supplies no
refinement-uniform sector weight or joint state limit. Other aggregate
certificates have an exhaustive, separate domain inventory.

The locality consumer instantiates the ambient Wilson influence bound on
the same physical preparation and Hamiltonian, with product link compression
preserving every interaction support. Its remote plaquette intervention
and complete local boundary detector are Gauss invariant; the untouched
canonical reference is retained. At source duration `1/10000`, the complete
output half-trace bounds for levels 0 through 3 are `2^-2`, `2^-7`, `2^-30`
and `2^-124`. These are within-member comparisons, not convergence of
states or records between levels. The electric evolution is onsite and
does not enter the propagation rate. The global Gauss projector is not
inserted as a local interaction. This supplies the finite locality premise
without a Lorentz or interacting continuum conclusion. In particular:

- scalar Q*=0 and Q*=6 have distinct states and charge responses;
- warm full208, cold full208, Gauss-projected78 and a pure ground trial
  are distinct preparations even when they use the same ambient H;
- periodic translation records do not become open-box records;
- smooth-sign, ordinary Wilson and total-degree/product cutoffs retain
  their own actions and regulators;
- internal spinor/flux or family regions do not become spatial balls.

The inventory rejects an unclassified certificate. None of those control
detector laws is admitted as a primary-field law merely by being listed.
The exact native208 global-ground/Gauss-ground equality has its own
spectral proof; transferring a finite preparation requires its error bound.
Empirically fitted parameters retain their ancestry in their separate
controls. The primary program uses no measured mass or coupling.

The sequence fixes all finite inputs; convergence of its joint volume,
spacing and link-cutoff refinement is a distinct output theorem. Fixing
this contract resolves construction item1. Physical state refinement,
relativistic/chiral fields, a particle dictionary and the Einstein regional
premises are not conclusions of a compatibility check.

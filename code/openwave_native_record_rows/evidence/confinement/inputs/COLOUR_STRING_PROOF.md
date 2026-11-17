# A positive Gauss-law source on finite open chains

## One exact two-cell attachment

Write C0 and C1 for the eight native colour generators in the two cells.
On the complete global colour-singlet subspace, `(C0+C1) psi=0` component
by component. Operators in different cells commute. Therefore
`C0² psi=C1² psi` and the following source identity is exact:

`(U-kappa/2)(C0²+C1²)+kappa C0² = U(C0²+C1²)`.

The global colour penalty vanishes in this domain. The adjusted local
coefficient is positive for the selected kappa=4. Every other strong term
is unchanged. The semileptonic current preserves colour, so the same
identity applies to the complete two-cell contact sector. This is an
operator transfer, rather than equality of selected fitted eigenvalues.

The executable checks the identity using integer native operators on the
entire 2,024-state three-quark occupation carrier and its complete
120-dimensional singlet basis. It also checks it on the exact complete
colour projector in the 17,000-state six-quark, zero spin/isospin projection
carrier. That projector has rank 240. Colour is a spin/isospin scalar, so
its zero action extends to all 2,520 magnetic singlet states. The full
coloured-complement difference is nonzero and explicitly retained as a
control. The two constructions are not identical on that complement.

## Open-chain Gauss law

For sites x=0,...,L-1, impose zero external colour flux and define the
axial-gauge electric generator across each cut by
`E_x=sum_(y<=x) C_y`. Then `E_x-E_(x-1)=C_x` and the right boundary flux
vanishes precisely on the global singlet. The declared electric action is
`kappa sum_x E_x²`.

This elimination has the usual exact tree realization. Attach a full
`L²(SU(3))` link coordinate to each edge and use its Haar measure. Fix the
left endpoint frame and recursively multiply link matrices to obtain site
frames. Acting with those frames on a global-singlet matter wavefunction
defines its dressed wavefunction. Haar integration makes this map an
isometry; local gauge invariance and the absence of loops make it onto the
zero-boundary-flux constraint kernel. Differentiating a link frame gives
the summed matter generator on one side of its cut, so the electric
Casimir becomes the stated prefix square. Native hopping carries its
fundamental link transport before this gauge elimination. The remaining
local colour Casimirs and colour-neutral spin/isospin currents are gauge
invariant. No link representation cutoff is imposed. A constraint-failure
projector completes any apparatus on the kinematic extension.

In the global three-quark singlet, the colour epsilon tensor gives
`T_i dot T_j=-2/3` for distinct quarks. For every subset X,

`C_X² = (2/3) n_X (3-n_X)`.

Thus a prefix containing one or two quarks has Casimir 4/3, while an empty
or complete prefix has zero. For ordered positions x1<=x2<=x3,

`V_string = (4 kappa/3)(x3-x1)`.

The program verifies the prefix identity on every occupation of each
executed carrier and retains exact static costs for separations zero to
sixteen. The identity itself holds for any finite open chain.

## Complete fixed-source spatial execution

The parameters of the selected two-cell source are fixed. Replace U by
`U-kappa/2` and the global colour penalty by the prefix electric action.
Extend the positive native hopping and the edge weak-current and
spin-gradient squares along the chain. For L=2 this is exactly the source
above. No coefficient is refitted at L=3,4,5.

The full three-quark colour singlet is `Sym³(C^(4L))`, with dimensions
120, 364, 816 and 1,540. The code constructs its explicit epsilon isometry
into all native occupation states, checks the native colour generators,
and checks hopping and spin-current intertwiners. Every singlet eigenvalue
is retained as a numerical diagnostic. Complete fixed-charge/spin blocks
of dimensions 14, 45, 104 and 200 have separate exact algebraic spectral
certificates from the native integer matrix elements. Their lowest bands
have spin and isospin one half.

The complete preparation program starts from the uniform full native
three-quark carrier, before selecting the colour singlet. Its successful
singlet branch admits the zero-boundary-flux Gauss realization. Coloured
preparation failures are retained in the native kinematic extension; they
are not asserted to be physical states of the zero-flux Gauss kernel.
On success, native occupation effects are compressed to the singlet and
then Gauss-dressed. Their square-root Kraus operators give a complete
physical probe with the same listed probabilities. Only the uniform-input
preparation has the stated zero-work identity; occupation-probe work is
not asserted to vanish.

The preparation is a uniform mixture on the complete isolated lowest
spectral cluster. Consecutive levels separated by less than 1/10000 source
units, or by an unresolved spectral error, belong to that cluster. This
is a declared measurement resolution, with both nearly degenerate edge
states retained. Its exact spectral projection on the uniform full native
carrier succeeds with probability rank/N; the complementary failure has
probability 1-rank/N. Every native occupation effect and every complementary
colour, spin, charge and Fock flag is specified. The projection preserves
mean energy on the uniform input. An idle neutral reference is retained.
The exact finite spectral certificate bounds the isolated-subspace
approximation; the tabulated eigensolver and occupation values are floating
diagnostics rather than independently certified machine transcendental
arithmetic.

The L=2 neutron energy is a calibration input, inherited through the exact
transfer. L=3,4,5 give approximately 0.9396642972, 0.939664299718 and
0.939664299718 GeV. Their native axial response is about 1.6658. The
calibrated two-cell axial response does not transfer to these larger
chains. Numerically close energies on the last two finite graphs do not
prove an infinite-volume limit.

## Uniform relative-confinement estimate

Let `b=A/9-(2/3)(U-kappa/2)>0`. The baryon-density and local-colour costs
on the three-quark singlet sum to
`6(U-kappa/2)+b sum_x n_x²`, with `sum_x n_x²>=3`. Every other displayed
term is positive. Hence the full source is bounded below by

`F + (4 kappa/3) span`,

where `F=3(mbar-delta/2)+6(U-kappa/2)+3b`.
A co-located spin/isospin-half baryon at an endpoint supplies a trial of
energy

`T=3mbar+delta/2+A+(3/4)(Js+Jw+eta+lambda)+3t`,

independent of L. The variational principle and positivity imply
`<span> <= (T-F)/(4kappa/3)` and
`P(span>=R) <= min(1,(T-F)/(4kappa R/3))` for every finite-chain ground
state in the admitted projection. This is a uniform bound on relative
separation, without fixing a centre, asserting a continuum limit, or
introducing pair creation and string breaking.

## Complete integer proof witness

The pinned `spectral_witnesses.json` contains every column and ordered
centre of a dyadic frame, with denominator `2^40`, for each full fixed
projection. These are proof coordinates, not stored spectral conclusions.
The checker rebuilds the exact native source and recomputes its entrywise
rounding enclosure, complete Gram error and complete residual using
integer arithmetic. It applies the unchanged six-quark bound
`epsilon = 2 rho + 4 K eta`, where `eta < 1/2` certifies that the entire
frame is invertible. Every eigenvalue interval follows from these
recomputed quantities. Integer ranges, dimensions and ordering are checked;
no supplied interval or error is trusted and no complementary column is
discarded. Native floating diagnostics still diagonalize the actual source.

## One-length comparison

The nuclear source provides an exact graph quadrupole interval q close to
one half. If one physical length a attaches both records, then
`Q=a² q` and `sigma=(4kappa u/3)/a`, where u is the electron GeV unit.
Eliminating a gives `sigma² Q=(4kappa u/3)² q`. The comparison inputs use
Q=0.285699 fm² from Puchalski, Komasa and Pachucki,
Phys. Rev. Lett. 125, 253001 (2020), arXiv:2010.06888v2, and the explicitly
illustrative tension benchmark sigma=1 GeV/fm. They are not fitting inputs.

At kappa=4 the quadrupole-matching length gives a tension near
0.00361 GeV/fm. Even allowing kappa to vary while preserving the two-cell
operator identity and nonnegative local colour requires kappa<=2U. The
exact outward quadrupole interval gives a maximum tension below
0.429 GeV/fm, incompatible with the supplied 1 GeV/fm benchmark. The
receipt verifies the inequality using rational squared quantities, without
relying on rounded square roots. This rejects the simultaneous benchmark
attachment in this specified positive-source transfer family. It is not
a general obstruction to OPH or to physical colour confinement.

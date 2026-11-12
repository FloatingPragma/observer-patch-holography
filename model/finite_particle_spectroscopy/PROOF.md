# Finite-time joint mass spectroscopy

## Unchanged source and a genuinely finite detector

Use the exact calibrated quadratic source, quasiparticle preparation and
five complete probes of `particle_joint_mass_identity`. The finite Fock
space has dimension 2^180 times 1653. Its one-excitation space has dimension
236. The parent constructs exact common projectors P_g for source charge,
color Casimir, helicity, momentum and energy. They include antiparticles
and whole color fibers. The vacuum and full many-excitation complement
are also retained. Energy and rest coefficients in this package retain
the parent's running-action and calibration conventions.

Let a label only charge, color Casimir, helicity, momentum and source
sector. Its projector is Q_a=sum_(g in a)P_g. No exact energy or family
measurement occurs. Let R_a be the least power of two at least one and at
least 4 max_(g in a)(E_g-1), in the declared GeV convention. Define

    phi_a(E)=1/4+(E-1)/R_a,       1/4 <= phi_a(E) <= 1/2.

A 32-dimensional pointer is prepared uniformly, conditionally evolved by
exp(2 pi i j phi_a(H)) at address j, and inverse-Fourier measured. This
uses the same source Hamiltonian, not a replacement family-label detector.
The scalar phase offset is part of the declared coherent phase reference.
The time sign can equivalently be implemented by reversing the Fourier
convention. The longest magnitude of coherent exposure is
2 pi (31)/R_a in GeV^-1. Five controlled power evolutions suffice in a
binary implementation, together with finite Fourier and readout services.
These are supplied services, not a native pulse compilation or a measured
laboratory clock duration.

The exact outcome amplitude and Kraus operators are

    A_y(phi) = (1/32) sum_(j=0)^31 exp(2 pi i j (phi-y/32)),
    K_(a,y) = sum_(g in a) A_y(phi_a(E_g)) P_g.

Complete the instrument with K_vac=P_vac and K_many=P_many. Then compose
it with each of the original five probes. Since the pointer Fourier
transform is unitary, sum_y |A_y(phi)|²=1. Orthogonality of P_g gives
sum_(a,y) K_(a,y)†K_(a,y)+P_vac+P_many=I. This proof holds for all source
inputs and after tensoring every branch with identity on any reference.
The complement is the entire original many-excitation subspace; it is not
deleted because the displayed preparation happens to give it zero weight.

Each K_(a,y) is scalar on a full color fiber and coherent between different
energies in the same a. Measuring the exact family label first would remove
these cross terms and give a different channel. A maximally entangled
pair of unresolved neutrino families supplies an explicit countercontrol:
the resulting reference-output half trace distance is approximately 1/2,
including every pointer outcome with its original weight.

## Full finite execution and numerical limits

The receipt evaluates every amplitude for all 112 common source groups,
all 32 pointer outputs and all original probe returns. The circulant matrix
whose columns are shifts of A is the full pointer unitary F† diag(exp(2 pi
i j phi)) F. Thus all pointer columns are supplied, including unused inputs,
through the stored kernel, and checked for unitarity. The actual source
charge, color, helicity and momentum commutators are independently replayed.

For the small neutrino energies the computation uses

    E-1 = m²/(sqrt(1+m²)+1),

with the exact rational mass square and 100 decimal digits. It never
subtracts binary64 E and p to infer a tiny mass. Parent rational energy and
mass intervals supply the exact separation bounds and resource planning.
The point amplitude execution has floating diagnostics, not a certified
interval error budget for the literal numeric matrices or physical controls.
Exact completeness and confidence concern the ideal declared instrument.

## One reusable pointer qubit gives the same complete instrument

The full M=2^n pointer can be replaced exactly by n sequential binary
services with one reusable quantum qubit and all n original classical
bits retained. Read the bits from least to most significant. After the
first j bits let r=sum_(l<j) b_l 2^l. Within each measured coarse class,
define U=exp(2 pi i phi_a(H)). Prepare the pointer in |+>, apply the
controlled source power U^(2^(n-1-j)) and the known phase
exp(-2 pi i r/2^(j+1)), then Hadamard-read its bit b_j. The source Kraus
operator is

    K_(b_j|r) = [I + (-1)^b_j exp(-2 pi i r/2^(j+1))
                         U^(2^(n-1-j))]/2.

For a complete word y=sum_j b_j 2^j, define z=exp(-2 pi i y/2^n) U.
Each factor is (I+z^(2^(n-1-j)))/2, since the higher bits contribute an
integer multiple of 2 pi and b_j contributes its sign. The polynomial
identity product_(l=0)^(n-1)(I+z^(2^l))=sum_(k=0)^(2^n-1) z^k therefore
proves that the complete word Kraus operator is exactly A_y(U).
This is an operator identity on every source input, including every
quantum reference, rather than only equality of measured probabilities.
It preserves every original output y by a bijective ordered bit record.
Each binary pair separately satisfies sum_b K_b†K_b=I.

After recording b, the known controlled X^b returns the pointer's measured
state to blank, while retaining b in the classical archive. No unknown
input or classical history is erased. Coherent source evolution, phase
feedback, binary readout and known-bit recycling are supplied services;
an autonomous blank factory is not derived. Vacuum and many-excitation
complement outcomes use the same original terminal branches.

The executable evaluates all 32 adaptive word amplitudes for each of the
112 exact source energies and compares them with the complete Fourier
kernel. The ideal all-input identity supplies the full reference domain;
the high-precision point comparison remains a floating diagnostic. The
executed equivalent detector uses one reusable quantum qubit, five retained
bits and five controlled source evolutions. The sum of their exposure
magnitudes is unchanged, (2^n-1) 2 pi/R_a.

Likewise an n-bit sufficient precision plan can use one reusable qubit,
n retained classical bits and n sequential source services, instead of n
simultaneous pointer qubits. The large planned n values and exposures
remain unexecuted. This reduces quantum storage, not coherent exposure,
classical record size or the assumed physical control accuracy.

## Confidence sets, including inconclusive records

Let M=32 and K=6. For circular lattice distance d=|M phi-y|_circle,
unitarity of the geometric sum and sin(pi x)>=2x for 0<=x<=1/2 give

    |A_y(phi)|² <= 1/(4 d²).

Outside d<=K, the points on each side have distances at least K+l,
l=0,1,... . Therefore

    Pr(d>K) <= (1/2) sum_(l>=0) 1/(K+l)²
             <= (1/2)(1/K²+1/K) = 7/72.

Since phi is in [1/4,1/2] and K/M<1/4, a good record has ordinary,
nonwrapped |phi-y/M|<=K/M. Intersect [y/M-K/M,y/M+K/M] with the declared
window [1/4,1/2]. An empty intersection is an explicit inconclusive outcome.
Otherwise its affine energy interval is mapped monotonically to
sqrt(E²-1). Rational outward square-root enclosures give the mass set.
The coverage is at least 65/72 for every eigenstate and hence for every
state and untouched reference, before any success conditioning. None of
the inconclusive or unlikely outputs is discarded. An implementation
error epsilon in half diamond norm adds at most epsilon to the failure
bound; that physical error has not been certified here.

The confidence sets concern the unknown readout energy of this fixed
calibrated source. They do not propagate uncertainty of the supplied
Yukawa, Higgs, neutrino or clock inputs. The displayed modal record is one
illustration of the complete law, not an experimental measurement or a
deterministic error bound.

## Precision planning is separate from the 32-state execution

For each of the 17 source target types, retain its calibrated mass and
ancestry. The planned nonzero-mass goal is one percent of that coefficient;
at a zero mass the goal is an absolute 1e-6 GeV. Exact rational enclosures
bound the energy changes for these mass goals and the nearest distinct
same-quantum-number energy gap. Select a power-of-two pointer M' so that
four times the confidence half-width K R_a/M' is smaller than both bounds.
On a good record, the whole reported mass interval then lies inside the
goal and excludes every other calibrated source mass with the same a.
For the source state with mass m, the measured center is at most one
half-width from E(m), and its reported interval extends by one more;
the factor four leaves a strict margin. The phase guard prevents aliases.

These planned large registers and exposure times are sufficient ideal
budgets. Only the 32-outcome law and its equivalent one-qubit adaptive realization
are executed here. In particular,
its neutrino-family outcomes remain unresolved; their exact mass-frame
distinctions cannot be promoted to measurements at infinite precision.
All charged-fermion coefficients are the inherited running-action inputs,
not the separate compiled free-lepton pole-like convention. No interacting
pole, hadronic quark mass or independent Standard Model spectrum is derived.

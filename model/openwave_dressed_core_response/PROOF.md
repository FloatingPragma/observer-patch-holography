# Dressed charge profiles and current response

The source is the positive-energy Pauli square root coupled to the complete
compact Maxwell plaquette of `openwave_dirac_maxwell_refinement`. It includes
the Pauli magnetic commutator and the positive covariant Wilson energy. The
native integer charge is q=3Q=-3; every electric probe here uses Q=-1. The
opposite quantum reference remains at vertex zero. Every integer loop flux,
both spin blocks and every native complement of the inherited reducing
one-particle sector belong to the source. This source supplies no transfer to
the signed Dirac sea, pair production or a relativistic scattering amplitude.

The electric and magnetic coefficients are 1/32 and 1/8, respectively. The
Pauli square has scale 1/512, giving free speed squared 1/128. The Wilson
coefficient is 1/16. The three masses and their electron calibration come
from the native parent construction. These coefficients use no charge-radius
or running-coupling data. The zero-Wilson member is a separate control with
all its states retained.

## Exact spin reduction and infinite-source ground certificate

At zero threaded external flux, the Pauli square commutes with sigma_z. Its
two orbital blocks are complex conjugates, and the real charge-position
operators intertwine them. The receipt retains the unpolarized mixture of
both blocks and all their outcomes. This removes neither a spin degeneracy
nor a current channel.

For m>=1, put A=D^2. Then 0<=A<=1/128 and

```
sqrt(m^2+A)-m = sum_(j=1)^6 binom(1/2,j) m^(1-2j) A^j + R,
||R|| <= m rho^7/(1-rho),  rho=1/(128m^2).
```

The polynomial is formed on padded integer-flux support before compression.
Every factor changes flux by at most twice the native charge. The integer
Pauli square is A times512. Its Gaussian-integer row sums are bounded by4;
all six integer products stay below2^53. The code constructs their exact
Gaussian rational sum and checks its equality with the independent operator
API. No floating matrix is silently promoted to exact arithmetic.

A complete complex dyadic eigenbasis certifies every finite eigenvalue.
Integer matrix products give its Gram error and full generator residual.
The actual rational errors must fit the fixed outward budget1/100. Both
real and imaginary entry-rounding errors are included. Integer-overflow
bounds are checked before every fixed-width product. A separate accepted
48-bit dyadic trial, stored in `witnesses.json`, has its full padded residual
verified exactly. These witnesses are certificates, not trusted solutions.
Changing a source matrix requires their residual inequalities to pass again.

Outside |n|<=C, the electric floor is at least (C+1-|q|)^2/8. The matter
kinetic energy and magnetic energy are positive. The boundary coupling is
bounded by b=1/(256m)+2r+1/16, where r is the Wilson coefficient. If a1 is the
certified second finite eigenvalue after the root remainder and c is the
outside floor, the Schur complement gives

```
E1(infinite) >= a1-b^2/(c-a1).
```

The displayed trial centre mu lies below this bound by more than1/10. Its
exact padded residual plus the root remainder identifies the unique ground
in each spin block. Spectral projection bounds its normalized vector error
by twice the residual divided by that gap. Every replay verifies that this
is at most10^-10. The complete spin doublet is retained. Squared-position
expectations differ from the accepted trial by at most2*10^-10.

## A measured radius and a physical charge probe

The four coordinates are (0,0),(1,0),(1,1),(0,1). A native occupation/flux
instrument measures x and its square. Write p=<x^2>; because x is zero or
one, the charge form factor is exactly

```
F(k)=<cos(kx)>=1-p(1-cos k).
```

A real scalar probe couples to the full charge and reference as
G_k=Q cos(kx)+R cos(0)=1-cos(kx). Its first energy derivative is therefore
Q[F(k)-1]. Reporting the separately known reference term gives the matter
coupling QF(k); the reference is retained in the instrument and its work.
The source is neutral as a whole, so its total uniform response is zero.
The radius is measured independently from the occupation records, rather
than inferred from a spectral threshold. For every real k,

```
|1-F(k)-k^2<x^2>/2| <= <x^2> k^4/24.
```

The regulated electron radius is about0.12118 cell. Its largest suppression
on the entire first Brillouin interval is 2p, about0.02937. The source therefore
has a nonconstant spatial coupling, but cannot resolve even a one-sixteenth
suppression within that interval. The zero-Wilson controls have smaller
radii and weaker resolved responses. These finite instruments supply no
order-one core onset, physical electron-radius fit or all-orders running
coupling. A photon/cavity gap is not identified with an inverse core radius.

## Current contact and certified frequency response

Let D=Qx. The complete pure-gradient probe is
H(theta)=exp(-i theta D)H exp(i theta D). Its current and contact are
J=i[H,D] and T=-[D,[D,H]]. Inserting the complete ground-state spectral
measure gives the static Ward identity

```
<T> = 2 sum_n (E_n-E0)|<n|D|0>|^2
     = 2 sum_(n>0) |<n|J|0>|^2/(E_n-E0).
```

The contact is essential; no channel is discarded to obtain cancellation.
The imaginary-frequency response including that contact is

```
Gamma(nu) = <T>-2 sum_n Delta_n |J_n0|^2/(Delta_n^2+nu^2)
          = 2 nu^2 <Dg,(H-E0)/[(H-E0)^2+nu^2] Dg>.
```

It is positive and increasing for nu>0 when the current is nonzero. The
finite tables retain every eigenvalue, dipole weight and current weight;
individual table entries are floating diagnostics. The quoted summed
response has a separate infinite-source certificate. For each frequency
an accepted dyadic vector solves (H-mu+i nu)x=Dv approximately. Its residual
is evaluated on all padded return paths using exact rational arithmetic.
The root remainder is added to this residual. Since the self-adjoint
resolvent norm is at most1/nu, this bounds the response error without a
finite spectral cutoff or a singular inverse on the ground-state subspace.

The kernel 2nu^2 t/(t^2+nu^2) has norm at most nu. The resolvent identity
bounds changes under a generator perturbation by twice its norm. Including
the ground vector error, the trial-centre error and the root remainder gives
an exact error below10^-9 on each reported response. Every replay verifies
this budget. The four regulated frequency intervals are pairwise ordered.
The small zero-Wilson heavy-family responses can be below the certified
absolute resolution; their floating values receive no stronger claim.

## Complete preparation, readout and work

The input is the matter/reference vacuum and bare loop n=0, with mean field
energy b=1/8. The rank-one success operator maps this input to the accepted
charged trial with amplitude1/sqrt(2). Its full positive square-root failure
complement retains the input with probability1/2. All orthogonal Fock and
reference inputs are in that complement. The charged spin is an unpolarized
mixture. The work of this preparation is

```
W_prepare = [m+<v,H_minus_mass v>-b]/2.
```

The bare loop itself is an externally prepared field state at energy b-E_gamma
above the true field vacuum. That input preparation is separate from the
charged preparation work.

For each wavevector the probe applies U=exp(-iG_k). Choosing X or Y uniformly
and retaining both signs gives four Kraus operators
(I+sU)/(2sqrt(2)) and (I-isU)/(2sqrt(2)). Their Gram sum is exactly the
identity. The receipt retains every pointer, spin, position and flux outcome,
plus preparation failure. The nonselective phase-pointer probe is
rho -> [rho+U rho U*]/2. Its original-input work is one quarter of the
conditional unitary energy change. The reported terminal spin, position and
flux projectors then dephase every pointer branch. Their additional source
work is retained separately. Both conjugate spin Hamiltonians are executed
with the same physical probe phase; conjugating the state does not conjugate
that intervention. Preparation, phase-pointer and terminal position/flux work
sum to the final complete mean source energy minus b, with failure included. No charge is absorbed by this phase probe. The reference and
Gauss constraint are unchanged. The executed point probabilities and work
are numerical diagnostics; the ground-radius and summed response intervals
have the independent exact certificates described above.

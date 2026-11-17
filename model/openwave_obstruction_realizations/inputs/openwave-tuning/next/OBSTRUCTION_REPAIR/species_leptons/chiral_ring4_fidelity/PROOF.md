# Unique finite-volume vacuum and its bounded current response

## The unchanged complete source

The Hamiltonian is the frozen four-cell, two-wall, one-physical-spin
U1 control, with sixteen CAR modes, all particle-number sectors and
the full infinite electric rotor. Its static reference, fractional-service
clock h=1/8, electric coefficient1/4 and boundary mass are unchanged.
The two mass cases are zero and the actual432 electron coefficient.
This is a finite-volume source theorem, not a three-dimensional lepton
mass prediction or a transfer to the full nonabelian gauge vacuum.

All other particle-number sectors have lower bounds above the cuts
used here. The neutral N8 sector has12870 complete Fock occupations.
Its electric energy splits exactly as

```
H_8=D+W(theta),
D=(−i d/dtheta+qbar)^2,
W(theta)=H_matter,8(theta)+B,
B=sum_x(q_x−qbar)^2/4>=0.
```

qbar is a constant diagonal matrix on the full Fock basis. W is a
bounded Hermitian matrix multiplication operator. D has compact
resolvent on the periodic rotor, so H has a discrete spectrum with
finite multiplicities and an attained ground energy. No commutation
between D and W is assumed.

## Shifted-rotor heat kernel and the full interacting trace

Golden–Thompson gives

```
Tr exp[−beta(D+W)] <= Tr exp(−beta D) exp(−beta W).
```

The unbounded version follows by finite spectral compression of D and
the trace-class limit; W is bounded and exp(−beta D) is trace class.
For each Fock word the diagonal rotor kernel contains

```
K_beta(q)=sum_(ell in Z) exp[−beta(ell+q)^2].
```

Poisson summation writes it as
`sqrt(pi/beta) sum_n exp(−pi^2 n^2/beta) exp(2pi i n q)`.
All paired Fourier coefficients are positive, hence
`K_beta(q)<=K_beta(0)=Z_rotor(beta)` for every real q. This includes
every occupation-dependent Gauss shift. With normalized Haar measure,

```
Tr exp(−beta H_8)
 <= Z_rotor(beta) average_theta tr_Fock exp[−beta W(theta)].
```

The kernel comparison is an operator inequality on the diagonal Fock
kernel. It remains valid when that kernel does not commute with W.
The centered rotor sum keeps ell=−16 through16 explicitly and bounds
its remaining tail by `exp(−256 beta)/(16 beta)`. Interval arithmetic
at40 decimal digits and outward binary conversion enclose that sum.

## Every free excited Fock level is retained

At each of1024 holonomies, all sixteen one-particle bands are evaluated.
Every eight-element occupation subset is included, giving12870 free
Fock energies E_j(theta). Positivity of B and the min-max principle
give `lambda_j(W)>=E_j` for every ordered eigenvalue. The frozen
Coulomb two-block calculation supplies the stronger pointwise ground
lower bound L_0(theta). Therefore

```
tr exp(−beta W)
 <= exp(−beta L_0)+sum_(j>=1) exp(−beta E_j).
```

Only one free ground contribution is replaced. All excited levels,
degenerate levels and high-energy occupations remain in the sum. The
implementation removes a minimizing computed occupation; removing any
other single occupation would make this upper bound larger. An energy
envelope1e−8 is applied in the conservative direction to every retained
term. The complete sixteen-band arrays and all12870 occupation indices
are stored, so every one of the13178880 evaluated Fock energies is
recoverable. The independent reviewer reconstructs the entire trace
by positive particle-hole elementary-symmetric sums, without deleting
or sampling excited occupations.

The actual potential W has derivative norm at most
`2||C_8||<=16 sqrt(2)`. The nearest-grid distance is pi/1024. Eigenvalue
comparison therefore gives the global integration bound

```
average_theta tr exp(−beta W(theta))
 <= exp[beta 16 sqrt(2) pi/1024]
    average_grid pointwise_trace_upper.
```

This relative quadrature factor controls the actual interacting
potential. It does not differentiate the comparison Slater frame or
assume smoothness through its Fermi-edge degeneracies.

The trace is scaled by a declared anchor energy to avoid large
intermediate exponentials. The anchor cancels in the eigen-count
formula and is not a source-energy subtraction. The beta values
1/4,1/2,3/4,1,5/4,3/2,2,3 and cut offsets1/16,1/8,1/4,1/2,1 are
verification filters; every test is retained.

## An actual eigen-count certificate

For any cut c,

```
number of eigenvalues <=c, counting multiplicity
 <= exp(beta c) Tr exp(−beta H_8).
```

A bound below2 implies at most one eigenvalue below the cut. A normalized
variational state with energy below c proves at least one. Together
they establish exactly one, including multiplicity. The retained
number-sector bounds exclude all other sectors below the same cut.

For zero bare mass, use beta=5/4 and the K4 Rayleigh upper bound plus1/8:

```
cut=−92.61918850609149,
eigen-count upper=1.91396856558887 <2.
```

The source thus has a unique actual vacuum, with spectral gap at least
1/8 above it. The electron-input member obeys
`||H_m−H_0||<=4|m|`; its heat trace is bounded by
`exp(4 beta |m|) Tr exp(−beta H_0)`. Its corresponding count upper is
1.91446652288262, also below2. No critical mass is selected.

Filters that fail the below2 gate are retained. Their failure is
insufficiency of that upper bound, not evidence for several vacua.
The source has no thermodynamic or continuum gap theorem here.

## Ground-state fidelity from the full infinite-flux residual

Let psi be a retained normalized preparation, mu its Rayleigh
expectation and sigma its complete spectral residual. The residual
contains the interior matrix residual and both exterior electric-flux
blocks. It is a residual of the infinite Hamiltonian. The frozen
calculation supplies a1e−6 Rayleigh envelope and a1e−8 residual envelope.

If the certified cut lies at least delta above mu, all excited
eigenvalues lie at least delta above mu. Consequently

```
1−|<Omega,psi>|^2 <= sigma^2/delta^2,
E_ground >= mu−sigma^2/delta.
```

The second relation is the Temple inequality. The stated Rayleigh
envelope is included in the reported energy interval. For zero mass
and K4 the admitted values give

```
E_ground in [−92.744190793400, −92.744188506091],
|<Omega,psi>|^2 >= 0.999997701532.
```

The electron-input K4 preparation has fidelity at least
0.999997701523. The coarser K2 residuals give only a zero fidelity lower
bound in this direct test. Their failed accuracy controls remain in
the output even though the source's vacuum uniqueness is established.
Sparse neutral Ritz ordering and the warm-seed control are not used
to prove this eigen-count or fidelity result.

## Transfer to the actual vacuum current

The retained particle and hole currents obey the complete CAR identity
and include the charge-reference translation. For a normalized boundary
spinor u they define one canonical fermion d. At a reference energy
below the charged spectrum, let R_+ and R_- be the full charged
resolvents. The zero-frequency response in direction u is the
expectation of

```
A_u=−d R_+ d*+d* R_- d.
```

If both charged gaps are at least g, then `||A_u||<=1/g`: its two
positive terms are bounded respectively by `d d*/g` and `d* d/g`,
whose supporting CAR projectors are complementary. The pure-state
trace distance from psi to the actual ground is at most sqrt(p), where
p is the admitted excited-state weight. Hence the Green-matrix
state-transfer error is at most `2 sqrt(p)/g`.

The resolvent identity bounds the additional reference-energy transfer
by `|Eref−E_ground|/g^2`, with the Rayleigh and Temple envelopes included.
Adding the frozen full-rotor prepared-current error yields a bound on
the actual vacuum Green matrix. At zero mass and K4:

```
prepared-current error <= 0.000037520541,
ground-state transfer <= 0.001007269351,
ground-energy transfer <= 0.000000252417,
actual vacuum Green error <= 0.001045042307.
```

The electron-input case gives at most0.001045116208. The admitted matrix
inverse perturbation bounds are about1.951. The local inverse beta
coefficient remains a local response observable, not an identified
propagating mass. This theorem supplies an actual finite-volume vacuum
response with explicit error; it does not solve the charged-family
ratios, their continuum poles, or neutrino mixing.

## Arithmetic and transfer scope

The rotor and holonomy factors use outward interval arithmetic. The
complete Fock sums use the stated1e−8 energy envelopes and positive-sum
roundoff factors. Their inputs retain the Coulomb comparison's explicit
one-body and moment arithmetic admissions. The independent reviewer
uses a distinct canonical trace identity and rechecks every count and
vacuum consumer. This is a finite numerical certificate with declared
arithmetic scope, not a claim of interval sparse eigensolvers.

No source coefficient, clock, reference, field state or interaction was
fitted to a measured mass. The proof applies to the declared one-spin,
one-dimensional compact-gauge source. It transfers neither a Euclidean
Ginsparg-Wilson measure nor a free coframe regulator into a physical
interacting four-dimensional theory. The independently scoped Wilson
additive-mass obstruction remains a property of the fractional Wilson
member; the wall vacuum theorem does not erase it.

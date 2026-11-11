# A neutral scalar pole on the actual local native ground

The source is the full thirty-mode native fermion algebra and five-state
scalar at unit up, down and lepton coefficients. Its Hamiltonian is
`H=V+N`, where N is the scalar occupation projector and V the literal
native interaction. The reference is the exact local native-current
singlet ground mixture, with the same retained dual reference and clock
`exp(-itH)` as `defect_scalar_response`. The calculation uses the complete
source through its proved invariant reduction; it does not assert that
a sixteen-dimensional block is the whole Fock carrier.

The exact source identities are

```
N^2=N,  {2N-I,V}=0,  H^2-H=V^2,  [N,V^2]=0.
```

The ground belongs to the selected V² spectral projector with eigenvalue
lambda. Its two energies are e=(1-Delta)/2 and (1+Delta)/2, with
`Delta=sqrt(1+4lambda)=1-2e`. The projector is derived from the complete
rational annihilator in `defect_scalar_response.ground_pair`; it is not
a numerically chosen pair of states. Since N commutes with V², applying
N or I-N cannot leave this pair. Every other full-source energy has zero
weight for these preparations. The physical ground multiplicity is
retained by its normalized projector; the identities apply to each member.

Its occupation and connected variance are exactly

```
n=<N>=(1-1/Delta)/2,
w=<N^2>-<N>^2=n(1-n)=lambda/(1+4lambda).
```

The connected scalar operator O=N-n has a single positive spectral line,
of energy Delta and weight w. Its retarded response and canonical spectral
coordinate are

```
D_N(z)=w[(z-Delta)^(-1)-(z+Delta)^(-1)]
      =2Delta w/(z^2-Delta^2),
q_N=(N-n)/sqrt(2Delta w),    D_q(z)=1/(z^2-Delta^2).
```

These are exact finite scalar response identities on the stated ground.
The elementary charged doublet field has a different multi-line response.
For the literal capped scalar fields the exact identity is
`sum_k(phi_k^dagger phi_k+phi_k phi_k^dagger)=4I-3N`.
Thus N is a neutral scalar fluctuation operator constructed from the
actual fields. It is not replaced by an unbounded continuum radial field.

## Complete observation law

Measure the projective instrument {N,I-N}, retaining both occupation
outcomes, and then measure the complete Hamiltonian energy. With the
initial ground energy subtracted, the nonzero probabilities are

| Occupation | Energy 0 | Energy Delta |
|---|---|---|
| Occupied | n² | n(1-n) |
| Blank | (1-n)² | n(1-n) |

Each row also includes the full outside-pair outcome with exact probability
zero. These six original records sum to one without conditioning.
Both occupation and energy operators preserve the native singlet and spin
symmetries; the unchanged reference is retained. The measurement supplies
work: its mean source-energy increase is `2Delta w`, so the positive
energy is not obtained from a free deterministic reset.

The pole energy is approximately 23.4611705885 source units. The receipt
contains exact algebraic coefficients and outward rational intervals.
There is no spatial propagation in this local calculation. A scalar pole
at a single cell does not establish a relativistic scalar rest mass or an
identification with the observed Higgs boson. A common GeV factor requires
a same-source mass/clock attachment. No W, Z, photon or gluon operator is
introduced by this scalar response calculation.

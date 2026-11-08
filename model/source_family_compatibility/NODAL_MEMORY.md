# Nonzero source memory near every kinetic node

The centered cubic full-CAR ground has a strictly positive composite-response
remainder near every zero of its kinetic symbol. A local operator inequality
and a source variational state supply the lower bound. At `U/c=8`, the bound
is uniform in the invariant thermodynamic reference and in a fixed momentum
neighbourhood. It also gives a positive dissipative memory at every nonzero
imaginary frequency. It supplies neither a zero-frequency atom nor a
particle mass.

## The actual star and its ground energy

Use the centered source, complete invariant ground construction and proper
cubic average of `MOMENTUM_RESPONSE.md`. Write

    K=T+U sum_x[n_xup n_xdown-N_x/2],
    Q=K+U|Lambda|/2=T+(U/2) sum_x(N_x-1)^2.

The scalar shift leaves every excitation clock unchanged. At the origin the
six signed hopping links have coefficient `c/2`. Their normalized neighbour
orbital has two modes `r_s=6^(-1/2)sum_y a_ys`; it introduces no extra mode.
Put `r=c sqrt(6)/2` and

    Tstar=r sum_s(a_0s* r_s+r_s* a_0s),
    Qstar=S_0.S_r+eta_0.eta_r,
    Astar=I/4-Qstar.

The eta signs of the origin and neighbour orbital are opposite. On their
full sixteen-dimensional CAR space the exact matrices satisfy

    [Tstar,Qstar]=0,
    spectrum(4 Astar +/- Tstar/r)={0 [10],2 [5],6 [1]}.

Consequently the source operator inequality is

    -4r Astar <= Tstar <= 4r Astar.                    (1)

The zero eigenspaces are retained. No occupation sector is projected out.
The coefficient four is sharp: the `Qstar=0` hopping block has eigenvalues
`+/-r`, while its `Astar` equals `1/4`.

A product of disjoint actual bond ground states is a variational state for
Q on any even cubic box. Each dimer has two fermions and energy
`(U-sqrt(U^2+4c^2))/2`; every bond between dimers has zero hopping
expectation because each factor has fixed particle number. Dividing by two
sites gives the source energy-density bound

    e <= -delta,    delta=(sqrt(U^2+4c^2)-U)/4 >0.       (2)

Boundary unmatched sites have vanishing density. The finite full ground
minimizes energy against this trial. Its translation average and proper
cubic average preserve the local energy density, so (2) holds on the
specified weak-star ground limit. This is a variational inequality on the
same Hamiltonian, not a replacement of the interacting ground by dimers.

The contact contribution to Q is positive, hence `e_T<=e<=-delta`.
Every kinetic bond is counted twice in `Tstar`, giving
`omega(Tstar)=2e_T`. Taking the negative inequality in (1) yields

    qstar=omega(Qstar) <=1/4-delta/(c sqrt(6)).          (3)

This strengthens the purely algebraic `qstar<=1/4` on the actual ground.
It does not require clustering, a ground gap or a factorization of local
correlations.

## A whole neighbourhood of nonreducing response

The exact source symbol and Gram from `MOMENTUM_RESPONSE.md` are

    h(k)^2=c^2 S(k) I,  S(k)=sum_j cos(k_j)^2,
    R*R(k)=[g0-6g2+(4g2-kappa^2 c^2)S(k)]I,
    g0=c^2(9/2-12qstar), g2=c^2 q2,
    -3/4<=q2<=1/4, |kappa|<=3.

Equations (3) and these exact local spectra give

    g0-6g2 >=2 sqrt(6)c delta,
    R*R(k) >=[2 sqrt(6)c delta-12c^2 S(k)]I.            (4)

Thus the residual is strictly positive at every kinetic node `S=0`.
For `U/c=8`, the rational inequalities `sqrt(17)>103/25` and
`sqrt(6)>12/5` imply `delta/c>3/50` and

    R*R(k_node) >(36/125)c^2 I,
    R*R(k) >=(18/125)c^2 I  if S(k)<=3/250.             (5)

The second statement concerns a set of positive momentum measure and all
coefficient packets supported there. It does not require assigning an
independent Hilbert-space fibre to an isolated measure-zero point.
The continuous finite-range Gram supplies the value at the node itself.
Both hidden-spin components and the complete magnetic parity cell are
retained. The original fourth-moment difference from the two-block
compression is `alpha^2 R*R`, `alpha=U/2`; (5) is therefore visible in
original response moments, without normalizing by a herald probability.

## Positive frequency-dependent memory on the same generator

Let W, R and the selfadjoint orthogonal compression `Lperp` be those of
`COMPOSITE_CAR.md`. All translation projections commute with this source
construction. The actual memory is

    Sigma(z)=R*(z-Lperp)^(-1)R.

A quantitative source bound for its first energy moment follows from the
literal odd operator `Y_x=[b_x,T]`. Each incident term is `(c/2)D_xy`,
where the full four-mode calculation gives `||D_xy||=sqrt(5)`. The sum of
these six term norms is at most `7c`. Each term has support on two sites.
With `B=U/2+12c`, the local commutator estimate in `UV_TAIL.md` therefore
bounds the sum of norms in `[K,Y_x]` by `14cB`. Its support lies in the
lattice ball of radius two. Two such odd operators anticommute beyond
distance four, whose ball has 129 sites. The full two-flavour Gram row
bound gives

    ||L Y|| <= sqrt(4*129)14cB <322cB.

Here L is the same signed retarded generator. Since `R=Y-V1 C`,
`Lperp R=(I-WW*)L Y-R C`, and the source bounds
`||R||<=sqrt(18)c`, `||C||<=3sqrt(3)c` imply

    ||Lperp R|| <= M,   M=322c(U/2+12c)+23c^2.          (6)

These inequalities use local source derivatives, not a norm of the
extensive Hamiltonian. They hold on all finitely supported response
vectors and extend by closure. At `U/c=8`, `M=5175c^2` is an explicit
conservative bound.

Let f be a unit coefficient packet supported in `S<=3/250`, and put
`r0=(18/125)c^2`. For `eta>0`, the spectral measure of `Rf` has mass
`n=||Rf||^2>=r0` and second moment at most `M^2`. Cauchy--Schwarz (or
scalar Jensen for the inverse) gives

    -Im <f,Sigma(i eta)f>
       =eta <Rf,(eta^2+Lperp^2)^(-1)Rf>
       >=eta n^2/(eta^2 n+M^2)
       >=eta r0^2/(eta^2 r0+M^2)>0.                   (7)

All frequencies and vectors refer to the same state and source clock.
The local inequalities prove a nonzero dynamical correction throughout
this momentum region; replacing the second block by a frequency-independent
Hermitian term cannot reproduce (7). The lower bound tends to zero with
eta. It neither supplies nor excludes a zero-frequency memory atom, a
surviving ballistic low-energy residue, or a continuum mass shell.

The executable verifies the full local star inequality, the sharp
coefficient, the dimer trial energy and the exact rational constants.
Its bounded scores for the local operators retain all outcomes. A local
state or native-service half-trace error is an additional measured-score
error; it is not absorbed by replacing the source ground or dropping
records.

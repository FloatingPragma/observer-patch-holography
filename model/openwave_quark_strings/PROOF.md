# Fractional native charge on colour-carrying strings

## Construction and preparation

Use the complete side-two periodic native source with 240 fermionic modes,
five scalar states per cell, 119 states per link, and the thirty-mode static
reference. Its Hamiltonian is H=H0+lambda T+4Nref with lambda=2^-34. The
specified input fills every source fermion mode, leaves every scalar and link
blank, and has the reference vacuum. This is a declared finite preparation,
not the interacting ground or a physical vacuum. All twelve local current
traces vanish, so the filled state is a gauge singlet.

For either of the three-dimensional weak-singlet families uc and dc, use
B=sum_i a_ref,i^* c_0,i/sqrt(3) with reference spin zero. On reference-vacuum
inputs B^*B=N_family,spin0/3<=I. Its failure complement and the input
reference/preparation failures complete the instrument. B succeeds with
probability one on the specified filled state. It creates the actual
reference-particle/source-hole singlet, without replacing the native current
by a label. The reference has its full Fock carrier; only this input occupies
its indicated one-particle sector.

For a path of ell edges, contract the actual link matrix indices with the
reference and hole indices, giving coefficients 3^(-(ell+1)/2). The normalized
link labels are (+,i,j) on a forward edge and (-,j,i) on a reverse edge. The
matter endpoint is the actual CAR vector c_i|filled>, including its occupation
sign. This identifies every reduced hole coordinate with a vector of the
original carrier. Its generator is -T^transpose because the filled state is
a singlet. Reference, internal-link and hole generator actions telescope at
each vertex. No links or complementary matter states are removed.

The producer checks all twelve exact constraints on all shortest paths and
also applies the original full CAR/scalar/link Gauss compiler on an embedded
path to each cell, including spectator cells. Omitting the reference or the
hole occupation signs fails nontrivial Gauss constraints.

## Actual charge, flux and centre phase

The native current is Q=T3+Y. The selected endpoints are weak singlets,
so their source-hole charges are +2/3 for uc and -1/3 for dc. The reference
has the opposite charge. Central hypercharge Gauss cancellation attaches
these charges to the actual string flux, as for any region excluding the
reference. Nontrivial colour flux is retained along the path.

Let C be the native colour Cartan induced by diag(1,1,-2)/3. The exact
fifteen-weight identity Q_i+C_i in Z lifts to the full Fock carrier:

    exp(2 pi i dGamma(Q)) exp(2 pi i dGamma(C)) = I.

On each selected hole family, the electromagnetic two-pi phase is
-1/2-i sqrt(3)/2. The colour centre supplies its inverse. A three-turn
electromagnetic rotation returns the state. A controlled endpoint central-Y
rotation and a coherent pointer measure this phase; on the weak-singlet
endpoint it is the Q phase. Choosing the X or Y quadrature with probability
one half gives four complete outcomes

    (1/8, 3/8, (2-sqrt(3))/8, (2+sqrt(3))/8).

At three turns these become (1/2,0,1/4,1/4). All failure/complement records of
the preparation and arrival experiment are retained with their own pointer
outputs. The central phase is an internal gauge-path holonomy. It is not a
proof of fractional spatial texture winding, a colour-flux-tube profile or
an observed quark identity. In particular the electromagnetic path does not
close after one turn; calling its partial phase an integer winding would be
incorrect. This distinction is preserved in the result record.

## Native hopping and complete time-dependent observations

Apply the original full PW119 hopping operator T to the prepared origin
singlet. The result lies exactly in the six orthogonal one-edge string
subspaces, each retaining both hole spin states. The direct CAR and compressed
coordinate calculation gives squared amplitude 1/16 on each route and total
norm squared 3/8, for either family. These are the actual kinetic words,
including their orientation, spin phase and reverse links.

Let F_e be the spin-pair projector for route e and u the prepared input.
Since H0 preserves all link occupations, F_e exp(-it H0)u=0. The reference
cost four is constant on the retained one-reference sector and factors out.
The inherited full-carrier bounds are ||H0||<=480 and ||T||<=360. Duhamel's
identity therefore gives

    ||F_e exp(-itH)u + i lambda t F_e T u||
        <= lambda t^2 ||T|| (2||H0||+lambda||T||)/2.

Consequently each original-input probability lies between
(lambda t)^2 max(0,1/4-epsilon)^2 and
(lambda t)^2(1/4+epsilon)^2, where
epsilon=360(960+360lambda)t/2. At t=2^-24 this interval is approximately
[6.91604e-37,8.15583e-37]. The exact rational bounds are emitted. There is
zero arrival at time zero and under the no-hopping control. The complete
arrival instrument consists of the six mutually orthogonal F_e and
I-sum_e F_e; quantum references and all pump/preparation failures persist.

These positive probabilities establish a particular microscopic charged
transport experiment from a supplied highly occupied state. They do not
establish a low-energy quark quasiparticle, its mass, confinement dynamics,
a hadron, or the absence of screening in other states. The existing equal
up/down native-energy obstruction is unaffected.

# Pure, exactly neutral and cubic-invariant weak-source family

This is a proposed next preparation, after the geometric representation gate.
No gravity or stationarity result is claimed in this directory yet.

Let P_minus(k) be the actual flat source vacuum projector and
P_plus=I−P_minus. For the cubic Gaussian envelope f_W(k), define four-column
mode matrices

    M_plus/minus=Σ_k f_W(k)² P_plus/minus(k),
    X_plus/minus(k)=f_W(k) P_plus/minus(k) M_plus/minus^(−1/2).

Their Gram matrices are exactly I and their cross Gram matrix vanishes.
The four-particle/four-hole Slater covariance is

    C_W=P_minus−X_minus X_minus^dagger+X_plus X_plus^dagger.

C_W is a projector of the same rank as P_minus. Thus every prepared state
has exactly the vacuum total charge, is pure, and is Gaussian; no external
preparation flag or Gaussian approximation to a mixed state is used.

On the cubic source, M_plus/minus commute with the spin representation of
every cubic rotation. The normalized complete four-mode subspaces are
therefore invariant. This must be checked numerically, together with all
Gram and projector identities, before using the preparation.

The actual excitation energy is

    E_W=tr(X_plus^dagger H X_plus)−tr(X_minus^dagger H X_minus)>0.

W changes the profile as well as the energy. In the resolved massless limit
E_W scales as1/W and local stress as1/W^4. It is not an arbitrary continuous
amplitude at a fixed profile. A weak-source sequence must retain its actual
energies, box-size errors, profile changes and geometries. Fixed-volume
stationarity is tested on the resulting state/geometry family, with unused
balls held out. The periodic mean-energy condition remains present.

At zero geometry response, direct regional entropies are exactly the
Gaussian covariance entropy of C_W restricted to the observer factor.
At nonzero geometry, the vacuum projector and the prepared modes must be
computed from the SAME perturbed Hamiltonian, then transformed by the same
declared observer map. Reusing the flat vacuum plus a linear entropy
counterterm would not discharge that test.

# Original spectral tails from local repair-source commutators

The same centered cubic occupation source supplies a volume-independent
exponential energy-tail bound for every normalized original fermionic
wavepacket. The bounds retain both original charge weights and can control
clock aliasing in `LOCAL_SPECTRAL_READOUT.md`. They do not establish a pole,
a positive weight in a shrinking energy window, or a physical particle.

## Local derivatives and all-volume moments

Use the source, invariant full ground and retarded generator of
`COMPOSITE_CAR.md`. A cubic nearest-neighbour graph has maximum degree
D <= 6. Write the centered onsite term, up to a scalar, as
W_x = U(n_up-1/2)(n_down-1/2). Its norm is U/4; each two-flavour bond has
norm c. Set

    A = (U+cD)/2,       B = U/2+2cD,
    C_n = A B^(n-1) (n-1)!      (n >= 1).

The first derivative [a_i,K] is a sum of one-site odd operators of total
norm at most A: the hopping annihilators and (U/2)b_i. For a term with
support size s, the sum of commutator norms is at most sB times its norm.
Only onsite terms and incident bonds contribute, because the Hamiltonian
is even. Each commutator increases support size by at most one. Expanding
the derivation from those one-site first terms therefore gives

    ||ad_K^n(a_i)|| <= C_n,
    support(ad_K^n(a_i)) subset ball(i,n).                 (1)

This estimate applies to the original Weyl links as well as the matched
real signed chart. Each Pauli link has one nonzero entry of magnitude c/2
per row; no protocol equivalence for the bare repair is inferred.

The cubic l1-ball cardinality is

    V_3(r) = (4r^3+6r^2+8r+3)/3.

Two odd derivatives anticommute when their centres are farther apart than
2n. Their anticommutator norm is at most 2C_n^2. With two modes per site,
each Gram row has at most 2V_3(2n) possibly nonzero entries. The Schur bound
then proves, for every unit l2 smearing f,

    integral E^(2n) d(mu_f^+ + mu_f^-)(E)
        <= 4 V_3(2n) C_n^2 =: M_(2n).                    (2)

Here both mu measures have positive physical energy K-E0 and retain their
original weights, whose sum is one. Separate weights need not be 1/2 for
a nonlocal packet. The invariant-ground identity V0*L^2V0=h^2+U^2/4
improves M_2 to (cD/2)^2+U^2/4. Equation (2) is uniform in volume and f.
It passes to the same infinite-ground representation by local convergence
and closure of the derivation domain. For an infinite l2 packet, first
take finite smearings and use the uniform bounds and spectral projections.

## Explicit tails and moment determinacy

Markov's inequality gives the original total charge tail

    (mu_f^+ + mu_f^-)((E_*,infinity))
        <= min(1, M_(2n)/E_*^(2n)).                      (3)

The executable searches finitely many integer n using rational arithmetic.
There is also a closed exponential bound. If B>0 and
m=floor(E_*/(2B)) >= 1, choose n=m+1 and use m! <= m^m:

    tail <= min(1, 4 V_3(2m+2) (A/E_*)^2 4^(-m)).        (4)

When B=0 the source derivative vanishes and all these tails vanish.
Cauchy--Schwarz and 2 sqrt(V_3(2n)) <= 23 n^2 imply, for 0<=s<1/B,

    integral exp(sE) d(mu_f^+ + mu_f^-)(E)
        <= 1 + 23 A s/(1-sB)^2.                         (5)

Indeed the n-th absolute moment divided by n! is at most
23 A B^(n-1) n; its power series sums to (5). Thus each separate charge
characteristic function is analytic in |Im t|<1/B. If all moments,
including the zeroth charge weight, agree, their power series agree near
zero, then their analytic functions agree throughout the strip, and Fourier
uniqueness gives identical charge measures. This is determinacy from the
entire sequence, not identification from finitely many moments.

For c=1/8,U=1,D=6, A=7/8,B=2. The exact eighteenth-moment bound at E_*=20
is less than 1.055e-5, below one sixteenth of the original local band weight
25/100352. It can be subtracted directly as an alias budget; it is not a
renormalized conditional probability. Under ballistic scaling B itself
grows: this fixed-source tail does not provide uniform physical-energy
tightness, and is consistent with `LOCAL_CHARGE_ESCAPE.md`.

`uv_tail.py` checks moments through order twelve on the complete 16-state
dimer and 64-state three-site source. Integer products use arbitrary
precision; exact rational positive-semidefinite elimination checks the
full original moment Gram. An independent spectral test checks the
original charge measures and exponential moment bound. These finite
controls verify the implementation without replacing the all-volume proof.

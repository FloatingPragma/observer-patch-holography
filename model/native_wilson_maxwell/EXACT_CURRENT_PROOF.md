# Independent exact current contraction

The actual cube one-body Hamiltonian and twelve Peierls derivative
matrices are assembled by Kronecker products in `literal_matter`, using
the executed oriented cube. This construction is independent of the
block-insertion producer in the matter calculation. Both matrices must
agree entry by entry when the common-source receipt consumes them.
All matrices are Hermitian over QQ(i). The derivative current is
`J_phase=partial H/partial theta`; its negative is the oriented charge
flow for D=head-tail. The coefficient is unchanged by reversing all
current signs, but the continuity equation retains that convention.

Write h=H_W/m. Its degree-twelve annihilator has roots
`plus/minus(sqrt(r) plus/minus 1/2)`, r=3,7,13. The h² factors are

```
q_r(z)=z²-(2r+1/2)z+(r-1/4)².
```

Construct rational Chinese-remainder projectors E_r(z) modulo
`Q=product_r q_r`. They satisfy `E_r²=E_r`, sum E_r=1 and
E_r E_s=0 for r!=s. The positive inverse square-root on the two roots
of q_r is the linear polynomial

```
f_r(z)=(3r+1/4-z)/(2sqrt(r)(r-1/4)).
```

Substitution at `(sqrt(r)+/-1/2)²` gives exactly
`1/(sqrt(r)+/-1/2)`. These denominators are positive. Let M_r be the
rational matrix h times the numerator of E_r f_r, leaving the factor
1/sqrt(r) outside. Exact matrix multiplication gives

```
M_r M_s=delta_rs r E_r(h²),
sign(h)=sum_r M_r/sqrt(r),
P_-=(I-sign(h))/2.
```

The spectral sector ranks are8,16,8, and tr M_r=0. Thus P_- is the
complete rank-sixteen negative projector. Each `tr(M_r J_phase,e)` is
exactly zero, proving all twelve filled-vacuum current means vanish.
No polynomial truncation, eigenvector solver or numerical spectral
threshold is used.

For a full finite fermion Fock vacuum filled on P_-, the canonical CAR
contraction is

```
<J_e J_f>-<J_e><J_f>=Tr(P_- J_e (I-P_-) J_f).
```

All32 fermion modes and all occupations remain in the source. This
formula evaluates a full-Fock matrix element without enumerating its
2^32-state matrix. It deletes no pair or multi-pair intermediate sector
from the coupled Hamiltonian.

The actual photon covariance is
`C_Q=P4/4+P6/(2sqrt6)`. After contracting its symmetric edge indices,
the terms linear in sign(h) cancel by cyclicity. Hence

```
sigma²/m²=(1/4) sum_ef C_Q,ef
  [Tr(J_e J_f)-sum_rs Tr(M_r J_e M_s J_f)/sqrt(rs)].
```

The executable uses the direct four-index contraction of the sparse
J matrices, distinct from the interval producer's rounded P*J support
algorithm. All contractions are rational before the explicit radicals
are inserted. The result is

```
sigma²/m²=965/1456+5sqrt6/36-sqrt14/48+sqrt21/336
           +sqrt91/546+7sqrt546/1872
         =1.0435157377429135586350084518253....
```

Outward rational radical isolation proves this value greater than
`(2043/2000)²`. The independent binomial-projector calculation encloses
the same value in its broader interval. Both checks use the literal
finite photon two-word vacuum moment. The full cutoff commutators,
nonlinear Peierls remainder and all subsequent many-body evolution are
handled by the common finite-source proof, not by a Gaussian transfer.

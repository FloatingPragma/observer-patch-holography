# Compact defects, native memory and interacting spectral readouts

## One spatial band and one nonlinear source

The source is the prescribed compact opposite-degree texture of
[compact_defect_pair](../compact_defect_pair/PROOF.md), with the actual
two-spin hopping, exterior mass `M>=8`, half-separation `a>=1` and an exterior
cube containing the constant outer layer. Its four-dimensional isolated band
has isometry `W`, compressed kinetic operator `A=W* h W` and pairing matrix
`S=W* C conjugate(W)`. The source identities are

```math
hW=WA,\qquad C\overline W=WS,\qquad
S^T=-S,\quad S^*S=I,\quad AS=-S\overline A.
```

The complete native internal representation has 45 annihilator slots: fifteen
states and three declared families. Set `E=W tensor I45`, `P=EE*` and
`N_P=dGamma(P)`. For each complete native memory record `r`, use the same
charged pairing kernel `K_r` as the memory construction, retaining the quark
design and all disabled records. Add the chosen interaction

```math
V=\lambda(N_P-90)^2,\qquad \lambda\ge0.
```

There are 180 localized annihilators and 360 redundant Nambu coordinates.
On the displayed radius-three cube the full field has 30,870 annihilators
and 61,740 Nambu rows. None of these dimensions is a particle species count.
The reference energy unit is the declared `v(P_pixel)` pulse unit; the pixel
capacity differs from the band projector `P` in this proof.

Every record's quadratic Nambu operator reduces the same localized sector:

```math
B_{r,W}=\begin{pmatrix}
A\otimes I_{45}&S\otimes\overline K_r\\
(S\otimes\overline K_r)^*&-\overline A\otimes I_{45}
\end{pmatrix}.
```

The particle-hole identity ensures that pairing has no cross term between
the localized and complementary annihilators. All Hamiltonian terms are even,
so their Fock factors commute and

```math
H_r^{\rm full}=H_{r,W}\otimes I+I\otimes H_{r,\perp},\qquad
H_{r,W}=H_{r,W}^{(2)}+\lambda(N_W-90)^2.
```

This is exact for each record, including noncommuting record kernels. It
retains the nonzero compact kinetic matrix `A`; setting it to zero changes
the operator. Tensoring with a parity-even complementary state gives a state
inclusion for the localized algebra. The construction is a specified source
interaction with a band projector. It is not a strictly onsite interaction
or an autonomous texture formation law.

Every spatially uniform native current `I tensor G` commutes with `P`, and
its second quantization commutes with `V`. Thus all twelve currents act on
the same nonlinear core. A fixed Higgs pairing does not commute with all
twelve currents. The actual charged Higgs tensors obey `QK+KQ^T=0` for
`Q=3(B-L)` with charges `(1,-1,-1,-3,3)` on `(Q,u^c,d^c,L,e^c)`, so the
charged Hamiltonian conserves this charge. Its full compact Fock range is
`[-144,108]`. Neutral Majorana pairing is outside that conservation claim.

The compact texture and a uniform north reference have the same pairing
law. Choose no localized quartic in the uniform reference. Since
`0<=N_P<=180`, the interaction norm is `8100 lambda`; the existing texture
bound, repeated on all 45 internal slots, gives

```math
\|H_{\rm compact}-H_{\rm uniform}\|
\le 90M\bigl((2a+3)^3-1\bigr)+8100\lambda.
```

The same bound controls the difference of their finite-volume ground
energies and is uniform in exterior volume. It includes the explicit
reference choice and is not a sum of constituent particle rest masses.

There is also a full-core variational energy bracket, without replacing its
vacuum by a four-mode slice. Let `m_j`, `j=1,...,21`, list all charged singular
values with color multiplicity. Each has two pairing blocks at each spatial
energy `e_i`. Their individual quadratic spectra are `{+/-m_j,+/-e_i}`;
the twelve neutral modes contribute ground energy `-3(e1+e2)`. Hence

```math
E_{\rm quad}=-2\sum_{j=1}^{21}
 [\max(m_j,e_1)+\max(m_j,e_2)]-3(e_1+e_2).
```

Since `V>=0`, `E_quad` is a lower bound for the interacting ground. A
quadratic product ground has mean `N=90`; each even pairing block has
number variance one and each odd blocked pair has zero variance. If `b`
blocks are even, that trial state has energy `E_quad+lambda b`. A fully
blocked state has exact number 90 and energy `-45(e1+e2)`. Therefore

```math
E_{\rm quad}\le E_{\rm full}\le
\min\{E_{\rm quad}+\lambda b,-45(e_1+e_2)\}.
```

The source calculation retains all 180 modes in this bracket; the inherited
`lambda=1/128` witness has 28 even quadratic blocks. These bounds do not
determine the full interacting vacuum or its excitation spectrum.

## Complete-record dynamics with a bound uniform in exterior volume

The compact spectrum is `{-e2,-e1,e1,e2}`, where
`0<e1<=e2<=r=b^2/(M-4)` for the certified source and `b^2=7/4` at `a=1`
(`b^2=3/2` otherwise). Consequently
`||dGamma(A tensor I45)||=45(e1+e2)<=90r`. A pairing operator on 180 modes
has norm at most `90||K_r||`, and the declared quark and memory kernels obey
`||K_r||<=1+g` for feedback gain `g`. Therefore every record satisfies

```math
\|H_{r,W}\|\le M_W:=90r+90(1+g)+8100\lambda.
```

For the inherited forward source and its first memory inverse, rational
mass and inverse-probability enclosures prove `g<1` and `p<1/12`.
The receipt's bound with those upper values therefore applies to that
actual conditional memory design, including its duty factor.

The native reset law is `T=(1-p)R+pS`, with `R` a fixed reset and `S`
stochastic. Its trace-norm contraction on zero-sum memory blocks is at most
`p<1`. The field-independent probe probabilities and the complete transition
law are exactly those of the
[memory realization](../mcpr_native_memory_bridge/PROOF.md). Applying its
Doeblin averaging argument with the displayed compact norm bound, step
`h=t/n`, arbitrary initial classical-quantum memory correlations and an
untouched finite reference gives

```math
\|\rho_{W,n}(t)-e^{-it\bar H_W}\rho_W(0)e^{it\bar H_W}\|_1
\le \delta_n:=
{4M_Wt\over n(1-p)}+
\left({8\over1-p}+4\right){M_W^2t^2\over n}.
```

For stationary product initial memory the first term vanishes. The mean
Hamiltonian contains the actual memory-derived `K_bar`, the same `A` and
the quartic `V`. The estimate sums all fine outcomes. Retaining the full
history register gives trace-norm error at most `sqrt(2 delta_n)`, capped
by two, by the purification/Jensen argument in
[interacting_defect_response](../interacting_defect_response/PROOF.md).
For correlated initial memory, the ideal law retains its conditional
history probabilities and evolves each initial conditional core state by
the mean unitary. The result is uniform in exterior volume on the same
180-mode algebra. It gives no uniform rare-conditioned-history estimate,
intervening field-measurement guarantee or lattice-spacing continuum limit.

The numerical control retains all 512 fine records per memory input, with
both enabled and disabled actions. Its sixteen-state probe has the compact
kinetic splitting and a noncommuting quartic. This probe checks the general
averaging calculation; it is not the full core or a reducing sector of
every family-dependent pulse.

## Exact interacting spectrum on a charged reducing sector

Choose a positive-energy spatial vector and its particle-hole partner for
each `e_i`. Charge-preserving singular bases for the neutral-Higgs kernel
give 84 charged pairing blocks and twelve unpaired neutrino modes. For one
charged-lepton singular value `m>0`, select a charged pair associated with
`e1` and one associated with `e2`. Fix one occupation in each of the other
82 charged blocks, and fill six neutral energy eigenmodes. Pair creation
and annihilation vanish on these spectator states. Their total number is
88 and their fixed kinetic energy is subtracted. The mean Hamiltonian and
quartic then restrict exactly to

```math
H_4=e_1(n_1-n_2)+e_2(n_3-n_4)
+m(a_1^\dagger a_2^\dagger+a_2a_1
   +a_3^\dagger a_4^\dagger+a_4a_3)
+\lambda(N_4-2)^2.
```

This sector is not asserted to contain the full compact core ground and is
not invariant under every unaveraged family pulse. It carries charges
`(-3,3,-3,3)` and retains both spatial pair energies.

Each pair's fermion parity is conserved. When both pairs are even, the
four eigenvalues are `0`, `4lambda`, and
`2lambda +/- 2sqrt(lambda^2+m^2)`. With both pairs odd the four eigenvalues
are `+/-e1 +/-e2`. With one pair odd they are
`lambda +/-e_i +/-m`, for `i=1,2`. These sectors exhaust all sixteen states;
the independent symbolic characteristic-polynomial check and numerical CAR
diagonalization agree. For nonnegative spatial energies the ground energy is

```math
E_0=\min\{2\lambda-2\sqrt{\lambda^2+m^2},
          -e_1-e_2,\lambda-m-\max(e_1,e_2)\}.
```

Thus a compact spatial splitting can change the ground's pair-parity
sector. For a unique blocked ground `|01;01>`, a chosen positive-energy
annihilator has the complete addition-plus-removal spectral measure

```math
\mu(\omega)=\tfrac12\delta\bigl(\omega-(e_1+\lambda-m)\bigr)
           +\tfrac12\delta\bigl(\omega-(e_1+\lambda+m)\bigr).
```

For the unique paired ground, write `r_m=sqrt(lambda^2+m^2)`. Its four
poles are `2r_m-lambda +/-m +/-e1`, each with residue
`(1-m/r_m)/4` for the `+m` sign and `(1+m/r_m)/4` for the `-m` sign.
Coincident poles have their weights added. All outcomes obey the unit
addition-plus-removal sum rule. These are excitation energies read from
the same Hamiltonian, with the supplied pulse clock.

At `M=64,a=1`, the finite positive spatial energies are approximately
`0.003905899328` and `0.003906437151`. Keep the inherited interaction choice
`lambda=1/128` and the forward pixel-memory mass design. Exact rational
mass enclosures are obtained by cubing its positive mass formula, bounding
the exponential and trigonometric factors, and enclosing the cube root by
rational bisection. The complete dyadic compact spectrum and its exterior
error give intervals for both `e_i`. These enclosures prove that all three
charged-lepton slices have the unique blocked ground, in finite and infinite
exterior volume. Their numerical positive-energy-orbital readouts are:

| Source mass coefficient | Two single-fermion poles | Residues |
| --- | --- | --- |
| `0.000002072030404` | `0.01171632730`, `0.01172047136` | `1/2`, `1/2` |
| `0.0004284302055` | `0.01128996912`, `0.01214682953` | `1/2`, `1/2` |
| `0.007205207119` | `0.004513192209`, `0.01892360645` | `1/2`, `1/2` |

The compact source's finite-to-infinite spatial energy error `delta_L`
implies a four-mode Hamiltonian error at most `2 delta_L`, ordered energy
errors at most `2 delta_L` and gap errors at most `4 delta_L`. In the
certified blocked regime each displayed pole has error at most `delta_L`.
The retained exact intervals include source-matrix rounding and all
exterior modes; the displayed decimals are numerical summaries.

Setting `e1=e2=0` changes every retained slice ground to the paired sector
and changes its poles and residues. This explicit control rules out
attaching the radial interacting readout unchanged to the compact band.
The operator join does not identify the pixel-dependent coefficients or
these slice poles with nature's particle masses. The mass-law architecture
has target-informed provenance, a physical unit/dispersion map is absent,
the full 180-mode interacting vacuum is not determined by these slices,
and their spectral measure has not been inserted into a common
electromagnetic feedback calculation. Interacting chiral continuum fields
and responsive Einstein geometry require further common-source results.

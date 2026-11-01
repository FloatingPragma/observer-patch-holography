# The complete ground-family spectrum and its internal spin

Every linear creation probe on every vector in the original cube's complete
five-dimensional ground space has support in the same 78 positive energy
lines. The line `lambda=sqrt(3)/16` has total hidden spin `5/2`; the other
77 lines have total hidden spin `3/2`. Keeping all five ground inputs and
all sixteen creation modes gives a 1256-dimensional reachable subspace of
the full 4368-dimensional `N=5` source. This closes a finite particle-search
loophole: choosing a different ground polarization or coherent ground
superposition does not reveal additional creation energies.

These are **total internal spin labels in the presence of a spin-two
ground**, not the intrinsic spin of a relativistic particle. The field
operator is an internal spinor. Neither this spinor nor the final total
spin is identified here with physical Lorentz spin or a measured species.

The source is unchanged: the open `2x2x2` Weyl/contact cube at `c=1/8,U=1`,
zero chemical potential, and ground energy `E_G=-sqrt(3)/4`. All arrays use
the exactly equivalent hidden-spin frame. The actual onsite Fock frame
transports the result back to the original physical source, including
intersite spin coherences. The previous [matrix response](../matrix_spectral_response/README.md)
and [ground theorem](../interacting_vacuum/README.md) supply the authenticated
78 residues and complete spin-two ground space. Their source bytes are
checked before the new result is admitted.

## Separating the two reachable spin sectors

Write `G_0` for the normalized highest-weight ground, `V_up,x=c_x,up^*G_0`
and `V_down,x=c_x,down^*G_0`. Let `T=S_plus` from the **entire**
`(N_up,N_down)=(4,1)` sector to `(5,0)`. The actual 56-by-560 integer
source matrix gives

```
T V_down=V_up,       T T^*=5 I_56,
H_(5,0) T=T H_(4,1),
V_up^* V_up=P_plus,  (H-E_G)V_up=lambda V_up,
P_plus=(I+A/sqrt(3))/2.
```

The `M=5/2` sector is pure spin `5/2`. At `M=3/2`, only spin `5/2` and
`3/2` are possible. Their exact projector and seed decomposition are

```
P_(S=5/2)=T^*T/5,
V_high=T^*V_up/5,     V_low=V_down-V_high,
T V_low=0,            V_high^*V_high=P_plus/5.
```

The high seed has only energy `lambda`. The retained matrix theorem gives
exactly the same residue `P_plus/5` at that line: its local weight is
`1/10`, its band is plus, and its rank is four. Positivity then leaves zero
low-spin residue at `lambda`. Every other retained residue belongs to
`V_low`, so is spin `3/2`. This identifies the **reachable parts** of the
energy eigenspaces; additional dark states at the same energies are not
excluded.

## All ground superpositions and all creator modes

Use the unnormalized ground basis `G_k=S_minus^k G_0`, `0<=k<=4`, with
squared norms

```
n_k=k! 4!/(4-k)! = 1,4,24,144,576.
```

The code constructs every vector in the full occupation basis, checks its
original Hamiltonian eigen-equation, and verifies the ladder coefficients.
The CAR give, for each site,

```
c_down^* G_k = S_minus^(k+1) V_up/5 + S_minus^k V_low,
c_up^* G_k   = (5-k) S_minus^k V_up/5 - k S_minus^(k-1) V_low.
```

The second term is zero at `k=0`. The implementation checks all eighty
creator columns with integer arithmetic, including their signs. The
spin-`5/2` ladder terminates after six weights and the spin-`3/2` ladder
after four. Since the original Hamiltonian commutes with these ladders,
no new energy is introduced. Linear combinations cover every complex
coherent ground vector and every linear CAR probe, rather than only
rotations of `G_0`.

For a spin `j` highest-weight vector, the squared norm multiplier after
`k` lowering steps is `k! (2j)!/(2j-k)!`. Dividing by `n_k` gives the
following coefficients multiplying the full spatial spectral measures of
`V_up` and `V_low`:

| Ground | Creator | High measure coefficient | Low measure coefficient |
|---|---|---:|---:|
| `G_k/sqrt(n_k)` | up | `(5-k)/5` | `k/4` |
| `G_k/sqrt(n_k)` | down | `(k+1)/5` | `(4-k)/4` |

These coefficients apply to off-diagonal spatial entries as well as
local probabilities: lowering acts on both sides of each spectral Gram
matrix and commutes with its energy projector. The two total-spin
summands are orthogonal.

## Complete reachable and dark dimensions

For `N=5`, the positive-weight dimensions are
`d_(5/2)=56`, `d_(3/2)=560`, and `d_(1/2)=1568`. Differences give the
multiplicities of the irreducible spin representations. The code checks
`[S_plus,S_minus]=2S_z` and the Hamiltonian intertwiner on **every** weight
space of the full `N=5` sector.

| Total hidden spin | Full multiplicity | Reachable multiplicity | Reachable dimension | Dark dimension |
|---|---:|---:|---:|---:|
| `5/2` | 56 | 4 | 24 | 312 |
| `3/2` | 504 | 308 | 1232 | 784 |
| `1/2` | 1008 | 0 | 0 | 2016 |
| Total | | | **1256** | **3112** |

At `lambda` there are four reachable spin-`5/2` copies. At each of the
other 77 lines there are four reachable spin-`3/2` copies. Distinct
energies and distinct spin weights are orthogonal, so these dimensions
add. The generated subspace is reducing for both the Hamiltonian and
spin action. Its orthogonal complement is dark for every single-creator
probe from the complete ground space and all subsequent time translates.
The dark-sector energies and full energy degeneracies are not computed.

Annihilation is also closed. With `B_x=c_x,up G_0`,

```
c_up G_k=S_minus^k B,
c_down G_k=k S_minus^(k-1) B.
```

The highest seed fills three negative-band orbitals and has energy
`-3lambda`, hence removal costs `lambda` relative to `E_G=-4lambda`.
Its four spatial modes generate sixteen spin-`3/2` states. The complete
`N=3` sector has dimension 560 and spin multiplicities 56 at `3/2` and
168 at `1/2`; the remaining 544 states comprise 52 spin-`3/2` copies and
168 spin-`1/2` copies. The code checks all eighty annihilator identities.

## A rotation-invariant reference and its full matrix readout

Take the normalized full-ground projector `rho_G=P_G/5`. Averaging the
preceding table gives high coefficient `3/5` and low coefficient `1/2`
for either spin. Cross-spin entries vanish in the hidden frame by spin
invariance. With `G_down` the previous highest-ground down-addition
matrix, the exact averaged addition and retarded responses are

```
G_add_bar(z) = [G_down(z)+P_plus/(z-lambda)]/2 tensor I_spin,
G_ret_bar(z) = [G_down(z)+P_plus/(z-lambda)
                          +P_minus/(z+lambda)]/2 tensor I_spin.
```

The removal pole is negative in the retarded convention; the positive
energy needed to remove a fermion is `lambda`. Per site and fixed spin:

- Addition at `lambda` has weight `3/10`; the other 77 weights are `w_i/2`.
- The total addition weight is `3/4`. Conditional on addition, the weights
  are `2/5` at `lambda` and `2w_i/3` elsewhere.
- Removal has weight `1/4` at retarded pole `-lambda`.

Every addition weight is strictly positive. `results.json` retains all
78 energy intervals and both the unconditional and conditional weight
intervals, with their internal spin labels. There are 79 distinct signed
retarded poles, each with matrix residue rank eight. Their minimal
16-column mixed-state GNS realization therefore has dimension 632.
This is a different object from the 1256-dimensional physical `N=5`
range that retains all five ground inputs separately.

Direct full-source integer Gram matrices independently check all 16-by-16
entries of the zeroth and first moments, including off-diagonal spatial
coherence:

```
M0_add=(3/4)I+(sqrt(3)/12)(A tensor I_spin),
M0_remove=(1/4)I-(sqrt(3)/12)(A tensor I_spin),
M0_ret=I,              M1_ret=(A tensor I_spin)/16+I/4.
```

For the original source frame `c_x=R_x^* d_x`, a hidden scalar spin block
`g_xy I` becomes `g_xy R_x^*R_y`. The exact source one-body frame and
Hamiltonian are checked. Intersite physical spin mixing is retained;
only the same-site block stays proportional to the spin identity.

## What this closes in the particle search

Every positive single-addition energy from any ground choice lies inside
the previously certified envelope. Ground choice, polarization, linear
probe choice, and conditioning on a nonzero addition outcome cannot
increase its maximum-to-minimum ratio. Thus the previous exact
common-positive-scale charged-lepton dictionary rejection applies to
this entire ground family. The invariant mixed reference also shows that
no previously retained line disappears from the full family readout.

This conclusion does not turn the finite gaps into relativistic masses,
exclude nonlinear composite probes or other reference states, or rule out
an infinite-volume/continuum particle sector. The native preparation of
the full mixed reference and complete finite-shot line resolution are not
executed by this package. No source parameter is fitted to measured data.

Run the source-bound replay with

```
python -m model.spin_resolved_spectrum --check model/spin_resolved_spectrum/results.json
```

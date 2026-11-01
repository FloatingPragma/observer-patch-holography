# A global interacting cube ground and positive charged excitations

The unchanged open `2×2×2` Weyl/contact source at `c=1/8,U=1` has **global ground energy `−√3/4` and exactly five ground states**, all in `N=4` and forming one hidden-spin `S=2` multiplet. This comparison covers all 65,536 Fock states. The empty Fock state is not this ground. No chemical potential, energy recentering, new mass term, or Hamiltonian replacement is introduced.

Above a specified pure ground `G`, the original local CAR probes have positive addition/removal energies. Write `λ=√3/16≈0.108253175473`. The following values are in the existing model's energy units:

| Local probe at the origin | Final number | Spectral weight | Positive gap above `G` |
|---|---:|---:|---:|
| `c_up G` | 3 | `1/2` | `λ` |
| `c_up† G` | 5 | `1/2` | `λ` |
| `c_down G` | 3 | `0` | absent |
| spin-`5/2` component of `c_down† G` | 5 | `1/10` | `λ` |
| eta-`5/2` component of `c_down† G` | 5 | `1/10` | `1+λ≈1.108253175473` |
| remaining component of `c_down† G` | 5 | `4/5` | positive measure, all gaps `>7/1024` |

The last three rows are an orthogonal decomposition of the down-addition measure. The two displayed `1/10` weights are exact symmetry-resolved contributions; without additionally resolving the remaining measure they are lower bounds on the total weights of any coincident atoms. The residual has strictly positive variance, so it is not a single unresolved sharp line. These are finite-volume excitation energies, **not relativistic rest masses or Standard Model assignments**. A native preparation and detector for this new ground have not been executed.

## Source and exact ground construction

`ground.py` consumes the original `cube_eta_protection.cube.operators()` full-Fock action and its previously constructed onsite unitary frame. That frame changes coordinates, not the Hamiltonian. The one-body hidden-frame action is the spin-independent signed cube adjacency `A/16`, with `A²=3I`. With even/odd sublattices,

```
A = [0 B; B^T 0],  BB^T=3I4,
K = [√3 I4; −B^T],  A K=−√3 K,  K†K=6I4.
```

Filling the four orthonormal negative-band orbitals `K/√6` with hidden spin up gives `G`. Its 62 nonzero occupation coefficients are explicitly retained as `(a+b√3)/36`. It is normalized, contains no doublons, and satisfies the original physical-frame eigenvalue equation after the existing Fock-frame adjoint. At the origin the onsite frame is the identity, so the table's probes are original physical CAR operators.

On `N=4`, positivity of the contact term and the one-body lower bound imply `H≥−4λ`. Equality requires four negative-band orbitals and zero contact. The negative band has eight spin/orbital modes, so its four-particle space has dimension70. Stacking all eight local pair-annihilation maps gives an exact integer matrix of shape224×70, equal to six times the normalized annihilation map. In the five hidden-spin weights its dimensions are `1,16,36,16,1`, and its ranks are `0,15,35,15,0`. Nonzero integer minors modulo101 give the rank lower bounds; five explicit spin-lowered vectors give matching upper bounds over the complex numbers. Thus the entire `N=4` ground kernel is exactly one `S=2` multiplet. Omitting odd-sublattice contacts increases the kernel to dimension16, a retained negative control.

## All-number comparison without a floating eigenvalue premise

For `N≤3`, `H≥−Nλ>−4λ`. For `N≥9`, the occupancy inequality `D≥N−8` gives

```
H ≥ −Nλ+(N−8) > 0.
```

Only `N=5,6,7,8` need additional bounds. The source commutes with hidden `SU(2)`. Every integer-spin representation contains weight0 and every half-integer representation contains weight1/2. Therefore every fixed-number eigenvalue occurs in the balanced hidden-spin block; this is not a projection that discards another multiplet's minimum. The four dimensions are `1568,3136,3920,4900`.

`global_ground.py` uses floating Cholesky only to propose an integer matrix `Q`: round a factor of `H+27/64 I` to20 dyadic bits. The proof accepts the proposal only after **exact integer** matrix multiplication verifies

```
H = Q Q^T / 2^40 −27/64 I + R,
||R||₂ ≤ ||R||∞ < 1/512.
```

All product/accumulation bounds are checked before `int64` arithmetic, and the residual row-sum bound is checked separately. No inverse, nonzero pivot, Cholesky correctness assertion, or floating eigensolver admission is needed: `Q Q^T` is positive semidefinite for any integer `Q`. Duplicate sparse entries are canonicalized before residual addition. The replay may produce a different valid proposal on another platform; its retained claim is the same exact rational ceiling, not a particular floating factor's hash. Large temporary factors are not repository artifacts.

Consequently `H_N>−217/512` for these sectors. The exact square comparison

```
√3/4 > 441/1024 > 217/512
```

separates them from `G`. Every `N≠4` lies strictly more than `7/1024` above the ground. This is a conservative **charged-sector** bound; it is not an asserted bound on neutral excitations within `N=4`. Together with the exact five-state kernel it proves the global ground statement. Every density matrix supported on that ground space is passive for this finite Hamiltonian and satisfies its positive-energy spectral condition. Uniqueness, an infinite-volume limit, and Lorentz-invariant vacuum structure remain open.

## Positive Lehmann measures and the genuine interaction response

`charged.py` evaluates the vectors and moments over `Q(√3)`, with no eigensolver. `c_up G` removes a negative orbital; `c_up†G` adds a positive orbital. Both have norm squared1/2, contact zero, and energy `−3λ`. The down-removal vector vanishes. The down-addition vector `v=c_down†G` has norm1 and contact expectation1/2, so this channel actually encounters the contact interaction.

Let `S±` be global hidden-spin ladders and `η±` the source eta ladders. The exact vectors

```
v_s = S−S+v/5,       v_eta = η+η−v/5,
v_r = v−v_s−v_eta
```

are mutually orthogonal, with squared norms `1/10,1/10,4/5`. The first two satisfy source energy equations at `−3λ` and `1−3λ`. The complete down-addition gap measure is therefore

```
μ_down = (1/10)δ_λ + (1/10)δ_(1+λ) + (4/5)μ_r.
```

All support is positive by the global ground comparison. Exact complete-measure moments are

```
E[ω]=1/2,  E[ω²]=131/256,  Var(ω)=67/256.
```

For the normalized residual,

```
E_r[ω]=1/2−√3/64,
E_r[ω²]=131/256−√3/64,
Var_r(ω)=1069/4096>0.
```

The contact-deleted comparator, applied to this same vector, has mean gap zero instead of1/2. The actual interacting result is not obtained by relabeling the previous empty-sector Weyl modes. Single addition/removal changes the conserved model number by `±1`, has odd CAR grading relative to this even ground, and does not supply an independently calibrated electric charge or physical spin.

## Same-process repair attachment

The already declared `cube_eta_protection.repair` has24 active columns supported entirely in `N=6`. Each of its130 complete outcomes is either an active native Kraus operator or the outside projector `Q`. Hence on `N3⊕N4⊕N5`, a6748-dimensional space, all active outcomes vanish and `Q` is identity. This includes cross-number coherences and arbitrary untouched references. The original Hamiltonian preserves number, so any finite interleaving of its evolution and these repairs preserves that statement and the normal field evolution.

This is **passive outside-service protection** supplied by the actual full-Fock completion. It is not active error correction, autonomous preparation of `G`, or invariance under arbitrary repair laws. The full outcome dictionary remains present; the outside record is certain for ideal states in this sum.

## Replay and remaining work

```
.venv/bin/python -W error -m pytest -q model/interacting_vacuum
.venv/bin/python -W error -m model.interacting_vacuum --write
.venv/bin/python -W error -m model.interacting_vacuum --check
```

The receipt binds the source, exact checks, this proof, and the focused paper fragment. Its claim is a globally certified finite interacting reference plus positive charged spectral data, not complete resolution of the down-addition spectrum. The new ground's native source/detector, a controlled growing-volume interacting vacuum, stable relativistic massive poles, physical charge/spin dictionaries, independently fixed units, and Einstein coarse-graining premises remain required. Earlier protected excited-reference spectra remain valid as excited-reference experiments; their negative-frequency obstruction does not apply to this newly proved global ground.

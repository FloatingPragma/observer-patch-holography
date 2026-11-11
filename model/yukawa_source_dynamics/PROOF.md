# Complete family response and canonical normalization

## Native completion and its full eigenspaces

Let `T` be the twelve-port first-moment response obtained by replaying every
accepted seam, placement coin and terminal selector timeout of
`native_family_dynamics.proposal_branches(32)`. The 66-state instrument has
unit column sums. Its first moments obey `loads*T66=T*loads`. The following
resolvent calculations concern this one-proposal marginal. They neither
restart a terminal timeout nor identify its attempt count with physical time.

Put `H=I-T`. Every accepted endpoint averaging map is a positive contraction,
so `H` is positive semidefinite. The protected core has two isolated ports;
the other ten ports form a connected repair graph. Thus `rank(H)=9`, and its
kernel is spanned by the two core indicators and the complement indicator.
The source occupation matrix is exactly

`D=diag(29/90,29/90,61/450,...,61/450)`.

Every endpoint branch preserves these two occupation sectors, so `[D,H]=0`.
All functions of this specified pair admit a common full-mode eigenbasis.
This statement does not classify the algebra of all native interventions.

Let the existing slow-band frame be `B`, with `B^T B=g I3`, and put
`P=BB^T/g`, `Q=I-P`, `F=B^T D B/g`, `A=B^T H B/g`. The compressed matrices do
not commute. Direct expansion gives the exact identity

`[F,A]=B^T(H Q D-D Q H)B/g`.

The full commutator vanishes; the displayed terms contain every intermediate
port excluded by the compression. Consequently the exact matrix units of
the compressed pair do not establish a noncommuting full source vertex.
The owning compressed-algebra certificate retains its algebraic domain.

## Effective matrix and residue metric

Choose any complete complement frame `Qb` with `B^T Qb=0`. Define

`C=Qb^T(H+J/12)Qb`, `V=B^T H Qb`, `G=Qb^T Qb`,

where `J` is the all-ones matrix. The added projector lifts the sole
complement null vector, the uniform mode. Its coupling to `B` vanishes, so
it contributes no virtual amplitude. The exact zero-frequency Schur matrix
and its residue metric are

`S=A-V C^-1 V^T/g`,

`Z=I3+V C^-1 G C^-1 V^T/g`.

Writing `L=B-Qb C^-1 V^T` proves both identities
`S=L^T H L/g` and `Z=L^T L/g`. Hence `Z` is positive and is the required
kinetic normalization for this linearized effective response. The canonical
coefficient is `Z^-1/2 S Z^-1/2`, with exact rank one. Its diagnostic
spectrum is `(0,0,0.03381748039929137)`. The unnormalized nonzero entry of
`S` is approximately `0.0371744984`. Neither number is equated with a
fermion mass or a pole of the full frequency-dependent resolvent.

The source kernel projector `P0` gives the exact zero-pole residue

`B^T P0 B/g=diag((5-sqrt(5))/20,3*(5+sqrt(5))/50,0)`.

Its rank is two and `S` annihilates its image. Direct diagonalization of the
complete twelve-port operator supplies nine positive poles, each with its
full three-by-three residue. Their sum, including the zero-pole residue,
is `I3`. The numerical entries are diagnostics; rank, normalization and
conservation statements are exact. A Schur resolvent check at nonreal and
negative arguments independently agrees with inversion on all twelve ports.

## Complete cyclic source amplitudes

The actual native half-turn product `R` obeys `R^3=I` and `tr(R)=0`.
All 128 native record probabilities are replayed, with the original
exclusive routing. Applying that routing to the family register is the
declared service inherited from `koide_source`; native record probabilities
alone do not force it. Their weights for `I,R,R^-1` are
`31/36,5/72,5/72`. For Kraus amplitudes
`K_r=sqrt(w_r) R^e(r)`, the complete Stinespring map therefore satisfies

`V^dag V=sum_r K_r^dag K_r=I3`.

Its three singular values are one. A common unitary change of record basis,
including the inherited Fourier analyzer, preserves this Gram matrix.
After tracing the complete record, the exact channel is

`E=(19/24) Id+(5/24) Twirl_C3`.

The twirl is an orthogonal projector of operator-space rank three. For every
integer `n>=1`,

`E^n=Twirl_C3+(19/24)^n (Id-Twirl_C3)`.

The principal generator is `log(19/24)*(Id-Twirl_C3)`. It has no Hamiltonian
commutator component, modulo an unobservable scalar. The same formula with
real `t>=0` is a completely positive continuous embedding. Thus choosing
neither a repetition count nor a continuous time produces a coherent
three-mass vertex in this family. Multiplying any fine branch by a scalar
phase cancels between `K_r` and `K_r^dag`; this exclusion covers every such
phase without a phase scan.

The selected default analyzer effect has eigenvalues

`(5/72-5 sqrt(2)/108,5/72-5 sqrt(2)/108,5/36+5 sqrt(2)/54)`.

All eight effects sum to identity. A feedback interaction using one analyzer
projector is a specified additional operation. Its repeated pair cannot be
interpreted as three distinct charged generations. The construction uses
no measured masses, fitted phases or chosen iteration count.

## Quantum records and relational phase

Retaining quantum records avoids the trace in the preceding channel. In the
actual three character sectors, the normalized record vectors are

`|e_k>=sum_r sqrt(w_r) omega^(k e(r)) |r>`, `omega=exp(2 pi i/3)`.

Their complete Gram matrix has diagonal one and off-diagonal `19/24`.
The phase of the invariant cyclic product
`<e0|e1><e1|e2><e2|e0>=(19/24)^3` is exactly zero. Common record unitaries
preserve it. For a uniform family input, the three Schmidt probabilities
are `31/36,5/72,5/72`. Canonical normalization of this full-rank Gram gives
an orthonormal record embedding, with no omitted quantum output.

The exact inherited weak-feedback word is also evaluated at its specified
`z=1`, half-angle tangent `1/64`. Every one of its eight Kraus records is
retained as a coherent row. The six columns are fixed by Nambu `X` and the
actual cyclic characters. Their Gram has unit diagonal. The oriented
three-family invariants in both Nambu sectors have exactly zero imaginary
part. This is an evaluation of that word, including its specified analyzer
projector; arbitrary coherent relational measurements are not excluded.

The completion response and cyclic routing contain no source attachment
that distinguishes physical up, down and charged-lepton poles. Their
normalization and full-outcome spectra are constructive finite results.
Additional matter dynamics can have different spectra and must be evaluated
on its own complete state and clock domain.

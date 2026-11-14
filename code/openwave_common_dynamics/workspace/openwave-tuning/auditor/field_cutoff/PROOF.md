# Cyclic registers and the integer field process

The source has matter transport with integer edge current, electric phases
`exp(-i alpha sum_e F_e^2/2)`, and magnetic factors
`exp(i beta_f (U_f+U_f*)/2)`. The face shift adds the actual integer column
of the declared face matrix. Every column entry has absolute value at most
one. A finite register uses residues modulo even `m`, with centered electric
labels `-m/2,...,m/2-1`. The comparison source uses `ell^2(Z^E)` with the
identical integer currents, face columns, coefficients, order and clock.

Let the initial state have field support `|F_e| <= E0`. Suppose all matter
transport in the observation interval changes each edge by at most `J` in
absolute value, including intermediate transport steps. For the one-particle
three-step source, `J=3` is a valid bound. Define `Lambda` as the sum of
absolute magnetic coefficients over every face and every update in the
interval. Put `d=m/2-E0-J` and assume `d>0`.

Under the centered embedding of the cyclic field basis into the integer
basis, the vector-norm difference between the two outputs is at most

`delta = min(2, 2 sum_{k=d}^infinity Lambda^k/k!)`.

This is uniform over the admitted initial subspace and any attached quantum
reference. Complete terminal measurements, including joint matter, field and
event records, have total-variation distance at most `min(1,delta)`, provided
the two terminal instruments agree under the same embedding. Native compiler
and apparatus errors are additional.

To prove the statement, expand every magnetic exponential in powers of its
bounded shift generator. Retain electric phases and matter transport exactly
between expansion terms. A term containing a total of `k<d` elementary face
shifts stays inside the common register window at every prefix. Its integer
and cyclic executions coincide, including their electric phases. Products
of the exponential series have total norm bound `Lambda^k/k!` at order `k`.
Each source's omitted tail is bounded by the corresponding exponential
tail; the triangle inequality gives the factor two. This proof does not use
a bounded electric Hamiltonian: its diagonal exponential is unitary on both
spaces and agrees on the interior. The source is a finite ordered product,
so no assertion about the domain of a sum of unbounded generators is needed.

For an initial state with probability `epsilon` outside the stated support,
the same argument applied to its two orthogonal components gives the bound
`delta + 2 sqrt(epsilon)`, capped at two. An unbounded family of states near
the moving register wall has no uniform convergence: a positive cyclic
shift takes `|m/2-1>` to `|-m/2>`, while the integer shift takes it to
`|m/2>`. Those outputs are orthogonal, with vector distance `sqrt(2)`.

The executable evaluates the exponential tail with exact rational
arithmetic. Choose an integer `k0 >= d` with `k0+1>Lambda`. Sum the finite
prefix exactly and bound the remaining positive series by
`(Lambda^k0/k0!)/(1-Lambda/(k0+1))`. The input coefficient is the exact
rational represented by the source's binary64 value, so decimal rounding
does not reduce this bound.

The complete sparse numerical control independently implements the actual
one-particle reflection, current transport, one-face magnetic kick and
electric phase. It compares registers of size 16, 32 and 64 through a common
signed-integer output alphabet. No amplitude threshold or success
normalization is used. The size-16 implementation is compared with the
producer's actual source. The size-64 member has its own analytic bound to
the integer process and is not treated as an exact infinite array.

The theorem is a finite-time transfer to an integer compact-rotor field.
It supplies neither a real-valued canonical field nor a Maxwell continuum.
The spatial refinement, prepared quantum state, electric and magnetic
normalizations, weak-angle approximation, interacting matter limit and
physical clock require their specified domains and error estimates. The
homogeneous 208-face source uses all 208 coefficients in `Lambda`; its
cutoff cost cannot be inferred from the executed one-face member. Changing
the register size also changes its native preparation and pulse program.

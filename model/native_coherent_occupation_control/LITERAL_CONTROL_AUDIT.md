# Exact evaluation and the local beam handoff

The declared source law uses exact generator entries in Q(sqrt(5)). The
owning `certified_reference.source` interprets the actual twelve integrated
control outputs as their exact binary64 rational values. It evaluates the
exact full six-level exponential and retains every proper-code complement.
Thus there is an analytic source semantics independent of a host `expm`
approximation. Its literal controls are not symbolic Clifford angles.

The native programme's `source_word` stores every angle as the exact
Fraction of the emitted binary64 value. Admission explicitly checks that
the rational is exactly binary64 representable. The 22-pulse beam contains
dyadic multiples of that binary approximation to pi. Replacing the payload
by symbolic pi would change the word. The coefficient conversion in
`integrated_controls` also has its own literal binary64 values.

`literal_angles.py` first isolates the angle issue: even with exact ideal
base generators, directed 120-bit sine/cosine and square-root enclosures
prove nonzero output from the one-occupation rail subspace to occupations
zero and two. Every off-sector entry has a real or imaginary interval
excluding zero. This rules out identifying this exact-angle word with a
number-preserving Clifford-field target.

`literal_controls.py` checks the full literal source semantics. It reads
all twelve actual integrated controls as exact rationals, combines them
with the exact source generators, and encloses each full six-level
exponential by a 96th-order matrix Taylor series. If B bounds the maximum
absolute row sum, the remaining norm is at most

    B^97 / [97! (1-B/98)].

Every multiplication uses directed 120-bit interval rounding. The exact
anti-Hermitian source law gives a unitary full processor response. All 80
proper-code complement rows across the 22-pulse word remain defined as
first-failure outputs. The success block is composed without normalization
or projection onto an occupation sector.

The resulting successful work operator has, for example,

    -35969104373888654701/2^120
       <= Re U[|00>,|01>]
       <= -35969104373888587730/2^120 < 0.

All eight real/imaginary off-sector entries have certified nonzero
intervals. An independent 100-digit full-source calculation lies inside
each enclosure. Exact evaluation therefore does not make this literal
beam preserve the fixed native occupation number. This is a finite-word
counterexample, not an obstruction to every possible source synthesis.

There is a concrete alternative control direction in the retained source
certificate: `real_12` exactly mixes processor basis vectors 1 and 2. On
the existing proper four-dimensional code, these are |01> and |10>, so
its exponential preserves total occupation for every real parameter.
The retained exact coefficient identity establishes the response direction.
The programme alphabet exposes only X/Y/Z/I/H/J and explicitly rejects an
additional kind. A direct `real_12` service therefore requires an explicit
control-domain and admission attachment. Reading that algebraic identity
does not silently enlarge the admitted source-word grammar or authorize
substituting exact algebraic coefficients for literal binary64 controls.

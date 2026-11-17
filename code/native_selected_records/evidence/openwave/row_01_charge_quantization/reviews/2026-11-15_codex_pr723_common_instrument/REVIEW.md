# Common positive instrument v2: independent same-group review

**Verdict:** the bounded source-interface compatibility component is supported.
This review is by another Codex worker in the same group. It does not satisfy
G6, establish a complete common controller, compile the coherent apparatus from
native operations, or award a physical row.

`SOURCE_BINDINGS.json` binds the exact source, producer audit, amendment,
preregistration and original thirty-seam graph. The reviewed source SHA-256 is
`7ad24094b5fdf527c9b4a34d562d71f6580896c7f8ad26a2a27d630c3509f019`;
the producer audit SHA-256 is
`4b2047513653d720eaeafcbf27aed18df4dc749fe152adc5f5aed784b4b505f2`.

## Independent executable check

Run `python3 test_common_instrument.py` from this directory. It uses only the
Python standard library, checks the input hashes before importing the reviewed
source, and emits `RESULTS.json`-equivalent JSON. The deterministic input is a
normalized complex superposition of all 289 endpoint pairs in `[-8,8]^2`,
with twelve zero spectator loads and seed `72314001`. We test one ticket
`[(0,1)]` and three tickets `[(0,1),(1,0),(12,13)]`.

The expected fine instrument is reconstructed using separate floor/ceiling
arithmetic. Every source-state key and amplitude agrees exactly. The complete
fine and erased instruments have norm one to floating-point precision. Applying
the inverse Householder transformation to the returned success state and every
failure state reconstructs the fine instrument, with maximum absolute error
`7.144028522217822e-17` and `3.5381555418772973e-17`, respectively. The two tests
retain 33 and 73 failure flags, including currents beyond `[-1,1]`. The exact
allowed graph is the thirty cell seams plus the reference seam; nonedges
`(0,2)` and `(0,12)` are rejected.

This check was first performed inline during review and then persisted here;
it is a verification of an existing registered component, not a new prospectively
registered scientific experiment. The saved run adds exact state-key comparison
and independently recomputes returned norms.

## Why the finite instrument is complete

For a fixed ticket, coin and current `j`, an output recovers its input as
`q = q' + j e_a - j e_b`. Thus each fine Kraus branch is a partial injection.
Its coefficient is `1/sqrt(2T)`, and every input has exactly one current for each
of `2T` ticket/coin choices. The fine branches therefore form an isometry on the
declared bounded domain. Averaging preserves `[-8,8]`. The Householder matrix is
unitary on the declared ticket/coin/current `[-1,0,1]` subspace and the identity
on higher current flags. Keeping all its output source kets preserves norm and
the full conditional instrument, including unsuccessful outcomes.

Version 2 resolves the concrete version-1 concerns raised in review: it restricts
addresses to the original graph and supplies the actual failure source states,
rather than only their aggregate probability. The success branch retains the
same anchor and normalization as the original row-7 analyzer.

## Scope that remains open

The source interface has two explicit observation modes. Classical committed
records and coherent erasure require different declared apparatus settings;
committed flags cannot later be silently treated as coherent. Ticket lists,
controller memory, feedback, reference ports, phases, eraser, analyzer and
observation timing remain part of the proposed experiment. Matching primitive
branches does not prove that a single native controller implements every menu.

The bounded fourteen-port source covers the declared finite embeddings. It is
not the unbounded row-1 theorem domain, the old 240-port comparator, or an
automatic physical identification of charge, mass and spatial spin. The dense
complex-state test is a finite implementation check; the partial-injection
argument supplies the general completeness statement within the bounded domain.
No source adoption, new platform run, outside-group acceptance or physical-row
promotion follows from this review.

`MUTABLE_INDEX_SCAN.md` records a separate read-only snapshot of stale current
tracking language observed during integration. Its findings were sent to the
integrator for correction; it is not an assertion that those mutable files
remain unchanged after this review was sealed.

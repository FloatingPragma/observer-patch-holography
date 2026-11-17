# Row23 local submission: native core-dependent effective coupling

**Requested review:** the literal finite criterion “coupling varies with
scale (onset at the core radius)” under the operational dictionary below.
Candidate `A-bound-core-response-v1`, class N+, unpinned. No value of nature
is used. Official status remains untested until the platform process accepts
and tests the proposed member. This packet is ready for local claim review,
not a report that OpenWave has awarded the row.

## What the task measures

The object is the lower native pair band at total momentum(pi,...,pi) on an
even periodic lattice. Its projector is determined from the complete native
transition, not a supplied wavefunction or field Hamiltonian. The observer
first records the pair's center and orientation, then retains every event
of the paired-source current response in that frozen intrinsic frame.

The effective coupling is the measured Fourier current-response amplitude
of that object divided by the independently measured point-source response
on the same graph, with a fixed total charge normalization of two. This is
an intrinsic current-probe effective charge. It is not a universal QED/QCD
beta function or SI electric charge. The spatial scale is the allowed probe
wavevector, and the core radius is separately measured from the same band's
position instrument. The center is a retained measurement outcome, not an
unobserved choice of origin.

## Source-to-output map

1. Public native binary preparation and every selected edge/coin/current
   reconstruct the full one- and two-load transfer matrices.
2. The translation projector and the source transfer identify the entire
   compact lower band. `PROJECTOR_PROOF.md` proves positivity, rank and exact
   spectral identification; `analyze.py` checks every matrix entry against
   the sealed native branch archive.
3. The complete diagonal current instrument on the normalized band functional
   equals the executed uniform-bond experiment outcome by outcome. This
   equality is sufficient for the stated finite measurement. It makes no
   assertion about a different off-diagonal instrument or coherent working
   preparation.
4. Actual paired currents and point-control currents yield
   `g_eff(k)/g_eff(0)=cos(kappa/2)` along `k=(kappa,...,kappa)`. The graph
   propagator and attempt rate cancel using the independently executed point
   response. Charge and clock normalization are fixed across scales.
5. The same object's `R²=1/4`, `R_parallel²=1/(4d)` and exact half-suppression
   scale give `|k|R_parallel=pi/3`. The radius is not fitted to the response.

## Required task outputs and fixed checks

The source packet executes dimensions1,2,3 and retains every branch. On the
L=12 member, the six measured normalized coefficients at
`kappa/pi=(1/6,1/3,1/2,2/3,5/6,1)` are

    ((sqrt(6)+sqrt(2))/4, sqrt(3)/2, sqrt(2)/2,
     1/2, (sqrt(6)-sqrt(2))/4, 0).

The point coefficient is exactly one; holding native events gives zero
response. The full free pair sector has higher source energy, while the
actual compact band's lower spectral projector is checked independently.
All zeros, comparison currents, source configurations and hold-control
records remain in the archive. The classical adjacent mixture equivalence
is retained as a property of this diagonal instrument, not used as proof
of the source spectrum.

The source run has1,114,640 literal native calls and complete records,
an independent all-record reader and byte-identical replay. This submission
adds a retrospective theorem and separately preregistered archive analysis;
it does not claim the original run preregistered the later projector proof.
Exact integer identities and rational radical bounds carry the conclusion;
floating tables only display it.

## Reproduction and custody

The source packet is `../20261114_codex_bound_core_response/`, manifest SHA256
`bf686b73a184fcfcfc16934cf30717c56b47a085435d4dc3cffece304b5a3c05`.
Its `PREREG.txt`, `PROOF.md`, `DESIGN.md`, `PIN.md`, `CELL.md`,
`CLAUSES.md`, `CONTROLS.md`, `METHOD_NOTE.md`, `headline.py`,
`verify_records.py`, `RAW_RECORDS.jsonl.gz`, receipt and replay remain frozen.
The current packet's preregistration binds that ancestry and its new files.

Run `python3 analyze.py` from this directory or use its absolute path. It
reads the frozen source archive, verifies all new projector identities and
emits `RECEIPT.json`. Re-executing the original source program is optional
for this new analysis; it is not needed to manufacture missing outcomes.

## Acceptance boundary

The local scientific claim is the finite operational effective-coupling
criterion above. A physical-radius continuum, universal beta coefficient,
physical unit conversion and arbitrary coherent preparation are outside
its scope, not added prerequisites for that finite claim. Review must
explicitly accept or reject the proposed effective-coupling dictionary;
the result should not be relabeled as universal charge running.

The proposed member still needs owner pin selection, an in-platform task
and review from the other group under OpenWave's process. A same-group
check does not pass G6. This packet makes those remaining decisions concrete
and reproducible rather than requesting approval of an unfinished design.

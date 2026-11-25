# A native written core, finite orbital preparation and load-dependent erasure

An empty native echo prepares the stationary written core used by the orbital
experiment. A sharp local load and a finite native evolution clock prepare its
lowest standing state with **96.05% acceptance**, including every refusal.
The public port records retain **98.14%** of admitted loads on the two core
endpoints through exposure 64 on both the 20- and 80-cell boards. The
zero-contrast control on the same admitted input gives **2.76% and 0.87%**
there at exposure 64. These are finite bound-load observations around an
addressed stationary memory. They supply no mobile nucleus or hydrogenic
shell series. [Executable](measure.py) · [Complete receipt](results.json.gz)

OpenWave's literal OW-31 criterion is “hydrogen-like discrete levels from
standing waves”. Its Interval Echo entry is untested at
[official head 0df73dc6](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L137).
The adopted candidate pin is 4fe10f9a. The local full criterion is partial.
The [two-load sector](../README.md) and [written-core census](../written_transitions/README.md)
exclude the prospective hydrogenic excited shells on their respective
sources. This assay concerns physical preparation and a source-removal
intervention on the written core.

## Formation and source selection

Begin with empty ports and blank collars. Apply one native 99-state echo at
the emitter seam, while its endpoints are known to hold 00. It preserves
the empty loads and produces the coherent zero-contrast memory

\[
|r_0\rangle=(|c=0,d=0\rangle+|c=1,d=0\rangle)/\sqrt2.
\]

The formation has zero native mean work on this vacuum. The seam address
and the decision to write it are preparation inputs. The memory is private
in the candidate specification. It is kept coherent; no intervening collar
measurement controls the experiment.

Supply one amount-one load sharply at the first endpoint. The native
one-load generator relative to the independently calibrated uniform blank
reference is H=L_cut/2−8P, with P the two endpoint projector. The blank
reference is 48 exchange units. This potential follows from the native
written contact, which annihilates all its vacuum and single-load columns.
The full 132-state one-cell memory-and-port audit reconstructs this reduction
directly from the native seam columns and degree references.

A supplied 4096-position clock gives the accepted operation

\[
A={1\over4096}\sum_{j=0}^{4095}
 e^{-ij(H-E_*)/4},\qquad E_*=-6.136.
\]

It uses forward native evolution, controlled timing and coherent clock
recombination. The accepted clock output and its orthogonal complement
are the two declared source flags. The refusal is kept in every denominator.
E_* is a calibration input from the independent written-core binding
measurement. No global spatial Householder compiler prepares this source.
The coherent clock and controlled free evolutions are supplied apparatus;
their autonomous implementation by native local services is unproved.

The admitted probabilities are 0.960514793 and 0.960511163. Their inferred
ground-state fidelities exceed 0.99999988. The source-selection ledger
retains the energy of both flag outcomes; their weighted sum equals the
sharp source's initial relative energy. Selecting a negative-energy line
does not represent unconditioned cooling.

## A native echo does not erase an occupied core

On the relevant local columns the exact native echo U obeys

\[
U|r_0,00\rangle=|\mathrm{blank},00\rangle,
\quad U|r_0,10\rangle=|r_0,10\rangle,
\quad U|r_0,01\rangle=|r_0,01\rangle.
\]

For a load wave ψ, one additional echo therefore gives

\[
U(|r_0\rangle\psi)=|\mathrm{blank}\rangle(1-P)\psi
                       +|r_0\rangle P\psi.
\]

The load is unchanged at the intervention. Only its exterior component
has an erased collar. The native block calculation gives an erased fraction
of **1.8565%** of admitted sources. Subsequent free evolution preserves
the blank and written blocks. Two immediate additional echoes, with no
read or free flight between them, restore the entire initial state because
U²=I. The terminal public port marginals of the zero- and two-echo arms
agree exactly in the calculation.

At exposure 64, one additional native echo leaves about 95.1% of admitted
loads on the core. A separate supplied memory-only swap maps r_0 to blank
at every load position. This comparator erases the full core and gives only
1.74% and 1.18% endpoint retention. It is a declared apparatus control, with
no native reachability claim for that memory-only swap. The native echo's
load-dependent action cannot be represented by a load-independent removal
of the orbital potential.

All work is retained. One native echo adds about 0.27457 exchange units
per original source start. The supplied full erase adds about 7.8842.
The unchanged and two-echo arms have zero net intervention work. The
changed-law arm retains the identical native-prepared source, switches
contrast to zero and records both the relative-energy quench and its
absolute work including the shift of the 6κ one-load reference. Release
means free evolution after these specified operations; it is no assertion
of work-free ionization.

## Public records and private diagnostics

Every public spatial bank contains all terminal load ports and source-tag
refusal. Its source-denominator is 2^22 original starts, seed 317031. Each bank is
drawn from integer weights on the 2^-40 grid, with source-tag refusal as
the exact complement, so last-bit roundoff of the computed probabilities
cannot change a seeded count. The L1 distance of every sampled law from
its computed probabilities is recorded; the largest is 6.9×10⁻¹⁰, inside
the a priori bound of 11·ports·2^-40 plus the 2×10⁻¹⁰ normalization budget.
Both boards, both flight contrasts, all four interventions and exposures
0, 1, 4, 16 and 64 are retained. The public count decoder reports the
endpoint fraction and RMS graph distance to the nearest endpoint. This
distance is zero on either endpoint; it differs from the midpoint radius
used in the core-scattering assay. At exposure 64 the unchanged source's
measured RMS distances are 0.13889 and 0.13929 links.

The 99% simultaneous Clopper–Pearson intervals cover all 80 endpoint
comparisons, conditional on each bank's admitted count. Every unchanged
and double-echo native interval has lower endpoint above 0.9810. Each
interval uses the shared forward-tail-checked count implementation, with
an independent replay of every endpoint in the tests. The finite propagation
uses Hermitian eigendecompositions and explicit norm/energy tolerances;
it is no interval-arithmetic enclosure or continuum-limit proof.

The receipt also retains an eleven-collar-label by port table and an
independent formation-label bank. These are **simulator diagnostics**.
They are not public collar measurements. The private labels never enter
the source control or the public spatial decoder. The erased fractions
above are native block calculations. A
[separate public material instrument](public_readback/README.md) reads the
terminal load port, routes that load to the first core endpoint without
crossing the core seam, and applies a native π contact. The load flips
endpoints for blank memory and stays put for the written zero-contrast
subspace. At exposure 64 it gives 74,603 flips in 4,028,408 admitted
20-cell sources and 74,590 in 4,028,686 admitted 80-cell sources after one
additional echo. The unchanged and two-echo sources give zero; the supplied
complete erase flips every admitted load. This material response infers the
collar property without a private-memory detector. Its routing and assay
occur after the free orbital flight.

On independent launches from the same finite accepted source, the supplied
transition clock and one-endpoint phase pulse give **331,295 and 332,148**
selected port records near frequency 0.47993. A common phase on both
endpoints gives **1,527 and 1,468**. At the binding-referenced hydrogenic
forecast 4.602, both banks give zero selected records; the exact probabilities
per original start are below 2.1×10⁻⁸. Transition-clock refusal and source-tag
refusal have separate bins. This spectroscopy retains the coherent
conjugated-evolution apparatus of the parent protocol, including inverse
evolution. The line settings use its calibrated spectrum. The finite source
clock and local source intervention do not derive the hydrogenic hierarchy.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -O -W error \
  -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.native_core_release.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -W error \
  -m pytest -q model/candidate_laws/interval_echo/openwave/row_31_orbitals/native_core_release
```

An explicit --output PATH creates a separate receipt exclusively. A .gz
suffix selects deterministic gzip compression. Serialization succeeds
before the output is opened. The replay compares every probability, count,
work account and source hash under ordinary and optimized Python. Seven tests
reconstruct the native echo, sum a finite clock directly, embed the full
memory-and-port source, audit the public records and count intervals,
check the changed-law and supplied-erasure discriminants, and redraw a bank
after one-bit changes to its tied probabilities.

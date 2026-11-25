# Transitions around a written material core

OpenWave records [OW-31 as untested](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L468); the full local criterion is partial. A load bound around one empty written seam has a directly readable transition at **0.479936** on 20 cells and **0.479928** on 80 cells. Independent binding measurements predict hydrogenic first-excitation frequencies near **4.6021**, where the native detector gives zero counts in 4,194,304 launches on each board. The same detector resolves a planted line there. This tests another natural atomic-source interpretation of Interval Echo, with no installed radial potential. [Protocol](PROTOCOL.md) · [Native source and complete records](results.json) · [Replay](measure.py)

## What the observer reads

The core is the existing stationary empty seam write used for the one-load response experiment. A supplied spatial Householder compiler maps a sharp endpoint load to the native ground state. Its ideal acceptance is one and refusal zero; a realization of this global compiler by native local gates is unestablished. The free native law holds the compiled state. A material degree-reference phase pulse changes the sign at one endpoint. A coherent transition clock selects the pulse's upward-frequency response, and a complete terminal port detector keeps every selected spatial outcome plus clock refusal.

Even frequency settings fit the line; odd settings test its shape. Held RMS errors are 0.0000402 and 0.0000286 in probability per original launch. The scan and its fit interval were chosen with the native spectrum known. It is a count-based local line measurement, without a prospective public prediction claim. The selected line gives **345,311 and 344,877 accepted outcomes** on the two boards. Applying a common phase to both endpoints instead gives **1,516 and 1,552**. The exact transition weights are approximately 0.08231 and 0.0003643. These different responses distinguish a local endpoint excitation from a common phase of the core.

A separate Ramsey comparison measures the ground binding relative to the blank uniform one-load reference. Both coherent branches contain one material load and the same empty written record; the unused written vacuum contributes zero energy. The inferred bindings are 6.136152 and 6.136133, with simultaneous 99% sampling half-width 0.000347. This reference fixes the hydrogenic forecast at three quarters of the binding. Its source clock, preparation and recombination are supplied metrology, with no laboratory energy or proper-time calibration.

## Why the absent target is decisive for this source

Removing the written edge gives a nonnegative graph Laplacian. Subtract eight at each of its two endpoints:

H = L_cut/2 − 8P.

On the endpoint span its matrix is diag(−6,−5.5). The variational principle gives at least two negative levels. On the orthogonal complement the matrix is nonnegative, so there can be at most two. Every negative level lies in [−8,−5.5]. This is an exact two-sided census on all three native boards, with every other level in [0,6]. The source has two localized standing states, without a hydrogenic series below the measured blank dissociation zero.

The observed native line joins approximately −6.13616 and −5.65623. A hydrogenic first excited state would sit near −1.534, inside the forbidden interval. A uniform change of energy origin leaves the binding-referenced transition comparison fixed. Finite Fourier resolution does not replace this spectral conclusion: including the binding allowance, its worst native target acceptance is below **1.642×10⁻⁶**, and its actual acceptance is about **2×10⁻⁸**. The clock period is 25.1327; the entire allowed detuning lies inside its alias-free half-period.

The detector sensitivity control moves the second level to one quarter of the ground energy by a supplied rank-one spectral term. It gives **345,160 and 344,949** target clicks. This installs the tested line and earns no atomic credit. Removing contrast while retaining the native initial vector and observer pulse gives **138 and 98** selected-line clicks. The unchanged vector is not stationary under that alternative law; every clock refusal remains in its denominator.

## Scope and replay

The core address is fixed by its written memory. Its nuclear interpretation, electric charge, autonomous state preparation, hydrogenic shells and physical atomic line scale are unresolved. The clock filter uses supplied coherent timing, inverse evolution and energy-selective recombination. It reads native frequencies; it does not create an independently established atom.

From the RER root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -O -B -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.written_transitions.measure
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_31_orbitals/written_transitions
~~~

The replay reads the deposited archive and raises on a mismatch, including in optimized Python. An explicit --output PATH creates a separate record exclusively. Six tests independently assemble the cut-edge generator, verify the native 99-state empty-write reduction and phase-reference pulse, sum the finite clock directly on a nonstationary input, decode complete counts and held predictions, and check negative probabilities and optimized mismatches. The 320-cell census uses the sparse graph and its two endpoint rows. No source spectrum or potential is fitted to a hydrogen line.

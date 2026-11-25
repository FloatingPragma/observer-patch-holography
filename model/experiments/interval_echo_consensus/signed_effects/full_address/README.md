# A complete material address register for one cell

One vacancy or excess at any of the twelve ports of a data cell can write its address into four fixed neutral pointer dimers. The same native pulse schedule works for either sign and for coherent inputs with an untouched charge reference. Every routed data port is restored. The final instrument resolves all twelve addresses, preserves the charge coherence within each address, and retains unused and invalid outcomes. Native controls, prepared registers and the final sharp read are supplied. [Schedule and complete records](results.json) · [Implementation](measure.py) · [Independent flat-word replay](test_deposited_word.py)

## Apparatus and exact geometry

Use the registered 20-cell, 240-port federation. The data cell is owner 0, ports 0 through 11. Its admitted state has exactly one vacancy or excess in unit-filled material. The sign is unknown to the controller. Four pointer/host services are:

| Bit | Owner | Service source | Neutral pointer (park, cross) | Host |
| --- | ---: | ---: | --- | ---: |
| 0 | 1 | 13 | (16,17) | 19 |
| 1 | 2 | 25 | (28,29) | 31 |
| 2 | 4 | 49 | (52,53) | 55 |
| 3 | 5 | 61 | (64,65) | 67 |

Each pointer starts at loads (2,0); every host starts at one. The source and cross ports of all services have degree five on this actual board, providing their native degree-reference terms. The service contacts join source to cross, park to cross, and source to host. Every contact is checked against the stored native incidence. The board file and seam list are hashed in the receipt.

For coherent charge inputs, port 239 carries the complementary load: two with a vacancy and zero with an excess. It is an untouched input reference throughout acquisition. A classical unknown sign needs no instruction from that reference. All other material starts at one and every seam memory is blank.

The protocol supplies a serial controller, calibrated native-duration waits and the completion order of local services. An implementation by local requests and acknowledgments is assumed and is not simulated. Exposure calibration, isolation and the clock implementation are apparatus. No shared physical proper time is inferred. The four final readings are local records; their ordered tuple is compared after acquisition, without an instantaneous common observer.

## Local native gate

The [fixed neutral pointer](../FIXED_POINTER.md) gives a 27-wait word on a qutrit control, a neutral pointer and a load-one host. On every source load \(n\in\{0,1,2\}\) and either pointer orientation it is

\[
G:\ |n\rangle|b\rangle|1\rangle
\longmapsto h\,|n\rangle|b\mathbin\oplus\mathbf1_{n\ne1}\rangle|1\rangle,
\qquad h=e^{-133\pi i/260}.
\]

The word uses only native contacts and degree references. Both odd-difference and extreme-load (0,2) exchange are refocused. The full local space has dimension \(3^4=81\); every output row of all six source/pointer code columns is checked, including their common scalar and absence of leakage. Exact symbolic tests derive the diagonal interaction and common phase. The pulse sequence is unchanged by the incoming sign, source position or earlier records. [Local identities and operator checks](test_address.py)

A native contact for \(\pi\) is exactly a tensor-factor swap on arbitrary ternary loads with blank seam memory. Equal loads are fixed. The unequal blocks have eigenvalues 8/9 or 32/33, whose \(\pi\) phases give the swap. The full 99-state native operator confirms both its blank-memory restoration and the swap action; the statement includes opposite loads.

## Routing is a conjugation

For each address whose binary label has bit \(j\) set, transport its load to service \(j\), execute \(G\), and reverse the transport. The twenty calls cover bit sets of sizes 6, 6, 4 and 4. The receipt contains every route and the complete flattened 736-wait pulse word.

Every forward path avoids all pointer and host ports and reference 239. It may pass through another data port or an idle service source. Those registers are restored by the reverse path. All route contacts are native seams. Forward path lengths sum to 98, with maximum seven; their reverses use another 98 swaps.

Let \(T\) be a route permutation taking data port \(a\) to service source \(s\). It fixes the pointer and host factors. Since \(G\) retains its control load and returns its host, the physical operation \(T^*G_{s,P,H}T\) is exactly the defect-presence copy controlled by data port \(a\). This is an operator identity on the code, including arbitrary entanglement with other data, route registers and untouched references. It is not inferred from equal classical marginals. The tensor-slot permutation is checked explicitly for every route.

Composition therefore gives

\[
U\,|\sigma,a\rangle|r\rangle
=h^{20}|\sigma,a\rangle|r\mathbin\oplus a\rangle,
\qquad \sigma\in\{-,+\},\quad 0\le a<12.
\]

All 384 sign/address/prior-register words are executed as physical 240-port arrays. Each returns every nonregister port, including the balancing reference, to its input value. The exact local operator identity and linearity extend this factorization to all coherent inputs. This avoids a \(3^{240}\) state vector without discarding branches. The execution is local matrix checking and reduced operator/permutation composition; it counts no literal native-event calls.

An independent implementation parses the archived flat pulse sequence directly, without the router or logical copy primitive. It reconstructs the ternary two-port gates, evaluates each four-port service on all 81 states and propagates all 240 physical port values and phases on the 384 coded inputs. It also checks compatibility with the four native collar deposits. The separate implementations agree with one common phase across addresses and signs. The small floating local and swap residuals are regression diagnostics, not physical error certificates or full-board numerical evolution.

## Complete instrument, repeat and controls

Read the eight pointer-port occupations and the four hosts sharply. A well-formed pointer is (2,0) or (0,2), and every host is one. Decode the four bits as an integer 0 through 15. Any malformed pointer or host gives the complete complement outcome 16. This aggregates the other possible fine occupation words rather than dropping them. On the admitted one-defect, zero-register input, words 12 through 15 and the invalid outcome have zero ideal probability. Outside that preparation, a four-bit parity code is not a count or a complete description of several defects; in particular an empty data cell aliases word zero.

The branch map for address \(a\) is

\[
\Phi_a(\rho)=P_a\rho P_a,\qquad
P_a=|-,a\rangle\langle-,a|+|+,a\rangle\langle+,a|.
\]

The twelve projectors sum to identity on the 24-dimensional signed address domain. Each has rank two and retains the relative charge phase with its balancing reference. The seventeen outcomes form a complete instrument, including the zero unused/complement branches. Reading a classical address does not copy an arbitrary quantum state.

The inverse native program resets the pointer register while leaving the corresponding conditional source state unchanged. Repeating with a reset register gives the same address, because \(P_bP_a=\delta_{ab}P_a\). Without an intervening read, inverse playback restores any source/reference superposition. With the pointer read and ignored first, the intended address dephasing remains. A native half contact on data ports 0 and 1 and its inverse distinguish these cases: return probability is one for coherent playback and one half after the pointer read.

Twenty-four deterministic synthetic banks, 4,096 starts each, cover every address and both signs with all seventeen outcome bins retained. These total 98,304 starts. Two actual schedule controls give an observable wrong answer with the same reader: omitting bit 3, or routing its address-8 query from data port 9, makes either sign at address 8 decode as zero. The data and other ports are restored in both controls. They test the measurement program; they establish no changed-law physical-emergence claim.

The complete pulse word has 396 contacts and 340 references. For preparation half-trace error \(\epsilon_p\) and absolute duration error at most \(\delta\) per wait, telescoping and complete-channel contractivity give output half-trace error at most

\[
\min(1,\epsilon_p+23948\delta),\qquad 23948=396\cdot33+340\cdot32.
\]

The bound includes arbitrary untouched references and every outcome, without conditioning on acceptance. Physical fidelity, clock calibration, arbitrary adaptive event grammars and external custody are not supplied.

## Classical record storage is a separate extension

The four output dimers coincide with the [native collar-storage stations](../COLLAR_STORAGE.md). A protected native basis word and echo deposit each orientation in a signed-contrast parity, with one common phase and endpoints (1,1). Its classical bit projectors commute with every term of the free federation Hamiltonian, even when neighbouring contacts change the endpoint loads. Thus the four-bit record of the measurement can persist during subsequent material motion. It records the measured address, rather than the later position of the source.

Retrieval supplies unit endpoint loads and protected inverse operations, followed by a sharp material read. A bare write and a half-exchange law give readout errors on the disturbed endpoint controls. Coherent free storage is distinct: the two memory lines can acquire different material-dependent phases and entangle with the environment. The conserved classical address and software hashes supply no trusted external custody or final observer read.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.full_address.measure
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/full_address
```

Replay preserves the receipt. Explicit `--write` exclusively creates a missing first record. Source hashes bind the software version; the bank is a construction replay, without a prospective experimental or OpenWave award claim. The full-cell address clause is realized on its declared signed one-defect domain. PR-03, PR-64 and PR-65 retain their broader physical, event-grammar and custody obligations.

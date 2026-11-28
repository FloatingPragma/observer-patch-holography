# Two independent orbital inputs on a written fan path

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The sender prepares an arbitrary superposition of the two archived outer-ring parity inputs with the charge at the source gateway. Every other fan is in the uniform vacuum of its admitted head words. The receiving observer accepts charge and the local orbital code at the last gateway. Separate observers at all other fans also read their uniform vacua, using the head-word sectors appropriate to charge arrival. Accept a trial only when all these records coincide. This gives a single linear amplitude from two input components to two output components.

Use both charge signs, the archived 14→24 bridge on the 20-cell board, and the first qualifying asymmetric three- and four-fan prefixes on the 320-cell board. The asymmetric rule and retained ports are fixed by the [terminal protocol](../fan_path/TERMINAL_CONTROL.md). The symmetric 20-cell routes of the same lengths serve as reflection controls. The source, receiver, projector and all material, ring and degree terms are those of the [finite path construction](../fan_path/fixture.py).

Check whether the determinant of the accepted amplitude is identically zero. Scan terms through order 72 modulo a prime. Certify the integer coefficients through the first nonzero residue by combining primes until the product exceeds twice the common absolute coefficient bound. If no residue appears, limit the inference to the tested order. On the reflection controls test the exact charge-word involution, the unchanged sender and the transformed detector. Finite-wait probabilities and a receiver that discards spectator records are separate measurements of the [path protocol](../fan_path/PROTOCOL.md).

## Exact transfer on both charge signs

| Native board | Accepted route | Leading determinant term | Certified $0<|t|<\delta$ |
| --- | --- | --- | --- |
| 20 cells | 14→24 | $-(1227664/893025)t^{16}$ | $\delta\ge 1.08\times10^{-18}$ |
| 320 cells | 98→241, 240→1011 | $-(39285248/2210236875)t^{20}$ | $\delta\ge 2.40\times10^{-23}$ |
| 320 cells | 98→241, 240→1011, 1010→1045 | $-(314281984/2844574858125)t^{24}$ | $\delta\ge 3.78\times10^{-28}$ |

The coefficients are exact in the archived source and receiver basis. They give strictly positive success probability for every normalized input within the certified punctured intervals. The small radii are lower bounds from a uniform analytic estimate. The certificate assigns no detection rate at the observation appointments. The result concerns joint acceptance and does not by itself establish entanglement when the spectator records are ignored, or energy protection of the same excitation. [Exact receipts](transfer_series.json) · [Finite-wait protocol](../fan_path/PROTOCOL.md)

For the symmetric 20-cell three- and four-fan routes, terminal reflection at gateway 116 fixes the sender and commutes with the complete generator. It interchanges the two code vectors with a minus sign. The exact permutation leaves each spectator word fixed on the accepted receiver support. Each word-resolved receiving amplitude therefore lies in the same one-dimensional fixed subspace. The joint accepted determinant is identically zero and tracing the spectator records cannot leave the receiver entangled with the sender's reference on these controls.

## Certificate and how to read it

Let $A=2H$ in the configuration basis. Its entries are integers. Divide the integer launch and receive columns $S,R$ by $\sqrt{g_S},\sqrt{g_R}$, where the squared factors count the joint word support. The program checks these columns against the fixture's actual preparation and decoder, including each spectator vacuum. Put $M_j=R^\mathsf{T}A^jS$. The coefficient of $(-it)^n$ in the determinant is $D_n/(2^n n!g_Sg_R)$, where
$$
D_n=\sum_{j=0}^n {n\choose j}
\left[
 (M_j)_{00}(M_{n-j})_{11}-(M_j)_{01}(M_{n-j})_{10}
\right].
$$

For the symmetric integer $A$, its maximum absolute row sum $B$ bounds its operator norm. This gives $|D_n|\le 2g_Sg_R(2B)^n$. The product of the recorded primes exceeds twice this bound at the reported order, so Chinese remaindering recovers the unique signed integer and certifies all earlier zeros. With leading coefficient $a_k$, the remaining terms for $|t|\le 1/B$ have magnitude at most $6(B|t|)^{k+1}/(k+1)!$. Taking $\delta$ no greater than $|a_k|(k+1)!/(12B^{k+1})$ keeps their sum below half the leading magnitude. The receipt gives $\delta$ as an exact fraction.

An observer can test the conclusion by preparing half of a maximally entangled pair in the sender code. The jointly accepted reference–receiver record has a negative partial-transpose eigenvalue $-|\det T(t)|/2$ before dividing by its success probability. This follows by the singular-value decomposition of the two-input amplitude $T(t)$. The full-rank result establishes entanglement in that accepted record. If the same experiment traces over spectators, its accepted output is a sum over distinct outcomes and requires a separate partial-transpose test.

Replay from the repository root.

~~~sh
python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_recoil/transfer_series.py
~~~

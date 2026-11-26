# Which windings can a receiver identify on the released guide?

Consider the three charged winding preparations labelled 0, + and − followed by the registered conversion to a neutral vacancy–excess pair. A receiver may read any part of the 14-port guide once during free flight, with arbitrary local analysis. The two pendant leaves and their written heads stay at the source. For this domain, **the two conjugate source labels give exactly the same distribution of every receiver outcome**, including an empty region or the arrival of only one constituent. This statement holds at every free-flight time. It gives an upper bound on how well the receiver could perform even with an ideal detector. [Finite proofs](RemoteFlavour.lean) · [Release and source audit](../audit.py)

## From two source symmetries to all receiver outcomes

1. For any sharp three-bit tip word, at least two of the bits agree. The transposition of those arms fixes the word. Each transposition interchanges the two source winding recipes up to an overall phase and preserves the charged-flight generator. The tip-word effect sums over the unread hub records. It therefore has the same probability for + and −. Lean proves the fixing transposition and the equality for every reflection-covariant pair of word distributions. Summing words that give the same two active launch bits preserves this equality. [Charged source symmetry check](../audit.py)
2. Let the two types of source-side pendant records be the leaf loads and the written bits. For each pendant, XOR its written bit with the indicator that the leaf load differs from one. The result equals the original written bit at launch. The integer contact audit checks preservation on every edge of the finite sector. Since these sectors are invariant under the generator, each power of the generator, and hence its finite-time evolution, preserves them too. Every residual record fixes both launch bits even on trials where a constituent stays on a leaf. [Contact and sector check](../audit.py)
3. Write $V_{exs}$ for the propagated amplitude at residual record $e$ and guide word $x$ from launch sector $s$. Let $g(e)$ be the sector read from that residual record. The contact result gives $V_{exs}=0$ when $s\ne g(e)$. For any input density $ρ$, tracing the residual record gives

   $$
   (ρ_{\rm guide})_{xy}
   =\sum_{e,s,t} V_{exs}ρ_{st}\overline{V_{eyt}}
   =\sum_e V_{ex,g(e)}ρ_{g(e),g(e)}\overline{V_{ey,g(e)}}.
   $$

   Thus every matrix element on the guide depends only on the four original sector weights. Lean proves this identity for any finite launch, residual and guide alphabets. Contracting it with an arbitrary receiver effect gives identical outcome probabilities when the input weights agree. In particular a receiver can retain interference between different arriving words within each sector. A joint measurement that coherently reads the old source heads as well as the guide accesses other matrix elements. [Remote detector and source-dependent measurement](../ERASER_PROTOCOL.md)

## The optimal chance of guessing the original label

Give the three windings equal prior probability. For any fixed complete measurement let $p(o)$ be the distribution for label 0, and $q(o)$ the common distribution for + and −. On outcome $o$, write $a(o)$ for the probability of guessing 0. The other two guesses contribute $q(o)(1-a(o))/3$ to total success, regardless of how the observer divides them. The optimal decision takes $a(o)=1$ when $p(o)\ge q(o)$ and zero otherwise. Lean proves its exact success and the upper bound

$$
P_{\rm best}=\frac{1}{3}\sum_o\max\{p(o),q(o)\}
=\frac{1+D_{\rm TV}(p,q)}{3}\le\frac{2}{3},
\qquad D_{\rm TV}(p,q)=\tfrac12\sum_o|p(o)-q(o)|.
$$

For the enrolled owner-17 aperture at time 2048, the archived model gives a complete-trial trace distance of $9.1638120\times10^{-8}$ between the zero and positive windings. This record allows arbitrary analysis of the six received pair words and includes separate outcomes for source refusal, a failed conversion arrival and a missed remote pair. The positive and negative eigenspaces of the difference of the two density matrices attain their trace distance as total variation. Even that ideal read would guess one label correctly with probability **0.333333364** per original launch. This is a calculation from the simulated outcome law. A larger receiver that keeps partially arriving loads can change this number. The two-conjugate equality and the $2/3$ ceiling continue to hold for the whole guide. Weak-source and charged-partner tags defining physical neutrino flavours remain unresolved on this construction. [Archived distributions and ideal receiver distances](../records.json)

From the repository root, run the six kernel checks using the Lean project environment.

~~~sh
cd Lean
lake env lean -j1 -DwarningAsError=true ../model/candidate_laws/interval_echo/openwave/row_10_oscillations/theory/RemoteFlavour.lean
~~~

The fixing-reflection lemma is axiom-free. The remaining five theorems use only Lean's standard logical axioms. The physical support of the registered generator is verified separately by the executable audit.

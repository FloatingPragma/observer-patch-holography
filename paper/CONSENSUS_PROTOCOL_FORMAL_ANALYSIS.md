# Formal Analysis of Reality as a Consensus Protocol: Byzantine Agreement, Repair Maps, and Quantum Error Correction

**Contribution to Observer-Patch Holography (OPH) — Paper 4 Extension**

> **Status notice:** The technical content of this document has been reformatted as a LaTeX appendix in `paper/appendix_B_bft_qecc_extensions.tex`, intended to be merged into `paper/reality_as_consensus_protocol.tex`. This Markdown file is retained for readability and issue-tracker cross-reference. All four issues raised in the PR #145 review have been addressed here and in the `.tex` file.

> **Honesty label key:**
> **[Established]** — follows from cited prior work or a complete argument given here.
> **[Conditional]** — true under explicitly stated additional assumptions not yet derived from OPH first principles.
> **[Conjecture / Proposed]** — a plausible open direction, not a settled result.

---

## Merge structure (for the maintainer)

The content below should be integrated into `paper/reality_as_consensus_protocol.tex` as follows:

1. **Short subsection** — paste the brief prose from Section 1 as a new subsection titled *Conditional BFT and QECC Extensions* under the existing *Discussion and Open Problems* section.
2. **Appendix B** — add `\input{appendix_B_bft_qecc_extensions}` immediately after the current quantum/algebraic-lift appendix (Appendix A), with section title *Conditional Distributed-Systems and QECC Extensions of the Consensus Formalism*.

---

## 1. Summary for Discussion and Open Problems (short subsection)

The consensus formalism of OPH has natural analogies to distributed Byzantine agreement. Observer patches correspond to protocol nodes, overlap repair corresponds to a quorum vote, and the repair fixed-point corresponds to consensus. Under explicit structural assumptions (quorum size ≥ 2f+1, partial synchrony, one-vote-per-view, certificate semantics, and DLS-style view-change), a QBFT-style interpretation of OPH repair satisfies safety and liveness — see Appendix B. The repair map admits a quantum-channel formulation via the Petz recovery map, but the CPTP property on all inputs requires either full-rank $\mathcal{N}(\sigma)$ or domain restriction; trace-preserving completion is not automatic when $\mathcal{N}(\sigma)$ has a non-trivial kernel. A QECC interpretation is possible under additional topological structure, but the code-distance/min-cut equality requires explicit conditions on logical-operator homology and boundary geometry. All of these extensions are conditional or conjectural and are not part of the core theorem package of Paper 4.

---

## 2. Theorem 1 — QBFT Safety Bound

### 2.0 What counts as a QBFT-style protocol (Definition 2.0)

**Definition 2.0 (QBFT-style protocol).** A consensus protocol is called *QBFT-style* in this analysis if it satisfies the following three structural properties. The theorem does not hold for protocols that lack any of them without a compensating change to the proof.

- **(P1) One-vote-per-view.** Each honest node casts at most one vote per view number. A node that has already voted in view $v$ ignores any later request to vote in view $v$.
- **(P2) Certificate semantics.** A decision requires a valid quorum certificate: $2f+1$ distinct, unforgeable, authenticated votes for the same value in the same view.
- **(P3) DLS-style view-change.** If no certificate is produced within a timeout, every honest node increments the view number by one and a new leader is selected by a fixed deterministic rule. At GST, timeouts fire correctly and the view-change terminates in bounded rounds.

The Istanbul BFT / QBFT protocol family (Moniz 2020; Saltini et al.) satisfies (P1)–(P3) and is the intended instance.

### 2.1 Setting and Explicit Assumptions

Let $\mathcal{O} = \{O_1, \ldots, O_n\}$ be a finite set of observers. Each observer $O_i$ holds a local patch state $\rho_i \in \mathcal{D}(\mathcal{H}_i)$. Agreement rounds are modelled as authenticated message-passing over a communication graph $G = (\mathcal{O}, E)$.

**[A1] Partial synchrony (DLS model).** Fixed but initially unknown bounds $\Delta$ (message delay) and $\Phi$ (processing rates). *Safety* holds unconditionally; *liveness* holds after the Global Stabilisation Time (GST).

**[A2] Byzantine fault model.** At most $f$ observers behave arbitrarily; the remaining $n - f$ are honest.

**[A3] Optimal fault bound.** $n \geq 3f + 1$ (necessary and sufficient).

**[A4] Strong quorum connectivity.** Every quorum $Q$ with $|Q| = 2f+1$ is *strongly connected* within $G$: for any $u, v \in Q$ there is a directed path in $G$ passing only through $Q$. This is strictly stronger than weak connectivity of the quorum-overlap graph and is required to propagate signed votes within a quorum.

**[A5] Message authentication.** All messages carry unforgeable digital signatures.

**[A6] OPH quorum overlap.** Any two quorums of size $2f+1$ satisfy $|Q_a \cap Q_b| \geq f+1$ (guaranteed by A3).

### 2.2 Theorem Statement **[Established, conditional on A1–A6 and Definition 2.0]**

**Theorem 1 (QBFT Safety Bound).** *Under assumptions A1–A6, any consensus protocol satisfying (P1)–(P3) of Definition 2.0 and run over the OPH observer graph satisfies:*

- *(Safety) No two honest observers finalise conflicting patch states.*
- *(Liveness) After GST, every honest observer finalises a state within $O(f \cdot \Delta)$ wall-clock time.*
- *(Optimality) The bound $f < n/3$ is tight.*

### 2.3 Proof Sketch

**Safety.** Suppose $O_a$ and $O_b$ finalise $s_a \neq s_b$ in the same view. By (P2), each required a certificate of $q = 2f+1$ votes (sets $Q_a, Q_b$). By (A3): $|Q_a \cap Q_b| \geq f+1$. By (A2), the intersection contains an honest $O^*$. By (A4), $O^*$'s signed vote is path-reachable within both quorums. By (P1), $O^*$ voted for at most one value — contradiction.

**Liveness and Optimality** follow from Dwork, Lynch, Stockmeyer (1988, Thm. 4.4) and Lamport, Shostak, Pease (1982), cited directly.

**Note on FLP.** Fischer, Lynch, Paterson (1985) is an *impossibility* result for fully asynchronous systems; it does not bear on achievability under partial synchrony (A1).

---

## 3. Theorem 2 — Convergence of the OPH Repair Map

### 3.1 Definition

**Definition 3.1 (OPH Repair Map — Petz form).** Let $\sigma \in \mathcal{D}(\mathcal{H})$ be a full-rank reference state and $\mathcal{N} : \mathcal{B}(\mathcal{H}) \to \mathcal{B}(\mathcal{K})$ a quantum channel. The *OPH repair map* is:

$$\mathcal{R}_{\sigma,\mathcal{N}}(\rho) := \sigma^{1/2}\, \mathcal{N}^\dagger\!\left(\mathcal{N}(\sigma)^{-1/2}\, \rho \,\mathcal{N}(\sigma)^{-1/2}\right) \sigma^{1/2},$$

where $\mathcal{N}^\dagger$ is the adjoint channel and inverses are taken on $\mathrm{supp}(\mathcal{N}(\sigma))$.

**Remark 3.2.** The Petz map and the trace-distance projection $\mathcal{P}_{\mathcal{S}}(\rho) := \arg\min_{\tau \in \mathcal{S}} \frac{1}{2}\|\rho - \tau\|_1$ are different objects. All subsequent properties refer to Definition 3.1.

### 3.2 CPTP Property **[Established, subject to domain restriction]**

**Proposition 3.3 (Petz map CPTP — domain-restricted statement).** *Let $\sigma$ have full support on $\mathcal{H}$. Then:*

*(a) $\mathcal{R}_{\sigma,\mathcal{N}}$ is completely positive.*

*(b) $\mathcal{R}_{\sigma,\mathcal{N}}$ is trace-preserving on $\mathrm{supp}(\mathcal{N}(\sigma))$, i.e., on inputs $\rho$ for which $\mathcal{N}(\sigma)^{-1/2}\rho\,\mathcal{N}(\sigma)^{-1/2}$ is well-defined.*

*(c) If additionally $\mathcal{N}(\sigma)$ has full rank on $\mathcal{K}$, then $\mathcal{R}_{\sigma,\mathcal{N}}$ is CPTP on all of $\mathcal{B}(\mathcal{K})$.*

**Domain restriction note.** If $\mathcal{N}(\sigma)$ is *not* full rank on $\mathcal{K}$, then either (i) the domain of $\mathcal{R}_{\sigma,\mathcal{N}}$ must be restricted to $\mathrm{supp}(\mathcal{N}(\sigma))$, or (ii) pseudoinverses must replace the inverses (generalised Petz map; cf. Junge et al. 2018), or (iii) a regularisation $\mathcal{N}(\sigma) \mapsto \mathcal{N}(\sigma) + \varepsilon\,\mathbf{1}$ must be introduced. Full-rank $\sigma$ alone does not prevent $\mathcal{N}(\sigma)$ from being rank-deficient: the channel $\mathcal{N}$ may map the support of $\sigma$ into a strict subspace of $\mathcal{K}$. In the OPH setting, whether $\mathcal{N}(\sigma)$ is full rank depends on the specific channel describing the overlap measurement and must be verified separately (open issue #62).

*Proof.* Complete positivity: composition of three CP operations, (i) sandwiching by $\mathcal{N}(\sigma)^{-1/2}(\cdot)\mathcal{N}(\sigma)^{-1/2}$ on $\mathrm{supp}(\mathcal{N}(\sigma))$, (ii) the adjoint channel $\mathcal{N}^\dagger$, (iii) sandwiching by $\sigma^{1/2}(\cdot)\sigma^{1/2}$. Trace preservation in the full-rank case: Petz (1986) and Fagnola–Umanità (2010). $\square$

### 3.3 Contraction **[Conditional]** and Spectral Gap **[Conjecture]**

**Conditional Proposition 3.4.** *If $\mathcal{N}$ is strictly contractive with gap $\lambda < 1$, then $\mathcal{R}_{\sigma,\mathcal{N}} \circ \mathcal{N}$ is contractive. Establishing $\lambda < 1$ for the OPH overlap channel requires analysis of the OPH Hamiltonian (open issue #62).*

**Conjecture 3.5.** *The transfer operator $\mathcal{T}$ has spectral gap $\delta > 0$, implying exponential convergence (open issue #63).*

### 3.4 Theorem 2 **[Conditional on Conjecture 3.5]**

**Theorem 2.** *Assuming Conjecture 3.5 with gap $\delta > 0$:*
$$\tfrac{1}{2}\left\|\mathcal{R}_{\sigma,\mathcal{N}}^{\circ t}(\rho) - \sigma\right\|_1 \leq C\, e^{-\delta t}.$$

---

## 4. Theorem 3 — QECC Correspondence

### 4.1 Notation

$N = \dim(\mathcal{H}) = 2^n$ for $n$ physical qubits. Standard notation: $[[n,k,d]]$ stabilizer code; $K = 2^k$; quantum Singleton bound: $k \leq n - 2(d-1)$.

### 4.2 Code Distance and Graph Min-Cut **[Conditional — tightened]**

The identity *code distance = graph min-cut* is specific to topological codes (surface/toric codes) where logical operators correspond to non-contractible homological cycles. It does not hold for a generic overlap graph.

**Conditional Claim 4.1 (tightened).** *Suppose the OPH observer network is equipped with a surface-code-type construction on a planar or toroidal graph $G_\text{OPH}$, and the following conditions are satisfied:*

- *(i) Logical $X$-type operators correspond to minimum-weight non-contractible cycles in the primal chain complex of $G_\text{OPH}$; logical $Z$-type operators correspond to minimum-weight non-contractible cycles in the dual complex.*
- *(ii) Boundary conditions are chosen such that no logical operator of weight strictly less than the min-cut of $G_\text{OPH}$ exists.*
- *(iii) For non-planar geometries, the relevant group is $H_1(G_\text{OPH}; \mathbb{F}_2)$ and the code distance equals the minimum over all non-trivial homology classes of the weight of a representative cycle.*

*Under (i)–(iii), the code distance $d$ of the OPH encoding equals the min-cut of $G_\text{OPH}$.*

*This claim is conditional on an explicit construction of the topological encoding map, the primal/dual complex of $G_\text{OPH}$, and verification of (i)–(iii), none of which have been provided (open issue #113).*

### 4.3 Communication Complexity **[Conjecture]**

**Conjecture 4.2.** *Per-round complexity is $O(n \cdot \mathrm{poly}(d))$ (open issue #72).*

### 4.4 Theorem 3 **[Conditional / Conjecture]**

**Theorem 3 (QECC Correspondence).** *Suppose conditions (i)–(iii) of Claim 4.1 hold. Then: (i) [Conditional] code distance = min-cut; (ii) [Established] Knill–Laflamme conditions hold for $t < d/2$; (iii) [Conjecture] per-round complexity is $O(n \cdot \mathrm{poly}(d))$.*

---

## 5. Theorem 4 — Asynchronous Convergence

### 5.1 Why fairness alone does not give probability-1 convergence

Standard strong fairness guarantees that every enabled action fires infinitely often along any fair schedule; it does not impose a probability space on the set of schedules. A statement of the form "converges with probability 1" requires a measure on schedules and does not follow from fairness alone. The FLP impossibility result (Fischer, Lynch, Paterson 1985) confirms that even strong fairness is insufficient for bounded-time consensus in a fully asynchronous system. The phrase "with probability 1" has been removed from Theorem 4a; the convergence is per-schedule and topological.

### 5.2 Additional assumptions for a quantitative bound

**[B1]** Finite known bound $\Delta$ on message delay after GST. **[B2]** Finite bound $\Phi$ on processing rates. **[B3]** $f < n/3$.

### 5.3 Two Statements

**Theorem 4a (Eventual Convergence) [Established under fairness only].** *In a fully asynchronous OPH observer network under standard strong fairness, iterated application of $\mathcal{R}_{\sigma,\mathcal{N}}$ converges to a consensus state $\sigma^*$ in the following sense: for every $\varepsilon > 0$ and every strongly fair schedule, there exists a step $T(\text{schedule}, \varepsilon) < \infty$ such that $\frac{1}{2}\|\mathcal{R}^{\circ t}(\rho) - \sigma^*\|_1 < \varepsilon$ for all $t \geq T$. This is a per-schedule topological statement. No probability measure on schedules is assumed or needed; no uniform finite bound on $T$ follows from fairness alone.*

**Theorem 4b (Quantitative Convergence) [Conditional on B1–B3].** *In a partially synchronous network satisfying B1–B3, after GST every honest observer reaches consensus within $T = O(f \cdot \Delta)$ wall-clock time (DLS framework, Theorem 4.4).*

---

## 6. Open Problems

- **#62** — Derive repair map from OPH dynamics; verify full-rank condition for $\mathcal{N}(\sigma)$.
- **#63** — Prove spectral gap from OPH Lyapunov functional (upgrades Conjecture 3.5).
- **#68** — Quantum observable-level confluence.
- **#69** — Continuum/refinement limit of consensus theorems.
- **#72** — Communication complexity (upgrades Conjecture 4.2).
- **#73** — Re-export repair map into OPH language.
- **#113** — Construct topological encoding map for Claim 4.1.

---

## 7. References

1. Fischer, M.J., Lynch, N.A., and Paterson, M.S. (1985). *JACM*, 32(2):374–382.
2. Lamport, L., Shostak, R., and Pease, M. (1982). *ACM TOPLAS*, 4(3):382–401.
3. Dwork, C., Lynch, N.A., and Stockmeyer, L. (1988). *JACM*, 35(2):288–323.
4. Castro, M. and Liskov, B. (1999). *OSDI 1999*, pp. 173–186.
5. Moniz, H. (2020). arXiv:2002.03613.
6. Saltini, R. et al. Consensys/qbft-formal-spec-and-verification, GitHub.
7. Petz, D. (1986). *CMP*, 105(1):123–131.
8. Fagnola, F. and Umanità, V. (2010). *IDAQP*, 13(3):459–486.
9. Junge, M. et al. (2018). *AHP*, 19(8):2505–2555.
10. Knill, E. and Laflamme, R. (1997). *PRA*, 55(2):900–911.
11. Gottesman, D. (1997). PhD thesis, Caltech.
12. Kitaev, A. (2003). *Annals of Physics*, 303(1):2–30.
13. Dennis, E., Kitaev, A., Landahl, A., and Preskill, J. (2002). *JMP*, 43(9):4452–4505.
14. Buhrman, H., Cleve, R., and Wigderson, A. (1998). *STOC 1998*, pp. 63–68.
15. Attiya, H. and Welch, J. (2004). *Distributed Computing*, Wiley.
16. Chandra, T.D. and Toueg, S. (1996). *JACM*, 43(2):225–267.
17. Aminof, B., Kupferman, O., and Rubin, S. (2014). *CONCUR 2014*.

---

*This document is an exploratory extension to OPH Paper 4, intended for merge into `reality_as_consensus_protocol.tex`. The author thanks Bernhard Mueller and the FloatingPragma team for their detailed review.*
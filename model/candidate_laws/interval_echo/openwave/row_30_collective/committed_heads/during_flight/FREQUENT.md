# Effective evolution under rapid local coin readings

Use the same strand, keep every coin word, and pool the terminal wall hits over the words. After a first coin measurement consider any prescribed sequence of additional ideal coin readings at the same head. The source and the detector stay fixed. The analytic calculation below is a follow-up to the two-read acquisition and supplies an arbitrary-cadence prediction.

## A sequence of intervals

Write $U_+(t)$ and $U_-(t)$ for the matter propagators in the conserved even- and odd-coin **parity** sectors of this sign-even head. The sharp binary label that the home observer reads is the complementary basis $(|+\rangle\pm|-\rangle)/\sqrt2$. Prior to the first read the head has parity $+$. After that read, summing over its labels leaves the two parity-diagonal entries equal to $\rho/2$ for the current matter state $\rho$. Free flight for a duration $s$ acts on these entries by $U_+$ and $U_-$. Reading the coin again and summing over outcomes averages the two entries and assigns the same half of their sum to each parity. This follows directly by adding the two orthogonal sharp projectors on either side of the density matrix. Hence the matter channel in each interval is exactly

$$\mathcal M_s(\rho)=\tfrac12 U_+(s)\rho U_+(s)^\dagger+
                       \tfrac12 U_-(s)\rho U_-(s)^\dagger.$$

Composition applies to arbitrary intervals by induction. An unmeasured final coin gives the same reduced matter state as reading it and pooling its labels. The same algebra applies when $\rho$ includes a separate unmeasured reference, by acting identically on that reference. It requires a sign-even input, a first parity-even head independent of the matter at that point, complete retention, a fixed schedule, and no outcome-conditioned intervention. Selecting the full head words gives the different conditional records in the acquisition.

For two readings, let $V_+$ and $V_-$ be the two propagators between the readings. If the first reading gives $a$ and the next $b$, their unnormalized matter amplitude immediately after the second read is $(V_+ +(-1)^{a+b}V_-)\psi/(2\sqrt 2)$; $\psi$ is the matter amplitude at the first read. The probability of a changed coin summed over both values of $a$ is $\|(V_+-V_-)\psi\|^2/4$. Thus immediate repetition gives identical readings for every source, and a later reading can distinguish the launches. The subsequent flight inside each conditioned branch begins with a sharp coin. Summing the final unmeasured parities gives the two terms in $\mathcal M$ again.

## A finite error bound

On this material the cut head costs $4\kappa I$ in the odd-coin sector and $4\kappa(I-P)$ in the even-coin sector. Here $P$ selects a physical connector with both endpoints occupied. Write $H=(H_++H_-)/2$ and $B=(H_+-H_-)/2=-2\kappa P$. Modulo a scalar, $B$ has norm $\kappa$ for $\kappa\geq0$. All terms elsewhere are the same in both sectors. Consequently $H$ has one half of the original interaction at the monitored connector, modulo a constant, while every other interaction keeps its original coefficient.

Let $F(\lambda)$ be conjugation by $\exp[-is(H+\lambda B)]$. Differentiating the exponential twice inserts two commutators with $B$, integrated over the time-ordered triangle of area $s^2/2$ and multiplied by two. Each intervening conjugation preserves the trace norm; each commutator has induced trace norm at most $2\kappa$ after dropping the scalar. Thus $\|F''(\lambda)\|_{1\to1}\leq4\kappa^2s^2$. Averaging Taylor's formulas from zero to $\lambda=\pm1$ cancels the linear term and bounds $\|\mathcal M_s-F(0)\|_{1\to1}$ by $2\kappa^2s^2$. Telescoping the contractive channels extends this bound to $2\kappa^2\sum_j s_j^2$ over any fixed sequence. All channels preserve trace. For each initial density matrix and every terminal detection event the absolute probability error is bounded by

$$\min\!\left(1,\,\kappa^2\sum_j s_j^2\right).$$

For $N$ equal intervals over duration $T$, the bound is $\min(1,\kappa^2 T^2/N)$. The frequent-read limit at fixed duration therefore gives unitary matter evolution under $H$. The bound is loose at the cadence of the two-read experiment. The finite channel calculation and the small-strand sharp-label check give its reported probabilities directly. No assumption about a Klein–Gordon dispersion enters the result.

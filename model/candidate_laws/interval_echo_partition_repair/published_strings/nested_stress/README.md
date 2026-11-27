# Can one stress account for the same memory on three proper balls?

On the published-string apparatus a native two-port memory can change its entropy while preserving the complete distributions of total load and seam energy. A comoving observer can enclose that entire memory in each of three concentric proper balls. Each survey then records the same entropy change. A local energy perturbation confined inside the smallest survey has a different, known dependence on the three radii. This obstructs the identification of the vacuum modular response with a single local stress on all preparations of this apparatus. It also gives a signed experimental requirement for any proposed response in the intervening annuli. [Native memory and instruments](../NATIVE_PAIR.md) · [Positive result for the homogeneous stock](../README.md)

## Transfer measured by the native seam

Use the same flat blank stock, finite region algebras, faithful product reference, and comoving proper-distance instrument as the published-string experiment. Keep the external sheet density and rod setting fixed. Let $U$ be the native exposure at $\theta=\pi/6$. Its initial four spectral tickets have weights $(5,1,1,1)/8$. Transfer probability $x$ from $U|00\rangle$ to $U|11\rangle$, leaving the other two words alone. Both input words are physical opposite-order loads $01,10$ on the first alphabet ($12,21$ on the second). This perturbation can be prepared by modifying the ticket probabilities locally before applying the same gate. We compare the regional states accessible on the memories, with the preparation selector discarded, as in the published protocol.

On those opposite-order loads the actual blank-seam generator is $(17I-\mathrm{SWAP})/2$. Both words have equal squared overlaps with the energy-8 and energy-9 eigenspaces. Evolution by $U$ preserves these overlaps. Their total loads agree. The full energy distribution stays $(1/4,3/8,3/8)$ at $0,8,9$ under the transfer. The complete load distribution stays $(1/8,3/4,1/8)$ at the three allowed totals. Individual port probabilities may change. In the vacuum eigenbasis the slope $s$ of the von Neumann entropy and the change in the modular expectation satisfy
\[
 s=\left.\frac{dS}{dx}\right|_{0}=\log 5,\qquad
 \Delta\langle-\log\sigma\rangle=x\log 5,\qquad
 0\le x\log 5-\Delta S=D(\sigma_x\Vert\sigma)\le(48/5)x^2 .
\]
The last bound follows from $\log z\le z-1$ applied to the four commuting weights. A forward complementary exposure followed by two sharp reads measures the spectral tickets. For $x=1/16$ their exact numerators out of 16 change from $(10,2,2,2)$ to $(9,2,2,3)$. For the first alphabet with left mounting the high-word counts were 10,153 and 9,220 in two independent simulated banks of 16,384. The archive replays the full 99-state generator restricted to the native columns, all four load/mount combinations, and every read bin. These counts concern the single-pair control. The fine-scale estimates below use the exact preparation probabilities. [Native transfer and replay](audit.py) · [Counts, distributions and finite radii](results.json)

## One preparation at every enrolled radius

Extend the experiment's existing geometric enrollment $1/(32\cdot4^j)$ over $j\ge0$. It places a comoving pair at proper distances $s(1\pm1/160)$ on a fixed ray for every rest enrollment $s$. Define
\[
 R_k=(512\cdot1024^k)^{-1},\qquad s_k=R_k/4,\qquad w_k=R_k^4 .
\]
Change each pair enrolled at $s_k$ by $x_k=\epsilon w_k$ for a common small real $\epsilon$. All other stocks and pairs retain their reference preparations. At dyadic $\epsilon$ each finite preparation uses a finite fair-ticket experiment. Increasing the enrollment keeps the same reduced preparation on every collection of pairs the two enrollments share. These compatible finite marginals define the limit. Entropy **differences** and first variations converge absolutely since $\sum w_k<\infty$ and the relative-entropy remainder sums at most $(48/5)\epsilon^2\sum w_k^2$.

Consider the three actual rest surveys of proper radii $R_k,4R_k,16R_k$. The outer end of pair $k$ is at $161R_k/640<R_k$. For $k>0$ the inner end of pair $k-1$ is at $159\cdot1024R_k/640>16R_k$. Monotonicity puts all perturbed pairs $i\ge k$ entirely inside each survey, and all coarser perturbed pairs outside. Unchanged memory pairs give zero entropy contrast even when cut. Additivity on the independent preparations gives the same first variation on all three balls,
\[
 L_k=\log5\sum_{i\ge k}R_i^4=\frac{\log5}{1-1024^{-4}}R_k^4.
\]
In a finite prefix ending at $m\ge k$, the sum stops at $m$ and still exceeds or equals $(\log5)R_k^4$. Thus every scale comparison concerns marginals of one common preparation, with a bound uniform in truncation. The unchanged rods report zero area variation at fixed proper volume on this variation. The original homogeneous-stock test couples a different family of simultaneous sheet and memory changes. [Material worldlines and membership read](../experiment.py) · [Published preparation protocol](../PROTOCOL.md)

## A necessary identity for a local stress

Consider a replacement first-order energy reading on the material rest slice given by one locally finite signed measure $\mu=\delta(T_{ab}u^au^b)\,d^3x$ (a signed measure also allows point or sheet support). Use the common comoving frame, the frozen flat vacuum, its calibrated proper distances, and its vacuum-subtracted flat diamond charge
\[
 B_t=2\pi\int_{B_t}\frac{t^2-r^2}{2t}\,d\mu,\qquad r=|x|.
\]
Ball selectors act by restriction of this single measure. Whenever the stress variation vanishes on $R<r<16R$, putting $E=\int_{r\le R}d\mu$, $M=\int_{r\le R}r^2d\mu$ gives $B_t=\pi(tE-M/t)$ for $t=R,4R,16R$. Atoms at $16R$ have zero kernel. This assumption holds for stress depending only on the local changed memories, with support confined to a small neighborhood of each changed pair disjoint from the intervening annuli. It is a physical locality hypothesis a stress readout can test. A state-independent additive scalar in the vacuum modular generator drops out of every normalized state difference.

The weights $(4,-17,4)$ annihilate both $t$ and $1/t$ at these radii. Write $e_t=L_k-B_t$. Then $4e_R-17e_{4R}+4e_{16R}=-9L_k$, and hence
\[
 \max_{q=1,4,16}\frac{|e_{qR_k}|}{(qR_k)^4}
 \ \ge\ \frac{9L_k}{(4+17\cdot4^4+4\cdot16^4)R_k^4}
 \ \ge\ \frac{9\log5}{266500}=5.4352499857\times10^{-5}.
\]
The constant $9/266500$ is optimal for unrestricted signed moments on this triple. In units $R=L=1$, let $c=9/266500$ and fit the two moments to $B_1=1+c$, $B_4=1-4^4c$. The same identity forces $B_{16}=1+16^4c$. The three normalized residuals are $(-c,c,-c)$. Positive and negative moment contributions are both allowed by this bound. For fixed finite $\epsilon$, regional modular changes are exactly $\epsilon L_k$. Regional entropy changes differ from them by $O(\epsilon^2R_k^8)$, so both readings have the same obstruction to fourth-order matching. One comoving observer's failure suffices for a claim quantified over all ten frames. [Formal signed-moment and error identities](../../../../../Lean/ObserverPatchHolography/EinsteinBranch/RegionalMoments.lean) · [Exact rational moment and sign checks](test_nested_stress.py)

## Which annular signal could supply the missing response?

Dropping the support assumption gives a measurable signed requirement. For $Q_R(r)$ the same combination of the three kernels, including $2\pi$, direct evaluation gives
\[
 \frac{Q_R(r)}{\pi}=
 \begin{cases}
 0 & 0\le r\le R\ \text{or}\ r\ge16R,\\
 4(r^2/R-R) & R<r\le4R,\\
 64R-r^2/(4R) & 4R<r<16R .
 \end{cases}
\]
The kernel is positive between the inner and outer radii and peaks at $60\pi R$. To match the positive entropy slope to $o(R_k^4)$ a responsive field would have to deliver $\int Q_{R_k}\,d\mu=-9L_k+o(R_k^4)$. In particular its negative variation in the annuli would be at least $9L_k/(60\pi R_k)-o(R_k^3)$. This is a necessary property of the **change** in energy density, with the ordinary coupling held fixed. The total energy can stay positive. A positive-only extra response at this preparation cannot supply the measured sign.

There is a stronger regularity test even for stress extending into the annuli. If $|\mu-\tau_0\,d^3x|(B_t)=o(t^3)$ at the central event, the kernel bound $\pi t$ gives
\[
 \lim_{t\downarrow0}B_t/t^4=(8\pi^2/15)\tau_0 .
\]
This condition includes a continuous density and a density having this event as a Lebesgue point. The survey entropy slope divided by $t^4$ has two subsequential limits, $C$ on $t=R_k$ and $C/16^4$ on $t=16R_k$, where $C=(\log5)/(1-1024^{-4})$. For any limiting regular stress reading, the limsup of the absolute normalized discrepancy is at least $C(1-16^{-4})/2$. A stress-only response fitting all the bands would need signed annular changes and failure of this regularity condition at the accumulation event.

The result applies to this flat reference, the accessible regional memory algebra and the preparation space allowing independent spectral transfers. The infinite preparation has structure at every scale. Its realization as a freely formed, smooth macroscopic source under the active Interval Echo law requires a dynamical argument. Any fixed finite truncation becomes unperturbed on sufficiently small balls inside its innermost pair. At each finite enrollment a smooth function of position can interpolate the chosen probabilities at the finitely many distinct pair centers, using disjoint neighborhoods. This gives no uniform bound on derivatives of all orders as the enrollment grows. The finite three-radius estimates exclude matching uniform over unrestricted preparations. The compatible limit excludes pointwise matching at its accumulation center if that limit belongs to the claimed sector. A smooth-sector claim needs an independent dynamical or operational criterion selecting the allowed preparations. Another route is a changed vacuum or regional matter algebra with different whole-pair modular gaps, together with independently measured stress and geometry on the affected sector. [Matched Einstein theorem's common-source requirements](../../../../../code/matched_einstein_limit/PROOF.md#the-common-source-interface)

A restriction on changes of vacuum follows directly from the existing native gate. Suppose the replacement faithful reference has $-\log\sigma=f(H,Q)$ with $H$ the native energy and $Q$ the total load. Since $H$ commutes with $Q$, the pulse $e^{-2i\theta H}$ leaves $\sigma$, and hence both of its endpoint marginals, invariant. The endpoint entropy slope about that reference is zero at every exposure. The measured slope about the published-string reference is $(\sqrt3/4)\log(5/3)=0.22119398\ldots$ per radian. The replay compares both slopes on the native seam using $\sigma\propto e^{-(H+Q)/8}$ as a control. Interpreting a Gibbs reference of the old charges as the physical vacuum would therefore require a different source of the positive entropy response. [Pulse and vacuum control](audit.py)

Run the checked native replay and the exact moment tests from the repository root with the normal project Python environment.

~~~sh
PYTHONPATH=code:. python -m model.candidate_laws.interval_echo_partition_repair.published_strings.nested_stress.audit
PYTHONPATH=code:. python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/nested_stress
~~~

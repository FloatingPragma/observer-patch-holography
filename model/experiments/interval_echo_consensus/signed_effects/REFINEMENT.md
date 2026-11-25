# A coherent observer refinement and its free-flight test

One coarse owner has twelve port addresses. Its four actual child owners each have twelve corresponding local addresses. On both signs of a single material defect, the protected refinement preserves the whole coarse matrix algebra, including address coherences. Every one of the 100 parent blocks in the two stored refinement steps implements the same isometry. Unrestricted native flight has a strictly positive first-order departure from this encoding for some input in **each** parent, regardless of the four independent local address charts. [Physical programs, integer leakage matrices](refinement.py) · [Results](refinement.json)

## What the parent and child observers read

Write \(a=0,\ldots,11\) for local port, \(c=0,\ldots,3\) for child, and \(s\) for vacancy or excess. The encoding and public observable map are
\[
V|a,s\rangle=\frac12\sum_c|c,a,s\rangle,\qquad
\iota(A)=I_4\otimes A,\qquad V^*\iota(A)V=A .
\]
The coordinate and child labels are supplied by the carrier. The protected sphere instrument starts with the unknown state on child 0. For each port \(a\), split from child 0 into child 1, then from 0 into 2 and from 1 into 3. Route the two addressed positions into contact by known native \(\pi\)-swaps inside the four-child region and undo the route. Each splitter uses \(\pi/2\). A reference held for \(3\pi/16\) on the new coordinate removes the relative \(i\) for both source signs: \(\exp[-i q_s 3\pi/16]=-i\), where \(q_s=-8,24\). Every terminal branch traverses two splitters. The remaining common phase is independent of \(a,c,s\). Only three partial exchanges and three phases per address are needed. Punctured reference sites are reached by reversible sharp exchanges. Across all 100 parents the 36 exchanges and 36 phases need another 354–402 routing swaps per parent. The checks apply the native gates to every input column on both signs. They therefore apply also to inputs entangled with an untouched reference.

Replicate a sign-preserving unitary \(U\) on the four children. An edge acts identically in each child; degree-six addresses obtain their phases by routing to a degree-five reference. Across child blocks the operations act on disjoint defect positions. The background contributions from these references are common global phases. The logical identity is \((I_4\otimes U)V=VU\). General coarse effects have the same pullback identity, and the four-rail compatible instruments from the signed-cell experiment can be replicated with the same outcome labels. Products of the replicated gates and addressed outcome Kraus operators intertwine by induction. If a terminal observer actually records the child index and then discards it, child coherence is lost but every subsequent replicated coarse record still has the same probability. Repeating the coordinate refinement \(L\) times gives a uniform \(4^L\)-child amplitude and the algebra embedding \(I_{4^L}\otimes A\). These identities concern corresponding protected preparations and control schedules, with an arbitrary external reference. The native boards realize the first two steps. The broader observer contract includes freely selected overlap histories, geometry and a refinement-natural current implementation. [Single-region controls](README.md) · [Earlier classical ancestor-pair experiment](../refined_projection/README.md)

## A finite quantum entropy maximizer and an invisible variation

For any signed coarse density \(\rho\), fix the entire coarse marginal by requiring \(\operatorname{Tr}_c\sigma=\rho\). Competitors include every fine signed density with those coarse expectations, arbitrary coherences and any correlations between child and address. The unique largest-entropy state is
\[
M(\rho)=I_4/4\otimes\rho,\qquad
D(\sigma\Vert M(\rho))=\log4+S(\rho)-S(\sigma)\ge0 .
\]
The support of any such \(\sigma\) lies inside child space times the support of \(\rho\): its positive diagonal mass on the orthogonal complement is zero. Thus the identity also holds for singular \(\rho\). It follows by inserting \(\log M(\rho)=-\log4+I_4\otimes\log\rho\) on that support. Equality of relative entropy to zero forces \(\sigma=M(\rho)\). More generally with any supplied faithful coarse reference \(\tau\),
\[
D(\sigma\Vert I_4/4\otimes\tau)=
D(\rho\Vert\tau)+D(\sigma\Vert M(\rho)).
\]
With any positive weights on all resolution marginals including the finest, iteration of this identity implements the same unique conditional minimizer at every fixed coarse marginal. The faithful uniform signed-position and child references can be prepared by scheduled uniform sharp coordinates. Choices of reference, physical parent map and measured constraints remain explicit inputs.

There is a native preparation, distinct from the coherent lift. On child \(c\), multiply every address by \((-1)^{k\cdot c}\), using the two-bit child label and a scheduled two-bit mask \(k\). A physical reference held for \(\pi/8\) on a port gives the same relative \(-1\) for either sign. Applying this phase on the twelve addresses selected by a child sign implements the word, with reversible routing where required. The four masks give mutually orthogonal child states. Use equal numbers of all four masks, discard the mask register, and the resulting physical ensemble is exactly \(M(\rho)\). The executable checks all four 48-rail words on all parent codes; applying the inverse encoding returns to the first child with certainty for mask zero and never for the other masks. Complete cycles of four scheduled starts therefore give recovery \(1/4\); the unmixed coherent code gives one. Access to the retained mask would change the accessible state.

The mixture
\[
\sigma_\lambda=[(1-\lambda)I_4/4+\lambda|u\rangle\langle u|]\otimes\rho,\qquad
u=(1,1,1,1)/2
\]
has exactly the same full coarse state and all the same protected coarse histories at every \(0\le\lambda\le1\). Its fine child eigenvalues are \((1+3\lambda)/4\) and three copies of \((1-\lambda)/4\). At \(\lambda=1/2\) the scheduled mask multiplicities \((5,1,1,1)\) give recovery \(5/8\), and the fine entropy derivative is \(-3\log5/4\). The coherent and uniform endpoints differ by \(\log4\) nats. A fixed information capacity and complete coarse agreement thus coexist with a signed fine entropy response in the executable model. The proper-volume, stress and small-ball stationarity conditions in the Einstein theorem require independent geometric measurements.

## Why the protected identification fails for autonomous flight

Allow each child to rechart its twelve-dimensional address space by any unitary \(U_c\). The lift on one parent is then \(V=(U_0,U_1,U_2,U_3)^T/2\). Let \(P\) be the full lifted parent projector on the fine sphere, and \(H_s\) the actual all-open signed-defect Hamiltonian. Its off-diagonal coefficient is \(-1/2\) on every physical seam for either sign. Each coarse owner has six fine seams into different parents. The six sending and receiving port positions are distinct because external contacts form a matching. The squared Hilbert–Schmidt norm of fine amplitude sent out of the parent by \(H_s V\) is therefore exactly \(6/16=3/8\), independent of all \(U_c\): each sending row of \(V\) has squared norm \(1/4\).

At each adjacent receiving parent the two external crossing seams enter at most two children. For any vector \(w\) supported in at most two of its four children,
\[
\|V^*w\|^2
=\tfrac14\left\|\sum_c U_c^* w_c\right\|^2
\le\tfrac12\sum_c\|w_c\|^2 .
\]
At least half the outgoing squared norm thus lies orthogonal to the entire receiving code. Summing over all neighbors and twelve orthogonal inputs gives exterior leakage at least \(3/16\) per sending parent. Some unit input in every parent consequently obeys
\[
\|(I-P)H_sV\psi\|\ge\sqrt{(3/16)/12}=1/8 .
\]
The coefficient is in native Hamiltonian time. Replacing or rescaling the coarse generator adds an in-code vector, leaving this transverse derivative untouched. The argument remains valid after any separate unitary change of chart on each child. Unequal child weights, time rescaling on the fine board in a scaling limit, other carriers or a different chosen coarse observable have separate tests.

For the original aligned port coordinates the verifier constructs the full integer numerator \(4BW-W(W^\mathsf{T}BW)\), where \(V=W/2\) and \(B=2H_s\) after subtraction of a scalar. This matrix is exactly \(16(I-P)H_sV\). There is nonzero transverse response on 120 of 240 coarse sharp columns and 444 of 960 on the next step for each sign. The minimum exterior squared trace per parent is \(3/16\) on both steps. Thus a protected internal observer alignment can be exact while autonomous signed propagation fails that same finite intertwining condition. The model's other propagating sectors and the OpenWave continuum criteria require their own comparison. [Independent native signed free generator](../../interval_echo_angular/charged_orientation/measure.py) · [Geometric parent certificates](../../mobile_memory/refinement.py)

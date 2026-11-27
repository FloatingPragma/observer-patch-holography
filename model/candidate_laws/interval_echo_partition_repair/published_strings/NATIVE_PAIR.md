# Endpoint entropy from the native Interval Echo seam

The two-qubit memory in the parallel-sheet experiment can be prepared on an isolated blank seam of unmodified Interval Echo. Three independent fair choices select a sharp two-port load, then free exchange creates the required coherence. A single port has precisely the entropy and modular change used by the small-ball experiment. A coherent read of both ports gives the predicted whole-pair entropy and relative entropy. The calculation uses actual columns of the original 99-state seam Hamiltonian. [Native write and local generator](../../contrast_memory_spherical_v1/law.py) · [State, readout and terminal banks](native_pair.py) · [Probability and entropy checks](test_native_pair.py)

## Which loads carry the memory?

Choose physical loads $q,q+1$ at each end. Both $q=0$ (vacancy relative to unit filling) and $q=1$ (excess) work. At the left port use the logical letter $a=n_L-q$. At the right port use $b=1-(n_R-q)$. Exchanging left and right gives the second mounting. In this order the two logical words $00,11$ are physical $q(q+1),(q+1)q$. Both have total load $2q+1$. The other two logical words have two equal loads and are spectators.

On blank memory the native write, reversal generator and unwrite reduce exactly to $H=(I-\mathrm{SWAP})/2+8(n_L-n_R)^2$. Binding the full operator to these four allowed columns gives zero leakage. The two equal loads have energy zero. On the other two,
\[
 H= (17I-\mathrm{SWAP})/2,\qquad
 e^{-2i\theta H}=e^{-17i\theta}(I\cos\theta+i\,\mathrm{SWAP}\sin\theta).
\]
Assign the initial logical $00$ word to five of eight fair tickets and each other word to one. The joint state after exposure $2\theta$ is $I/8+|\psi_\theta\rangle\langle\psi_\theta|/2$, with $|\psi_\theta\rangle=\cos\theta|00\rangle+i\sin\theta|11\rangle$. One fixed logical phase converts it to the real rotation convention of the endpoint model. The state and the charge-one interference measurement are block diagonal in total load, including the total-three realization at $q=1$. Preparing it needs no superposition of distinct conserved total loads. Eigenweights of the full state are $(5,1,1,1)/8$ at every exposure. This uses the original native interaction, three fair preparation commands and a selected exposure.

## Separate endpoint and joint reads

At either logical endpoint the probability of one is $p(\theta)=1/4+\sin^2\theta/2$. Use a known comparison port of load $q$, perform the native echo write, read the absolute contrast without resolving the coin, then unwrite. The two projectors pulled back to the endpoint are precisely the sharp digits $0,1$ on $q,q+1$. Both orientations of the physical comparison give this instrument. At the reversed port complement its reported digit. The reference angle is $\theta_*=\pi/6$, hence $p_*=3/8$ and $p'(\theta_*)=\sqrt{3}/4$. With $h$ the binary entropy the cut slope is
\[
 {d\over d\theta}h(p(\theta_*))={\sqrt3\over4}\log(5/3)
 =0.22119398354\ \text{nats},
\]
or $0.11059699177$ nats per native exchange time. The finite change is exactly $h(p)-h(p_*)$ and its modular charge is $(p-p_*)\log(5/3)$. Their difference is the binary relative entropy. On a fully included pair the entropy change vanishes; its modular charge and relative entropy both equal $\tfrac12\log5\,\sin^2(\theta-\theta_*)$. The nonzero phase between logical $00$ and $11$ matters for this whole-pair conclusion.

The prepared spectral basis can be measured by letting the same seam evolve forward for the complement of its exposure modulo $2\pi$, then performing two sharp port reads. All three original energies are integers, so the complement implements the inverse on this entire code. At $\theta_*$ a coherent inverse gives high-word probability $5/8$. A sharp endpoint read inserted before the inverse removes the relevant coherence, giving $7/16$ instead. With 16,384 independent simulated trials per bank, the four charge/mount combinations gave high-word counts $10186,10195,10279,10169$ after the coherent inverse and $7224,7187,7289,7202$ after the intermediate read. All other outcomes are included in the archive. A single terminal sharp read on each endpoint, without the inverse, would give the entropy of the dephased readouts. The two instruments have different outcome laws. The exact native written-sign reversal on the mixed-load sector also reads exchange parity. Along with a quarter-period exchange followed by a load read and a direct load read this gives the three Bloch components within the conserved sector. [All banks, held signed increments and native energies](native_pair.json) · [Native written-parity pullback](../../../experiments/material_inertia/sourced_survey/README.md)

## Which coupling is external?

The native energy probabilities are $1/4$ at zero and $3/8$ at each of 8 and 9 throughout this rotation. Its actual mean is $51/8$ in exchange units, with derivative zero. The endpoint's vacuum modular gap $\log(5/3)$ is a property of this preparation. Converting the sheet's independently measured regional work to an additional exposure by dividing by $[\log(5/3)/(2\pi)]\,\sqrt3/4$ gives entropy slope $2\pi$ per unit applied work. That pulse setting requires the proper radial weight, the observer direction and the sheet stress specified by the extension. The seam alone supplies neither the world-sheet action nor the source-dependent rod deflection.

The joint operation is local when the two loads share the enrolled seam. A mounted apparatus may subsequently move sharp or coherent loads using addressed full swaps, and park them with interfering exchange switched off between observations. An inverse read of a separated pair requires collecting both loads coherently before either sharp read. The published-string construction supplies the enrollment, placement, schedules and geometry used to interpret these quantum memories on overlapping balls. The native result here covers the preparation, free gate and endpoint instruments on the source seam. [Exposure and register primitives](../../interval_echo/SPEC.md) · [Source-to-ball interpretation](README.md)

~~~sh
python -m model.candidate_laws.interval_echo_partition_repair.published_strings.native_pair
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/test_native_pair.py
~~~

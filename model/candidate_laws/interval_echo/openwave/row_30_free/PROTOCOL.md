# Free doubled-load spectroscopy

## Source and flight

The sphere has 20 echosahedra, 240 load ports and 630 seams. A basis state carries one load of two units at a port, with the other ports vacant. All seam memories begin blank. The native 99-state term preserves this sector. The code computes its two-port block from the full contrast write at κ = 8, sums the actual native terms over the seam list, and includes the private contrast reference on each degree-five port. Degree-six ports have six physical seams. The readout code verifies the outgoing detector effect on the doubled-left, doubled-right and empty inputs of the same native echo. [Native law](../../../contrast_memory_spherical_v1/law.py) · [Carrier](../../../../executable_model/carrier.py)

Let $u$ be the normalized uniform mode. Choose five eigenvalues of the unweighted incidence Laplacian nearest to 0.12, 0.32, 2.9, 5 and 7. In each eigenspace, project a port of the first cell and normalize the result to obtain a real mode $v$. Use the port with the largest projector diagonal in that cell, breaking ties by port index. The geometry fixes this choice before the free dynamics are evaluated. The five selected eigenvalues are distinct. The values 5 and 7 are selection targets, and the actual eigenvalues for those targets are 6 and 6.846204.

The designer supplies each real profile and both coherent reference superpositions $(u+v)/√2$ and $(u-v)/√2$. The code compiles each superposition from a sharp double on port zero. It distributes its weight along a fixed breadth-first spanning tree with native two-port exchange pulses. An addressed degree-five reference fixes each output phase. At a gateway, a π exchange takes that amplitude to a native degree-five neighbor and a second π exchange restores the port ordering. Each of the ten preparations uses 359 seam contacts and 240 addressed degree phases. These are controlled sources. Preparations involve the whole sphere. The flight itself activates all seams together and uses no route, alternating gates, repeated measurements or protected trajectory. The stationary uniform component supplies a phase reference within the same charge sector.

## Spatial observable

On a seam $ab$, the native echo writes outgoing loads $(1,1)$ whenever the double occupies $a$ or $b$. The stored contrasts of these two alternatives are orthogonal. A quiet outcome includes the possibility that the double occupies another port. The resulting probability is $s_{ab}=p_a+p_b$. In a native triangle $(a,b,c)$, three separate terminal seam measurements therefore reconstruct $p_a=(s_{ab}+s_{ac}-s_{bc})/2$, without a vacant witness rail. All three measurements use fresh preparations.

For each port use its lexicographically first incident face. A spatial scan on repeated preparations measures both signs of the source before flight. Each sign uses independent binomial counts on each seam in the carrier's sorted seam list. Shared seam estimates can serve several triangles. With $N=240$, their difference reconstructs the signed real mode $v_a=√N(p_a^+-p_a^-)/2$. Its squared wave number is the normalized Dirichlet energy

```text
λ = Σ_{ab ∈ seams} (v_a − v_b)² / Σ_a v_a².
```

This graph-distance convention uses the registered seam incidence. A whole-board spatial coefficient is assembled across repeated finite-patch measurements. The protocol assumes access to those comparison records and controlled preparations. It supplies no material clock for the whole board.

## Temporal observable and held data

Propagate a separate plus-source preparation through the assembled sparse native Hamiltonian with a matrix-exponential action. Subtracting the uniform energy removes a common phase from the calculation. The selection of the temporal frequency uses only the receiver counts, the supplied exposure intervals and a fixed search range. It receives neither the spatial coefficient nor the restricted generator spectrum as input.

For each selected mode the receiver measures a triangle in cell zero at 513 exposures spaced by 0.375, through exposure 192. The interval is a setting of the free Hamiltonian. Each edge and exposure uses an independent terminal measurement on a fresh preparation, with $2^{25}$ simulated shots, seed 430029. All shots enter each edge probability. Reconstruct the receiver port with the triangle identity above. Fit an arbitrary constant, sine and cosine amplitude on the 257 even-indexed exposures over frequencies 0.003 to 4.15 radians per unit. The even sample spacing of 0.75 gives a Nyquist angular frequency above 4.18. Evaluate the fitted signal on the 256 interleaved odd exposures without refitting. The receipt includes both simulated counts and the fits using noiseless Born probabilities. It tabulates seam counts in the order fixed by the carrier hash. Temporal counts follow $(ab,ac,bc)$ in the named triangle.

## What the curvature tests

In the doubled sector the assembled Hamiltonian equals $H=24κI+L/2$. The coefficient follows from a native exchange of amplitude one half and six native or degree-reference contributions of $4κ$ per port. A contrast between the uniform mode and a mode of spatial eigenvalue $λ$ therefore beats at $β=λ/2$. A common frequency offset is invisible to these population measurements.

Suppose a positive-frequency Klein–Gordon branch with rest frequency $m$ and one speed $c$ described these same modes and the adopted squared wave number. It would give

```text
β = sqrt(m² + c² λ) − m
c² = (2mβ + β²) / λ = m + λ/4.
```

Distinct positive $λ$ require different squared speeds for every finite $m$. This is an exact obstruction for the specified sector and spatial convention. For finite precision, calibrating $c$ to the first mode at an assumed $m$ gives a useful sensitivity check. For assumed rest frequencies 1, 10, 192 and 1,000, the largest held relative beat errors from the sampled data are 46.34%, 12.76%, 0.886% and 0.194%. Choosing a large $m$ makes the exact square-root curve approach a linear relative beat over a bounded range. The counts alone cannot exclude an arbitrarily large $m$ at free $c$. The absolute energy 192 comes from the operator and has no independent rest-clock calibration in this experiment.

Other excitation sectors, operational definitions of spatial scale, patch clock units and continuum limits require their own measurements. This run measures the free sector on the finite spherical carrier. [Experiment](experiment.py) · [Data](results.json)

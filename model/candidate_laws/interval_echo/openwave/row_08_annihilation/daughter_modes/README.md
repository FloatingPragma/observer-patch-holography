# Which excitations reach the two wings?

After the opposite walls disappear from the source ring, the two receiving wings can each carry **one excitation of a gapped material wave**. In the full energy read at exposure 32, 2,064,278 of 2,066,989 ordered-ring records carry exactly one mode in each wing. There were 67,108,864 launches. In 99.55% of those one-per-wing records, exactly one wing carries the middle-frequency mode and the other carries an outer-frequency mode. A second read distinguishes a coherent pair of frequency assignments from a statistical mixture. [Protocol](PROTOCOL.md) · [Complete counts, controls and frequencies](results.json)

## Deriving the frequencies

Keep the ordered ring endpoint fixed at its measured load, zero or two. Write each dimer's (2,0) and (0,2) configurations as a two-state system. The native exchange matrix element gives a transverse field h = 1/2. Half the difference between unlike and equal endpoint energies gives J = (4κ + 1/2)/2 = 16.25. Reversing alternate orientation labels turns the guide into a ferromagnetic transverse-field Ising chain with a fixed spin at the inner boundary and a free outer end. Reversing the measured ring phase gives the same frequencies.

For n dimers add the fixed boundary as a site with transverse field zero. The two Majorana operators at each site can be written as the product of X on all earlier spins times Z or Y on that site. Under this change of basis each transverse field joins the two Majoranas at its site and each Ising interaction joins successive sites. Omitting the isolated boundary Majorana leaves an alternating open path with 2n+1 vertices. Its rectangular coupling matrix has the two diagonals J and h. Its squared singular values are the eigenvalues of the n-dimensional tridiagonal matrix with diagonal J²+h² and off-diagonal −Jh, where alternating signs fix the wave-number convention. The zero mode accounts for the two possible fixed boundary sectors. Choosing one boundary leaves n independent two-level occupancies.

The sine eigenvectors give standing wave numbers kⱼ = πj/(n+1) and one-mode excitation costs

$$
\omega_j = 2\sqrt{(J-h)^2+4Jh\sin^2(k_j/2)},\qquad j=1,\dots,n.
$$

Each three-dimer wing therefore has three one-mode costs, 31.80076, 32.51538 and 33.21463 exchange units. The eight subset sums are distinct. These values assign all the observed wing energies to occupations without using the collision data to fit a frequency. With the attachment cut, the square bidiagonal coupling matrix has entries h on its diagonal and J on its adjacent diagonal. The lowest mode is then a soft boundary mode. The native and zero-contrast connected arms and the cut arm each use their own frequency dictionary for this energy analysis. Finite spin matrices from endpoint loads independently check the complete spectra on both boundary phases for lengths 1, 2, 3 and 5. [Mode calculation](measure.py) · [Direct finite-chain checks](test_modes.py)

A uniform extension of the same guarded native wing has one positive-frequency branch, gap 31.5 and maximum group speed one dimer per exchange unit. The maximum follows by differentiating the dispersion and occurs at cos(k) = min(J/h,h/J). These are properties of that material guide. They do not identify its excitations with transverse photons or show that the incoming electric charges disappear.

## Reading a phase between the daughters

After an ordered ring read, each wing admits only the two lowest one-mode energies of its measured boundary phase. The modes stay coherent on admission. The outer-dimer swap fixes the relative sign of their real amplitudes by a ground-to-mode matrix element. A supplied coherent analyser compares either the energy labels (Z) or their equal superpositions (X). Both wings use the same choice. All 19 final outcomes include both reported ring phases, each local rejection and a source miss. The input flight is the same time-independent, guarded material collision as in the spectral acquisition. The local filters, central flag, exchange-time appointment and local superposition reads are supplied measurement services.

For a product of logical qubits with Bloch vectors u and v, Cauchy–Schwarz gives |uₓvₓ|+|u_zv_z|≤1. Convexity gives the same threshold for any separable mixture upon this fixed admission. The ordered ring may localize correlations in its two wings. This test characterizes their state *after* that published central observation.

| Native wall release | Accepted X/X trials | X/X correlation | Accepted Z/Z trials | Z/Z correlation |
| --- | ---: | ---: | ---: | ---: |
| Coherent wings | 139,485 | +0.99104 | 138,422 | −0.98854 |
| Energy dephased | 138,182 | −0.00230 | 138,642 | −0.98849 |

Each setting used 8,388,608 launches. About 1.65% of all launches pass the central and local filters. A simultaneous 99% lower bound for |⟨XX⟩|+|⟨ZZ⟩| is **1.9574** on the coherent native pair. The dephased control retains the energy correlation; its X correlation vanishes in the model. The native sham yields only 17 and 14 admitted trials on the two coherent settings. Neither the cut nor the zero-contrast wall run resolves this witness at the common native receiver settings. The intervals use bounded ±1 sample means conditional on the accepted count and a union bound across all 24 settings. The reconstructed native wing excitations carry energy and bipartite entanglement; an independent radiation or particle–antiparticle identification is a further obligation.

From the RER root, replay with `OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_08_annihilation.daughter_modes.measure`.

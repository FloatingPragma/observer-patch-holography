# A compact lowest mode from native seam geometry

This supporting theorem specifies a candidate construction. It supplies no energy measurement, apparatus compiler, physical mass identification or external review. The target eigenvector is derived from graph symmetry; dynamic preparation must be supplied by a separately tested process.

Take a connected N-port core with unweighted Laplacian L0 and spectral gap at least mu>0. Multiply every core seam weight by c. Add two ports p,q, each joined to the same core port 0 with weight w, and join p to q with weight v>0. The weighted native one-particle transfer is T=I-L/(2W), where W is total seam weight. This identity follows by counting the two fair endpoint outcomes on each addressed seam; no Hamiltonian is installed.

The normalized contrast f=(e_p-e_q)/sqrt(2) has Lf=(w+2v)f. It is a strict lowest nonconstant mode under the sufficient bound below. Thus the lowest mode is spatially supported on two ports, even when the core is large.

Choose w=2Nv. Require c*mu > (w+2v)+2w(2N+3)(N-1)/N. Every mean-zero vector orthogonal to f has equal leaf coordinates a. Write its core vector as x=m*1+z, with sum z=0. Mean zero gives a=-Nm/2. Its squared norm is beta*m^2+||z||^2, beta=N(N+2)/2, and its energy form is c*z^T L0 z+2w(alpha*m+z0)^2, alpha=(N+2)/2.

For epsilon=1/[2(N+2)], (u+v)^2 >= (1-epsilon)u^2-(1/epsilon-1)v^2. Also z0^2 <= (N-1)||z||^2/N. Substitution bounds the m part below by (w+3v)*beta*m^2 and the z part below by [c*mu-2w(2N+3)(N-1)/N]||z||^2. Both coefficients strictly exceed w+2v, proving the assertion. The full graph is connected, so its zero mode is exactly constant. This is a variational proof, not a conclusion from a numerical eigenvalue ordering.

For B original twelve-port cells joined by all matching-port adjacent-cell seams, the core is the Cartesian product of the cell graph and a B-vertex path. Its gap is the smaller factor gap. The cell has gap 5-sqrt(5)>2; for B>=2 the path gap is 4sin^2(pi/(2B))>=4/B^2. A rational usable bound is mu=min(2,4/B^2). With v=1 and w=2N, choosing c as a power of two above the sufficient threshold gives examples B=1,2,3 with c=1024,8192,32768.

The declared protocol may address the same leaf triangle with other observer programs, but that alone does not prove a common controller or transfer energy/spin interpretations. Source-arc unitary dynamics, state acquisition, operational phase/energy reading, independent-copy additivity and full outcomes require their own implementation and checks.

Verification script input: original thirty-seam cell, B=1,2,3, fixed formulas above. Exact rational f identities and conservative lower bounds must pass; numerical spectra are an additional arithmetic check. A failure of the spectral inequality, endpoint transfer identity or support identity rejects this construction. This calculation is an analytic construction check, not a prospective empirical prediction.

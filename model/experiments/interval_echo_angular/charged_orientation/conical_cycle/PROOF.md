# Why the two paths have opposite orbital phase

Write |j⟩ for the ternary word having the unique vacancy or excess on face port j in a unit background. Blank native contacts preserve these three words as a space, and each private reference is diagonal. A terminal material address therefore identifies the signed load in every trial. Let u = (1,1,1)/√3. Summing the three native face contacts gives

\[
H_\triangle=(2\kappa+3/2)I-(3/2)uu^T.
\]

Each word sees two unequal contacts. Each such contact moves the defect to the neighboring port with coefficient −1/2. These facts give the diagonal and off-diagonal entries. The relative response at the exceptional port to its private reference is \(f=(1+q)^2-1\), equal to −1 for q = −1 and 3 for q = +1. The addresses, native write and contrast supply this coefficient. The experiment solves the exact product of scheduled pulses. It uses no time-average approximation. [Source reduction](measure.py) · [Independent word-space check](test_cycle.py)

## A single conical degeneracy

Remove common phases from a single step. Write \(x=b+\cos\theta,\ y=\sin\theta\) and \(v_j=0.35 f(x\cos\phi_j+y\sin\phi_j)\). With \(D_{jj}=e^{-ihv_j}\), the step is

\[
F=D^{1/2}\{I+(e^{3ih/2}-1)uu^T\}D^{1/2}.
\]

It is symmetric and unitary, so its real and imaginary parts are commuting real symmetric matrices. Every simple eigenline admits a local real unit eigenvector. An accumulated adiabatic geometric phase on such a line is its orientation sign. Both arms may carry arbitrary dynamical phase.

The step is similar to a rank-one perturbation of the diagonal matrix D, with all components of the perturbing vectors nonzero. For an eigenvalue different from the three diagonal entries, its eigenvector is fixed up to scale by the rank-one equation. For an eigenvalue equal to a diagonal entry of multiplicity m, that equation forces the rank-one scalar to vanish. The remaining vectors have support on its m equal ports and satisfy one nonzero linear constraint, giving dimension m−1. Unitarity makes algebraic and geometric multiplicities equal. A repeated eigenvalue therefore requires three equal diagonal entries. In the disk containing the tested controls, h times the largest possible difference of the v values is at most 1.26, below 2π. Equality of the three phases requires x = y = 0. There the two-dimensional space perpendicular to u is exactly degenerate. Along these bounded controls all eigenphases can be ordered on a single cut after removing their common phase. Indeed the total phase excursion from the identity is at most h(3.15+1.5) < π by the sum of the spectral angle bounds of the three unitary factors.

In the ordered even and odd basis (2,−1,−1)/√6 and (0,1,−1)/√2 of the doublet, first-order degenerate perturbation of the palindrome gives the compression of the private generator. Apart from its nonzero scalar 0.35 f/2 this is

\[
\begin{pmatrix}x&-y\\-y&-x\end{pmatrix}.
\]

On a circle its normalized eigenvector turns through half the polar angle. Each of the two split real eigenlines therefore reverses its orientation once around the degeneracy. The separated u line extends over the center. The absence of any other degeneracy extends these orientation signs by continuous deformation to the finite enclosing and displaced loops. The supplied odd preparation is an exact endpoint eigenvector. Along the positive x axis it belongs to one of the two split lines by continuity from the center and the absence of crossings. A loop enclosing the center once gives a negative adiabatic sign for this source. A loop displaced outside the center or enclosing it twice gives positive sign. The smooth step families have simple eigenphases with a positive minimum separation on each compact path, as required by the discrete adiabatic theorem. The recorded probabilities come from independent finite products. [Tanaka, discrete adiabatic theorem](https://arxiv.org/abs/1107.2989) · [Berry, geometric phase of real eigenlines](https://doi.org/10.1098/rspa.1984.0023)

## What an observer reads at finite speed

Let S exchange the two side ports. In each arm the opposite exposure sequence obeys \(F(-\theta)=SF(\theta)S\) exactly, including the common phases. Thus both arms have the same instantaneous spectrum. For initial odd state a, the two terminal states L and R obey R = −SL at any speed. Recombination on a balanced coherent arm gives amplitudes \((I-S)L/2\) at the bright exit and \((I+S)L/2\) at the dark exit. These are orthogonal complementary mirror-parity subspaces. Their summed probabilities are one without selecting a band or a return event. In the adiabatic limit the matched dynamical phases cancel. For the finite experiment departures from perfect bright or dark report transitions into a different mirror parity. [Lean proof for an arbitrary pulse list](../../../../../Lean/SourceLaw/OrbitalMirrorCycle.lean)

Under actual spatial template turns, an isolated native π contact acts as an ordinary port interchange because \(e^{-i\pi\kappa(n_a-n_b)^2}=1\) at κ = 8. The pulse compiler implements the complete proper 120° turn about this face, whose cube fixes all twelve addresses, and two perpendicular template half-turns with commuting address maps. Their complete unfiltered lifts give identical arms in both spatial controls. Winding the private potential changes a material orbital by encircling a conical intersection in control space. A law-selected spinor would require an independent physical lift of rotations and a consistent exchange on the charged excitation. [Native contact theorem](../mobile_core_exchange/PROOF.md)

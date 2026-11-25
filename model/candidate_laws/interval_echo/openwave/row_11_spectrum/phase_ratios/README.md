# Can the phase-2/9 line ratios be exact finite formation-energy ratios?

The marked-edge reader sets an angle of **2/9 radians** from an exclusive-event probability. Its nominal intensity ratios are approximately 1 ∶ 206.770316 ∶ 3477.472837. Each exact ratio between distinct lines is transcendental. Consequently, a finite autonomous Hamiltonian with algebraic entries and an algebraic host energy cannot produce any one of those exact pair ratios as a ratio of nonzero host-referenced eigenenergies. Even if the common energy zero and scale can take arbitrary real values, the exact *triple* still cannot be matched. The test concerns the unrounded formula, with the count-to-angle calibration declared by the reader. [Native reader and its observed calibration](../../../../../experiments/interval_echo_lepton_register/README.md)

## The algebraic certificate

Let $K$ be the field of complex algebraic numbers, $z=\exp(2i/9)$, and $\omega=\exp(2\pi i/3)$. The latter is algebraic because $\omega^2+\omega+1=0$. Since $2i/9$ is nonzero algebraic, the Hermite–Lindemann theorem makes $z$ transcendental over $\mathbb Q$ and hence over $K$. A putative polynomial relation over $K$ would have only finitely many coefficients, generating a finite algebraic extension of $\mathbb Q$, so algebraicity would descend by transitivity. [Eberl, formal Hermite–Lindemann–Weierstraß theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)

The source gives $P_k=[\sqrt 2+2\cos(2/9+2\pi k/3)]^2/36$. Define
$$
 F_k(X)=(\omega^k X^2+\sqrt 2 X+\omega^{-k})^2.
$$
Then $F_k(z)=36z^2P_k$. Each $F_k$ has nonzero leading coefficient, so $F_k(z)\ne0$. For distinct $i,j$ the determinant of their constant and leading coefficients is
$$
 \omega^{2(j-i)}-\omega^{-2(j-i)}\ne0.
$$
Hence $F_i-aF_j$ is a nonzero polynomial over $K$ for every $a\in K$. If $P_i/P_j=a$, this polynomial would vanish at $z$. The contradiction proves transcendence for all three unordered pairs and their reciprocals. The Lean lemma proves this rational-function reduction given a transcendental coordinate and nonproportional polynomials. The external transcendence theorem and the cyclotomic coefficient checks are explicit inputs. [Lean reduction](../../../../../../Lean/SourceLaw/ChargedRegisterPole.lean) · [Exact coefficient certificate](certificate.py)

For the calibration-independent comparison, in increasing intensity order the indices are 1, 2, 0. Apply the same argument to $(F_0-F_1)/(F_2-F_1)$. Its constant/leading determinant has magnitude $3|\omega-\omega^2|$, so the exact gap quotient $(P_0-P_1)/(P_2-P_1)$ is also transcendental. The certificate computes every coefficient and all four determinants in $\mathbb Z[T]/(T^8-T^4+1)$, where $T=\exp(2\pi i/24)$, $\omega=T^8$ and $\sqrt2=T^3+T^{-3}$. The regression test independently compares each polynomial with the executable reader formula, including the zero-phase degeneracy. [Exact results](results.json) · [Replay tests](test_phase_ratios.py)

## Application to the search

A finite matrix with entries in $K$ has a monic characteristic polynomial in $K[X]$, so all its eigenvalues lie in $K$. This covers any finite autonomous sector admitting an algebraic matrix representation, including the registered blank formation Hamiltonian with rational contrast costs and exchange entries. Subtracting an algebraic host energy and dividing two nonzero results stays in $K$. A common energy unit cancels even if the physical conversion to keV is unknown. For the whole triple, *any* common affine calibration of three distinct algebraic eigenvalues preserves their algebraic quotient of differences. Multiplying or shifting all three together therefore cannot reproduce the exact register triple. [Registered sector generator](../measure.py)

This proof gives no lower bound on numerical disagreement. Algebraic numbers are dense, so empirical masses and their finite uncertainties still require an actual spectral measurement. The root and coefficient argument also works for other nonzero algebraic angles in radians, including a supplied nonzero finite frequency estimate. A driven measurement can use an algebraic exposure whose propagator contains transcendental phases. A register extension can put such values directly into its mass operator. Infinite-volume limits, resonance parameters, logarithms of Floquet eigenvalues, nonlinear energy conventions and effective fits require separate analysis. For example, the three distinct positive algebraic energies $1,4,(6+\sqrt{39})^2$ have exact Koide quotient $2/3$. Writing their positive square roots as $1,2,u$, the Koide identity reduces to $u^2-12u-3=0$. The obstruction isolates the specific phase prescription. [Guided effective mass operator](../../electron_field/lepton_family/dirac_read/README.md)

From the repository root run the following commands.

```sh
python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.phase_ratios.certificate
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/phase_ratios
```

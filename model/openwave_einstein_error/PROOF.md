# Quantitative finite Einstein tensor reconstruction

The existing Lean theorem `tomographyDecoder_charge` reconstructs the nine
trace-free coordinates of a symmetric four-dimensional tensor from nine
specified null-null charges. This package reuses that design and inverse;
it does not infer a physical metric from the nine labels.

Let eta=diag(-1,1,1,1), c=Tr(eta B)/4 and R=B-c eta. Suppose every supplied
charge differs from k_r^T B k_r by at most epsilon. Because the decoder is
linear, each coordinate error is at most epsilon times the absolute row sum
of its actual inverse matrix. The nine row sums are exactly

```
3/4, 1/2, 1/2, 1/2, 5/4, 2+sqrt(3), 2+sqrt(3), 5/4, (5+sqrt(3))/2.
```

The remaining spatial diagonal has row sum5/4. Therefore every reconstructed
tensor entry has error at most `(2+sqrt(3))*epsilon`. This constant is sharp:
the receipt gives the sign vector attaining the sixth decoder-row absolute
sum. Summing the entry budgets in each matrix row yields the operator norm
bound `(23/4+2sqrt(3))*epsilon`. This is a valid upper bound, not a claim of
optimal operator-norm constant. The executable constructs the design from
the null vectors and independently inverts it in exact Q(sqrt(3)).

`EinsteinResidualBound.lean` proves the general finite weighted-error and
decoder-difference inequalities. The explicit coefficients are checked by
exact arithmetic against the existing finite tomography construction.
No approximate rank test or numerical matrix inversion is used.

## Ward residuals and the correct connectivity

Write B=c eta+R. If the charge residuals are bounded by epsilon, every entry
of R is bounded by C epsilon, C=2+sqrt(3). If the divergence of B is bounded
by w, a forward derivative of spacing a gives the edge estimate

```
|c(x+e_j)-c(x)| <= a*w+8*C*epsilon.
```

There are four derivative terms, each with two endpoint errors. The metric's
diagonal entries have magnitude one. For centered derivatives, the same
argument instead gives

```
|c(x+e_j)-c(x-e_j)| <= 2*a*w+8*C*epsilon.
```

Summation along a path gives the corresponding finite connectivity bound.
The derivative graph matters: c(x)=(-1)^x has zero centered derivative but
is not constant across nearest-neighbour edges. The package retains this
counterexample. A centered-difference proof cannot silently substitute the
forward-step connectivity premise of the existing theorem.

The actual finite native-profile curvature from `openwave_source_ward`
supplies four exact reconstruction examples, including a boundary halo.
These are conditional finite metric records. Their source-to-physical-metric
assignment, physical stress identification and continuum errors remain open.
Null records alone never determine c: any constant multiple of eta has
zero null-null charges. A vacuum/trace input is still required to determine
the cosmological term.

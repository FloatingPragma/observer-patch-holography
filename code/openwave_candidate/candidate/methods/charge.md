# Electric charge on the five-level code

Let the five occupation bits be `(c1,c2,c3,s1,s2)`, with
`n3=c1+c2+c3` and `n2=s1+s2`. The declared hypercharge and electric
generators give

```text
Y = -n3/3 + n2/2,
T3 = (s1-s2)/2,
Q = T3+Y = s1-n3/3,
k = 6Y mod 6 = (3n2-2n3) mod 6.
```

Thus `Q` lies in one-third integer units, and no code weight has
`Q=+1/2` or `Q=-1/2`. The colour representations are
`Lambda^n3(C^3)`, so their singlets occur at `n3=0,3`. Each of the eight
colour-singlet weights has integral `Q`. More generally any tensor-product
weight has `Q_total=sum(s1)-sum(n3)/3`; total colour triality zero makes
this integer. Triality zero is a necessary condition for a colour singlet;
the calculation does not assert its sufficiency.

The kernel element has phase
`n3/3+n2/2+Y=n2` turns, which is integral. The independent combinatorial
count `sum binom(3,n3)binom(2,n2)` within each class agrees with the full
occupation enumeration and yields class multiplicities `(4,6,6,4,6,6)`.

| Equation or check | Implementation |
| --- | --- |
| All 32 occupations, `Y`, `T3`, `Q`, `k` | [charge.py:23](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/charge.py#L23), `enumerate_code` |
| Charge lattice, singlets and kernel | [charge.py:43](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/charge.py#L43), `verify_charge_lattice` |
| Independent exterior-dimension count | [charge.py:59](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/charge.py#L59) |
| All modular halves and charge-balanced pairs | [charge.py:75](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/charge.py#L75), `run` |
| Broken normalization, weak weight and missing outcome | [charge.py:100](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/charge.py#L100) |

## Counterexample to the proposed modular shortcut

`2j=k mod 6` has no solution for odd `k`. This is a statement about two
equal centre classes. An electron-like colour/weak singlet has `Q=-1`, `Y=-1` and
`k=0`, for which `j=0,3` are both solutions. The odd-class argument does
not establish the electron's inability to split into two half-charge
objects. The absence of half-charge weights is instead a consequence of
the specific electric generator and code stated above.

The script also retains a pair of coloured weights with charges `-2/3`
and `-1/3`, whose sum and centre class equal the electron-like weight.
This is an algebraically allowed decomposition, not an executed decay.
Its energetic availability and the isolation of the daughters require
dynamics. Charge arithmetic alone cannot establish their impossibility.

**Scope**: this calculation proves the exact electric-weight lattice of the
declared exterior code and identifies the boundary of the proposed
no-halving argument. It supplies no topological winding integral, charge
isolation energy, physical quark assignment or confinement calculation.

The source is the public RER construction at commit
`3bcacf1276ddb3ead35be2b9adcd63b6643e2d35`, specifically
[`model/defect_particle_classes/PROOF.md`](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/defect_particle_classes/PROOF.md).
The portable implementation evaluates the equations independently, with
exact rational arithmetic and no imported result receipt.


## Parameters, omissions and criteria

The five-level decomposition, exterior occupation rule, hypercharges,
weak generator and electric normalization are architectural inputs to
this script. No measured masses, charge magnitude or fitting target are
used. The RER derivation of their native ancestry is cited rather than
silently reproduced by this small calculation.

OpenWave row 1 concerns quantized particle charge and its resistance to
fractionation. Row 12 additionally requires confined fractional charges.
This artifact addresses the algebraic subquestions; `full_criterion_met`
is false for both rows. An OpenWave maintainer's cell assignment is not
claimed.

Run after installing the candidate package:

```sh
python -m oph_openwave.charge
python -m unittest discover -s tests -p test_charge.py
```

Every run rejects three broken inputs: halving the electric normalization,
changing the weak hypercharge, and deleting a code weight. The exact
scope qualification and the `k=0` counterexample are part of the output.
The [independent adversarial review](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/evidence/charge_spin_review.json) confirms the exact census with a separate Laurent-polynomial calculation and rejects a halved-normalization mutation. It refutes the electron-like singlet shortcut and qualifies the physical row claims; [the audit script](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/audits/charge_spin_review.py) is executable.

# Complete local class conjugation

For integer loads `a,b` in `0..5`, coin `c` in `{0,1}` and classes `k+,k-`
in `Z6`, define

```text
a' = floor((a+b+c)/2),   b' = a+b-a',   j=a-a'
(k+',k-') = (k++j,k--j) mod 6 when accepted, otherwise (k+,k-)
accepted iff j != 0 and
  cost[k++j]+cost[k--j]-cost[k+]-cost[k-] <= 5 j(a-b-j)
cost = (0,5,8,10,8,5).
C(a,b,k+,k-) = (5-a,5-b,-k+,-k-) mod 6 on the class coordinates.
```

All 1,296 local states and both coin records are retained. The complete
record is `(coin,current,accepted)`, including no-operation and rejected
twist branches. Complementing the coin and conjugating the state gives
`F(C(s),1-c)=C(F(s,c))`, with current reversed and acceptance unchanged.
Loads and the two-face class sum are conserved by every attempted update.

For the fair-coin kernel this gives `K(Cs,Ct)=K(s,t)`. The odd observable
`f(s)=1_(k+=1)-1_(k+=5)` consequently has zero expectation after one step
from each conjugate preparation pair. Kernel equivariance supplies the
premise of the existing finite-time induction theorem. The theorem also
requires a conjugation-symmetric preparation and an equivariant schedule;
it is not a statement about an arbitrary biased ensemble.

| Equation or observable | Candidate implementation | Pinned public source |
|---|---|---|
| Balanced repair, current and acceptance | [`step`, conjugation.py:27](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/conjugation.py#L27) | [literal repair and cost law](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/defect_integer_lift/model.py#L112) |
| Complete branch intertwining | [`branch_mismatches`, conjugation.py:45](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/conjugation.py#L45) | [owning branch certificate](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/charge_conjugation_dynamics/symmetry.py#L37) |
| Kernel defect and odd expectation | [`kernel_defect` and `paired_odd_mean`, conjugation.py:68](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/conjugation.py#L68) | [finite-kernel theorem](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/Lean/Dynamics/ChargeConjugation.lean) |

Candidate links refer to the OPH package, without assigning an OpenWave
model number or implying an upstream publication.

## Computed results and controls

All 2,592 complete branches satisfy the identity exactly. The fair-coin
kernel has maximum total-variation equivariance defect zero over all
1,296 states. Every conjugate preparation pair has zero mean for the stated
odd class observable after one attempt. The calculations use integers and
fractions; there is no sampling error, truncation or numerical tolerance.

The asymmetric cost control `(0,5,8,10,8,6)` breaks 20 branch identities.
The uncomplemented-coin control changes actual target states on odd load
sums, in addition to failing the retained coin-record relation. Giving coin
one probability `3/4` produces maximum kernel defect exactly `1/2`, with
the source-state witness `(0,1,0,0)`. A conjugate pair based at `(0,1,0,5)`
then has odd-class mean `-1/4`. These controls demonstrate failures on the
same finite state space as the claimed identity.

## Boundary and omissions

The cost table, acceptance coefficient and classical integer-load dynamics
are declared construction inputs. This local proof does not establish
physical charge conjugation for every quantum source used elsewhere in
OPH. It computes no pair annihilation, released radiation, decay spectrum,
baryon number, baryon-to-photon ratio or cosmic time. Classes zero and three
are self-conjugate under this involution, so the six class labels alone
cannot distinguish every proposed particle/antiparticle interpretation.
No federation run or Lean compilation is performed by this portable module.

**Scope**: exact complete-record conjugation and its fair-kernel consequence
on the declared local integer-load law; this is a static symmetry result
and does not satisfy OpenWave's annihilation-to-radiation criterion.

## Reproduction and review

```sh
python -m oph_openwave.conjugation
python -m unittest discover -s tests -p 'test_conjugation.py'
```

Minimal inspection set: `oph_openwave/conjugation.py`,
`tests/test_conjugation.py`, this note and the common run receipt. The
producing implementation was also compared directly with the pinned RER
`local_step` on all 2,592 states/coins, including each complete record.
The [independent adversarial audit](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/evidence/source_parity_audit.json) confirms every branch against the pinned native law with its own source-parity checks. The [audit script](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/audits/source_parity.py) supplies the computation. Class symmetry is confirmed; physical annihilation is not established.

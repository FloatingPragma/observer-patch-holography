# The entire half-filled dimer charged family

For every repulsive contact strength `U>=0` and hopping `t=c/2>0`, use the
same complete two-site field and the half-filled ground of `K=H-(U/2)N`.
An onsite unitary `d_1=i sigma_axis c_1` puts any of its three orientations
in the spin-independent Hubbard form. Independently assembled symbolic CAR
matrices give the characteristic polynomials of all five number sectors;
no parameter scan, fitted root or visibility cutoff enters the certificate.

Writing `d=sqrt(U^2+16t^2)`, the unique ground energy is `-(U+d)/2`.
Uniqueness follows across every number sector: the vacuum and full state
lie `(U+d)/2>0` above it, the lowest N=1 and N=3 levels lie `d/2-t>=t>0`
above it, and the N=2 triplet lies `(d-U)/2>0` above it since
`d^2-U^2=16t^2>0`. The other N=2 singlets have gaps `(U+d)/2` and `d`.
These inequalities include `U=0`; no noninteracting degeneracy is hidden
in the ground assumption.
Each charge-one sector has exactly two doublets with ground-relative gaps
`g_minus=d/2-t` and `g_plus=d/2+t`. Since `d>=4t`,

```
g_minus >= t > 0,       3*g_minus-g_plus = d-4t >= 0,
1 < g_plus/g_minus <= 3.
```

Equality at the upper end requires `U=0`; there the upper local single-mode
probe residue vanishes. Retaining that dark doublet makes the bound
conservative. The lower endpoint is an infimum as `U/t` diverges.
No tuning within this whole family can reproduce the electron–muon
hierarchy by one common positive linear gap-to-mass conversion. More
basically, two distinct same-charge energies cannot supply three distinct
same-charge lepton generations. The pair states have charge two and cannot
be relabeled as a missing charge-one family.

This is a finite particle-search exclusion, not an exclusion of OPH.
The reference, action and dictionary are explicit hypotheses. Moving the
chemical potential while keeping the half-filled ground can make the
two-gap ratio arbitrarily large; an executable countercontrol retains that
escape from the ratio bound. It still leaves only two doublets per charge.
More precisely, writing `delta=mu-U/2`, this same N=2 state is the unique
ground exactly while `|delta|<g_minus`. Addition gaps shift by `-delta`
and removal gaps by `+delta`. The N=0/N=4 gaps stay positive in this open
interval because `d-U<=4t`; internal N=2 gaps do not change. At either
endpoint another number sector becomes degenerate with the ground, so the
unique-ground contract ends there. Approaching an endpoint from inside
makes one same-charge ratio diverge without increasing its two-level count.
Changing the carrier, action, reference or mass observable requires its own
test. No independent physical comparison or measured mass is added here.

Run `python -m model.dimer_family.family` and
`python -m pytest -q model/dimer_family` from the repository root.

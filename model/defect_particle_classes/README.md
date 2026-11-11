# Defect classes and particle classes

A protected defect on a carrier face holds a flux class `k` in
`Z6 = Z3 × Z2`, a colour-centre charge and a weak-centre parity. Every
Standard Model field carries the same label, `k = 6Y mod 6`. This package
computes the label on all 32 occupation patterns of the five-level code and
on the Standard Model fields, and checks the consequences in exact
arithmetic.

| Class | Colour | Weak | Fields |
| --- | --- | --- | --- |
| 1 | triplet | doublet | `Q = (u_L, d_L)` |
| 5 | antitriplet | doublet | antiquarks of `Q` |
| 4 | triplet | singlet | `u_R`, `d_R` |
| 2 | antitriplet | singlet | their conjugates `u^c`, `d^c` |
| 3 | singlet | doublet | `L`, its conjugate, the Higgs doublet |
| 0 | singlet | singlet | `e_R`, `nu_R`, their conjugates `e^c`, `nu^c`; gluons, `W`, `B` |

The class is the centre class. It does not carry the generation, the
weak isospin `T3` or the rest mass, and a colour singlet and a colour octet
share a class.

The checks cover the class rule and its multiplicities `(4,6,6,4,6,6)`, the
`Z6` kernel, the Aharonov-Bohm phase `exp(2 pi i f k / 6)` of each class
around each flux, the 16 reaction types of the register and the Standard
Model vertices with the same classes, class neutrality on closed carriers,
the grading fixed by the rotation by `2 pi` (every pattern is a fermion,
and defects, links and the Higgs register are bosons), and the plaquette
cost of one flux quantum of each class. The shortest generators give
`8 pi^2 / tr_5 H^dag H = 12/5, 3, 4, 3, 12/5` for classes 1 to 5; the
weak-centre value 4 is the coefficient of the electroweak transmutation
exponent.

From the repository root:

```sh
python -m model.defect_particle_classes --write
python -m model.defect_particle_classes --check
python -m pytest -q -W error model/defect_particle_classes
```

[PROOF.md](PROOF.md) states each check and separates exact computation from
stated premises. [results.json](results.json) holds the tables, the 61 named
checks and the input hashes. Run counts of the level-10 universe run are
imported from the public simulator and labelled as imported.

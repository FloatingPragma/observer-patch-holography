# The integer lift of the defect register

**Result.** The lift is exact and satisfies Stokes on every carrier: the twenty face values of a carrier
sum to zero on every stored cut. It is not a conserved label of the member inside a class. Flat faces of
value `+6` and `-6` occur on 1.04 percent of the flat faces at level 6, where the spin-statistics analysis of
the classes expected none.

The universe runs of the public simulator carry a `Z6` register on the faces
of every carrier and on the plaquettes of the screen. A repair with integer
current `j` proposes the twist `j mod 6`. This package keeps the current
unreduced: every applied twist adds `j` to an integer link on its seam, and
a face takes the signed sum of the links around it. The class of the face is
that integer modulo 6. The integer is the only reading of the register
beyond `Z6`.

| Statement | Label |
| --- | --- |
| The value modulo 6 equals the stored class on every stored cut | exact |
| The twenty face values of a carrier sum to zero, and so do the plaquettes of the screen | exact |
| The law reads the class only, so no step of a run depends on the value | exact |
| A face of absolute value 4 or more does not move with its value | exact |
| 99.7 percent of the defect faces hold a value that is `6Y` of a Standard Model multiplet | counted in the runs |
| Class-4 faces hold `-2` on 92 percent and `+4` on 8 percent | counted in the runs |
| The value is the hypercharge of the member inside a class | reading, not established |

The reading `m = 6Y` separates `u_R` (`+4`) from `d_R` (`-2`). It does not
separate `u_L` from `d_L`, which share `+1`, and `T3` stays outside it. A
face is the flux that such a quantum sees. A pure flux content is a boson:
the word `(+4,-2,-2)` is the flux that a `u d d` content sees. It is not
a neutron.

The runs decide against the reading as an exact law. One composition
carries several integer words, flat faces hold `+6` and `-6`, values of
size 5, 7, 8, 9 and 10 occur, and a value of size 4 or more does not travel
with its defect. The lift records how often the applied current wound
around a face.

The package holds a pure Python reference implementation of law C, the
register and the lift on a small closed federation. It equals the public
engine attempt for attempt at that size. The counts of the runs at levels 4
to 7 are imported and labelled as imported.

From the repository root:

```sh
python -m model.defect_integer_lift --write
python -m model.defect_integer_lift --check
python -m pytest -q -W error model/defect_integer_lift
```

[PROOF.md](PROOF.md) states each result with its label.
[results.json](results.json) holds the tables, the 52 named checks and the
input hashes. [imported_census.json](imported_census.json) holds the run
counts with the simulator commit and the file hashes.

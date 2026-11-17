# One scale and finite-refinement corrections

The new family's one anchor is the selected **limiting** frequency2/sqrt5,
assigned511keV. The common conversion is511sqrt5/2keV per reconstructed-frequency
unit. It is fixed once, before interpreting any of the seven refinements.
No refinement is separately reset to511keV. CALIBRATION.json uses only the
signed poles decoded from public Ramsey records to report both positive bands.

At finite epsilon the single source frequency is slightly above its limiting
value, so the finite pair is correspondingly above1022keV. This is a controlled
refinement correction, not a separated-copy excess: at the same epsilon the
pair frequency is exactly twice the measured single frequency. The limiting
pair is1022keV. At epsilon1/64, the measured single energy is about511.00416keV
and the pair about1022.00832keV. The rigorous single-energy error bound is about
0.00581keV. All seven corrections are retained in the JSON table.

The second positive source band tends to511sqrt(3/2)=625.8446292811keV without
another fitted scale. The complete population records contain this band with
small weight. It is an internal unanchored frequency test, not the assertion
that nature has an electron excitation at this energy. Exact finite outcome
populations and their numerical reader are checked here. A finite-shot
experimental error bar or power calculation is not claimed; the reader's
condition number is retained, and a few noiseless moments cannot be passed
off as a short physical sampling campaign.

The first fixed-word candidate used a different, explicitly declared phase
convention and remains in its original receipt. Its fixed-word alias was a
failure of that interpretation. The refined family is a new candidate and
its common scale is not applied retrospectively to the old logarithmic
E6/E12 packet or to the fixed-word packet.

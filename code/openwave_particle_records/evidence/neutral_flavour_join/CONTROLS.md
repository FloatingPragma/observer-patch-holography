# Controls and identities

All controls act on the same source positions, tag register and final reader.

| Control | Fixed test | Result |
| --- | --- | --- |
| Erase transported tag update | k=1, exactly the primary native coin tape; load/current records unchanged | Maximum final-reader probability change0.96200501569481 |
| Reverse four-cycle route | k=1, every native wait/move and final pointer retained | Maximum reader probability change0.7084371833674205 |
| Hold internal native source | k=1, first stage never moves | Internal cap timeout retained; outgoing payload still reaches receiver |
| Disable initial contact | Same slot program with zero signed source pair | Payload remains at source; terminal outcome is no arrival |
| Read family after each completed cycle | Complete factorized measurement channel | Maximum k=3 return loss0.4724799919636907; monotonicity is not assumed |
| Force holding coins at first stage | Eight attempts | Explicit internal timeout, partial quantum state retained |
| Planted source family | All six input/orbit labels through actual compiled frame and native pointer | All matching conditional probabilities exceed1-10^-9; constant-bin and omitted-orbit readers fail |

The k=3 identity follows from the actual source-derived three-cycle and its
chosen clock boundary. It is not a generic property of every native law:
the erased-tag source has no k=1 appearance, and the holding source cannot
complete a cycle. Complete-outcome normalization, source conservation,
fixed total load and the two-orbit/padding accounting are instrument or
architecture identities. They do not by themselves establish PMNS values.

The first registered protocol placed erased-tag controls at k=3, where both
maps would be identity. That protocol was never executed. Its exact files
and root registration are preserved in attempts/00_registered_unrun. The
corrected k=1 protocol was registered before the source run. The repeated-read
control's ideal channel and implementation envelope are separated in
APPARATUS_BOUND.md.

# Row 4 · Persistence of a prepared bound pair

A single echosahedron on a 20-cell sphere reads whether the two released loads occupy one of its 30 internal seams. Each terminal observation uses a fresh preparation. All 630 physical seams evolve during each flight. The source supplies the initial adjacent pair at ports 156 and 157 and blank seam memories. At native exposure 37 the patch detected 4,028 pairs in 4,096 trials. With the same initial state, detector and exchange dynamics, setting the contrast coupling to zero gave 824 detections. Every quiet read is included. A second calculation bounds the probability of a local pair below by 94% at every exposure from 0 through 37. The species, charge and spin remain unassigned. [Evidence archive](row-04-evidence.zip)

## Reproduce the result

Extract the archive into an empty directory. The archived runtime used Python 3.12.6, NumPy 2.3.5 and SciPy 1.18.1. The fixture and count records are in the ZIP. From the extraction directory run the following commands.

```sh
LC_ALL=C shasum -a 256 -c SHA256SUMS
export OPENBLAS_NUM_THREADS=1
python -B -W error row_04_check.py --board FEDERATION_LEVEL0.npz
python -B -W error row_04_window.py --board FEDERATION_LEVEL0.npz
python -B -W error row_04_check.py --board FEDERATION_LEVEL0.npz --remove-contrast
```

The first simulation compares every terminal outcome and probability with the receipt. The second independently bounds the probability across the entire window. The last command should exit 1 with the message that the freely released pair missed the local persistence gate. The archive preserves the earlier reviewed rehearsal, including its historical PR text. The current proposed wording and filing materials are in the submission folder.

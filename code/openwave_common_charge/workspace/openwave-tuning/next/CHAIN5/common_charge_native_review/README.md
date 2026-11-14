# Common-source charge native opposing review

Independent full literal ROW1 review of the frozen g=1 common charge/topology bridge, plus a separately authored terminal-observable execution annex. Read `REVIEW.json` for the review verdict and `PROOF.md` for exact scope.

`audit.py` verifies all frozen inputs, reconstructs exact Gauss/current and regional cuts, replays every archived native pulse and first-failure effect, recomputes every C/E record and bound, audits native reset on all 78 inputs, and runs the new terminal experiment. It imports neither the producer's executor nor its topology/readout implementations. Source definitions and the prior frozen raw-six-level executor are declared dependencies.

Run `./REPRODUCE.sh` from this directory, or `bash openwave-tuning/next/CHAIN5/common_charge_native_review/REPRODUCE.sh` from the repository root. One numerical thread, small matrices only. No huge full-field word is traversed.

The local `inputs` tree contains the producer packet and new dependencies. Shared frozen inputs remain in `../dirac_kg_opposing_review/inputs`, whose manifest and executor are pinned in `INPUTS.json`. Every individual reused byte is checked before and after execution. `MANIFEST.json` pins all review outputs.

The new annex computes only terminal regional charge/flux and common-clock records by exact transfer through the unchanged physical g=1 generator. Its complete native apparatus has its own export/reset failures. It does not identify its native history with the original huge field program or claim to execute that program. The reviewer authored this annex; its acceptance needs an independent audit.

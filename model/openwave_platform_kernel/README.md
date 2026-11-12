# Native observer waves in OpenWave

A headless Taichi kernel evolves the complete native 180-coordinate
kinetic/pairing source. It executes all three charged-family preparations
on periodic cubes, preserving their complete arrival and failure records.
An independent Fourier solution, a full-carrier derivative comparison and
a changed-sign control test the actual action. The finite integration
bound and numerical roundoff are distinguished.

Use Python 3.12 with OpenWave and RER scientific dependencies installed.
Given a checkout of OpenWave at
ef55175d04737e73c92f33f9f9cf3d17c593da54, run from the RER root:

    python -m model.openwave_platform_kernel --write --openwave-root /path/to/openwave
    python -m model.openwave_platform_kernel --check --openwave-root /path/to/openwave

The command copies the portable engine into the checkout's untracked
openwave/xperiments/oph_observer directory. Existing unrelated content at
that path is rejected. No upstream tracked file or scientific status is
changed. The engine can also run directly with --source source_packet.json
and --output result.json.

The carrier, analyzer and empirical scale are the stated OPH comparison
member. This is a runnable platform integration, without an admission,
physical particle-identification or continuum claim.

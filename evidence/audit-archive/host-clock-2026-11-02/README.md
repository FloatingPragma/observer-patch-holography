# Archived independent host-clock capture

This directory preserves an independent capture found in the original working
checkout before synchronization on 2026-11-02. It is separate from the
[canonical host-clock attachment](../../../model/native_clock/HOST_CLOCK.md).
It adds no physical identification, calibration claim or completion status.

The original checkout HEAD was
`fc2aab8c0e87279c26f45828ece2d6be79913f51`.
The [original report](HOST_CLOCK.md), [raw records](host_clock_records.json)
and [original analysis](host_clock_results.json) are preserved byte for byte.
Their hashes and the backup-manifest hash are in [validation.json](validation.json).
The capture itself does not record an authenticated acquisition date or machine
identity; the directory date is the archival date.

All 24 declared capture-source hashes match that original commit. An independent
replay loaded the exact historical analyzer and checked both complete streams:
264,532 native attempts and 40,000 movements per run, including zero-current
attempts, seeded choices, workload sinks and all block endpoints. It reproduced
the entire original analysis exactly. The first baseline block is 31,185,042 ns;
the canonical capture records a different execution at 32,356,500 ns.
The three held-out baseline blocks meet the original engineering tolerance,
and all four workload-control blocks reject that calibration.

This check verifies recorded-event and analysis consistency. It does not
reproduce or independently authenticate the measured timestamps. No new capture
was performed, and neither canonical timing file was overwritten. Commands in
the preserved report describe its historical source context; they should not be
read as commands to replace the current canonical capture.

The `published-linux-capture/` subdirectory preserves the exact raw tape and
receipt from published commit `412165b0380914681ff96908bd6624d0f6a7607a`.
Its baseline calibration block is 45.033045 ms. The current executable
attachment uses the separately preserved, compatibility-checked historical
Mac capture with a 32.3565 ms block; the original local communication
capture above has a 31.185042 ms block. These are distinct recorded runs.
Neither selecting a retained tape nor replaying its analyzer acquires new
timestamps or establishes a physical energy calibration.

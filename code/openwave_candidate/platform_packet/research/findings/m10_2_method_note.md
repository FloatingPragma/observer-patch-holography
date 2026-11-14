# M10.2 proposed method note: native current subgate

Preparation only; task ID has not been allocated upstream. The independent local
arithmetic control ran and its configuration is already known. No official cell
change is requested from that run.

## Equations and code map

For integer d, abs(d) modulo 2 has the parity of d. Therefore the numerator
of j = [d - (2 coin - 1)(abs(d) mod 2)] / 2 is even, so j is integer.
The updates a' = a-j and b' = b+j conserve total load. With the recorded
oriented current increment E' = E+j, a'+E' = a+E and b'-E' = b-E.
Summing the two coin currents cancels their parity corrections and gives d.

`../scripts/m10_2_native_charge.py` enumerates the entire stated local domain
and tests these equations, then corrupts a record increment as a rejecting control.
The corresponding physical criterion also requires electron/topological charge
identification; the arithmetic does not supply that identification.

## Joined candidate reproduction and receipts

The author bundle's `openwave-tuning/integration/assess.py` dispatches all 31
literal gates from `openwave-tuning/platform_packet/packet.json` against one
candidate. It writes a source-bound assessment with `--output assessment.json`;
`--verify assessment.json` checks those source hashes. The bundle's
`SOURCE_MANIFEST.json` independently records each exported source digest.
The exporter includes the actual joined candidate and geometry, not merely this
independent arithmetic control. Run records stay at the author's pinned commit.

Before submitting this note, replace the proposed task designation with its
allocated roadmap row, pin the author source commit and attach run environment,
command/stdout hashes, elapsed time, all failed controls and independent audit.
The standalone local control and existing integration probes cannot be recast as
fresh pre-registered cell claims. No method note substitutes for the per-row
missing physical observables recorded in the packet.

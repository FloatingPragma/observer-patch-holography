# Prospective trajectory supplement

Motivation: root review correctly distinguished the complete source-column
census from execution of a persistent controller. This supplement adds the latter
without changing any original result. All previous code, freeze and receipts stay.

Code: trajectories.py; independently written reader: trajectory_check.py.
Source is an opaque object exposing only (coin,current). Observer holds the public
initial one-unit pointer, complete records, discovered bond and policy memory.
Source state is never supplied to the controller. Each observer position follows
incidence reconstruction. First nonzero search current sets the bond. Actual
native events displace the unit to a fixed outside target, and the persistent
controller returns it using the same records before work sampling begins.

Inputs: all12 positions, signs-1/+1. Seeds72302000+(sign+1)*100+p. Search cap256;
16 warm-up epochs; native recorded displacement cap32; return cap32;32 separate
work pulses, each preceded by8 hold epochs. All768 work outcomes retained.
Held-law controls use seeds72302500+(sign+1)*100+p and search cap16.
The finite-seed work frequency gate is max absolute deviation<1/10 from the
independently counted population probabilities(1/6,1/2,1/3). This is a disclosed
simulation diagnostic; exact scientific work-law evidence remains the census.

Paired trajectories: independent copies initially(0,+1),(7,-1), with local seeds
72302901+exclusive and72302911+exclusive. After discovery and8 warm-up epochs,
retain2048 joint epochs. Independent mode advances each once; exclusive mode
uses a separate blank native reference coin (seed72302921+exclusive) to select
one. Both complete local tapes, all selection coins and every epoch boundary
are stored. No physical clock identification is added. Exact pair-kernel
statement remains a consequence of the complete executed source columns.

Gates:24 discovery/recovery successes;24 held-law discovery failures retained;
768 work draws retained; stated frequency deviation; all independent-controller
records reconstruct; both pair tapes have2048 retained epochs. Every failure
will remain in output and cause its gate to fail. No conditional normalization
or omission of failed trials. No physical energy, mass or spin gate is added.

Commands: python3 trajectories.py --write; python3 trajectories.py;
python3 trajectory_check.py. No installation, fetch or package execution.

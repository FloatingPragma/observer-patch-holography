"""Replay the frozen distributed clock, local winding, passive tail and equipped controls."""
import argparse
import json

from tools.strict_json import load
from . import clock
from .equipped import race
from .equivalence import passive_outcomes, verify_fourier
from .frontier import exact_tails
from .kernel import require
from .local_clock import replay as local_clock
from .receipts import HERE, archive, receipt, verify_freeze


def replay_distributed(spec):
    records = HERE / "records"
    timing = archive(records / "timing.npz", clock.acquire(spec, "timing"))
    rest = archive(records / "rest.npz", clock.acquire(spec, "rest"))
    forecast = clock.predict(spec, timing, rest)
    targets = archive(records / "finite_forecast.npz", clock.finite_forecasts(spec, rest))
    reference = forecast[f"{max(spec['sizes'])}_out"]["reference"]
    channels = race.forecast(spec["race"], reference)
    receipt(records / "forecast.json", {"clock": forecast, "channels": channels})
    moving = archive(records / "moving.npz", clock.acquire(spec, "moving"))
    clocks = clock.open_pulse(spec, forecast, targets, moving, rest)
    seed = bytes.fromhex(load(records / "release.json")["sender_entropy"])
    require(len(seed) == 32, "invalid released entropy")
    raw, relay = race.acquire(spec["race"], channels, seed)
    archive(records / "channel.npz", raw)
    tails = exact_tails(spec["tail_waits"], spec["tail_speeds"])
    receipt(records / "result.json", {"clock": clocks, "channels": relay, "tails": tails})
    return clocks, relay, tails


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--passive-outcomes", action="store_true",
                        help="emit exact inside/outside and both polarization weights at every frozen threshold")
    args = parser.parse_args()
    verify_freeze()
    spec = load(HERE / "protocol.json")
    local = load(HERE / "local_clock/protocol.json")
    fibers = verify_fourier(spec, local)
    clocks, relay, tails = replay_distributed(spec)
    winding = local_clock.replay()
    outcomes = passive_outcomes(spec["tail_waits"], spec["tail_speeds"], tails)
    print(json.dumps({"free_fourier": fibers,
                      "distributed": {k: row["decision"] for k, row in clocks.items()},
                      "winding": winding,
                      "equipped": {k: [row["accepted"], row["ahead"]] for k, row in relay.items()},
                      "passive": outcomes if args.passive_outcomes else
                          {t: len(row["thresholds"]) for t, row in outcomes.items()}}, indent=2))


if __name__ == "__main__":
    main()

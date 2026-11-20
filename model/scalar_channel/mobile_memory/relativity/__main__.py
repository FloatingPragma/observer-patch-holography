"""Replay the unchanged-law audits and the separate incidence construction."""

import json

from . import clock_read, feedback, incidence, ring_clock, unchanged_clock, unchanged_dispersion


def run():
    return {
        "mobile_memory_reference": {
            "dispersion": unchanged_dispersion.run(),
            "clock": unchanged_clock.run(),
        },
        "incidence_transport_variant": {
            "dispersion_and_motion": incidence.run(),
            "two_read_clock": clock_read.run(),
            "circulating_clock": ring_clock.run(),
            "record_controlled_cage": feedback.run(),
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))

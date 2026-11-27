"""Read an exact material pulse redshift from one published observer survey."""

import argparse
from fractions import Fraction
import json
from pathlib import Path


def pulse_reading(survey, source, emissions, service="strict"):
    first, second = emissions
    if not survey["first_emission"] <= first < second <= survey["last_emission"]:
        raise ValueError(
            "Both emissions must be in the recorded window and in increasing order"
        )
    label = [str(source[0]), int(source[1])]
    i = survey["source_labels"].index(label)
    segments = survey[service]["records"][i]
    received = []
    for e in emissions:
        waits = next(row[1] for row in reversed(segments) if row[0] <= e)
        delay = waits[survey["reference_port"]]
        received.append(None if delay is None else e + delay)
    shift = None
    if None not in received:
        shift = str(Fraction(received[1] - received[0], second - first) - 1)
    return dict(
        source=label,
        receiver=survey["receiver"],
        port=survey["reference_port"],
        emitted=list(emissions),
        received=received,
        redshift=shift,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("survey", type=Path)
    cli.add_argument(
        "--source", nargs=2, type=int, required=True, metavar=("OBSERVER", "PORT")
    )
    cli.add_argument(
        "--emit", nargs=2, type=int, required=True, metavar=("FIRST", "SECOND")
    )
    cli.add_argument("--service", choices=["strict", "transparent"], default="strict")
    args = cli.parse_args()
    try:
        reading = pulse_reading(
            json.loads(args.survey.read_bytes()), args.source, args.emit, args.service
        )
    except ValueError as error:
        cli.error(str(error))
    print(json.dumps(reading))

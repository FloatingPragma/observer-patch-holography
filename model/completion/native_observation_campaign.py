"""Bind finite spatial, detector, intervention and modular results to their domains."""
import argparse
import json
from pathlib import Path

from tools import strict_json
from .finite_observable_campaign import ROOT, validate_family

PATH = Path('model/completion/native_observation_campaign.json')
CAMPAIGN = 'OPH-NATIVE-OBSERVATION-CLOSURE-2026-11-05'
COMPONENTS = {
    'composite_spatial_response': ('FA-153', ['CP-05', 'CP-07']),
    'native_composite_probe': ('FA-154', ['CP-05', 'CP-07']),
    'interacting_observation_process': ('FA-155', ['CP-06']),
    'native_regional_modular': ('FA-156', ['CP-02', 'CP-09']),
}
ISSUES = [1138, 1141, 1142, 1143, 1145]
OWNERS = dict(zip(['CP-02', 'CP-05', 'CP-06', 'CP-07', 'CP-09'], ISSUES))


def validate(packet, root=ROOT):
    return validate_family(packet, root, path=PATH, campaign=CAMPAIGN,
                           components=COMPONENTS, issues=ISSUES, owners=OWNERS,
                           schema='oph.native_observation_campaign.v1',
                           check='model/completion/test_native_observation_campaign.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()

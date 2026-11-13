"""Exact coefficients shared by the declared spatial source consumers."""
import json
from pathlib import Path
from fractions import Fraction as F


def degree_coefficients():
    d=json.loads((Path(__file__).parent/'inputs.json').read_text())
    return F(d['degree_colour_counterterm']),F(d['degree_number_counterterm'])

"""Conditional bare electric tension on the same declared staggered clock.

No interacting inter-charge potential or physical velocity map is proved.
"""
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent


def main():
    inputs=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    electron=F(inputs['electron_anchor_GeV'])
    light=F(inputs['432']['declared_light_midpoint'])
    E0=electron/light;t=F(1,8);speed=2*t
    hbarc=F('0.1973269804')  # Existing row22 unit-conversion decimal.
    spacing=hbarc/(E0*speed)
    rows=[]
    for kappa in [F(1,2),F(1),F(2)]:
        bare=4*kappa/3
        converted=E0*bare/spacing
        rows.append(dict(kappa=str(kappa),bare_source_cost_per_cell=str(bare),
            conditional_bare_tension_GeV_per_fm_exact=str(converted),
            conditional_bare_tension_GeV_per_fm=float(converted)))
    out=dict(schema='oph.obstruction_repair.staggered_conditional_tension.v1',
        energy_unit_GeV_exact=str(E0),free_staggered_slope_exact=str(speed),
        hbarc_GeV_fm_exact=str(hbarc),conditional_spacing_fm_exact=str(spacing),
        conditional_spacing_fm=float(spacing),rows=rows,
        condition='Identify the free staggered slope2t with physical c and use the same432analyzer energy anchor; neither identification is a proved commonDirac spatial transfer.',
        extracted_interacting_string_tension=False,string_breaking_channels_measured=False,
        physical_confinement_or_no_go_claimed=False,
        conclusion='Changing the mass normalization alone does not solve the confinement scale. This is the bare tree electric cost under an explicit hypothetical unit map, not the interacting Cornell potential.')
    (HERE/'CONDITIONAL_TENSION_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['conditional_spacing_fm','rows','extracted_interacting_string_tension']},indent=2))


if __name__=='__main__':main()

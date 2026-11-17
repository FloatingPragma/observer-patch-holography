"""Preserve the actual small electron coefficient on a proposed common clock.

This is a new number-source join, not an identification of the original
positive-Laplacian chain with the full Dirac/Higgs field. A same-clock
declaration is necessary and is not manufactured by converting units.
"""
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent


def run():
    inputs=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    unit=F(inputs['electron_anchor_GeV']);target=F('0.93827208943')
    ratios=json.loads((HERE/'inputs/3_results.json').read_text())['exact_free_mass_ratio_intervals']
    a=list(map(F,ratios['up']));b=list(map(F,ratios['down']))
    records=[]
    for name in ['432','3600']:
        elo,ehi=map(F,inputs[name]['light_interval'])
        assert 0<elo<=ehi
        # For every three-quark pattern at least one occupied cell has odd
        # quark number. Its local spin is half-integer and S_x²>=3/4. The
        # native spin coefficient remains J=1/4. All other terms are positive.
        lo=3*unit+F(3,16)*unit/ehi
        hi=3*unit+F(9,16)*unit/elo
        # The up/down free-gap join changes only the small scalar number term.
        vectorlike_lo=(2*a[0]+b[0])*unit+F(3,16)*unit/ehi
        assert vectorlike_lo>target
        records.append(dict(analyzer=name,electron_source_interval=list(map(str,[elo,ehi])),
            common_energy_unit_GeV_interval=list(map(str,[unit/ehi,unit/elo])),
            common_energy_unit_GeV_interval_decimal=list(map(float,[unit/ehi,unit/elo])),
            equal_light_baryon_interval_GeV=list(map(str,[lo,hi])),
            equal_light_baryon_interval_GeV_decimal=list(map(float,[lo,hi])),
            vectorlike_proton_lower_GeV=str(vectorlike_lo),
            vectorlike_proton_lower_GeV_decimal=float(vectorlike_lo),
            lower_exceeds_comparison=True,
            unchanged_coefficients=dict(t='1/8',U='1',J='1/4',kappa='1'),
            source_change='Replace original normalized m_light=1 by absolute analyzer e_light while retaining interactions in its proposed raw clock. This changes interaction/mass ratios; it is not a pure units conversion.'))
    # A common energy-reference change H->H+cI leaves E_particle-E_vacuum
    # unchanged. Global rescaling H->sH changes e and hadron energy together,
    # so the electron-anchored ratio also remains unchanged.
    q=F(7,13);e=F(1,10000);h=F(1,2);vac=F(3,17)
    assert (h+q)-(vac+q)==h-vac
    assert (3*h-3*vac)/(3*e)==(h-vac)/e
    out=dict(schema='oph.obstruction_repair.raw_mass_join.v1',inputs=inputs,
        coefficient_rows=records,reference_shift_cannot_change_excitation_gap=True,
        common_clock_rescaling_cannot_change_mass_ratio=True,
        division_by_electron_is_not_applied_to_interaction_twice=True,
        full_Dirac_to_positive_Laplacian_interaction_transfer=False,
        physical_proton_mass_solved=False,
        next_required_change='The fixed local colourful-spin cost J=1/4 alone exceeds the proton comparison on both actual raw-electron normalizations. A valid source change must derive a different spin interaction or a dressed vacuum subtraction that changes the excitation gap, with the complete altered source and electron sector retained; an additive reference shift cannot do it.')
    (HERE/'RAW_MASS_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(records,indent=2))


if __name__=='__main__':run()

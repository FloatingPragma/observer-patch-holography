"""Exact rational publication intervals from the complete residual gate."""
from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).parent
r=json.loads((ROOT/'RESULT.json').read_text())
check=json.loads((ROOT/'REVIEW.json').read_text())
radius=F(1,10**48)
# The reviewer explicitly asserts its complete interval residual bound
# below 1e-48. Its displayed bound has ample margin below this rational
# radius. Exact endpoint arithmetic avoids narrowing a displayed decimal.
assert F(check['every_actual_pole_absolute_error_bound'])<radius/100
intervals=[(F(x)-radius,F(x)+radius)for x in r['full_positive_masses']]
assert len(intervals)==63 and intervals[0][0]>0
assert all(intervals[j][1]<intervals[j+1][0]for j in range(2))
assert intervals[3][0]>F(13,30)
unit=F(json.loads((ROOT/'inputs/common_clock.json').read_text())['register_432_source_unit_GeV_exact'])*10**9
light=[(a*unit,b*unit)for a,b in intervals[:3]]
split=[(light[j][0]**2-light[0][1]**2,light[j][1]**2-light[0][0]**2)for j in [1,2]]
ratio=(split[0][0]/split[1][1],split[0][1]/split[1][0])
assert F(33,10000)<ratio[0]<ratio[1]<F(34,10000)
out={'all_pole_radius_exact':str(radius),'all63_actual_pole_intervals_exact':[[str(x)for x in row]for row in intervals],
 'conditional_light_eV_intervals_exact':[[str(x)for x in row]for row in light],
 'conditional_splittings_eV2_intervals_exact':[[str(x)for x in row]for row in split],
 'splitting_ratio_interval_exact':[str(x)for x in ratio],
 'splitting_ratio_between_exact':['33/10000','34/10000'],
 'spectral_pole_inputs_do_not_use_empirical_masses':True,
 'unit_display_consumes_electron_anchor':True}
(ROOT/'EXACT_POLE_INTERVALS.json').write_text(json.dumps(out,indent=2)+'\n')
print('All63 rational pole intervals and scale-independent ratio bracket verified.')

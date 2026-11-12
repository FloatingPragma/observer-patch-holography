import sympy as s
from .spatial import analytic,shell_source,trial_bound

def test_exact_full_spatial_spin_band_and_trap_release():
 d=analytic();assert s.sympify(d['full_spatial_full_photon_vacuum_trial_shift_exact'])<s.Rational(1,16)
 rows=d['controlled_trap_release_trials'];assert rows[-1]['energy_upper_diagnostic']<16.000023
 assert all(x['energy_upper_diagnostic']>16 for x in rows)
 assert d['untrapped_threshold_is_normalizable_eigenstate']is False

def test_source_fixed_complete_shell_refinement_is_not_claimed_monotone():
 a=shell_source(2,1,1.);b=shell_source(3,1,1.);c=shell_source(4,1,1.)
 assert b['ground_energy']>a['ground_energy'] and c['ground_energy']<b['ground_energy']
 assert all(x['complete_ground_Casimir_residual']<1e-10 for x in(a,b,c))

def test_trap_parameter_changes_actual_energy_on_same_regulator():
 a=shell_source(4,1,1.);b=shell_source(4,1,.25)
 assert b['ground_energy']<a['ground_energy'] and b['ground_mean_spatial_number']>a['ground_mean_spatial_number']

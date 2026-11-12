import unittest

import sympy as s

from .flavour import (example_matrices, jarlskog, loop_readout, native_constraints,
                      neutrino_phase_rank, quotient, relation_and_readout, rotation, zero)


class FlavourTests(unittest.TestCase):
    def test_exact_orbit_and_native_constraint_ranks(self):
        q=quotient();n=native_constraints()
        self.assertEqual((q['charged_exact_orbit_rank'],q['Weinberg_exact_orbit_rank'],
                          q['rank_two_Weinberg_exact_orbit_rank']),(41,44,44))
        self.assertEqual((q['charged_invariant_parameters'],q['Weinberg_generic_invariant_parameters'],
                          q['rank_two_Weinberg_invariant_parameters']),(13,22,20))
        self.assertEqual(n['extra_linear_family_constraints_from_native_gauge_Ward_equations'],0)
        self.assertEqual(n['charged_family_embedding_complex_rank'],27)
        self.assertEqual(n['family_generator_native_current_commutators'],540)

    def test_actual_field_redefinitions_preserve_invariants(self):
        yu,yd,_,_,_=example_matrices()
        q=rotation(s.Rational(5,13),s.Rational(7,25),s.Rational(3,5),1)
        u=s.diag(s.I,1,-1);d=s.Matrix([[0,1,0],[0,0,1],[1,0,0]])
        transformed=[q.H*yu*u,q.H*yd*d]
        for before,after in zip([yu,yd],transformed):
            for k in [1,2,3]:
                self.assertEqual(s.simplify(s.trace((before*before.H)**k)-s.trace((after*after.H)**k)),0)
        h1,h2=yu*yu.H,yd*yd.H
        g1,g2=[a*a.H for a in transformed]
        self.assertTrue(zero((g1*g2-g2*g1)-q.H*(h1*h2-h2*h1)*q))
        self.assertFalse(zero(s.diag(2,1,1).H*s.diag(2,1,1)-s.eye(3)))

    def test_CP_requires_coherent_loop_and_retains_failures(self):
        _,_,_,_,v=example_matrices();a=loop_readout(v);b=loop_readout(v.conjugate())
        pa={k:s.sympify(x)for k,x in a['probabilities'].items()};pb={k:s.sympify(x)for k,x in b['probabilities'].items()}
        self.assertEqual(sum(pa.values()),1)
        self.assertEqual(pa['path0_failure'],pb['path0_failure'])
        self.assertEqual(pa['path1_failure'],pb['path1_failure'])
        self.assertEqual(pa['joint_Y_plus'],pb['joint_Y_minus'])
        self.assertNotEqual(pa['joint_Y_plus'],pa['joint_Y_minus'])
        self.assertEqual(a['complete_native_current_record_count'],18)
        self.assertEqual(sum(s.sympify(row['probability'])for row in a['complete_native_current_records']),1)
        for row,conjugate in zip(a['complete_native_current_records'],b['complete_native_current_records']):
            opposite=next(x for x in b['complete_native_current_records']
                          if x['output_family_pair']==row['output_family_pair'] and x['Y']!=row['Y'])
            self.assertEqual(s.sympify(row['probability']),s.sympify(opposite['probability']))
        self.assertEqual(jarlskog(rotation(0,s.Rational(5,13),s.Rational(8,17))),0)
        r=relation_and_readout()
        self.assertTrue(r['native_common_feature_trial_relation_exact'])
        self.assertFalse(r['relation_for_entire_native_gauge_admissible_family'])

    def test_rank_two_neutrino_phase_and_oscillation_distinction(self):
        n=neutrino_phase_rank()
        self.assertEqual(n['rank_two_Majorana_phases'],1)
        self.assertTrue(n['oscillation_generator_U_m_squared_U_dagger_unchanged_by_Majorana_column_phases'])
        v=n['remaining_relative_phase_changes_Lepton_number_violating_invariant']
        self.assertNotEqual(s.sympify(v['before']),s.sympify(v['after']))


if __name__=='__main__':unittest.main()

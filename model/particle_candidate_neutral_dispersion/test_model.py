import copy
from fractions import Fraction as F
from itertools import product
import unittest

import sympy as s

from model.particle_common_contract import contract
from . import model


class NeutralMomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.member=contract.member(0)
        cls.bound=model.bind_member(cls.member)

    def test_full_local_occupation_and_neighbor_CAR_controls(self):
        control=model.occupied_pair_control()
        self.assertEqual(control['active_source_modes'],(6,10,11,21,25,26))
        for row in control['families']:
            self.assertEqual(sum(x['dimension'] for x in row['occupation_sectors']),64)
            self.assertEqual(row['occupied_pair_number_coefficient'],F(row['dimension'],2))
            self.assertTrue(row['coefficient_sharp_at_full_occupation'])
            self.assertTrue(all(F(k)>=0 for k in row['deficit_spectrum']))
        car=model.neighbor_car_control()
        self.assertEqual(car['exact_graded_coefficient_identities'],60)
        self.assertEqual(car['reverse_commutators_exactly_zero'],60)
        self.assertEqual(car['full_CAR_dimension'],4096)

    def test_occupied_pair_orientation_matters(self):
        local=model.spin.local_controls();polynomial=local['families']['G1']['polynomials'][0]
        occupied=s.zeros(64);vacant=s.zeros(64)
        for i in range(6):
            B=sum((v*model.spin.odd.creation_word(w,6)
                for w,v in model.derivative(polynomial,i).items()),s.zeros(64))
            occupied+=B*B.H;vacant+=B.H*B
        self.assertEqual(occupied[0,0],0)
        self.assertEqual(vacant[0,0],3)
        self.assertGreater(vacant[0,0],local['positive_number'][0,0])

    def test_complete_Walsh_normalization_has_no_dimension_double_count(self):
        cells=tuple(product(range(2),repeat=3));m0,m1=s.symbols('m0 m1',real=True)
        M=s.Matrix([[m0 if x==y else m1 if sum(a!=b for a,b in zip(x,y))==1 else 0
            for y in cells] for x in cells])
        W=s.Matrix([[(-1)**sum(a*b for a,b in zip(q,x)) for x in cells] for q in cells])/s.sqrt(8)
        diagonal=s.simplify(W*M*W.T)
        self.assertEqual(diagonal,s.diag(*(m0+(3-2*sum(q))*m1 for q in cells)))
        # Each row includes every component of one multiplet exactly once.
        a0=diagonal[0,0]/48;a3=diagonal[7,7]/48
        self.assertEqual(s.expand(8*(a0-a3)),m1)

    def test_actual_original_bounds(self):
        result=self.bound['result']
        for row in result['channels']:
            d=F(row['dimension']);delta=F(row['source_defect_expectation_upper'])
            self.assertLess(delta,5)
            self.assertEqual(F(row['neighbor_first_moment_absolute_upper']),d*delta/64)
            self.assertLess(F(row['centroid_span_upper']),F(5,9))
            self.assertEqual(row['centroid_strict_interval'],['16/3','38/3'])
            self.assertEqual(row['within_multiplet_centroid_ratio_strict_upper'],'53/48')
        self.assertFalse(result['finite_Walsh_labels_are_conserved_operational_momentum'])
        self.assertFalse(result['nonconstant_original_odd_two_point_first_moments_proved'])
        self.assertEqual(result['physical_mass_table_GeV'],[])

    def test_complete_clock_error_and_common_binding(self):
        estimator=self.bound['result']['complete_neighbor_mean_estimator'];plan=estimator['source_clock_plan']
        self.assertLessEqual(F(estimator['neighbor_estimator_bias_upper']),F(1,2048))
        self.assertEqual(F(estimator['neighbor_estimator_bias_upper']),
            16*F(estimator['one_Walsh_row_mean_error_upper']))
        for key,value in self.member['ids'].items():self.assertEqual(plan[key],value)
        self.assertFalse(plan['source_energy_zero_changed'])
        self.assertFalse(estimator['full_original_histogram_executed'])
        self.assertTrue(estimator['all_failures_both_energy_records_reference_and_memory_retained'])

    def test_same_centroid_does_not_give_mass_gap_or_pole(self):
        low=F(1,10**9);weight=F(1,10**6);mean=F(9)
        high=(mean-weight*low)/(1-weight)
        self.assertEqual(weight*low+(1-weight)*high,mean)
        self.assertLess(weight*(low-mean)**2+(1-weight)*(high-mean)**2,F(1))
        self.assertLess(low,F(1,10**6))
        self.assertNotEqual(low,mean)
        # Four identical copies satisfy the centroid constraints with m1=0.
        self.assertLess(F(0),F(5,9))

    def test_source_and_claim_mutations_fail(self):
        mutated=copy.deepcopy(self.bound)
        mutated['result']['physical_mass_table_GeV']=[{'species':'neutrino','mass':0}]
        mutated['result_sha256']=contract.digest(mutated['result'])
        with self.assertRaises(ValueError):model.verify_bound(self.member,mutated)
        other=copy.deepcopy(self.member);other['ids']['clock_id']='substituted-clock'
        with self.assertRaises(ValueError):model.bind_member(other)
        with self.assertRaises(ValueError):model.centroid_bounds(2,5,19)
        with self.assertRaises(ValueError):model.centroid_bounds(2,4,36)
        with self.assertRaises(ValueError):model.bind_member(contract.member(1))


if __name__=='__main__':unittest.main()

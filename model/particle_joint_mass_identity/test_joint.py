from fractions import Fraction as F
import unittest

from . import joint


class JointTests(unittest.TestCase):
    def test_actual_native_joint_operators_and_antiparticles(self):
        f=joint.fermion_joint();b=joint.boson_joint()
        self.assertLess(max(f['errors'].values()),3e-10)
        self.assertLess(max(b['errors'].values()),2e-12)
        self.assertEqual(sum(r['projector_rank'] for r in f['groups']),180)
        self.assertEqual(sum(r['projector_rank'] for r in b['groups']),56)
        self.assertEqual({r['charge'] for r in f['groups']},{'-1','1','-2/3','2/3','-1/3','1/3','0'})
        self.assertEqual({r['projector_rank'] for r in f['groups']},{1,3})
        self.assertEqual({r['projector_rank'] for r in b['groups']},{1,8})
        self.assertTrue(any(r['candidate_species']==['gluon'] and r['color_C2']=='3' for r in b['groups']))

    def test_full_domain_completion_and_colour_reference(self):
        r=joint.complete_joint_records()
        self.assertEqual(len(r['complete_records']),570)
        self.assertEqual((r['positive_records'],r['zero_records']),(113,457))
        self.assertEqual(sum(F(x['probability_exact']) for x in r['complete_records']),1)
        self.assertEqual(sum(int(x['rank_exact']) for x in r['complete_projectors']),int(r['source_Hilbert_dimension_exact']))
        self.assertEqual(sum(F(x) for x in r['all_fermions_filled_control']['pointer_probabilities']),1)
        self.assertGreater(F(r['negative_control_color_reference_dephasing']['unhalved_complete_output_trace_distance_exact']),1)

    def test_energy_resolution_never_uses_float_degeneracy(self):
        r=joint.complete_joint_records()
        self.assertEqual(r['fine_energy_binning']['class_count'],112)
        self.assertLess(r['finite_energy_binning']['class_count'],112)
        unresolved=[x for x in r['finite_energy_binning']['classes'] if len(x['candidate_species'])>1]
        self.assertEqual(len(unresolved),4)
        self.assertTrue(all(x['candidate_species']==['neutrino_mass_eigenstate_1','neutrino_mass_eigenstate_2','neutrino_mass_eigenstate_3'] for x in unresolved))
        nu=next(x for x in r['joint_groups'] if x['candidate_species']==['neutrino_mass_eigenstate_2'])
        lo,hi=map(F,nu['energy_GeV_interval'])
        self.assertGreater(lo,1)
        self.assertLess(hi-1,F(1,10**21))
        self.assertEqual(float(lo),1.)

    def test_exact_spin_refinement_preserves_original_mass_records(self):
        r=joint.compiled_lepton_joint()
        self.assertEqual(sum(r['joint_projector_ranks']),16)
        self.assertEqual(len(r['complete_records_per_family']),60)
        self.assertEqual(sum(F(x['probability_exact']) for x in r['complete_records_per_family']),1)
        self.assertTrue(r['original36_mass_records_recovered_by_forgetting_spin'])
        self.assertTrue(any(x['source_energy_in_same_m']==2 and x['native_electric_charge']==0 and F(x['probability_exact'])>0 for x in r['complete_records_per_family']))
        self.assertFalse(r['helicity_at_rest_claimed'])


if __name__=='__main__':unittest.main()

"""Independent algebra, boundary, charge and complete-outcome checks."""
import unittest

import numpy as np
from scipy import linalg as la
import sympy as s

from . import source, readout, packet
from model.common_source.chiral_candidate import kinetic


class CompleteSpectrum(unittest.TestCase):
    def test_actual_native_momentum_pair_charge_and_dispersion(self):
        f = source.fermions()
        h, delta = kinetic.momentum_pair_source((1, 2, 2), kernel=f['K'], speed=1.)
        actual = la.eigvalsh(np.block([[h, delta], [delta.conj().T, -h.T]]))
        expected = []
        for channel, masses in f['masses'].items():
            count = 4 if channel == 'neutrino' else (8 if channel == 'charged_lepton' else 24)
            for mass in masses:
                expected.extend([float(s.sqrt(9+mass*mass))]*count)
        np.testing.assert_allclose(actual, np.sort(np.r_[-np.array(expected), expected]), atol=3e-12, rtol=2e-13)
        self.assertLess(max(f['checks']['numerical_errors'].values()), 2e-10)
        self.assertEqual(len(f['modes']), 180)
        self.assertEqual(len(f['targets']), 12)

    def test_full_car_source_and_complete_effect_on_extreme_inputs(self):
        # Both occupied and empty high modes exercise the full 180-mode domain.
        rng = np.random.default_rng(7301)
        masks = [0, 1, (1 << 180)-1, 1 << 179, ((1 << 180)-1) ^ (1 << 179)]
        masks += [sum(int(x) << i for i,x in enumerate(rng.integers(0, 2, size=180))) for _ in range(9)]
        def compose(mask, i, create_i, j, create_j):
            first = readout.car_action(mask, j, create_j)
            if first is None:
                return {}
            second = readout.car_action(first[0], i, create_i)
            return {} if second is None else {second[0]: first[1]*second[1]}
        for mask in masks:
            effect = 0
            for j in range(180):
                effect += sum(readout.car_action(mask, j, create) is not None for create in (True,False))
            self.assertEqual(effect, 180)
            for i,j in ((0,0),(179,179),(0,179),(83,84),(17,131)):
                terms = compose(mask,i,False,j,True)
                for key,value in compose(mask,j,True,i,False).items():
                    terms[key] = terms.get(key,0)+value
                terms = {key:value for key,value in terms.items() if value}
                self.assertEqual(terms, {mask:1} if i==j else {})
        rows = readout.fermion_records()
        self.assertEqual(rows['positive_record_count'], 42)
        self.assertEqual(rows['zero_record_count'], 42)
        self.assertEqual(sum(s.sympify(x['probability']['exact']) for x in rows['records']), 1)

    def test_global_cap_boundary_and_all_source_columns(self):
        states, creators, no_op, error, complete_error = readout.boson_operators()
        self.assertEqual(len(states), 1653)
        self.assertLess(error, 2e-12)
        self.assertLess(complete_error, 2e-12)
        self.assertTrue(all(a[:,57:].nnz == 0 for a in creators))
        np.testing.assert_allclose(no_op[:1], [1/29], atol=2e-16)
        np.testing.assert_array_equal(no_op[1:57], np.zeros(56))
        np.testing.assert_allclose(no_op[57:], np.full(1596,28/29), atol=2e-16)
        # Literal vacuum q|0> spectral weights, with no pole fitted from a table.
        for mode, creator in zip(source.bosons()['modes'], creators, strict=True):
            omega = mode['energy_GeV']['numeric']
            weight = float((creator[:,0].T@creator[:,0]).toarray()[0,0])/(2*omega)
            self.assertAlmostEqual(2*omega*weight, 1., places=14)

    def test_every_boson_pointer_energy_record_from_native_cap_operators(self):
        b, records = source.bosons(), readout.boson_records()
        states, creators, _, _, _ = readout.boson_operators()
        classes = ['W','Z','Higgs','massless']
        categories = [m['species'] if m['species'] in classes[:3] else 'massless' for m in b['modes']]
        sectors = records['complete_energy_projectors']
        index_sector = []
        for state in states:
            counts = {label:sum(categories[j]==label for j in state) for label in classes}
            index_sector.append(next(x['index'] for x in sectors if x['occupation_by_energy_class']==counts))
        for pointer in range(24):
            group = b['groups'][pointer//2]
            actual = np.zeros(15)
            for mode in group['modes']:
                a = creators[mode] if pointer%2 == 0 else creators[mode].T
                column = a[:,0].tocoo()
                for row, value in zip(column.row, column.data, strict=True):
                    actual[index_sector[row]] += abs(value)**2/58
            expected = [x['probability']['numeric'] for x in records['records'] if x['pointer_index']==pointer]
            np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=1e-14)
        self.assertEqual(sum(row['projector_rank'] for row in sectors), 1653)
        self.assertEqual(records['joint_record_count'], 375)
        self.assertEqual(records['zero_record_count'], 362)

    def test_neutrino_and_antitau_complete_analyzers(self):
        record = readout.neutrino_flavour_records()
        self.assertEqual(record['mass_eigenfields'], 3)
        self.assertEqual(record['added_source_modes'], 0)
        self.assertTrue(record['anti_tau_neutrino_analyzer_included'])
        for flavour in ('electron','muon','tau'):
            for orientation in ('neutrino','antineutrino'):
                rows = [x for x in record['records'] if x['flavour']==flavour and x['orientation']==orientation]
                self.assertEqual(s.simplify(sum(s.sympify(x['probability']['exact']) for x in rows)), s.Rational(1,6))
                self.assertEqual({x['mass_eigenstate'] for x in rows}, {1,2,3})
        self.assertEqual(record['total_record_count_including_zero_outside'], 96)

    def test_complete_types_ancestry_and_scope(self):
        result = packet.calculate()
        self.assertEqual(len(result['elementary_particles']), 17)
        self.assertEqual(sum(x['positive_modes_per_momentum'] for x in result['elementary_particles'][:12]), 90)
        self.assertEqual(sum(x['physical_modes_per_momentum'] for x in result['elementary_particles'][12:]), 28)
        self.assertTrue(all(x['physical_source_derived_mass_GeV'] is None for x in result['elementary_particles']))
        self.assertTrue(all(not x['native_source_prediction'] for x in result['elementary_particles']))
        self.assertFalse(result['acceptance_flags']['physical_target_acceptance'])
        self.assertFalse(result['acceptance_flags']['neutrino_low_energy_to_MZ_running_and_matching_proved'])
        self.assertEqual({x['particle'] for x in result['composite_controls']}, {'pi+','pi0','pi-','proton','neutron','rho isospin triplet'})


if __name__ == '__main__':
    unittest.main()

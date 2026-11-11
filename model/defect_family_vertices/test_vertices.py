"""Independent complete scalar spectra and finite scalar-pulse controls."""
from fractions import Fraction
import unittest
import numpy as np
from scipy.linalg import eigh, expm

from .certificate import certificate
from .source import ground_pair, scalar_measure, scalar_population, scalar_spectrum
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory


def numeric(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def midpoint(interval):
    return float(sum(map(Fraction, interval))/2)


class FamilyVertices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = certificate()

    def test_complete_scalar_spectra_independent_generalized_eigensolver(self):
        source = ground_pair()['block']
        gs, hs = numeric(source['gram']), numeric(source['H'])
        energies, vectors = eigh(gs@hs, gs)
        state = vectors[:, 0]
        for sign in (1, -1):
            data = scalar_measure(sign)
            target = data['target']
            g, h = numeric(target['gram']), numeric(target['H'])
            target_energies, target_vectors = eigh(g@h, g)
            amplitudes = target_vectors.conj().T@g@numeric(data['map'])@state/2
            groups = []
            for energy, probability in zip(target_energies, abs(amplitudes)**2):
                if groups and abs(energy-groups[-1][0]) < 1e-8:
                    groups[-1][1] += float(probability)
                else:
                    groups.append([float(energy), float(probability)])
            rows = scalar_spectrum()['records']
            self.assertEqual(len(groups), len(rows))
            for (energy, probability), row in zip(groups, rows):
                self.assertAlmostEqual(energy, midpoint(row['target_energy_interval']), delta=1e-8)
                self.assertAlmostEqual(probability, midpoint(row['probe_probability_interval']), delta=1e-10)
            moment = sum(4*p*(energy-energies[0]) for energy, p in groups)
            expected = midpoint(self.receipt['complete_scalar_readout']['scalar_retarded_first_moment_each_component']['interval'])/2
            self.assertAlmostEqual(moment, expected, delta=1e-9)

    def test_literal_two_mode_full_reference_finite_pulse(self):
        # This independent 60-state control includes every scalar/reference state.
        a = np.zeros((4, 4), complex)
        b = np.zeros((4, 4), complex)
        for mask in range(4):
            if mask & 1:
                a[mask ^ 1, mask] = 1
            if mask & 2:
                b[mask ^ 2, mask] = (-1)**(mask & 1)
        pair = a@b
        fields = [numeric(p) for p in higgs_memory()[0]]
        raising = []
        for k in range(2):
            r = np.zeros((3, 3), complex)
            r[k+1, 0] = 1
            raising.append(r)
        small_x = sum(np.kron(p.conj().T, r)+np.kron(p, r.conj().T)
                      for p, r in zip(fields, raising))
        x = np.kron(np.eye(4), small_x)
        big_b = np.kron(np.kron(pair, np.eye(5)), raising[0].conj().T)
        pointer = 1j*(big_b-big_b.conj().T)/8
        self.assertLessEqual(np.linalg.norm(pointer, 2), 1)
        rng = np.random.default_rng(8419)
        z = rng.normal(size=(20, 20))+1j*rng.normal(size=(20, 20))
        rho = z@z.conj().T
        rho /= np.trace(rho)
        ref = np.diag([1., 0., 0.])
        extended = np.kron(rho, ref)
        channel = np.kron(pair, fields[0].conj().T)
        channel += channel.conj().T
        expected = float(np.trace(rho@channel).real/8)
        step = 1e-4
        values = []
        for theta in (-step, 0, step):
            u = expm(-1j*theta*x)
            evolved = u@extended@u.conj().T
            mean = float(np.trace(evolved@pointer).real)
            probabilities = ((1+mean)/2, (1-mean)/2)
            self.assertAlmostEqual(sum(probabilities), 1)
            self.assertTrue(all(0 <= p <= 1 for p in probabilities))
            values.append(mean)
        self.assertAlmostEqual((values[2]-values[0])/(2*step), expected, delta=1e-9)
        particle = np.kron(np.eye(4), np.diag([0, 1, 1, 0, 0]))
        antiparticle = np.kron(np.eye(4), np.diag([0, 0, 0, 1, 1]))
        vp = float(np.trace(rho@(particle@channel+channel@particle)).real)
        va = float(np.trace(rho@(antiparticle@channel+channel@antiparticle)).real)
        for theta in (0.41, np.pi/(2*np.sqrt(2))):
            u = expm(-1j*theta*x)
            mean = float(np.trace(u@extended@u.conj().T@pointer).real)
            formula = (np.sin(np.sqrt(2)*theta)*vp/np.sqrt(2)
                       +np.cos(np.sqrt(2)*theta)*np.sin(theta)*va)/8
            self.assertAlmostEqual(mean, formula, delta=1e-12)

    def test_complete_no_operation_spectrum(self):
        source = ground_pair()['block']
        gram, h = numeric(source['gram']), numeric(source['H'])
        energies, vectors = eigh(gram@h, gram)
        image = np.sqrt(3)*numeric(ground_pair()['n'])@vectors[:, 0]/2
        probabilities = abs(vectors.conj().T@gram@image)**2
        n = midpoint(self.receipt['complete_scalar_readout']['scalar_population'][1]['interval'])*4
        low = sum(p for energy, p in zip(energies, probabilities) if abs(energy-energies[0]) < 1e-8)
        high = sum(p for energy, p in zip(energies, probabilities) if abs(energy-(1-energies[0])) < 1e-8)
        self.assertAlmostEqual(low, 3*n*n/4, delta=1e-11)
        self.assertAlmostEqual(high, 3*n*(1-n)/4, delta=1e-11)
        self.assertAlmostEqual(sum(probabilities)-low-high, 0, delta=1e-11)

    def test_actual_finite_pulse_complete_probabilities(self):
        channels = self.receipt['actual_scalar_pulse_responses']['channels']
        for row in channels.values():
            probabilities = list(map(midpoint, row['exact_finite_pulse']['all_pointer_probability_intervals']))
            self.assertAlmostEqual(sum(probabilities), 1, delta=1e-15)
            self.assertTrue(all(0 < p < 1 for p in probabilities))
            vp = midpoint(row['scalar_particle_channel_energy']['interval'])
            va = midpoint(row['scalar_antiparticle_channel_energy']['interval'])
            total = midpoint(row['full_channel_interaction_energy']['interval'])
            self.assertAlmostEqual(vp+va, total, delta=1e-12)
            self.assertAlmostEqual(probabilities[0]-probabilities[1], vp/(row['pointer_divisor']*np.sqrt(2)), delta=1e-15)

    def test_contact_and_neutral_complete_counts(self):
        contacts = self.receipt['native_contact_vertices']
        self.assertEqual(contacts['all_fermion_pairs_checked'], 900)
        self.assertEqual(len(contacts['nonzero_components']), 56)
        neutral = self.receipt['neutral_channels']
        self.assertEqual(neutral['minimal_native_neutral_orientation_contact']['antisymmetric_contact_rank'], 28)
        retained = neutral['actual_gauge_singlet_neutral_composite']['complete_records']
        self.assertEqual(len(retained['recorded_settings']), 4)
        self.assertTrue(all(len(row['complete_target_energy_records']) == 9
                            for row in retained['recorded_settings']))
        self.assertEqual(neutral['declared_full32_sterile_spectator_extension']['sterile_pair_contact_checks'], 64)


if __name__ == '__main__':
    unittest.main()

"""Run charged Wilson-loop preparation and absorption of a finite Maxwell quantum."""
import argparse
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
import os
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from scipy.linalg import eigh, expm, qr
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
PIN = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"
# Original field.cochains(n), hashed as edge tuples and sorted (row, col, value) triples.
EXPECTED = {1: "a37ad965f700a7d4b8a9006fb72aaa2ed0eb4c886cae41b65e2174da6c7d8acb",
            4: "b0b3de060afbc99663588061d1c0579f1322c6d4fa619239cfc0d6880cf060a0",
            6: "80b3dcc0c3ee0384e67a76a0acfb520fcf403af51df367a722cce20ccae4a0a8",
            8: "11f233a61cfa7b7aafded6689d36db8ae3812daf3af76105a51ed216ebb3702b"}
G = 1 / 64
CAP = 3
SHOTS = 2**42
PHASES = np.arange(4) * np.pi/2


def check_directed_charge_transfer():
    """A pi/2 addressed Wilson hop implements the archived Peierls phase on a unit singular channel."""
    pauli = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1,-1])]
    beta = np.diag([1,1,-1,-1])
    q_fixture = np.array([[0.,1.,0.],[1.,0.,np.sqrt(2.)],[0.,np.sqrt(2.),0.]])
    phase = expm(1j*G*q_fixture)
    checks = []
    for matrix in pauli:
        alpha = np.block([[np.zeros((2,2)),matrix],[matrix,np.zeros((2,2))]])
        directed = (-(beta+1j*alpha)/2).conj().T
        left, singular, right_h = np.linalg.svd(directed)
        np.testing.assert_allclose(singular,[1,1,0,0],atol=1e-12)
        forward = np.kron(directed,phase)
        hop = np.block([[np.zeros_like(forward),forward.conj().T],[forward,np.zeros_like(forward)]])
        ground = np.eye(3)
        for spin in (0,1):
            at_tail = np.kron(right_h[spin].conj()[:,None],ground)
            at_head = np.kron(left[:,spin,None],phase)
            launched = np.vstack((at_tail,np.zeros_like(at_tail)))
            targeted = -1j*np.vstack((np.zeros_like(at_head),at_head))
            residual = np.max(abs(expm(-1j*np.pi*hop/2) @ launched-targeted))
            assert residual < 2e-12
            checks.append(float(residual))
    return checks


def inc_matrix(entries, shape):
    i, j, values = zip(*entries, strict=True)
    return coo_matrix((values, (i, j)), shape=shape).tocsr()


class Box:
    def __init__(self, n):
        self.n = n
        nodes = list(product(range(n+1), repeat=3))
        index = {x: i for i, x in enumerate(nodes)}
        edges, edge_data, d, c = [], [], [], []
        lookup = {}
        def shifted(x, axis):
            return tuple(v + (i == axis) for i, v in enumerate(x))
        for x in nodes:
            for a in range(3):
                if x[a] == n:
                    continue
                y = shifted(x, a)
                i = len(edges)
                edges.append((x, y))
                lookup[x, y] = (i, 1)
                lookup[y, x] = (i, -1)
                edge_data.append((index[x], index[y], a))
                d.extend([(i, index[x], -1), (i, index[y], 1)])
        faces = 0
        for x in nodes:
            for a, b in ((0, 1), (0, 2), (1, 2)):
                if x[a] == n or x[b] == n:
                    continue
                xy = shifted(shifted(x, a), b)
                loop = [x, shifted(x, a), xy, shifted(x, b), x]
                for t, h in zip(loop[:-1], loop[1:], strict=True):
                    edge, sign = lookup[t, h]
                    c.append((faces, edge, sign))
                faces += 1
        digest = sha256()
        for entries in (edge_data, sorted(d), sorted(c)):
            digest.update(json.dumps(entries, separators=(',', ':')).encode())
        self.hash = digest.hexdigest()
        assert self.hash == EXPECTED[n]
        self.edges, self.lookup = edges, lookup
        self.gradient = inc_matrix(d, (len(edges), len(nodes)))
        curl = inc_matrix(c, (faces, len(edges)))
        assert (curl @ self.gradient).nnz == 0
        stiffness = (curl.T @ curl).toarray()
        val, vec = eigh(stiffness.astype(float))
        physical = val > 1e-9
        assert physical.sum() == len(edges)-len(nodes)+1
        self.freq, self.frame = np.sqrt(val[physical]), vec[:, physical]
        self.q = self.frame / np.sqrt(2*self.freq)
        self.gradient_error = float(np.max(abs(self.gradient.T @ self.q)))
        assert self.gradient_error < 1e-11
        self.minfreq_prediction = 2*np.sqrt(2)*np.sin(np.pi/(2*(n+1)))
        assert abs(self.minfreq_prediction-self.freq[0]) < 2e-12

    def rotate(self, point):
        x, y, z = point
        return self.n-y, x, z

    def paths(self, z):
        m = self.n//2
        path = [(m,m,z), (m,m+1,z), (m,m+1,z+1), (m,m,z+1),
                (m,m-1,z+1), (m,m-1,z), (m,m,z)]
        paths = []
        for _ in range(4):
            paths.append([self.lookup[a,b] for a,b in zip(path[:-1], path[1:], strict=True)])
            path = list(map(self.rotate, path))
        return paths


def occupation_basis(dim):
    words = [word for q in range(CAP+1) for word in combinations_with_replacement(range(dim), q)]
    index = {word: i for i, word in enumerate(words)}
    levels = np.array([len(w) for w in words])
    rows, cols, scales, labels = [], [], [], []
    for i, word in enumerate(words):
        for j in sorted(set(word)):
            reduced = list(word)
            reduced.remove(j)
            rows.append(index[tuple(reduced)])
            cols.append(i)
            scales.append(np.sqrt(word.count(j)))
            labels.append(j)
    return (levels, np.array(rows), np.array(cols), np.array(scales), np.array(labels))


def coordinate_bank(directions):
    """The vacuum and every local pulse close on the span of these creation vectors."""
    edges = list(directions)
    vectors = np.stack([directions[e] for e in edges], axis=1)
    frame, triangular, _ = qr(vectors, mode="economic", pivoting=True)
    size = np.count_nonzero(abs(np.diag(triangular)) > 1e-12)
    frame = frame[:, :size]
    coefficients = frame.conj().T @ vectors
    np.testing.assert_allclose(frame @ coefficients, vectors, atol=3e-13)
    levels, rows, cols, scales, labels = occupation_basis(size)
    q_op = {}
    for i, e in enumerate(edges):
        low = coo_matrix((scales*coefficients[labels,i].conj(), (rows, cols)),
                         shape=(len(levels),len(levels))).tocsr()
        low.eliminate_zeros()
        q_op[e] = low + low.getH()
    return frame, levels, q_op


def drive(state, path, coordinates, reverse=False):
    steps = reversed(path) if reverse else iter(path)
    factor = -1j*G if reverse else 1j*G
    for e, direction in steps:
        state = expm_multiply((factor*direction)*coordinates[e],state,traceA=0)
    return state


def circuit_bank(box, z):
    """Each ordered product is exact on the active Fock space, including the cap."""
    paths = box.paths(z)
    edges = sorted({e for path in paths for e, _ in path})
    frame, levels, q_op = coordinate_bank({e: box.q[e] for e in edges})
    vac = np.zeros(len(levels), complex)
    vac[0] = 1
    rays = [drive(vac, path, q_op) for path in paths]
    rays = np.array(rays)
    np.testing.assert_allclose(np.sum(abs(rays)**2, axis=1), 1, atol=4e-13)
    photons = rays[:, levels == 1] @ frame.T
    test = box.frame @ photons.T
    turned = np.zeros_like(test)
    for i,(tail,head) in enumerate(box.edges):
        j,sign = box.lookup[box.rotate(tail), box.rotate(head)]
        turned[j] = sign*test[i]
    np.testing.assert_allclose(turned,np.roll(test,-1,axis=1),atol=2e-13)
    return paths, levels, rays, photons


def source_record(m, levels, rays, photons):
    weights = np.exp(1j*np.outer(m-np.arange(4), PHASES))/4
    output = weights @ rays
    probs = np.stack([np.sum(abs(output[:,levels==q])**2,axis=1) for q in range(CAP+1)], axis=1)
    assert abs(probs.sum()-1) < 3e-13
    accepted = probs[0,1]
    photon = (weights[0] @ photons)/np.sqrt(accepted)
    assert abs(np.vdot(photon,photon)-1) < 2e-12
    return probs, photon


def receive(field_state, omega, bank, t):
    arriving = np.exp(-1j*t*omega)*field_state
    amps = bank.conj() @ arriving
    analyzer = np.exp(-1j*np.outer(np.arange(4), PHASES))/4
    vacuum = abs(analyzer @ amps)**2
    return np.r_[vacuum, 1-vacuum.sum()]


def acquisition(rng, emitted, received):
    counts = rng.multinomial(SHOTS, emitted.ravel()).reshape(emitted.shape)
    distal = rng.multinomial(counts[0,1], received/received.sum())
    assert distal.sum() == counts[0,1] and counts.sum() == SHOTS
    return {"source_pointer_by_field_number_counts": counts.tolist(),
            "receiver_absorbed_pointer_0_1_2_3_and_missed_counts": distal.tolist()}


def pointer_only(box, source_paths, detector_paths, t, rng):
    """Trace the field out of both independent path-pointer instruments."""
    evolved = np.exp(-1j*t*box.freq)
    vectors = {e: box.q[e]*evolved for path in source_paths for e,_ in path}
    vectors.update({e: box.q[e] for path in detector_paths for e,_ in path})
    _, levels, q_op = coordinate_bank(vectors)
    vacuum = np.zeros(len(levels),complex)
    vacuum[0] = 1
    bank = np.stack([drive(vacuum, path, q_op) for path in source_paths])
    fourier = np.exp(-1j*np.outer(np.arange(4), PHASES))/4
    source_states = np.einsum('sk,k,kf->fs',fourier,np.exp(1j*PHASES),bank)
    # Vacuum as a fifth input directly calibrates the driven receiver's background.
    inputs = np.column_stack((source_states,vacuum))
    images = np.stack([drive(inputs,path,q_op,reverse=True) for path in detector_paths])
    outputs = np.einsum('rk,kfs->srf',fourier,images,optimize=True)
    probabilities = (abs(outputs)**2).sum(axis=2)
    joint_plus, calibration = probabilities[:4], probabilities[4]
    joint_minus = np.roll(joint_plus,-2,axis=0)
    np.testing.assert_allclose([joint_plus.sum(),calibration.sum()],1,atol=2e-12)
    np.testing.assert_allclose(joint_plus.sum(axis=0),joint_minus.sum(axis=0),atol=2e-14)
    plus = joint_plus[0]/joint_plus[0].sum()
    minus = joint_minus[0]/joint_minus[0].sum()
    # Independent global count is virtual in this arm: sum its complete probabilities.
    vac_given_plus = abs(outputs[0,:,0])**2/joint_plus[0].sum()
    vac_given_minus = abs(outputs[2,:,0])**2/joint_minus[0].sum()
    return {"time":float(t),"full_active_dimension":len(levels),
            "local_joint_for_positive_winding":joint_plus.tolist(),
            "local_joint_for_negative_winding":joint_minus.tolist(),
            "joint_counts_plus":rng.multinomial(SHOTS,joint_plus.ravel()/joint_plus.sum()).reshape(4,4).tolist(),
            "joint_counts_minus":rng.multinomial(SHOTS,joint_minus.ravel()/joint_minus.sum()).reshape(4,4).tolist(),
            "receiver_given_local_herald_plus":plus.tolist(),"receiver_given_local_herald_minus":minus.tolist(),
            "receiver_vacuum_input":calibration.tolist(),
            "receiver_inputs_winding_TV_given_local_herald":float(abs(plus-minus).sum()/2),
            "receiver_vacuum_vs_positive_TV_given_local_herald":float(abs(calibration-plus).sum()/2),
            "joint_absorption_and_receiver_pointers_given_local_source_only_plus":vac_given_plus.tolist(),
            "joint_absorption_and_receiver_pointers_given_local_source_only_minus":vac_given_minus.tolist(),
            "unconditioned_receiver_winding_TV":float(abs(joint_plus.sum(axis=0)-joint_minus.sum(axis=0)).sum()/2)}


def clock_receiver(box, source_paths, detector_paths, t, delay, rng):
    """Coherently probe the field at two local times and trace over all field occupations."""
    early_phase = np.exp(-1j*delay*box.freq)
    source_phase = np.exp(-1j*t*box.freq)*early_phase
    vectors = {e:box.q[e]*source_phase for path in source_paths for e,_ in path}
    late_ids = {e for path in detector_paths for e,_ in path}
    vectors.update({e:box.q[e] for e in late_ids})
    vectors.update({~e:box.q[e]*early_phase for e in late_ids})
    _, levels, q_op = coordinate_bank(vectors)
    vacuum = np.zeros(len(levels),complex)
    vacuum[0] = 1
    bank = np.stack([drive(vacuum,path,q_op) for path in source_paths])
    fourier = np.exp(-1j*np.outer(np.arange(4),PHASES))/4
    plus_source = np.einsum("sk,k,kf->fs",fourier,np.exp(1j*PHASES),bank)
    inputs = np.column_stack((plus_source,vacuum))
    out = []
    for early in (True,False):
        paths = [[(~e if early else e,sgn) for e,sgn in path] for path in detector_paths]
        images = np.stack([drive(inputs,path,q_op,reverse=True) for path in paths])
        out.append(np.einsum("rk,kfs->srf",fourier,images,optimize=True))
    tables = np.stack([(abs((out[0]+sign*1j*out[1])/2)**2).sum(axis=2)
                       for sign in (1,-1)],axis=2)
    positive, vacuum_receiver = tables[:4],tables[4]
    negative = np.roll(positive,-2,axis=0)
    np.testing.assert_allclose([positive.sum(),vacuum_receiver.sum()],1,atol=2e-12)
    np.testing.assert_allclose(positive.sum(axis=0),negative.sum(axis=0),atol=2e-14)
    signs = np.array([[0,0],[1,-1],[0,0],[-1,1]])
    shots = 2**49
    distributions, counts, scores, variances = [], [], [], []
    for joint in (positive,negative):
        dist = joint[0]/joint[0].sum()
        distributions.append(dist)
        row = rng.multinomial(shots,joint.ravel()/joint.sum()).reshape(4,4,2)
        counts.append(row.tolist())
        heralded = row[0]/row[0].sum()
        contrast = (signs*heralded).sum()
        scores.append(float(contrast))
        variances.append(float(((signs**2*heralded).sum()-contrast**2)/(row[0].sum()-1)))
    tv = abs(distributions[0]-distributions[1]).sum()/2
    return {"early_time":float(t),"separation_in_field_time":float(delay),"active_fock_dimension":len(levels),
            "joint_pointer_probabilities_positive":positive.tolist(),"joint_pointer_probabilities_negative":negative.tolist(),
            "receiver_vacuum_probabilities":vacuum_receiver.tolist(),
            "receiver_positive_given_local_herald":distributions[0].tolist(),
            "receiver_negative_given_local_herald":distributions[1].tolist(),
            "receiver_winding_TV_given_local_herald":float(tv),
            "unconditioned_receiver_winding_TV":float(abs(positive.sum(axis=0)-negative.sum(axis=0)).sum()/2),
            "signed_contrast_from_probabilities":[float((d*signs).sum()) for d in distributions],
            "source_angle_clock_counts_by_winding":counts, "enrollments_per_winding":shots,
            "signed_contrast_from_counts":scores,
            "counts_positive_minus_negative":scores[0]-scores[1],
            "standard_error_of_count_difference":float(np.sqrt(sum(variances)))}


def measure(n, rng):
    box = Box(n)
    source_paths, levels, rays, source = circuit_bank(box, 0)
    detector_paths, _, _, detector = circuit_bank(box,n-1)
    used_source = {e for path in source_paths for e,_ in path}
    used_detector = {e for path in detector_paths for e,_ in path}
    assert used_source.isdisjoint(used_detector)
    assert {v for e in used_source for v in box.edges[e]}.isdisjoint(
           {v for e in used_detector for v in box.edges[e]})
    times = np.arange(8*n+1)/4
    entries, states = [], []
    for m in (1,-1):
        emitted, photon = source_record(m,levels,rays,source)
        states.append(photon)
        rates = np.array([receive(photon,box.freq,detector,t) for t in times])
        target = m % 4
        wrong = np.max(np.delete(rates[:,:4],target,axis=1))
        assert wrong < 5e-26
        peak = int(rates[:,target].argmax())
        # Do not constrain peak arrival to a preset window or to cross the finite volume.
        assert rates[peak,target] > rates[0,target]
        source_0 = np.zeros_like(rays)
        source_0[:,0] = 1
        weights = np.exp(1j*np.outer(m-np.arange(4), PHASES))/4
        zero_emitted = np.sum(abs((weights @ source_0)[:,levels==1])**2)
        assert zero_emitted == 0
        np.testing.assert_allclose(receive(photon,box.freq,np.zeros_like(detector),n),[0,0,0,0,1])
        weights_low = abs(photon[abs(box.freq-box.freq[0])<1e-9])**2
        occupancy = np.arange(CAP+1)
        conditional_numbers = emitted[0]/emitted[0].sum()
        mean = float(conditional_numbers @ occupancy)
        factorial = occupancy*(occupancy-1)
        g2 = float(conditional_numbers @ factorial/mean**2)
        control = np.array([1-mean+mean**2/3, mean-mean**2/2, 0, mean**2/6])
        assert np.all(control>=0)
        np.testing.assert_allclose([control.sum(), control@occupancy,control@factorial], [1,mean,mean**2])
        assert 0 < g2 < 1
        entries.append({"source_circular_order": m,
                        "occupation_given_local_pointer_zero":conditional_numbers.tolist(),
                        "factorial_number_correlation_given_local_pointer_zero":g2,
                        "mean_matched_number_counter_control_with_correlation_one":control.tolist(), "source_pointer_by_field_number": emitted.tolist(),
                        "lowest_frequency_emitted_fraction":float(weights_low.sum()),
                        "receiver_times":times.tolist(),
                        "receiver_probabilities_given_launch":rates.tolist(),
                        "matching_pointer":target,"maximum_opposite_or_other_pointer":float(wrong),
                        "peak_time":float(times[peak]),"peak_to_zero_wait_ratio":float(rates[peak,target]/rates[0,target]),
                        "joint_peak_probability":float(emitted[0,1]*rates[peak,target]),
                        "joint_side_wait_probability":float(emitted[0,1]*rates[4*n,target]),
                        "zero_source_coupling_one_quantum_probability":float(zero_emitted),
                        "zero_receiver_coupling_absorption_probability":0.,
                        "frozen_flight_matching_pointer_probability":float(rates[0,target]),
                        "independent_count_settings":[{"time":int(t),**acquisition(rng,emitted,rates[4*t])}
                                                       for t in (0,n)]})
    overlap = abs(np.vdot(*states))
    assert overlap < 2e-12
    return {"side":n,"field_modes":len(box.freq),"edge_count":len(box.edges),
            "cochain_digest":box.hash,"lowest_frequency":float(box.freq[0]),
            "lowest_frequency_multiplicity":int(sum(abs(box.freq-box.freq[0]) < 1e-9)),
            "lowest_open_box_frequency_prediction":float(box.minfreq_prediction),
            "maximum_normalized_field_response_to_single_nodal_gradient":box.gradient_error,
            "absolute_overlap_of_opposite_circular_one_quantum_states":float(overlap),
            "source_paths":source_paths, "receiver_paths":detector_paths,"runs":entries,
            "two_local_pointers_field_unobserved":[pointer_only(box,source_paths,detector_paths,t,rng)
                                                   for t in sorted({0,n,entries[0]["peak_time"]})],
            "two_time_local_receiver":[clock_receiver(box,source_paths,detector_paths,t,dt,rng)
                                      for t,dt in sorted({(0,0), (n,np.pi/(2*box.freq[0])),
                                                         (entries[0]["peak_time"],np.pi/(2*box.freq[0]))})]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    unit = Box(1)
    np.testing.assert_allclose(unit.freq**2,[4,4,4,6,6])
    rng = np.random.default_rng(311786)
    out = {"archive_revision": PIN,"field_source_sha256":"d879dc9bf4351f0b2f04dd5fc400cfbb536df75ade09e734778a79e461f1d6bd",
           "one_cube_gradient_curl_hash":unit.hash, "addressed_Wilson_hop_errors": check_directed_charge_transfer(), "total_field_occupation_cap":CAP,
           "charge_coupling":G, "enrollments_per_independent_count_setting":SHOTS,
           "count_seed":311786, "boards":[measure(n,rng) for n in (4,6,8)]}
    with args.output.open("x") as f:
        f.write(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()

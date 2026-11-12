"""Portable headless Taichi solver; imports no RER module or empirical target."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import taichi as ti


def decode(rows):
    return np.array(rows, dtype=float)[..., 0] + 1j * np.array(rows, dtype=float)[..., 1]


def encode(a):
    return np.stack((np.asarray(a).real, np.asarray(a).imag), axis=-1).tolist()


def matrices_from_packet(source):
    matrices = []
    for entries in source['native_H0_and_three_kinetic_sparse_matrices']:
        m = np.zeros((180, 180), complex)
        for i, j, real, imag in entries:
            m[i, j] = real+1j*imag
        matrices.append(m)
    return matrices


@ti.data_oriented
class NativeKernel:
    def __init__(self, side, matrices, spacing=1.):
        self.side = side
        self.spacing = spacing
        self.components = len(matrices[0])
        self.width = max(int(np.count_nonzero(abs(m) > 1e-16, axis=1).max()) for m in matrices)
        shape = (side, side, side, self.components)
        self.state = ti.Vector.field(2, ti.f64, shape=shape)
        self.temp = ti.Vector.field(2, ti.f64, shape=shape)
        self.k = [ti.Vector.field(2, ti.f64, shape=shape) for _ in range(4)]
        self.cols = ti.field(ti.i32, shape=(4, self.components, self.width))
        self.values = ti.Vector.field(2, ti.f64, shape=(4, self.components, self.width))
        cols = np.full((4, self.components, self.width), -1, np.int32)
        values = np.zeros((*cols.shape, 2))
        for q, matrix in enumerate(matrices):
            for i in range(self.components):
                ids = np.flatnonzero(abs(matrix[i]) > 1e-16)
                cols[q, i, :len(ids)] = ids
                values[q, i, :len(ids), 0] = matrix[i, ids].real
                values[q, i, :len(ids), 1] = matrix[i, ids].imag
        self.cols.from_numpy(cols)
        self.values.from_numpy(values)

    @ti.func
    def mul(self, a, b):
        return ti.Vector([a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]])

    @ti.kernel
    def rhs(self, src: ti.template(), dst: ti.template()):
        for x, y, z, i in src:
            out = ti.Vector([0., 0.])
            for k in range(self.width):
                j = self.cols[0, i, k]
                if j >= 0:
                    v = self.mul(self.values[0, i, k], src[x, y, z, j])
                    out += ti.Vector([v[1], -v[0]])
            for axis in ti.static(range(3)):
                for k in range(self.width):
                    j = self.cols[axis+1, i, k]
                    if j >= 0:
                        plus = ti.Vector([x, y, z])
                        minus = ti.Vector([x, y, z])
                        plus[axis] = (plus[axis]+1) % self.side
                        minus[axis] = (minus[axis]+self.side-1) % self.side
                        diff = src[plus[0], plus[1], plus[2], j]-src[minus[0], minus[1], minus[2], j]
                        out -= self.mul(self.values[axis+1, i, k], diff)/(2*self.spacing)
            dst[x, y, z, i] = out

    @ti.kernel
    def intermediate(self, derivative: ti.template(), amount: ti.f64):
        for I in ti.grouped(self.state):
            self.temp[I] = self.state[I]+amount*derivative[I]

    @ti.kernel
    def finish(self, dt: ti.f64):
        for I in ti.grouped(self.state):
            self.state[I] += dt/6*(self.k[0][I]+2*self.k[1][I]+2*self.k[2][I]+self.k[3][I])

    def evolve(self, initial, dt, steps):
        self.state.from_numpy(np.stack((initial.real, initial.imag), axis=-1))
        for _ in range(steps):
            self.rhs(self.state, self.k[0])
            self.intermediate(self.k[0], dt/2)
            self.rhs(self.temp, self.k[1])
            self.intermediate(self.k[1], dt/2)
            self.rhs(self.temp, self.k[2])
            self.intermediate(self.k[2], dt)
            self.rhs(self.temp, self.k[3])
            self.finish(dt)
        a = self.state.to_numpy()
        return a[..., 0]+1j*a[..., 1]

    def derivative(self, initial):
        self.state.from_numpy(np.stack((initial.real, initial.imag), axis=-1))
        self.rhs(self.state, self.temp)
        a = self.temp.to_numpy()
        return a[..., 0]+1j*a[..., 1]


def numpy_derivative(initial, matrices, spacing):
    h = np.einsum('ij,...j->...i', matrices[0], initial)
    for axis in range(3):
        diff = np.roll(initial, -1, axis=axis)-np.roll(initial, 1, axis=axis)
        h += np.einsum('ij,...j->...i', matrices[axis+1], -1j*diff/(2*spacing))
    return -1j*h


def branch_packet(side, frame, matrices, mass, time):
    """Actual positive-band addition and complete source-coordinate detector."""
    amplitude = np.zeros((side, side, side, 180), complex)
    for address in np.ndindex((side,)*3):
        p = np.sin(2*np.pi*np.array(address)/side)
        H = matrices[0]+sum(x*m for x, m in zip(p, matrices[1:]))
        h = frame.conj().T@H@frame
        E = np.sqrt(np.dot(p, p)+mass*mass)
        v = (np.eye(4)+h/E)[:, 0]/2/side**1.5
        v = (np.cos(time*E)*np.eye(4)-1j*np.sin(time*E)*h/E)@v
        amplitude[address] = frame@v
    return np.fft.ifftn(amplitude, axes=(0, 1, 2), norm='ortho')


def run(source, *, sides=(4, 6), steps_list=(10, 20), time=1.):
    ti.init(arch=ti.cpu, default_fp=ti.f64, cpu_max_num_threads=4, offline_cache=False)
    matrices = matrices_from_packet(source)
    frames = [decode(x) for x in source['native_charged_frames']]
    masses = source['native_mass_coefficients']
    assert matrices[0].shape == (180, 180)
    assert all(np.linalg.norm(m-m.conj().T) < 1e-13 for m in matrices)
    bound = source['symbol_norm_bound']
    rows = []
    for side in sides:
        kernel = NativeKernel(side, matrices)
        rng = np.random.default_rng(49173+side)
        probe = rng.normal(size=(side, side, side, 180))+1j*rng.normal(size=(side, side, side, 180))
        probe /= np.linalg.norm(probe)
        derivative_error = float(np.linalg.norm(kernel.derivative(probe)-numpy_derivative(probe, matrices, 1.)))
        assert derivative_error < 3e-14
        mutation = matrices.copy()
        mutation[2] = -mutation[2]
        mutation_difference = float(np.linalg.norm(numpy_derivative(probe, mutation, 1.)-numpy_derivative(probe, matrices, 1.)))
        assert mutation_difference > .1
        cases = []
        for family, frame in enumerate(frames):
            initial = branch_packet(side, frame, matrices, masses[family], 0.)
            exact = branch_packet(side, frame, matrices, masses[family], time)
            assert abs(np.linalg.norm(initial)**2-.5) < 1e-12
            for steps in steps_list:
                dt = time/steps
                assert dt*bound <= 2
                approximate = kernel.evolve(initial, dt, steps)
                error_bound = time**5*bound**5/(120*steps**4)
                error = float(np.linalg.norm(approximate-exact))
                assert error <= error_bound+2e-12
                raw_norm = float(np.linalg.norm(approximate)**2)
                corrected = approximate*np.sqrt(.5/raw_norm)
                exact_records = abs(exact)**2
                approximation_records = abs(corrected)**2
                tv = float(np.sum(abs(exact_records-approximation_records))/2)
                assert tv < 2*error_bound+2e-12
                cases.append({'family':family, 'steps':steps, 'event_parameter':time,
                    'vector_error':error, 'analytic_RK4_operator_error_bound':error_bound,
                    'raw_success_norm':raw_norm, 'normalization_drift':raw_norm-.5,
                    'complete_detector_TV':tv, 'original_preparation_failure':.5,
                    'exact_source_coordinate_records':exact_records.ravel().tolist(),
                    'normalized_numerical_records':approximation_records.ravel().tolist(),
                    'probability_total_including_failure':float(exact_records.sum()+.5)})
        rows.append({'side':side, 'native_CAR_modes':90*side**3, 'Nambu_components':180*side**3,
                     'full_carrier_derivative_error':derivative_error,
                     'kinetic_sign_mutation_difference':mutation_difference, 'cases':cases})
    return {'schema':'oph.openwave_native_kernel_run.v1', 'runs':rows,
            'solver':'Taichi CPU float64 full180-row native stencil, fourth-order Runge-Kutta',
            'physical_interacting_spectrum_or_clock_identified':False,
            'roundoff_scope':'Floating readbacks are diagnostics; the analytic bound covers exact-arithmetic integration truncation only.'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    source = json.loads(args.source.read_text())
    result = run(source)
    result['source_packet_sha256'] = hashlib.sha256(args.source.read_bytes()).hexdigest()
    try:
        import openwave
        result['OpenWave_version'] = openwave.__version__
    except ImportError:
        result['OpenWave_version'] = None
    result['taichi_version'] = list(ti.__version__)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('PASS_FULL_NATIVE_KERNEL', len(result['runs']))


if __name__ == '__main__':
    main()

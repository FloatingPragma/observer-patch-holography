"""Select one immutable rational frame, then prove its exact residual bounds."""
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
import frame
import providers

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'vector_field'))
import cell_complex as cc


def run():
    before = providers.verify()
    if (HERE/'FRAME.json').exists() and (HERE/'CERTIFICATE.json').exists():
        # An immutable choice is consumed, never silently selected again.
        d, stored = frame.data()
        checked = frame.certificate_from_data(d)
        assert checked == stored['certificate']
        print(json.dumps({'reused_frozen_frame':True, 'H_error':float(frame.F(checked['cap1_H_error_including_vacuum']))}), flush=True)
        return
    if (HERE/'FRAME.json').exists():
        selected = json.loads((HERE/'FRAME.json').read_text())
        certificate = frame.certificate_from_data(selected)
        assert before == providers.verify()
        result = {'schema':'oph.quadratic.size3.immutable_frame.v1', 'candidate_id':'native-local12-dirac-noncompact-quadratic-v4-proposed',
            'certificate':certificate, 'frame_sha256':providers.sha(HERE/'FRAME.json'),
            'geometry_sha256':providers.sha(HERE/'GEOMETRY.json'), 'source_sha256':before,
            'checks':{'actual_size3_integer_geometry':True, 'all646positive_modes_selected':True,
                'integer_outward_frame_and_spectrum_bound':True, 'three_harmonic_modes_retained_separately':True,
                'provider_bytes_unchanged':True}, 'official_rows_changed':[]}
        (HERE/'CERTIFICATE.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'resumed_same_frame':True, 'H_error':float(frame.F(certificate['cap1_H_error_including_vacuum']))}), flush=True)
        return
    graph, B, A, periods, winding, labels = cc.build(3)
    integer_A = np.asarray(A, np.int64)
    geometry = {'size':3, 'ports':324, 'edges':list(map(list, graph.seams)),
        'face_entries':[[int(e), int(f), int(integer_A[e, f])] for e, f in zip(*np.nonzero(integer_A))],
        'period_times3_entries':[[i, j, int(3*periods[i, j])] for i in range(3) for j in range(972) if periods[i, j]],
        'winding_columns':[[int(x) for x in winding[:, j]] for j in range(3)],
        'face_labels':[list(x) for x in labels], 'cell_coordinates':list(map(list, graph.cell_coordinates))}
    (HERE/'GEOMETRY.json').write_text(json.dumps(geometry, separators=(',', ':'))+'\n')
    K = integer_A@integer_A.T
    values, vectors = np.linalg.eigh(K.astype(float))
    active = values > 1e-9; values, vectors = values[active], vectors[:, active]
    assert len(values) == 646
    # Sign convention is recorded; degenerate bases are fixed by the saved data.
    for j in range(646):
        i = int(np.argmax(abs(vectors[:, j])))
        if vectors[i, j] < 0: vectors[:, j] *= -1
    selected = {'denominator':frame.DENOMINATOR,
        'frame_integers':np.rint(vectors*frame.DENOMINATOR).astype(np.int64).tolist(),
        'curvature_numerators':np.rint(values*frame.DENOMINATOR).astype(np.int64).tolist(),
        'mode_order':'ascending frozen curvature approximation; mode j is column j-1, physical labels1..646',
        'selection_only':'numpy.eigh proposes integers; exact residual proof is the certificate',
        'numpy_version':np.__version__}
    (HERE/'FRAME.json').write_text(json.dumps(selected, separators=(',', ':'))+'\n')
    certificate = frame.certificate_from_data(selected)
    assert before == providers.verify()
    result = {'schema':'oph.quadratic.size3.immutable_frame.v1', 'candidate_id':'native-local12-dirac-noncompact-quadratic-v4-proposed',
        'certificate':certificate, 'frame_sha256':providers.sha(HERE/'FRAME.json'),
        'geometry_sha256':providers.sha(HERE/'GEOMETRY.json'), 'source_sha256':before,
        'checks':{'actual_size3_integer_geometry':True, 'all646positive_modes_selected':True,
            'integer_outward_frame_and_spectrum_bound':True, 'three_harmonic_modes_retained_separately':True,
            'provider_bytes_unchanged':True}, 'official_rows_changed':[]}
    (HERE/'CERTIFICATE.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks':result['checks'], 'H_error':float(frame.F(certificate['cap1_H_error_including_vacuum'])),
                      'minimum_curvature':float(frame.F(certificate['minimum_curvature']))}), flush=True)


if __name__ == '__main__': run()

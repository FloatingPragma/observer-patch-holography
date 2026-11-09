"""Literal primitive update and protected-sector boundary for the illustration."""
from model.observer import native
from model.defect_protection import protection


def certificate():
    loads = tuple(int(i in (0, 11)) for i in range(12))
    rows = []
    for weight, state in native.tick(native.State(loads)):
        rows.append({'weight': str(weight), 'input_seam_loads': [1, 0],
                     'output_seam_loads': list(state.loads[:2]),
                     'record': list(state.prefix[0]),
                     'next_seam': list(native.selected_seam(state))})
    if len(rows) != 2 or {tuple(r['record']) for r in rows} != {(0, 1), (1, 0)}:
        raise ArithmeticError('literal primitive observer changed')
    checked = 0
    for law in ('C', 'E'):
        for branch in protection.branches(law):
            mapping = dict(protection.lifted_pairs(branch))
            for a, b in mapping.items():
                flip = lambda i: (i+66) % 132
                if mapping.get(flip(a)) != flip(b):
                    raise ArithmeticError('native memory sector changes under repair')
                checked += 1
    zero = (0,)*12
    for seam in ((0, 1), (0, 5)):
        for coin in (0, 1):
            if native.native(zero, seam, coin)[0] != zero:
                raise ArithmeticError('load conservation at zero failed')
    return {'scope': 'literal two-event observer and separate ordered-cover loop',
            'example_initial_occupied_ports': [0, 11],
            'initial_state_supplied': True,
            'primitive_branches': rows,
            'controller_reads': 'phase and previous committed current only',
            'loop_C': protection.exchange_loop('C'),
            'loop_E': protection.exchange_loop('E'),
            'loop_seam_word_is_supplied': True,
            'loop_is_the_literal_two_event_policy': False,
            'charge_commutation_columns_checked': checked,
            'all_complete_branches_preserve_each_Q_eigenspace': True,
            'primitive_creates_nonzero_load_from_zero': False,
            'spatial_topological_defect_creation_derived': False,
            'CAR_excitation_creates_background_winding': False}

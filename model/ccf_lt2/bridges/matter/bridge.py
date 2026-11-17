"""Exact charge sectors, moments and clustering of the complete H_S source."""
from collections import defaultdict
import itertools
import json
from pathlib import Path
import sympy as sy
from ... import core, source, spatial, interactions as inter

def execute():


    def coord(i):
        return (i % 48)//3

    charge = core.native_charge_coordinates()
    weights = {
      'Q3': lambda i: charge[coord(i)],
      'B3': lambda i: 1 if coord(i) in inter.UP+inter.DOWN else
                    -1 if coord(i) in inter.ANTI_UP+inter.ANTI_DOWN else 0,
      'L': lambda i: 1 if coord(i) in inter.LEPTON else
                   -1 if coord(i) in inter.ANTI_LEPTON else 0,
    }

    def defects(poly):
        bad = {k: [] for k in list(weights)+['parity']}
        for word, shift in poly:
            for label, weight in weights.items():
                delta = sum((1 if c else -1)*weight(i) for i, c in word)
                if delta:
                    bad[label].append((word, shift, delta))
            if len(word) % 2:
                bad['parity'].append((word, shift, 1))
        return bad

    def certificate(src):
        poly = src['matter']
        bad = defects(poly)
        if any(bad.values()):
            raise AssertionError('Full-source conservation fails: '+repr({k:v[:1] for k,v in bad.items() if v}))
        assert core.gadjoint(poly) == poly
        shards = source.energy_shards(src)
        for p in shards:
            assert not any(defects(p).values())
            assert core.gadjoint(p) == p
        assert core.gadd(*shards) == poly
        colour_changes = family_changes = 0
        for word, _ in poly:
            colour = [(i,c) for i,c in word if coord(i) in inter.COLOURFUL]
            colour_changes += sum(1 if c else -1 for i,c in colour) != 0
            family_changes += any(sum(1 if c else -1 for i,c in colour if i%3 == f) != 0
                                  for f in range(3))
        assert colour_changes > 0 and family_changes > 0
        return {'monomials':len(poly), 'digest':core.polynomial_digest(poly),
                'conserved':['Q3','B3','L','fermion_parity'],
                'all_shards_conserve':len(shards), 'hermitian':True,
                'colourful_number_changing_words':colour_changes,
                'colourful_family_changing_words':family_changes}

    edge = source.assemble(spatial.edge())
    cert_edge = certificate(edge)
    # A changed input on the same U(1) charge and even-parity alphabet.
    mutation = core.mono(tuple((core.index(0,0,c,0),True) for c in (2,6,8,10)))
    mutation = core.add(mutation, core.adjoint(mutation))
    bad = defects(core.lift(mutation, 1))
    assert not bad['Q3'] and not bad['parity'] and bad['B3'] and bad['L']

    graph = spatial.Complex('disjoint_identical_edges',
      ((0,0,0),(1,0,0),(3,0,0),(4,0,0)), ((0,1),(2,3)), ())
    disjoint = source.assemble(graph)
    cert_disjoint = certificate(disjoint)
    embedded = {}
    for component in range(2):
        translated = { (tuple((i+192*component,c) for i,c in word),
                        (shift[0],0) if component == 0 else (0,shift[0])):coefficient
                       for (word,shift),coefficient in edge['matter'].items() }
        embedded = core.gadd(embedded, translated)
    assert embedded == disjoint['matter']
    assert disjoint['field']['electric'] == edge['field']['electric']
    assert disjoint['field']['magnetic'] == edge['field']['magnetic']
    assert graph.faces == () and graph.edges == ((0,1),(2,3))

    def labels(mask):
        occupied = [i for i in range(96*graph.nvertices) if mask & (1<<i)]
        return {**{name:sum(weight(i) for i in occupied) for name,weight in weights.items()},
                'parity':len(occupied)%2}

    def apply_state(src, state):
        out = defaultdict(lambda:sy.Integer(0))
        for (mask,flux), a in state.items():
            for target,b in source.apply_source_column(src,mask,flux,cap=6).items():
                out[target] += a*b
        return core.clean(out)

    def moment_record(state, out):
        norm = sy.expand(sum(sy.conjugate(a)*a for a in state.values()))
        energy = sy.expand(sum(sy.conjugate(a)*out.get(key,0) for key,a in state.items())/norm)
        second = sy.expand(sum(sy.conjugate(a)*a for a in out.values())/norm)
        variance = sy.expand(second-energy**2)
        assert norm.is_positive is True and energy.is_real is True
        assert variance.is_nonnegative is True
        original = {tuple(sorted(labels(key[0]).items())) for key in state}
        assert len(original) == 1
        assert all(tuple(sorted(labels(key[0]).items())) in original for key in out)
        return {'seed_norm_squared':str(norm), 'mean_energy':str(energy),
                'second_moment':str(second),'variance':str(variance),
                'output_count':len(out), 'charge_labels':dict(next(iter(original))),
                'all_output_labels_preserved':True}

    def dump_column(state):
        return [{'occupation':str(mask),'flux':list(flux),'coefficient':str(c)}
                for (mask,flux),c in sorted(state.items())]

    # The epsilon seed is a colour singlet; its spin/energy are not named a neutron.
    seed = defaultdict(lambda:sy.Integer(0))
    for perm in itertools.permutations(range(3)):
        sign = (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        word = ((core.index(0,0,inter.UP[perm[0]],0),True),
                (core.index(0,0,inter.DOWN[perm[1]],0),True),
                (core.index(0,1,inter.DOWN[perm[2]],0),True))
        applied = core.apply_word(word, 0)
        if applied is not None:
            seed[(applied[0],(0,0))] += sign*applied[1]
    seed = core.clean(seed)
    packet=json.loads((Path(__file__).resolve().parents[2]/'representation_snapshot.json').read_text())
    inter._colour_gate()
    generators=[sy.Matrix(row) for row in packet['colour_generators']]
    for generator in generators[:8]:
        colour_out = defaultdict(lambda:sy.Integer(0))
        for s in range(2):
            for f in range(3):
                for i in range(15):
                    for j in range(15):
                        c = generator[i,j]
                        if c == 0:
                            continue
                        word = ((core.index(0,s,i,f),True),(core.index(0,s,j,f),False))
                        for (mask,flux),a in seed.items():
                            acted = core.apply_word(word,mask)
                            if acted is not None:
                                colour_out[(acted[0],flux)] += c*a*acted[1]
        assert not core.clean(colour_out), 'epsilon seed is not a native colour singlet'

    states = {'native_empty':{(0,(0,0)):sy.Integer(1)},
      'neutral_first_component':{(1<<core.index(0,0,10,0),(0,0)):sy.Integer(1)},
      'neutral_second_component':{(1<<core.index(2,0,10,0),(0,0)):sy.Integer(1)},
      'epsilon_udd_preparation':seed}
    outputs = {name:apply_state(disjoint,state) for name,state in states.items()}
    moments = {name:moment_record(states[name],out) for name,out in outputs.items()}
    assert moments['neutral_first_component']['mean_energy'] == moments['neutral_second_component']['mean_energy']
    delta = sy.expand(sy.sympify(moments['neutral_first_component']['mean_energy'])-
                      sy.sympify(moments['native_empty']['mean_energy']))
    assert delta != 0

    empty_edge = source.apply_source_column(edge,0,(0,),cap=6)
    empty_embedded = defaultdict(lambda:sy.Integer(0))
    for component in range(2):
        for (mask,flux),coefficient in empty_edge.items():
            key = (mask << (192*component),(flux[0],0) if component == 0 else (0,flux[0]))
            empty_embedded[key] += coefficient
    assert core.clean(empty_embedded) == outputs['native_empty']
    witnesses = []
    for component in range(2):
        entries = []
        for (mask,flux),coefficient in sorted(outputs['native_empty'].items()):
            if not mask:
                continue
            outside = mask >> 192 if component == 0 else mask & ((1<<192)-1)
            if not outside and flux[1-component] == 0:
                entries.append({'occupation':str(mask),'flux':list(flux),'coefficient':str(coefficient)})
        assert entries
        witnesses.append({'component':component,'nonzero_local_outputs':len(entries),'first_entry':entries[0]})

    result = {'all_passed':True,'scope':'Full finite H_S; no physical mass, binding or decay claim',
      'edge':cert_edge,'disconnected':cert_disjoint,
      'source_cluster_additivity':True,'scope_is_matter_field_source_only':True,
      'control_Q_neutral_even_B_L_violating':{k:len(v) for k,v in bad.items()},
      'native_colour_singlet_checks':8,'moments':moments,
      'neutral_preparation_energy_contrast':str(delta),
      'neutral_contrast_is_positive':bool(delta.is_positive),
      'complete_disconnected_empty_column_identity':True,
      'nontrivial_component_witnesses':witnesses}
    columns = {'scope':result['scope'],'hard_wall_cap':6,
      'graph':{'positions':graph.positions,'edges':graph.edges,'faces':graph.faces},
      'states':{name:{'input':dump_column(state),'full_HS_action':dump_column(outputs[name])}
                for name,state in states.items()},
      'one_edge_empty_action':dump_column(empty_edge)}
    return result, columns

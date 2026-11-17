"""Source-compiled, locally gauge-covariant staggered matter transport.

The actual native currents are lifted on all32 auxiliary exterior levels.
Their15-dimensional invariant sector acts on separate physical CAR modes.
Conjugating the native exchange echo gives linked hopping on those modes;
every local four-column factor and every fermionic routing sign is retained.
The complete Fock word is specified by these factors, not by truncation to
one particle. Its numerical error bound uses full local operators.

Links, cubic geometry, inverse coherent controls, spin/exterior services and
proper-code loading remain declared apparatus. Links are background fields,
not dynamical gauge registers. The free limit retains two Dirac tastes per
internal component: this is not a chiral Standard Model or a mass derivation.
The nearest standard construction is Susskind, Phys.Rev.D16,3031(1977),
https://doi.org/10.1103/PhysRevD.16.3031. LocalGaugeJetAction supplies the
separate infinitesimal algebra; the new result is its finite source-linked
CAR realization on this same native internal representation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_matter_attachment as matter
from . import native_staggered_field as field
from .chiral_candidate import gaussian_source as gaussian
from ..observer.source_scalar import compiler


ROOT = Path(__file__).resolve().parents[2]
LOCAL_ADMISSION = F(1, 10**10)
CURRENT_ADMISSION = F(1, 10**8)
REPLAY_ADMISSION = 2e-9


def _freeze(value):
    value = np.asarray(value, complex)
    return np.frombuffer(value.tobytes(), complex).reshape(value.shape)


def _fraction(value):
    if isinstance(value, bool) or type(value) not in (int, F):
        raise ValueError('an original exact Fraction or integer is required')
    return F(value)


def _path(path):
    if type(path) not in (list, tuple) or len(path) > 64:
        raise ValueError('a finite chronological current path is required')
    result = []
    for entry in path:
        if type(entry) not in (list, tuple) or len(entry) != 2:
            raise ValueError('each path entry is a current index and exact angle')
        index, angle = entry
        if type(index) is not int or not 0 <= index < 12:
            raise ValueError('one of the twelve canonical native currents required')
        result.append((index, _fraction(angle)))
    return tuple(result)


def inverse_path(path):
    return tuple((index, -angle) for index, angle in reversed(_path(path)))


@lru_cache(maxsize=1)
def exact_currents():
    """Physical one-body generators computed from the literal native lift."""
    r = matter.representation(); e = r['encoder']
    result = tuple(matter.clean(e.T*matter.exterior_derivative(matter.spin_lift(g))*e)
                   for g in matter.canonical_currents())
    if any(not matter.zero(g+g.H) for g in result):
        raise ArithmeticError('native physical-current anti-Hermiticity failed')
    return result


@lru_cache(maxsize=256)
def _current_path(path):
    canonical = matter.canonical_currents(); word = []
    for index, angle in path:
        word += matter.onebody_word(matter.spin_lift(canonical[index]), float(angle))
    full = matter.execute_word(word, 5)
    ideal = np.eye(32, dtype=complex)
    for index, angle in path:
        g = matter.exterior_derivative(matter.spin_lift(canonical[index]))
        ideal = expm(float(angle)*np.asarray(g, complex))@ideal
    e = np.asarray(matter.representation()['encoder'], complex)
    leakage = np.linalg.norm((np.eye(32)-e@e.T)@full@e, 2)
    error = float(np.linalg.norm(full-ideal, 2))
    if max(error, leakage) > float(CURRENT_ADMISSION):
        raise ArithmeticError('full auxiliary current execution exceeds its admission')
    return {'full32': _freeze(full), 'matter15': _freeze(e.T@full@e),
            'ideal_full32': _freeze(ideal), 'ideal_matter15': _freeze(e.T@ideal@e),
            'word': tuple(word), 'path': path, 'resources': matter.resources(word),
            'source_error_upper_conditional': str(CURRENT_ADMISSION),
            'operator_error_diagnostic': error,
            'selected_to_auxiliary_complement_error_diagnostic': float(leakage),
            'all_auxiliary_levels_retained': 32}


def current_path(path):
    return _current_path(_path(path))


def _lift_two(one):
    """Exact number-conserving Gaussian lift, including vacuum/double fill."""
    out = np.zeros((4, 4), complex)
    out[0, 0] = 1; out[3, 3] = np.linalg.det(one)
    out[np.ix_([1, 2], [1, 2])] = one
    return out


@lru_cache(maxsize=128)
def local_factor(kind, real=F(0), imag=F(0), duration=F(0)):
    """Execute all four Fock columns through the actual native current map."""
    real, imag, duration = map(_fraction, (real, imag, duration))
    z = complex(float(real), float(imag))
    if kind == 'hopping':
        if real and imag:
            raise ValueError('canonical source hoppings have a real or imaginary coefficient')
        magnitude = abs(real) + abs(imag)
        signed = z*(1 if duration >= 0 else -1)
        phase = math.atan2((-signed).imag, (-signed).real) if magnitude else 0.
        echo_time = abs(duration)*magnitude/field.native_parameters()['exchange']
        word = field._phase_word(-phase)+field.echo_word(echo_time)+field._phase_word(phase)
        one = expm(-1j*float(duration)*np.array([[0,z],[z.conjugate(),0]], complex))
    elif kind == 'onsite':
        if imag: raise ValueError('real onsite coefficient required')
        word = field._phase_word(-float(duration*real))
        one = np.diag([np.exp(-1j*float(duration*real)), 1])
    elif kind == 'fswap':
        if real or imag or duration: raise ValueError('fermion swap has no free parameters')
        word = gaussian.local_word('fswap')
        one = np.array([[0,1],[1,0]], complex)
    else:
        raise ValueError('hopping, onsite or fswap source factor required')
    actual = matter.execute_word(word, 2)
    target = _lift_two(one)
    error = float(np.linalg.norm(actual-target, 2))
    if error > float(LOCAL_ADMISSION):
        raise ArithmeticError('complete local Fock source factor exceeds its admission')
    return {'kind': kind, 'parameters': tuple(map(str, (real,imag,duration))),
            'word': tuple(word), 'actual4': _freeze(actual), 'target4': _freeze(target),
            'one_body': _freeze(one), 'operator_error_diagnostic': error,
            'resources': matter.resources(word)}


class FockWord:
    """Complete many-mode word, stored as four-column local source factors."""
    def __init__(self, modes):
        if type(modes) is not int or not 2 <= modes <= 512:
            raise ValueError('two through512 physical fermion modes required')
        self.modes = modes; self.catalogue = []; self.order = []; self._keys = {}

    def append(self, factor, sites):
        if (type(sites) is not tuple or len(sites) != 2 or sites[0] == sites[1]
                or any(type(i) is not int or not 0 <= i < self.modes for i in sites)):
            raise ValueError('two distinct physical mode addresses required')
        if factor['kind'] != 'onsite' and sites[1] != sites[0]+1:
            raise ValueError('off-diagonal physical factors require adjacent CAR routing')
        # Literal occupation wire0 of each nine-bit native carrier; every
        # other decoded bit remains a spectator. This validates wiring only,
        # not a complete many-carrier preparation/C/E/service history.
        field.decoded_occupation_word(compiler.remap(factor['word'],sites),modes=self.modes)
        key = factor['kind'], factor['parameters']
        if key not in self._keys:
            self._keys[key] = len(self.catalogue); self.catalogue.append(factor)
        self.order.append((self._keys[key], sites))

    def term(self, kind, i, j, coefficient, duration):
        coefficient = s.sympify(coefficient)
        if coefficient.free_symbols or coefficient.has(s.Float):
            raise ValueError('an exact source coefficient required')
        real, imag = s.re(coefficient), s.im(coefficient)
        if not real.is_Rational or not imag.is_Rational:
            raise ValueError('canonical coefficients must lie in Q+iQ')
        factor = local_factor(kind, F(real), F(imag), _fraction(duration))
        if kind == 'onsite':
            self.append(factor, (i, (i+1) % self.modes)); return
        if not 0 <= i < j < self.modes: raise ValueError('ordered CAR endpoints required')
        route = [(k-1,k) for k in range(j,i+1,-1)]
        swap = local_factor('fswap')
        for pair in route: self.append(swap,pair)
        self.append(factor,(i,i+1))
        for pair in reversed(route): self.append(swap,pair)

    def current(self, path, offset=0, families=3):
        if type(families) is not int or not 1 <= families <= 3 or offset+15*families > self.modes:
            raise ValueError('complete fifteen-component family copies required')
        for index, duration in _path(path):
            h = s.I*exact_currents()[index]
            pieces=[]
            for i in range(15):
                if h[i,i]:
                    pieces.append(h[i,i]*matter.unit(15,i,i))
                    for f in range(families): self.term('onsite',offset+families*i+f,0,h[i,i],duration)
                for j in range(i+1,15):
                    if h[i,j]:
                        pieces.append(h[i,j]*matter.unit(15,i,j)+h[j,i]*matter.unit(15,j,i))
                        for f in range(families):
                            self.term('hopping',offset+families*i+f,offset+families*j+f,h[i,j],duration)
            if any(not matter.zero(a*b-b*a) for a,b in itertools.combinations(pieces,2)):
                raise ArithmeticError('canonical current terms require an additional splitting law')

    def ideal_one_body(self):
        out=np.eye(self.modes,dtype=complex)
        for block,sites in self.order:
            out[list(sites)]=self.catalogue[block]['one_body']@out[list(sites)]
        return out

    def complete_matrix(self):
        if self.modes > 8: raise ValueError('dense audit limited to eight physical modes')
        out=np.eye(1<<self.modes,dtype=complex)
        for block,sites in self.order:
            out=field.attachment._little_apply(self.catalogue[block]['actual4'],sites,self.modes,out)
        return out

    def error_bound(self):
        """All-Fock telescoping, conditional on every full local admission.

        Each embedding preserves its complete local operator norm. With m
        factors of error<=epsilon, prod(1+epsilon)-1 <=m epsilon/(1-m epsilon).
        No particle-number cutoff or one-body-to-Fock inference is used.
        """
        total=len(self.order)*LOCAL_ADMISSION
        if total >= 1: raise ArithmeticError('full-Fock source budget exhausted')
        return total/(1-total)

    def record(self):
        return {'physical_CAR_modes':self.modes,'complete_Fock_dimension':str(1<<self.modes),
                'complete_local_Fock_columns':4,'factor_order':self.order,
                'local_operator_admission':str(LOCAL_ADMISSION),
                'unrestricted_Fock_operator_error_upper_conditional':str(self.error_bound()),
                'source_error_is_outward_or_hardware_certificate':False,
                'all_occupation_sectors_and_JW_routing_retained':True,
                'native_nine_bit_carrier_occupation_wiring_checked':True,
                'source_local_factors':[{'kind':r['kind'],'parameters':r['parameters'],
                    'word':r['word'],'resources':r['resources']} for r in self.catalogue]}


def compile_bond(path, duration=F(1,8), *, families=3, eta=1):
    """Exact linked kinetic bond: h_xy=-i eta J U, with all CAR sectors."""
    path=_path(path);duration=_fraction(duration)
    if type(families) is not int or not 1<=families<=3 or type(eta) is not int or eta not in (-1,1):
        raise ValueError('one through three family copies and staggered sign required')
    d=15*families;word=FockWord(2*d)
    word.current(inverse_path(path),families=families)
    j=s.Rational(field.native_parameters()['exchange'])
    for i in range(d):word.term('hopping',i,d+i,-s.I*eta*j,duration)
    word.current(path,families=families)
    return word


def background(shape=(2,2,2), *, families=3, spacing=F(1)):
    """Fixed nonflat native-current paths and independent local frames."""
    coordinates=field._coordinates(shape);spacing=_fraction(spacing)
    if spacing<=0 or type(families) is not int or not 1<=families<=3:
        raise ValueError('positive spacing and one through three families required')
    indices={x:i for i,x in enumerate(coordinates)};d=15*families
    paths=(((0,F(1,7)),(8,F(1,11))),((1,F(1,5)),(9,F(1,13))),((6,F(1,9)),(11,F(1,17))))
    links=[];gauges=[];gauge_paths=[]
    for i,x in enumerate(coordinates):
        gp=((2*(i%3),F(1,i+9)),(8+i%3,F(1,i+13)),(11,F(i+1,101)))
        gauge_paths.append(gp);gauges.append(np.kron(current_path(gp)['ideal_matter15'],np.eye(families)))
        for axis in range(3):
            y=list(x);y[axis]+=1
            if y[axis]>=shape[axis]:continue
            target=indices[tuple(y)]
            links.append({'left':i,'right':target,'axis':axis,'eta':field._eta(x,axis),
                          'path':paths[axis],'matrix':np.kron(current_path(paths[axis])['ideal_matter15'],np.eye(families))})
    return {'shape':shape,'coordinates':coordinates,'families':families,'spacing':spacing,
            'links':links,'gauges':tuple(gauges),'gauge_paths':tuple(gauge_paths),
            'internal_dimension':d}


def lattice_hamiltonian(data, *, transformed=False):
    d=data['internal_dimension'];n=len(data['coordinates']);h=np.zeros((n*d,n*d),complex)
    j=float(field.native_parameters()['exchange']/data['spacing'])
    for row in data['links']:
        left,right=row['left'],row['right'];u=row['matrix']
        if transformed:u=data['gauges'][left]@u@data['gauges'][right].conj().T
        block=-1j*j*row['eta']*u
        h[d*left:d*(left+1),d*right:d*(right+1)]+=block
        h[d*right:d*(right+1),d*left:d*(left+1)]+=block.conj().T
    return h


def build(*, include_continuum=True):
    """Finite same-current transport witness; no pole or mass fit is run."""
    pins=source_bindings();data=background();h=lattice_hamiltonian(data)
    changed=lattice_hamiltonian(data,transformed=True)
    d=data['internal_dimension'];n=len(data['coordinates'])
    whole=np.zeros_like(h)
    for i,g in enumerate(data['gauges']):whole[d*i:d*(i+1),d*i:d*(i+1)]=g
    covariance=float(np.linalg.norm(changed-whole@h@whole.conj().T,2))
    lookup={(r['left'],r['axis']):r for r in data['links']}
    a=lookup[0,0];b=lookup[a['right'],1];c=lookup[0,1];e=lookup[c['right'],0]
    plaquette=a['matrix']@b['matrix']@e['matrix'].conj().T@c['matrix'].conj().T
    nonflat=float(np.linalg.norm(plaquette-np.eye(d),2))
    if covariance>1e-10 or nonflat<1e-3:
        raise ArithmeticError('independent local covariance or nonflat link control failed')
    path=a['path'];word=compile_bond(path);u=np.kron(current_path(path)['ideal_matter15'],np.eye(3))
    J=float(field.native_parameters()['exchange']);bond=np.zeros((2*d,2*d),complex)
    bond[:d,d:]=-1j*J*u;bond[d:,:d]=bond[:d,d:].conj().T
    ideal_word=word.ideal_one_body();bond_error=float(np.linalg.norm(ideal_word-expm(-.125j*bond),2))
    if bond_error>1e-10:
        raise ArithmeticError('native echo/current factors do not realize their linked kinetic bond')
    # The representation is the actual source-derived M15; no species table
    # selects or drops a kinetic component. Three neutral family copies are
    # explicit, and the eight-corner doubling is retained on every copy.
    clifford=field.clifford_matrices()
    if any(g.shape!=(8,8) for g in clifford):
        raise ArithmeticError('complete native corner representation changed')
    arrays={'lattice_hamiltonian':_freeze(h),'gauged_hamiltonian':_freeze(changed),
            'whole_site_gauge':_freeze(whole),'plaquette':_freeze(plaquette),
            'bond_one_body_generator':_freeze(bond),'bond_ideal_one_body_word':_freeze(ideal_word)}
    requests=[]
    paths=[((i,F(1,8)),) for i in range(12)]+[r['path'] for r in data['links']]+list(data['gauge_paths'])
    for i,path in enumerate(dict.fromkeys(paths)):
        row=current_path(path);arrays[f'current_{i}_full32']=row['full32']
        requests.append({'array':f'current_{i}_full32','path':[[j,str(t)] for j,t in path],
                         'word':row['word'],'resources':row['resources'],
                         'operator_error_diagnostic':row['operator_error_diagnostic'],
                         'selected_to_auxiliary_complement_error_diagnostic':row['selected_to_auxiliary_complement_error_diagnostic']})
    for i,block in enumerate(word.catalogue):
        arrays[f'bond_factor_{i}_actual4']=block['actual4']
        arrays[f'bond_factor_{i}_target4']=block['target4']
    report={'schema':'oph.native_gauge_matter_transport.v1',
        'source_parameters':{'exchange_J':str(field.native_parameters()['exchange']),
            'link_paths_are_declared_background_controls':True,'measured_inputs':[]},
        'same_current_carrier':{'native_ports':12,'full_auxiliary_levels':32,
            'selected_internal_components':15,'family_copies':3,
            'physical_CAR_is_separate_from_auxiliary_exterior_parity':True,
            'canonical_current_requests':requests},
        'local_covariance':{'site_count':n,'spatial_dimension':3,
            'oriented_bond_law':'h_xy=-i*J*eta(x,axis)*rho(U_xy)/a; h_yx=h_xy^dagger',
            'link_transformation':'U_xy -> G_x U_xy G_y^dagger',
            'independent_site_transformations':n,'operator_error_diagnostic':covariance,
            'nonflat_plaquette_operator_defect':nonflat,
            'time_dependent_frames_or_temporal_connection_constructed':False},
        'native_bond_word':word.record(),
        'native_bond_one_body_identity_error_diagnostic':bond_error,
        'complete_species':{'corner_components':8,'internal_components':15,'family_copies':3,
            'free_reduced_momentum_spinor_dimension':360,
            'Dirac_tastes_per_internal_component':2,'Dirac_fields_counting_internal_and_family':90,
            'both_chiralities_and_both_tastes_retained':True,'chiral_SM_projection_performed':False},
        'scope':{'same_native_current_linked_CAR_transport_constructed':True,
            'full_Fock_word_defined_by_all_local_source_factors':True,
            'dense_entire_Fock_operator_or_native_many_carrier_CE_history_executed':False,
            'local_factor_numerics_are_outward_or_hardware_certificates':False,
            'ideal_continuous_current_and_proper_code_services_remain_declared':True,
            'cubic_geometry_link_paths_inverse_controls_and_family_copy_are_architecture':True,
            'native_onsite_scalar_gaps_identified_with_SM_Higgs_masses':False,
            'gauge_field_dynamics_Gauss_constraint_and_gauge_measure_derived':False,
            'chiral_interacting_continuum_or_measured_spectrum_derived':False},
        'theorem_connection':{'Lean':'Screen/LocalGaugeJetAction.lean: covariant_first_jet and curvature_second_jet',
            'added_result':'finite same-current linked hopping and a native full-Fock factor construction',
            'nearest_lattice_source':'https://doi.org/10.1103/PhysRevD.16.3031'}}
    if include_continuum:
        from . import native_gauge_continuum as continuum
        report['continuum'],extra=continuum.build()
        arrays.update({'continuum__'+key:value for key,value in extra.items()})
    if source_bindings()!=pins:raise RuntimeError('source changed during transport construction')
    return report,arrays


def source_bindings():
    from . import native_gauge_continuum as continuum
    names=set(matter.source_bindings())|set(field.source_bindings())|{
        'model/common_source/native_gauge_matter_transport.py',
        'model/common_source/test_native_gauge_matter_transport.py',
        'Lean/Screen/LocalGaugeJetAction.lean'}
    names.update(continuum.source_bindings())
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _json_copy(value):
    return strict_json.loads(json.dumps(value,allow_nan=False))


def _compare(actual,expected,path='receipt'):
    """Typed exact structure; only floating diagnostics have a numerical gate."""
    if type(actual) is not type(expected):raise ValueError(path+' has a different type')
    if isinstance(expected,dict):
        if set(actual)!=set(expected):raise ValueError(path+' has a different key set')
        for key in expected:_compare(actual[key],expected[key],path+'.'+key)
    elif isinstance(expected,list):
        if len(actual)!=len(expected):raise ValueError(path+' has a different length')
        for j,(a,b) in enumerate(zip(actual,expected,strict=True)):_compare(a,b,f'{path}[{j}]')
    elif type(expected) is float:
        if not math.isfinite(actual) or not math.isfinite(expected) or abs(actual-expected)>REPLAY_ADMISSION*max(1,abs(expected)):
            raise ValueError(path+' differs beyond the diagnostic replay gate')
    elif actual!=expected:raise ValueError(path+' differs')


def verify_arrays(report,arrays):
    """Admit retained arrays on their own complete domains before source replay."""
    from . import native_gauge_continuum as continuum
    shapes={'lattice_hamiltonian':(360,360),'gauged_hamiltonian':(360,360),
        'whole_site_gauge':(360,360),'plaquette':(45,45),
        'bond_one_body_generator':(90,90),'bond_ideal_one_body_word':(90,90)}
    requests=report['same_current_carrier']['canonical_current_requests']
    for row in requests:
        if row['array'] in shapes:raise ValueError('duplicate current array declaration')
        shapes[row['array']]=(32,32)
    factors=report['native_bond_word']['source_local_factors']
    word=report['native_bond_word'];total=len(word['factor_order'])*LOCAL_ADMISSION
    if (type(word['physical_CAR_modes']) is not int or word['physical_CAR_modes']!=90
            or word['complete_Fock_dimension']!=str(1<<90)
            or type(word['complete_local_Fock_columns']) is not int or word['complete_local_Fock_columns']!=4
            or word['local_operator_admission']!=str(LOCAL_ADMISSION) or total>=1
            or word['unrestricted_Fock_operator_error_upper_conditional']!=str(total/(1-total))):
        raise ValueError('complete Fock domain or exact source budget differs')
    for i in range(len(factors)):
        shapes[f'bond_factor_{i}_actual4']=(4,4);shapes[f'bond_factor_{i}_target4']=(4,4)
    if 'continuum' in report:shapes.update({'continuum__'+key:shape for key,shape in continuum.SCHEMA.items()})
    if type(arrays) is not dict or set(arrays)!=set(shapes):
        raise ValueError('complete exact transport/continuum array menu required')
    for key,shape in shapes.items():
        a=arrays[key]
        if not isinstance(a,np.ndarray) or a.shape!=shape or a.dtype!=np.dtype(complex) or not np.isfinite(a).all():
            raise ValueError('invalid complete array: '+key)
    canonical=matter.canonical_currents();encoder=np.asarray(matter.representation()['encoder'],complex)
    for row in requests:
        ideal=np.eye(32,dtype=complex)
        for index,angle in _path([(j,F(t)) for j,t in row['path']]):
            g=matter.exterior_derivative(matter.spin_lift(canonical[index]))
            ideal=expm(float(angle)*np.asarray(g,complex))@ideal
        actual=arrays[row['array']]
        error=float(np.linalg.norm(actual-ideal,2))
        leakage=float(np.linalg.norm((np.eye(32)-encoder@encoder.T)@actual@encoder,2))
        if max(error,leakage)>float(CURRENT_ADMISSION):
            raise ValueError('retained full current fails its unchanged source admission')
        _compare(float(row['operator_error_diagnostic']),error,'current own-array error')
        _compare(float(row['selected_to_auxiliary_complement_error_diagnostic']),leakage,'current own-array leakage')
    one=np.eye(90,dtype=complex);targets=[]
    for i,row in enumerate(factors):
        real,imag,time=map(F,row['parameters']);z=complex(float(real),float(imag));kind=row['kind']
        if kind=='hopping':u=expm(-1j*float(time)*np.array([[0,z],[z.conjugate(),0]],complex))
        elif kind=='onsite':u=np.diag([np.exp(-1j*float(time*real)),1])
        elif kind=='fswap':u=np.array([[0,1],[1,0]],complex)
        else:raise ValueError('unknown retained physical CAR factor')
        target=_lift_two(u);targets.append(u)
        if np.linalg.norm(arrays[f'bond_factor_{i}_target4']-target,2)>1e-12:
            raise ValueError('retained target is not the complete Gaussian Fock factor')
        if np.linalg.norm(arrays[f'bond_factor_{i}_actual4']-target,2)>float(LOCAL_ADMISSION):
            raise ValueError('retained full local source factor exceeds its unchanged admission')
    for block,sites in report['native_bond_word']['factor_order']:
        if type(block) is not int or not 0<=block<len(targets) or len(sites)!=2 or any(type(j) is not int or not 0<=j<90 for j in sites):
            raise ValueError('invalid complete factor order')
        one[list(sites)]=targets[block]@one[list(sites)]
    if np.linalg.norm(one-arrays['bond_ideal_one_body_word'],2)>1e-10:
        raise ValueError('retained factor order and complete one-body word disagree')
    h=arrays['bond_one_body_generator'];error=float(np.linalg.norm(one-expm(-.125j*h),2))
    if np.linalg.norm(h-h.conj().T,2)>1e-12 or error>1e-10:
        raise ValueError('retained linked bond fails its Hermitian generator identity')
    _compare(float(report['native_bond_one_body_identity_error_diagnostic']),error,'bond own-array identity')
    h,g,changed=map(arrays.__getitem__,('lattice_hamiltonian','whole_site_gauge','gauged_hamiltonian'))
    cov=float(np.linalg.norm(changed-g@h@g.conj().T,2))
    nonflat=float(np.linalg.norm(arrays['plaquette']-np.eye(45),2))
    if max(np.linalg.norm(h-h.conj().T,2),np.linalg.norm(changed-changed.conj().T,2),cov)>1e-10 or nonflat<1e-3:
        raise ValueError('retained local covariance/nonflatness fails')
    _compare(float(report['local_covariance']['operator_error_diagnostic']),cov,'lattice own-array covariance')
    _compare(float(report['local_covariance']['nonflat_plaquette_operator_defect']),nonflat,'plaquette own-array diagnostic')
    if 'continuum' in report:
        continuum.verify_arrays(report['continuum'],{key:arrays['continuum__'+key] for key in continuum.SCHEMA})


def check(path,arrays_path=None):
    path=Path(path);pins=source_bindings();retained=strict_json.loads(path.read_text(encoding='utf-8'))
    if type(retained) is not dict or retained.get('source_bindings')!=pins:raise ValueError('source bindings differ')
    arrays_path=path.with_name(retained['artifact']['path']) if arrays_path is None else Path(arrays_path)
    if hashlib.sha256(arrays_path.read_bytes()).hexdigest()!=retained['artifact']['sha256']:
        raise ValueError('array artifact hash differs')
    with np.load(arrays_path,allow_pickle=False) as archive:
        if len(archive.files)!=len(set(archive.files)):raise ValueError('duplicate array keys')
        arrays={key:archive[key] for key in archive.files}
    scientific={k:v for k,v in retained.items() if k not in ('source_bindings','artifact')}
    verify_arrays(scientific,arrays)
    fresh,new_arrays=build();fresh=_json_copy(fresh)
    _compare(scientific,fresh)
    if set(arrays)!=set(new_arrays):raise ValueError('fresh array menu differs')
    for key in arrays:
        if arrays[key].shape!=new_arrays[key].shape or arrays[key].dtype!=new_arrays[key].dtype or np.max(abs(arrays[key]-new_arrays[key]))>REPLAY_ADMISSION:
            raise ValueError('fresh complete array differs: '+key)
    if source_bindings()!=pins:raise RuntimeError('source changed during replay')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('native_gauge_matter_transport_results.json'))
    parser.add_argument('--arrays',type=Path)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    if args.check:
        check(args.check,args.arrays);print('native gauge matter transport/continuum receipt: PASS');return
    pins=source_bindings();report,arrays=build();report=_json_copy(report)
    verify_arrays(report,arrays)
    if source_bindings()!=pins:raise RuntimeError('source changed during generation')
    target=args.arrays or args.out.with_suffix('.npz')
    np.savez_compressed(target,**arrays)
    report.update(source_bindings=pins,artifact={'path':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
    args.out.write_text(json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('native gauge matter transport/continuum receipt generated')


if __name__=='__main__':main()

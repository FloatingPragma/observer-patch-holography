#!/usr/bin/env python3
"""Certified finite H_B excitation energies; no relativistic mass assignment.

This tool independently replays the exact256-state dyadic eigensystem proof.
Only its four fixed input files and NumPy's archive/integer-array facilities
are used; no eigensolver or broad package import is required.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import types
import numpy as np

PINS={
 'model/native_ordered_ground/source.py':'e0ca9a40181637bc9263b55ff41dc9b6aa895f723a6d497d182a011b80097dc2',
 'model/native_ordered_ground/certificate.py':'307fd37c11b4312495700b5c8cd411d2b20001cd1b7765576be21abf2aa8db38',
 'model/native_ordered_ground/witness.npz':'5df39bfaf63d5af40f80425f74011b367f9d9c5fd8951eb0419ced77da5396b4',
 'model/native_ordered_ground/results.json':'8412f3f752614b8265aa8ead9c14df169c44ab5ec66fe2f2759df625a6e76db2',
}
BITS=48
SIZE=256


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def load_json(path):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:raise ValueError('duplicate JSON key')
            result[key]=value
        return result
    def bad(_):raise ValueError('nonfinite JSON constant')
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=bad)


def checked_inputs(root):
    root=Path(root).resolve();result={}
    for name,digest in PINS.items():
        data=(root/name).read_bytes()
        if hashlib.sha256(data).hexdigest()!=digest:
            raise ValueError('reviewed source pin changed: '+name)
        result[name]=data
    return result


def sqrt_upper(value):
    value=F(value)
    if value<0:raise ValueError('nonnegative exact radicand required')
    scale=1<<100;n=isqrt(value.numerator*scale*scale//value.denominator)
    lo=F(n,scale)
    return lo if lo*lo==value else F(n+1,scale)


def decimal_interval(lo,hi,places=12):
    if type(places) is not int or not 0<=places<=30 or lo>hi:
        raise ValueError('ordered exact interval and bounded precision required')
    scale=10**places
    lower=(lo*scale).numerator//(lo*scale).denominator
    upper=-((-hi*scale).numerator//(-hi*scale).denominator)
    def render(n):
        sign='-' if n<0 else '';whole,part=divmod(abs(n),scale)
        return sign+str(whole)+('.'+str(part).zfill(places) if places else '')
    return [render(lower),render(upper)]


def verify_witness(witness,h,scale):
    """Independent exact reproduction of FA98's polar/Weyl certificate."""
    if set(witness)!={'energies','vectors'}:raise ValueError('complete witness required')
    if witness['energies'].shape!=(SIZE,) or witness['vectors'].shape!=(SIZE,SIZE):
        raise ValueError('full256-state eigensystem required')
    if any(a.dtype!=np.dtype('int64') for a in witness.values()):
        raise ValueError('literal int64 witness required')
    if h.shape!=(SIZE,SIZE) or h.dtype!=np.dtype('int64') or scale!=8 or not np.array_equal(h,h.T):
        raise ValueError('exact declared source matrix required')
    if max(F(sum(abs(int(x)) for x in row),scale) for row in h)>14:
        raise ValueError('source matrix exceeds declared norm bound')
    e=witness['energies'].astype(object);v=witness['vectors'].astype(object);d=1<<BITS
    if any(a>b for a,b in zip(e,e[1:])):raise ValueError('ordered energies required')
    if any(x<=0 for x in v[:,0]):raise ValueError('positive Marshall ground required')
    gram=v.T@v-d*d*np.eye(SIZE,dtype=object)
    defect=F(max(sum(abs(int(x)) for x in row) for row in gram),d*d)
    if defect>=F(1,10**6):raise ValueError('complete Gram defect too large')
    residual=h.astype(object)@v*d-scale*v*e[None,:]
    eps=sqrt_upper(F(sum(int(x)**2 for x in residual.flat),scale**2*d**4))
    norm_bound=F(14)
    if any(abs(F(int(x),d))>norm_bound for x in e):raise ValueError('proposal outside source norm bound')
    radius=eps+2*norm_bound*defect;e0,e1=F(int(e[0]),d),F(int(e[1]),d)
    if e1-e0<=2*radius:raise ValueError('positive isolated first gap required')
    norm=F(sum(int(x)**2 for x in v[:,0]),d*d)
    r0=sqrt_upper(F(sum(int(x)**2 for x in residual[:,0]),scale**2*d**4)/norm)
    error=2*r0/(e1-radius-e0)
    if error>=F(1,10**8):raise ValueError('ground-vector enclosure too loose')
    raw=sum(int(v[m,0])*int(v[m^(1<<site),0]) for m in range(SIZE) for site in range(8))
    mag=F(raw,16*d*d)/norm
    result={'dimension':SIZE,'vector_bits':BITS,'Gram_defect':str(defect),
      'complete_polar_Weyl_radius':str(radius),
      'ground_energy_interval':decimal_interval(e0-r0,e0+r0,13),
      'first_excited_interval':decimal_interval(e1-radius,e1+radius),
      'finite_field_gap_interval':decimal_interval(e1-radius-e0-r0,e1+radius-e0+r0),
      'normalized_ground_vector_error':str(error),
      'staggered_magnetization_per_site_interval':decimal_interval(mag-error,mag+error),
      'all_256_eigenvalues_certified':True,'untrusted_eigensolver_is_not_the_proof':True}
    return result,[F(int(x),d) for x in e],radius,(e0-r0,e0+r0)


def shared_ground_ratio(level,reference,ground,ordered_after_reference=False):
    """Exact box extrema with the SAME ground energy in numerator/denominator."""
    intervals=[tuple(map(F,pair)) for pair in (level,reference,ground)]
    if any(len(pair)!=2 or pair[0]>pair[1] for pair in intervals):
        raise ValueError('ordered endpoint pairs required')
    level,reference,ground=intervals
    if reference[0]<=ground[1] or level[0]<=ground[1]:
        raise ValueError('positive separated excitation intervals required')
    corners=[(e-g)/(r-g) for e in level for r in reference for g in ground]
    lo,hi=min(corners),max(corners)
    if ordered_after_reference:lo=max(F(1),lo)
    if lo>hi:raise ValueError('inconsistent ordered excitation intervals')
    return lo,hi


def extract(energy,radius,ground):
    if len(energy)!=SIZE or any(a>b for a,b in zip(energy,energy[1:])):
        raise ValueError('complete ordered energy centers required')
    e0lo,e0hi=ground;e1lo,e1hi=energy[1]-radius,energy[1]+radius
    if radius<0 or e0lo>e0hi or e1lo<=e0hi:raise ValueError('invalid spectral enclosures')
    rows=[];clusters=[]
    for i,e in enumerate(energy):
        lo,hi=ground if i==0 else(e-radius,e+radius)
        gaplo,gaphi=(F(0),F(0)) if i==0 else(lo-e0hi,hi-e0lo)
        if i==0:qlo=qhi=F(0)
        elif i==1:qlo=qhi=F(1)
        else:
            # Keep the common ground variable. With positive denominators,
            # the linear-fractional extrema occur among the eight box corners.
            qlo,qhi=shared_ground_ratio((lo,hi),(e1lo,e1hi),ground,True)
        rows.append({'ordered_state_index':i,'energy_interval':decimal_interval(lo,hi),
                     'gap_interval':decimal_interval(gaplo,gaphi),
                     'gap_over_first_gap_interval':decimal_interval(qlo,qhi),
                     'exact_energy_interval':[str(lo),str(hi)],
                     'exact_gap_interval':[str(gaplo),str(gaphi)],
                     'exact_ratio_interval':[str(qlo),str(qhi)]})
        if not clusters or lo>clusters[-1]['hi']:
            clusters.append({'lo':lo,'hi':hi,'indices':[i]})
        else:
            clusters[-1]['hi']=max(hi,clusters[-1]['hi']);clusters[-1]['indices'].append(i)
    if clusters[0]['indices']!=[0]:raise ValueError('ground cluster is not isolated')
    summaries=[]
    for n,c in enumerate(clusters):
        i,j=c['indices'][0],c['indices'][-1]
        summaries.append({'cluster':n,'ordered_state_indices':[i,j],
          'certified_state_count':j-i+1,'exact_degeneracy_claimed':False,
          'energy_interval':decimal_interval(c['lo'],c['hi']),
          'gap_interval':[rows[i]['gap_interval'][0],rows[j]['gap_interval'][1]],
          'gap_over_first_gap_interval':[rows[i]['gap_over_first_gap_interval'][0],rows[j]['gap_over_first_gap_interval'][1]]})
    if sum(c['certified_state_count'] for c in summaries)!=SIZE:raise ArithmeticError('incomplete cluster cover')
    return rows,summaries


def report(root):
    root=Path(root).resolve();tool_bytes=Path(__file__).read_bytes();inputs=checked_inputs(root)
    name='model/native_ordered_ground/source.py'
    source=types.ModuleType('_verified_ordered_spin_source')
    exec(compile(inputs[name],str(root/name),'exec'),source.__dict__)
    h,scale=source.matrix()
    import io
    with np.load(io.BytesIO(inputs['model/native_ordered_ground/witness.npz']),allow_pickle=False) as archive:
        witness={k:archive[k] for k in archive.files}
    certificate,energy,radius,ground=verify_witness(witness,h,scale)
    recorded=load_json(root/'model/native_ordered_ground/results.json')['finite_spectrum']
    if canonical(certificate)!=canonical(recorded):raise ValueError('independent finite certificate does not reproduce FA98')
    rows,clusters=extract(energy,radius,ground)
    if checked_inputs(root)!=inputs:raise ValueError('source changed during energy extraction')
    if Path(__file__).read_bytes()!=tool_bytes:raise ValueError('readout tool changed during extraction')
    return {'schema':'oph.ordered_spin_energy_readout.v1',
      'scope':'CERTIFIED_FINITE_SELECTION_ENERGIES_NOT_RELATIVISTIC_MASSES',
      'Hamiltonian':'H_B=J sum_12edges(S_x dot S_y-1/4)-B sum_x epsilon_x S_x^1',
      'graph':'open2x2x2cube','sites':8,'Hilbert_dimension':SIZE,'J':'1','B_over_J':'1/4',
      'source_energy_unit':'J','relative_reference':'Delta_1=E_1-E_0 of the same H_B',
      'source_bindings':dict(PINS),'tool_sha256':hashlib.sha256(tool_bytes).hexdigest(),
      'independent_exact_certificate':certificate,'all_ordered_levels':rows,
      'isolated_clusters':clusters,'cluster_count':len(clusters),'total_certified_states':SIZE,
      'exact_degeneracies_identified':False,'all_levels_detector_visible_claimed':False,
      'zero_field_source_ground_spectrum':False,'native_energy_measurement_executed_here':False,
      'relativistic_masses_identified':False,'physical_species_identified':False,
      'note':'Certified finite selection-Hamiltonian energy differences. B is absent from the zero-field limiting dynamics. Cluster counts do not assert exact degeneracy; detector weights and relativistic mass shells are not inferred.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path.cwd())
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',type=Path)
    group.add_argument('--check',type=Path)
    args=parser.parse_args();value=report(args.root)
    if args.write:args.write.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
    elif args.check:
        if canonical(load_json(args.check))!=canonical(value):raise ValueError('finite energy readout does not replay')
    else:print(json.dumps(value,sort_keys=True,indent=2))

if __name__=='__main__':main()

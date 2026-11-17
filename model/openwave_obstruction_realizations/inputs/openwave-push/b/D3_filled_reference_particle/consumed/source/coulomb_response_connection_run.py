"""Certified finite static charge energies on the actual LOCAL12 graph."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE/'inputs'
SPATIAL = INPUT/'reverse-engineering-reality/model/native_spatial_refinement/results.json'
WALK = INPUT/'openwave-tuning/auditor/native_walk_bridge/COMMON_WALK_RECEIPT.json'
LINKS = tuple(tuple(map(tuple,row)) for row in json.loads(SPATIAL.read_text())['intercell_oriented_port_pairs'])
INTERNAL = tuple(tuple(edge) for row in json.loads(WALK.read_text())['descriptor']['internal_matchings'] for edge in row)
DENOMINATOR = 3*(1 << 58)


def symbols(k):
    k=np.asarray(k,float)
    H=np.broadcast_to(6*np.eye(12,dtype=complex),k.shape[:-1]+(12,12)).copy()
    for a,b in INTERNAL:
        H[...,a,b]-=1; H[...,b,a]-=1
    for axis,row in enumerate(LINKS):
        phase=np.exp(1j*k[...,axis])
        for a,b in row:
            H[...,a,b]-=phase; H[...,b,a]-=phase.conjugate()
    return H


def laplacian(values):
    # The same expression works over floating or exact int64 coordinates.
    out=6*values.copy()
    for a,b in INTERNAL:
        out[...,a]-=values[...,b]; out[...,b]-=values[...,a]
    for axis,row in enumerate(LINKS):
        for a,b in row:
            out[...,a]-=np.roll(values[...,b],-1,axis=axis)
            out[...,b]-=np.roll(values[...,a],1,axis=axis)
    return out


def green_column(size):
    if type(size) is not int or size < 2 or size & (size-1):
        raise ValueError('declared dyadic cubic member required')
    k=2*np.pi*np.arange(size)/size
    yz=np.array(np.meshgrid(k,k,indexing='ij')).transpose(1,2,0)
    fourier=np.empty((size,size,size,12),complex)
    for x in range(size):
        points=np.empty((size,size,3)); points[...,0]=k[x]; points[...,1:]=yz
        H=symbols(points)
        rhs=np.zeros((size,size,12,1),complex); rhs[...,0,0]=1
        if x==0:
            H[0,0]+=np.ones((12,12))/12
            rhs[0,0,:,0]-=1/12
        fourier[x]=np.linalg.solve(H,rhs)[...,0]
    spatial=np.fft.ifftn(fourier,axes=(0,1,2))
    assert np.max(abs(spatial.imag)) < 1e-12
    return spatial.real


def certify(column):
    size=column.shape[0]; vertices=12*size**3
    assert DENOMINATOR % vertices == 0
    integer=np.rint(column*DENOMINATOR).astype(np.int64)
    # Each operation is exact int64 arithmetic, with an explicit worst-case
    # intermediate bound rather than an unchecked overflow assumption.
    maximum=max(abs(int(integer.min())),abs(int(integer.max())))
    assert 12*maximum+DENOMINATOR < 1 << 63
    residual=laplacian(integer)
    residual+=DENOMINATOR//vertices
    residual[0,0,0,0]-=DENOMINATOR
    assert sum(map(int,residual.ravel()))==0
    square=sum(int(x)**2 for x in residual.ravel())
    root=math.isqrt(square)
    if root*root < square: root+=1
    # Internal diameter3, each cell crossing uses at most3+1edges.
    diameter=6*size+3
    gap=Fraction(1,vertices*diameter)
    residual_upper=Fraction(root,DENOMINATOR)
    energy_upper=residual_upper/gap  # a=1/2 and sqrt2 <=2.
    return integer, {'denominator':str(DENOMINATOR),
        'integer_residual_squared':str(square),'residual_norm_upper_exact':str(residual_upper),
        'graph_diameter_upper':diameter,'positive_laplacian_gap_lower_exact':str(gap),
        'pair_energy_absolute_error_upper_exact':str(energy_upper),
        'pair_energy_absolute_error_upper':float(energy_upper),
        'integer_vector_sha256':hashlib.sha256(integer.astype('<i8').tobytes()).hexdigest(),
        'zero_mean_projection_needed_for_vector_bound':True,
        'energy_differences_invariant_under_that_projection':True}


def direct_control(column):
    size=2; vertices=96
    # Independent direct graph assembly from oriented internal and cell links.
    L=np.zeros((vertices,vertices),int)
    for x in np.ndindex((size,size,size)):
        cell=np.ravel_multi_index(x,(size,size,size))
        for a,b in INTERNAL:
            i,j=12*cell+a,12*cell+b
            L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
        for axis,row in enumerate(LINKS):
            y=list(x);y[axis]=(y[axis]+1)%size
            nxt=np.ravel_multi_index(tuple(y),(size,size,size))
            for a,b in row:
                i,j=12*cell+a,12*nxt+b
                L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    assert np.all(np.diag(L)==6)
    P=np.ones((vertices,vertices))/vertices
    G=np.linalg.inv(L+P)-P
    residual=float(np.max(abs(G[:,0]-column.ravel())))
    assert residual<2e-14
    frozen=json.loads((INPUT/'openwave-tuning/common_dynamics/quadratic_field/SOURCE.json').read_text())
    errors=[]
    for record in frozen['static_charge_cases']:
        a,b=record['vertices'];q=np.eye(vertices)[a]-np.eye(vertices)[b]
        energy=float(q@G@q)/4
        errors.append(abs(energy-record['electric_minimum_energy']))
    assert max(errors)<2e-13
    # Independent check of the stated internal diameter used in the exact bound.
    neighbours=[set() for _ in range(12)]
    for a,b in INTERNAL:neighbours[a].add(b);neighbours[b].add(a)
    diam=0
    for source in range(12):
        seen={source};front={source};distance=0
        while len(seen)<12:
            front=set().union(*(neighbours[x] for x in front))-seen
            assert front;seen|=front;distance+=1
        diam=max(diam,distance)
    assert diam==3
    return {'direct_96vertex_green_error':residual,
            'frozen_quadratic_source_energy_errors':errors,'internal_diameter':diam}


def fit(distances,energies,error):
    x=1/np.asarray(distances,float); y=np.asarray(energies,float)
    design=np.column_stack((np.ones_like(x),x))
    coefficient=np.linalg.lstsq(design,y,rcond=None)[0]
    residual=y-design@coefficient
    denominator=float(np.sum((y-y.mean())**2))
    r2=1-float(residual@residual)/denominator
    # Exact-energy uncertainty gives an explicit interval for the fitted slope.
    centered=x-x.mean(); slope_radius=error*float(np.abs(centered).sum()/np.dot(centered,centered))
    return {'intercept':float(coefficient[0]),'inverse_distance_coefficient':float(coefficient[1]),
            'coefficient_certificate_radius':slope_radius,'R_squared':r2,
            'fit_is_diagnostic_not_a_supplied_source_parameter':True}


def main():
    pins=json.loads((HERE/'INPUTS.json').read_text())['files_sha256']
    for path,h in pins.items():assert hashlib.sha256((INPUT/path).read_bytes()).hexdigest()==h
    cases=[]; independent=None
    for size in (2,4,8,16,32):
        column=green_column(size);integer,certificate=certify(column)
        if size==2:independent=direct_control(column)
        ranges=range(2,max(4,size//4+1)) if size>=16 else range(1,size//2+1)
        ranges=[r for r in ranges if r<size//2]
        axes=[]
        for axis in range(3):
            values=[]
            for r in ranges:
                index=[0,0,0,0];index[axis]=r
                # Two opposite unit charges at the same local port, translated.
                energy=Fraction(int(integer[0,0,0,0])-int(integer[tuple(index)]),2*DENOMINATOR)
                values.append({'separation_cells':r,'energy_exact_trial':str(energy),'energy':float(energy)})
            row={'axis':axis,'opposite_charge_energies':values}
            if len(values)>=3:
                row['inverse_distance_fit']=fit(ranges,[v['energy'] for v in values],certificate['pair_energy_absolute_error_upper'])
                assert row['inverse_distance_fit']['inverse_distance_coefficient']<0
                row['same_sign_cross_term_has_opposite_slope']=True
            axes.append(row)
        cases.append({'side':size,'vertices':12*size**3,'seams':36*size**3,
                      'certificate':certificate,'axes':axes})
        np.savez_compressed(HERE/f'green_{size}.npz',integer=integer,denominator=np.array(str(DENOMINATOR)))
        print('Certified static member',size,'energy bound',certificate['pair_energy_absolute_error_upper'],flush=True)
    result={'schema':'oph.actual_local12.static_charge_energy.v1',
        'source':'native-local12-dirac-noncompact-quadratic-v4-proposed; static charge-sector generator',
        'electric_coefficient':'1/2','independent_size2_control':independent,'cases':cases,
        'predicted_leading_inverse_distance_coefficient':-1/(12*math.pi),
        'leading_coefficient_ancestry':'actual harmonic conductance3/2 and declared electric coefficient1/2',
        'coupling_fitted_to_observed_electromagnetism':False,
        'native_full_field_preparation_or_energy_instrument_executed':False,
        'freely_moving_force_executed':False,'spatial_infinite_volume_error_proved':False,
        'official_rows_changed':[],'physical_rows_validated':[],'source_pins':pins}
    for path,h in pins.items():assert hashlib.sha256((INPUT/path).read_bytes()).hexdigest()==h
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()

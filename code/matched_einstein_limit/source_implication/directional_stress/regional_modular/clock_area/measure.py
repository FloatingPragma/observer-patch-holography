"""Test one shared clock and flux-area conversion on six native sharp images."""
import argparse
from itertools import permutations
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from ....action_response.four_ticks import parse_interval
from ....jacobi_screen import certificate as optical, measure as jacobi
from .. import measure as regional

c = regional.c
HERE = Path(__file__).resolve().parent
NOISE = iv.mpf([-1, 1])


def image_records():
    formation = regional.bridge.native.iso.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = regional.bridge.beam.screens(board, route)
    weights = regional.bridge.native.prior.weights()[[7]]
    tagged = []
    for charge in (-1, 1):
        by_direction = [[], []]
        for odd in (False, True):
            h = regional.bridge.native.generator(board, route, charge, odd)
            source = np.array(optical.amplitude_jets(h, *rings, weights), dtype=object)
            for keep, amps in zip(by_direction, (source, source.transpose(0, 1, 2, 4, 3, 6, 5)), strict=True):
                keep.append(optical.probability_jets(amps)[:, :, 0])
        tagged.append(np.asarray(by_direction, dtype=object))
    return tagged


def moment_model(operators, multiplicity, states):
    denominator = multiplicity*optical.prior.FRAC**2
    base = np.asarray(operators[:,0,0], dtype=float)/denominator + 1j*np.asarray(operators[:,0,1], dtype=float)/denominator
    population = np.array([np.array(x, dtype=complex)*p for x,p in zip(states, (.8657, .1343))])
    outgoing, incoming = optical.translated_basis(optical.beam.probabilities(base, population))
    axes = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    centers = np.empty((2,3,3,3),dtype=object)
    coeffs = np.empty((2,15,3,3,3),dtype=object)
    norms = np.empty_like(centers)
    for k in range(3):
        error = c.ri(sp.Rational(100*10*64**k, 10**19))
        denom = 2**k*denominator*optical.beam.UNIT**2
        for i, a in enumerate(outgoing.T):
            for j, b in enumerate(incoming.T):
                coeff = np.outer(a,b)
                assert sum(abs(x) for x in coeff.ravel()) <= 100*optical.beam.UNIT**2
                re, im = (operators[:,k]*coeff[None,None,:,:,None,None]).sum((2,3)).swapaxes(0,1)
                for s, (r, z) in enumerate(zip(re, im, strict=True)):
                    block = (sp.Matrix(r) + sp.I*sp.Matrix(z))/denom
                    centers[s,k,i,j] = c.ri(c.inner(block,states[s])) + NOISE*error
                    moments = [c.inner(block, axis) for axis in axes]
                    for h, value in enumerate(moments):
                        coeffs[s,h,k,i,j] = c.ri(value) + NOISE*2*error
                    norms[s,k,i,j] = iv.sqrt(c.ri(sum(x*x for x in moments))) + 2*error
    return centers, coeffs, norms


def factor(x, value):
    # Object ndarray on the left applies the interval factor entrywise.
    return x*value


def values(model, p, radius, derivatives):
    centers, coeffs, norms = model
    blocks = centers+factor(norms, NOISE*radius)
    m = blocks[0]+factor(blocks[1]-blocks[0],p)
    j = jacobi.matrix_jets(m)
    det, expansion, acceleration = jacobi.optical_matrices(j)
    assert 0 not in det
    result=np.array([expansion.trace(), -acceleration.trace()], dtype=object)
    if not derivatives:
        return result
    variables = [factor(x,1-p) for x in coeffs[0]]+[factor(x,p) for x in coeffs[1]]+[blocks[1]-blocks[0]]
    gradient=[]
    for x in variables:
        db, dc = jacobi.deformation_response(j, jacobi.response_jets(m, x))
        gradient.append([db.trace(), -dc.trace()])
    return np.array(gradient,dtype=object)


def maxabs(interval):
    return abs(interval).b


def hull(values):
    return iv.mpf([min(v.a for v in values),max(v.b for v in values)])


def odd_pair(models, ref, controls, splits=10):
    _, raw_radius, p = ref
    radius=c.ri(raw_radius)
    spans=[[],[]]
    changes={key:[[[],[]],[[],[]]] for key in controls}
    width = p.b-p.a
    for part in range(splits):
        lo = p.a + width*part/splits
        hi = p.a + width*(part+1)/splits
        segment = iv.mpf([lo.a,hi.b])
        middle = segment.mid
        center=values(models[0],middle,c.ri(0),False)-values(models[1],middle,c.ri(0),False)
        forward=values(models[0],segment,radius,True)
        slope=forward-values(models[1],segment,radius,True)
        for k in range(2):
            ru = sum(iv.sqrt(sum(maxabs(v)**2 for v in group))
                     for group in (slope[:15,k],slope[15:30,k]))*radius
            ru += maxabs(segment-middle)*maxabs(slope[-1,k])
            spans[k].append(center[k]+NOISE*ru)
        # Control directions are coordinates of the unnormalized parity blocks,
        # while the interval derivatives above act on conditional density matrices.
        for menu,gradient in enumerate((forward,slope)):
            cov=np.concatenate([gradient[:15]/(1-segment), gradient[15:30]/segment])
            for k in range(2):
                norm=iv.sqrt(sum(maxabs(x)**2 for x in cov[:,k]))
                for key,(center,error) in controls.items():
                    changes[key][menu][k].append(sum(x*y for x,y in zip(center,cov[:,k],strict=True))
                                                  + NOISE*norm*error)
    return list(map(hull,spans)), {key:[[hull(x) for x in row] for row in v]
                                   for key,v in changes.items()}


def determinant(m):
    return sum((-1)**sum(v[i]>v[j] for i in range(len(v)) for j in range(i+1,len(v)))
               * np.prod([row[j] for row,j in zip(m,v,strict=True)])
               for v in permutations(range(len(m))))


def common_readout(opt):
    jet_keys=("source_expansion_slope","focusing","expansion","source_focusing_slope")
    rows=[]
    base={}
    for side in ("forward","reverse"):
        vacancy=opt["charges"]["-1"][side]
        base[side]=[]
        for other in (opt["charges"]["1"][side],opt["fair_untagged_images"][side]):
            row=[parse_interval(vacancy[k])-parse_interval(other[k]) for k in jet_keys]
            base[side].append(str(-row[1]/row[2]))
            rows.append(row)
    det=determinant(rows)
    assert det.a>5e6
    return dict(required_clock_at_reference=base, equal_work_four_row_determinant=str(det))


def run():
    iv.dps=35
    raw=image_records()
    ref=regional.source()
    src=json.loads((regional.HERE/"energy_odd.json").read_text())
    opt=json.loads((regional.HERE/"optical_certificate.json").read_text())
    controls={key:([c.ri(sp.Rational(v)) for v in src[key]["center_coordinates"]],
                   c.ri(sp.Rational(src[key]["direction_error_with_serialization_upper"])))
              for key in ("equal_work_control","fixed_work_control","fixed_work_second_control")}
    current=list(map(parse_interval,src["reference_energy_current_by_charge"]))
    current.append((current[0]+current[1])/2)
    out={}
    matrix=[]
    slopes=[]
    for name,data,n,q in zip(("vacancy","excess","fair"),[*raw,raw[0]+raw[1]],(1,1,2),current,strict=True):
        models=[moment_model(x,n,ref[0]) for x in data]
        odd,change=odd_pair(models,ref,controls)
        slopes.append(change)
        out[name]=dict(expansion_odd=str(odd[0]), focusing_odd=str(odd[1]), energy_current=str(q))
        matrix.append([odd[0],q,odd[1]])
    det=determinant(matrix)
    assert det.a>0, str(det)
    # Matter probabilities at fixed charge mix and a common area obey the linear
    # mixture identity on either direction, irrespective of the even stress.
    def mix(v):
        return v[2]-(v[0]+v[1])/2

    def forward(k):
        return [parse_interval(x["forward"][k]) for x in
                (opt["charges"]["-1"],opt["charges"]["1"],opt["fair_untagged_images"])]
    fixed=[[v[2],v[0],v[1]] for v in matrix]
    fixed.append([mix(forward("focusing")),mix(forward("expansion")),c.ri(0)])
    clock={}
    for key in controls:
        gradients=[]
        for order in range(2):
            delta=[x[key][1][order] for x in slopes]+[mix([x[key][0][order] for x in slopes])]
            gradients.append(determinant([[*row,w] for row,w in zip(fixed,delta,strict=True)]))
        assert 0 not in gradients[0], (key,str(gradients[0]))
        clock[key]=-gradients[1]/gradients[0]
    assert clock["fixed_work_control"].b < clock["equal_work_control"].a, {k:str(v) for k,v in clock.items()}
    return dict(exchange_dwell=4,rows=out,absolute_odd_determinant=str(det),
                common_curvature=common_readout(opt),
                differential_mixture_required_clock={k:str(v) for k,v in clock.items()})


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    data=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        (HERE/"certificate.json").write_text(data)
    else:
        print(data,end="")

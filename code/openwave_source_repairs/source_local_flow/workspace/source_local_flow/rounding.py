"""Independent rational enclosures for floating trigonometry and frames."""
from fractions import Fraction as F
import math
U=F(1,2**53)

def trig_interval(x,cosine=False):
    x=F(x);term=F(1)if cosine else x;s=term
    for k in range(1,48):
        term *= -x*x/F((2*k-1)*(2*k)if cosine else(2*k)*(2*k+1));s+=term
    k=48;tail=abs(term*x*x/F((2*k-1)*(2*k)if cosine else(2*k)*(2*k+1)))
    assert abs(x)<=4
    return s-tail,s+tail

def trig_error(x):
    co=trig_interval(x,True);si=trig_interval(x,False)
    # L1 bound exceeds the complex absolute error.
    return max(abs(F(math.cos(x))-v)for v in co)+max(abs(F(math.sin(x))-v)for v in si)

def sqrt_interval(n):
    scale=10**80;k=math.isqrt(n*scale*scale)
    return F(k,scale),F(k+1,scale)

def frame_error():
    lo,hi=sqrt_interval(12);a=F(1/math.sqrt(12))
    coefficient=max(abs(a-1/lo),abs(a-1/hi))
    # The integer Gram identities are verified independently. Each second
    # orthonormal vector has squared integer norm12; l1 norm <=12. Two
    # factor frames and their product receive this conservative enclosure.
    return 100*coefficient+16*U

def guard(program,steps,pair,energy,time):
    trig=max([trig_error(-2*h)for _,h in program]+[trig_error(-energy*time/steps)])
    assert trig<8*U
    depth=len(program);local=512*U+8*trig
    norm=F(2)if pair else F(3,2)
    one=norm*depth*local/(1-depth*local)
    stored=norm*steps*depth*local/(1-steps*depth*local)
    frame=frame_error()
    cells=240*240*4 if pair else 240*2
    dot=16*cells*U/(1-16*cells*U)*4
    clock=steps*abs(F(float(energy*time/steps))-F(str(energy))*F(str(time))/steps)
    return {'directed_trig_L1_error_upper':str(trig),'input_frame_error_upper':str(frame),
            'one_step_state_arithmetic_upper':str(one+2*frame+clock),
            'stored_full_propagation_arithmetic_upper':str(stored+frame),
            'final_inner_product_readout_upper':str(dot),'clock_rounding_upper':str(clock),
            'scope':'Exact ideal local word versus one-step centre, stored whole-product centre and final reference Gram are separately enclosed. No unit-gain assumption is made for rounded matrices.'}

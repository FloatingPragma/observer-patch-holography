"""Exact full-Fock free comparator for the admitted five-block packets.

Finite Gaussian kernel derivatives evaluate the infinite-Fock-space
propagator between the actual finite input and detector vectors. No
oscillator cutoff is applied to the intervening free evolution.
"""
from functools import lru_cache
import itertools
import mpmath as mp


def compare(incoming,outgoing,words,clocks,cosines,exact_endpoint):
    iv=mp.iv
    V=5;D=iv.mpf(2);coupling=iv.mpf(1)/8
    omega=[iv.sqrt(D*D-4*D*coupling*c) for c in cosines]
    h=[D-2*coupling*c for c in cosines]
    support=[i for i,w in enumerate(words)
             if w.count(1)%2==1 and w.count(2)%2==1 and w.count(3)%2==0]
    subsets={i:tuple(tuple(x for x,label in enumerate(words[i]) if label==color)
                     for color in (1,2,3)) for i in support}
    one_color=sorted({sites for row in subsets.values() for sites in row})
    output=[]
    for clock in clocks:
        if not clock:
            output.append({'clock_Delta_times_t':str(clock),
                           'same_admitted_free_hit_probability_interval_exact':['0','0'],
                           'same_admitted_free_complement_probability_interval_exact':['1','1'],
                           'same_admitted_free_hit_display_midpoint':0.0})
            continue
        time=iv.mpf(clock.numerator)/clock.denominator
        Astar=[iv.mpc(iv.cos(w*time),hh*iv.sin(w*time)/w)
               for w,hh in zip(omega,h)]
        B=[iv.mpc(0,2*coupling*c*iv.sin(w*time)/w)
           for w,c in zip(omega,cosines)]
        X=[1/a for a in Astar]
        Z=[b/a for a,b in zip(Astar,B)]
        Y=[-iv.mpc(b.real,-b.imag)/a for a,b in zip(Astar,B)]
        # Each Astar remains in the upper half-plane for 0<T<=1;
        # the continuous square-root branch is therefore explicit.
        prefactor=iv.mpc(iv.cos(D*V*time/2),iv.sin(D*V*time/2))
        for a in Astar:
            assert a.imag.a>0
            modulus=iv.sqrt(a.real**2+a.imag**2)
            root=iv.mpc(iv.sqrt((modulus+a.real)/2),
                        iv.sqrt((modulus-a.real)/2))
            prefactor/=root
        def spatial(values):
            return [[sum((values[j]*cosines[(j*(x-y))%V]
                          for j in range(V)),iv.mpc(0))/V
                     for y in range(V)] for x in range(V)]
        XX,ZZ,YY=map(spatial,(X,Z,Y))
        kernel=[[iv.mpc(0) for _ in range(2*V)] for _ in range(2*V)]
        for x,y in itertools.product(range(V),repeat=2):
            kernel[x][y]=ZZ[x][y]
            kernel[x+V][y+V]=YY[x][y]
            kernel[x][y+V]=XX[x][y]
            kernel[y+V][x]=XX[x][y]
        @lru_cache(None)
        def hafnian(indices):
            if not indices:return iv.mpc(1)
            if len(indices)%2:return iv.mpc(0)
            return sum((kernel[indices[0]][indices[j]]
                        *hafnian(indices[1:j]+indices[j+1:])
                        for j in range(1,len(indices))),iv.mpc(0))
        entries={(left,right):prefactor*hafnian(left+tuple(V+j for j in right))
                 for left,right in itertools.product(one_color,repeat=2)}
        amplitude=iv.mpc(0)
        for i,j in itertools.product(support,repeat=2):
            entry=iv.mpc(1)
            for color in range(3):
                entry*=entries[subsets[i][color],subsets[j][color]]
            amplitude+=iv.mpc(outgoing[i].real,-outgoing[i].imag)*entry*incoming[j]
        probability=amplitude.real**2+amplitude.imag**2
        low=max(0,exact_endpoint(probability))
        high=min(1,exact_endpoint(probability,True))
        controls=[]
        for left,right in (((),()),((0,1),()),((0,),(0,)),((1,),(0,)),((0,1,2),(0,))):
            entry=entries[left,right]
            controls.append({'outgoing_sites':list(left),'incoming_sites':list(right),
                             'real_interval_exact':[str(exact_endpoint(entry.real)),str(exact_endpoint(entry.real,True))],
                             'imaginary_interval_exact':[str(exact_endpoint(entry.imag)),str(exact_endpoint(entry.imag,True))]})
        output.append({'clock_Delta_times_t':str(clock),
                       'same_admitted_free_amplitude_real_interval':str(amplitude.real),
                       'same_admitted_free_amplitude_imaginary_interval':str(amplitude.imag),
                       'same_admitted_free_hit_probability_interval_exact':[str(low),str(high)],
                       'same_admitted_free_complement_probability_interval_exact':[str(1-high),str(1-low)],
                       'same_admitted_free_hit_display_midpoint':float((low+high)/2),
                       'independent_one_color_kernel_controls':controls})
    return output

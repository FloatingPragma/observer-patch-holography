"""Native selector family; exact prefixes and coherent-label amplitudes."""
from fractions import Fraction
from itertools import product
import math
import numpy as np
import native_source as common


def selector(count,cap,graph):
    width=(count-1).bit_length();base=2**width;good=tuple(range(count));bad=tuple(range(count,base))
    words=[]
    for word in range(base):
        bits=[(word>>(width-1-j))&1 for j in range(width)];records=[]
        for bit in bits:
            out,record=common.native([0]*graph['N'],tuple(graph['blank']),bit)
            assert not any(out) and record[3]==0;records.append(record)
        words.append({'word':word,'bits':bits,'native_records':records,'probability':str(Fraction(1,base)),
                      'predicate':'accept' if word<count else 'retry'})
    # Attempt index is a coarse flag; all rejected words and the accepted word
    # stay as coherent memory, never replaced by a classical uniform mixture.
    bell=np.array([1.,0.,0.,1.],complex)/math.sqrt(2)
    target=np.outer(bell,bell.conj());attempts=[];success=Fraction(0);max_error=0.
    for t in range(1,cap+1):
        prefixes=list(product(bad,repeat=t-1))
        if not prefixes:break
        p=Fraction(len(prefixes)*count,base**t);success+=p
        amps=np.ones((len(prefixes)*count,1),complex)/math.sqrt(base**t)
        joint=amps*bell[None,:]
        retained=joint.T@joint.conj()
        error=float(np.max(abs(retained-float(p)*target)));max_error=max(max_error,error)
        assert error<1e-12
        attempts.append({'attempt':t,'prior_rejected_word_prefixes':prefixes,
            'accepted_words_coherent':good,'each_joint_word_amplitude_squared':str(Fraction(1,base**t)),
            'coarse_outcome_probability':str(p),'reference_crossblock_error':error})
    timeouts=list(product(bad,repeat=cap)) if bad else []
    timeout=Fraction(len(timeouts),base**cap)
    assert success+timeout==1
    expected_bits=sum(Fraction(width)*Fraction(len(bad),base)**n for n in range(cap))
    conditional_bits=sum(Fraction(a['coarse_outcome_probability'])*a['attempt']*width for a in attempts)/success
    return {'count':count,'width':width,'base':base,'cap':cap,'attempt_word_catalog':words,
        'accepted_attempts':attempts,'timeout_word_prefixes':timeouts,'timeout_probability':str(timeout),
        'success_probability':str(success),'expected_bits_per_entered_selector':str(expected_bits),
        'mean_bits_conditioned_on_success':str(conditional_bits),'maximum_bits':width*(cap if bad else 1),
        'coherent_reference_max_error':max_error,
        'conditional_choice_amplitude_squared':str(Fraction(1,count))}


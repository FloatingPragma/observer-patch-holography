import unittest
from itertools import product
import sympy as s

from .model import *
from .routing import certificate as routing_certificate


class MobileFluxTests(unittest.TestCase):
    def test_literal_domain_and_gauge_closure(self):
        base=observer().base
        self.assertEqual(base.dimension,86400)
        for word in base.link_words:
            self.assertIn(word^STAR,base.link_words)
            self.assertEqual(base.law.syndrome(word).bit_count(),2*charge(word))
            for edge in base.law.interior:
                target=base.law.transport(word,edge)[0]
                self.assertIn(target,base.link_words)
                self.assertEqual(charge(target),charge(word))

    def test_full_exchange_covariance_and_charge_transport(self):
        # Complete link columns, with both trivial and nontrivial F-frame signs.
        for left,right in product(observer().base.link_words,repeat=2):
            for occupied in ((),(mode(5,3,0,0),),(mode(5,3,0,0),MODES+mode(5,3,0,0))):
                state=PairState({((left,right),occupied):1})
                swapped=swap_links(state)
                self.assertEqual(swap_links(swapped),state)
                self.assertEqual(exchange(exchange(state),True),state)
                target=swapped.amplitudes[0][0][0]
                self.assertEqual((charge(target[0]),charge(target[1])),(charge(right),charge(left)))
                for cell in (0,1):
                    self.assertEqual(gamma(swapped,cell),swap_links(gamma(state,cell)))

    def test_physical_full_word_and_orbit_swap_agree(self):
        for word in observer().base.original_link_words:
            seed=positive(PairState({((word,0),(mode(5,3,0,0),)):1}))
            self.assertEqual(swap_links(seed),swap_links(seed,False))
            self.assertEqual(exchange(seed),exchange(seed,keep_gauge_bits=False))

    def test_vacuum_preparation_mobility_and_complete_measurement(self):
        raw=PairState({((0,0),()):1});vacuum=positive(raw)
        self.assertEqual(tick(vacuum),vacuum)
        self.assertEqual(insert_pair(insert_pair(vacuum,0),0),vacuum)
        seed=insert_pair(vacuum,0);out=exchange(seed)
        self.assertEqual(out.norm_squared(),1)
        self.assertEqual([(q,v.norm_squared()) for q,v in charge_measurement(out)],
                         [((0,0),0),((0,1),s.Rational(1,2)),((1,0),s.Rational(1,2)),((1,1),0)])
        self.assertEqual(exchange(out),swap_links(seed))
        self.assertEqual([v.norm_squared() for _,v in constraint_branches(raw)],[s.Rational(1,4)]*4)
        local=FockState({(0,()):1})
        self.assertEqual(sum(z.norm_squared() for _,z in observer().face_measurement(local)),1)

    def test_actual_onsite_wave_inverse_and_constraint(self):
        seed=positive(PairState({((1,0),(MODES+mode(7,3,0,0),)):1}))
        out=tick(seed)
        self.assertEqual(out.norm_squared(),1)
        self.assertEqual(tick(out,True),seed)
        for cell in (0,1):self.assertEqual(gamma(out,cell),out)

    def test_protected_object_routing(self):
        receipt=routing_certificate()
        self.assertEqual(receipt['exact_source_route_intertwiner_columns'],20)
        self.assertFalse(receipt['SM_identity_or_GeV_mass_derived'])


if __name__=='__main__':unittest.main()

#!/usr/bin/env python3
import unittest
from h0_ledger import *

def frame(vertices,components,finite):
    return dict(metadata=dict(vertices=str(vertices),components=str(components),merges=str(vertices-components)),
                pairs=Counter({(Decimal(0),Decimal(x)):n for x,n in finite.items()})+Counter({(Decimal(0),Decimal(-1)):components}))

class H0LedgerTests(unittest.TestCase):
    def test_disconnected_then_quotient_then_merge(self):
        # Three original vertices, collapse a connected pair after .5.
        fs=[frame(3,2,{'.1':1}),frame(2,1,{'1.9':1})]
        baseline=Counter({(Decimal(0),Decimal('.1')):1,(Decimal(0),Decimal('1.9')):1,(Decimal(0),Decimal(-1)):1})
        r=ledger(fs,[['0','.5'],['.5','3']],baseline,projected(baseline))
        self.assertEqual(r['status'],'PASS');self.assertEqual([s['new_mergers'] for s in r['steps']],[1,1])

    def test_do_not_add_old_local_mst_bars(self):
        fs=[frame(5,2,{'.2':3}),frame(3,1,{'.1':1,'1.1':1})]
        baseline=Counter({(Decimal(0),Decimal('.2')):3,(Decimal(0),Decimal('1.1')):1,(Decimal(0),Decimal(-1)):1})
        self.assertEqual(ledger(fs,[['0','1'],['1','2']],baseline,projected(baseline))['status'],'PASS')

    def test_conservation_failure_rejects(self):
        fs=[frame(3,2,{'.1':1}),frame(2,1,{'.2':1})]
        baseline=Counter({(Decimal(0),Decimal('.1')):1,(Decimal(0),Decimal('1.9')):1,(Decimal(0),Decimal(-1)):1})
        self.assertEqual(ledger(fs,[['0','.5'],['.5','3']],baseline,projected(baseline))['status'],'H0_NOT_VERIFIED')

    def test_zero_distance_is_known_difference(self):
        baseline=Counter({(Decimal(0),Decimal(0)):1,(Decimal(0),Decimal(-1)):1})
        self.assertEqual(ledger([],[],baseline,projected(baseline))['status'],'KNOWN_DIFFERENCE')

    def test_keeps_survivor_multiplicity(self):
        f=frame(5,3,{'.3':2});baseline=f['pairs']
        self.assertEqual(ledger([f],[['0','1']],baseline,projected(baseline))['status'],'PASS')

    def test_reference_mismatch_not_tolerated(self):
        f=frame(3,3,{});ref=Counter({(0.,1.):1,(0.,float('inf')):2})
        self.assertEqual(ledger([f],[['0','1']],f['pairs'],ref)['status'],'H0_NOT_VERIFIED')

if __name__=='__main__':unittest.main()

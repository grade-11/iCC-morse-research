#!/usr/bin/env python3
import json
import math
from pathlib import Path
import tempfile
import unittest
import selector

class SelectorTests(unittest.TestCase):
    def test_strict_pair_density_and_unit_geometry(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'p.csv';p.write_text('1,0,0\n-1,0,0\n0,1,0\n')
            x=selector.raw_features(p)
            self.assertEqual(x['n'],3)
            self.assertEqual(x['density']['14'],0)
            self.assertEqual(x['density']['15'],2/3)
            self.assertEqual(x['max_unit_norm_error'],0)
            self.assertEqual(x['minimum_distance'],math.sqrt(2))

    def test_noQ_and_final_only_cutoff(self):
        cmd=selector.command('/engine','/points','noQ')
        self.assertIn('ph',cmd);self.assertNotIn('--eps-breaks',cmd)
        for s,ticks in selector.SCHEDULES.items():
            if ticks is None:continue
            self.assertEqual(ticks[-1],20)
            self.assertEqual(len(set(ticks)),len(ticks))
            self.assertTrue(all(0<k<20 for k in ticks[:-1]))
            cmd=selector.command('/engine','/points',s)
            values=cmd[cmd.index('--eps-breaks')+1:cmd.index('--pv-cap-scale')]
            self.assertEqual(list(map(float,values)),[k/10 for k in ticks])

    def test_domain_fallback_and_excluded_schedule(self):
        x=dict(n=100,density={str(k):.2 for k in range(8,19)},minimum_distance=.1,max_unit_norm_error=0)
        model=dict(domain=dict(n=[64,200],density08=[.1,.3]),coefficients=[0]*10,eligible_schedules=['q12','noQ'],fallback='noQ')
        self.assertEqual(selector.choose(x,model)['schedule'],'noQ')
        x['n']=201
        self.assertEqual(selector.choose(x,model)['reason'],'predeclared_fallback')
        self.assertEqual(selector.choose(x,model)['schedule'],'noQ')
        x['n']=100;model['eligible_schedules']=['q12']
        self.assertEqual(selector.choose(x,model)['schedule'],'q12')
        decision=selector.choose(x,model)
        self.assertEqual(decision['G_seconds']['q12'],decision['predictions_s']['noQ']-decision['predictions_s']['q12'])

if __name__=='__main__':unittest.main()

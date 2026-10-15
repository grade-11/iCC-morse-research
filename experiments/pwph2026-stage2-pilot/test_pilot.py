#!/usr/bin/env python3
"""Checks for general schedule parsing and evidence failure/resume behavior."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from pilot import parse_compact, execute, write

class PilotTests(unittest.TestCase):
    def test_many_windows_and_multiplicities(self):
        text=''
        for a,b in [(0,1.2),(1.2,1.4),(1.4,1.6),(1.6,1.8),(1.8,2)]:
            text+=f'window lower bound, window upper bound\n{a}, {b}\ndimension, birth weight, death weight\n3, .9, -1\n'
        text+='[adaptive] next=-1\nrun time = 10\n'
        d,ignored,w=parse_compact(text,'piecewise')
        self.assertEqual(ignored,4);self.assertEqual(sum(d[2].values()),1);self.assertEqual(len(w),5)
        text=text.replace('3, .9, -1\n[adaptive]','3, .9, -1\n3, .9, -1\n[adaptive]')
        self.assertEqual(sum(parse_compact(text,'piecewise')[0][2].values()),2)

    def test_reject_incomplete_and_malformed(self):
        for text in ['run time = 0','window lower bound, window upper bound\n0, 1.2\nrun time = 1',
                     'window lower bound, window upper bound\n0, 2\n3, .5, nan\nrun time = 1',
                     'window lower bound, window upper bound\n0, 2\n3, .5, 1\nBOGUS\nrun time = 1']:
            with self.assertRaises(ValueError):parse_compact(text,'piecewise')

    def test_failure_timeout_memory_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            run=Path(tmp);write(run/'initial-context.json',{'memory_limit_bytes':256*1024**2})
            out,c=execute(run,'error',[sys.executable,'-c','import sys;print("kept");sys.exit(7)'])
            self.assertEqual(c['exit_status'],7);self.assertEqual(c['status'],'ERROR')
            self.assertEqual((out/'stdout.log').read_text(),'kept\n')
            self.assertEqual(execute(run,'error',[sys.executable,'-c','import sys;print("kept");sys.exit(7)'])[0],out)
            out,c=execute(run,'timeout',[sys.executable,'-c','import time;time.sleep(2)'],timeout=.1)
            self.assertEqual(c['status'],'TIMEOUT');self.assertNotEqual(c['exit_status'],0)
            write(run/'initial-context.json',{'memory_limit_bytes':1})
            out,c=execute(run,'memory',[sys.executable,'-c','import time;time.sleep(2)'])
            self.assertEqual(c['status'],'MEMORY_LIMIT')
            self.assertEqual(len((run/'index.jsonl').read_text().splitlines()),3)

if __name__=='__main__':unittest.main()

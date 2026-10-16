#!/usr/bin/env python3
"""Check that process peak RSS survives execv; not a PH benchmark."""
import json
import os
import resource
import sys

buffer=bytearray(64*1024**2)
print(json.dumps(dict(before_exec_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)),flush=True)
os.execv(sys.executable,[sys.executable,'-c',"import json,resource;print(json.dumps(dict(after_exec_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))"])

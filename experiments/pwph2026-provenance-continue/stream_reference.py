#!/usr/bin/env python3
"""Explicit F2 append-column laboratory reference, NOT native iCC continuation.
Pre-enumeration to final cutoff is charged and disclosed; no timing comparison
with iCC. Quotient is explicitly unsupported. No G/controller training uses it.
"""
import math,itertools,json,hashlib
from collections import Counter
from pathlib import Path

def cells(points,maxdim,cutoff):
    n=len(points);dist={}
    for i,j in itertools.combinations(range(n),2):
        total=0.
        for a,b in zip(points[i],points[j]):delta=float(a)-float(b);total+=delta*delta
        dist[i,j]=math.sqrt(total)
    out=[]
    for d in range(maxdim+1):
        for s in itertools.combinations(range(n),d+1):
            w=max((dist[a,b] for a,b in itertools.combinations(s,2)),default=0.)
            if w<cutoff:out.append((w,d,sum(math.comb(v,k+1) for k,v in enumerate(s)),s))
    return sorted(out)

class StreamF2:
    def __init__(self,ordered,maxdim):
        self.stream=ordered;self.maxdim=maxdim;self.cursor=0;self.ids={};self.reduced={};self.owner={};self.alive={};self.finite=[];self.history=[];self.last=0.
    def advance(self,eps):
        if eps<=self.last:raise ValueError('Checkpoint must increase')
        before=self.cursor
        while self.cursor<len(self.stream) and self.stream[self.cursor][0]<eps:
            w,d,bindex,s=self.stream[self.cursor];idx=self.cursor;self.ids[s]=idx
            column={self.ids[s[:k]+s[k+1:]] for k in range(len(s))} if d else set()
            while column and max(column) in self.owner:column^=self.reduced[self.owner[max(column)]]
            if column:
                pivot=max(column);self.owner[pivot]=idx;self.reduced[idx]=column;birth,dim=self.alive.pop(pivot);self.finite.append((dim,birth,w))
            else:self.alive[idx]=(w,d)
            self.cursor+=1
        self.last=eps
        self.history.append(dict(eps=eps,processed_new_columns=self.cursor-before,total_processed_columns=self.cursor,retained_reduced_columns=len(self.reduced),live_at_checkpoint=len(self.alive)))
        return self.diagrams()
    def diagrams(self,include_zero_h0=True):
        bars={d:Counter() for d in range(self.maxdim)}
        for d,b,e in self.finite:
            if d in bars and (b<e or (d==0 and include_zero_h0)):bars[d][b,e]+=1
        for b,d in self.alive.values():
            if d in bars:bars[d][b,-1.]+=1
        return bars
    def quotient(self,*args):raise NotImplementedError('No quotient/label-remapping or native Morse-state migration supported; refuse explicitly')

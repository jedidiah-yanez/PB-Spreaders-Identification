#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
    Copyright (c) 2021, Jedidiah Yanez-Sierra, Cinvestav-Guadalajara
    All rights reserved.
    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:
    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following
      disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
      following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the name of Cinvestav-Guadalajara nor the names of its contributors may be used to endorse or
      promote products derived from this software without specific prior written permission.
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import time 
import numpy as np
from math import log
from igraph import VertexClustering


def HybridRank(g, attr=None, directed=True):
    if directed:
        mode = 'IN'
    else:
        mode = 'ALL'
    
    g.vs['ECC'] = g.eccentricity(mode=mode)
    if attr==None:
        attr='ICC'
        g.vs[attr] = 0
        g.vs['KS'] = g.coreness(mode=mode)
        
        for v in g.vs:
            for w in v.neighbors(mode=mode):
                v[attr] += w['KS']
        del g.vs['KS']
    
    g.vs['HC'] = 0
    for v in g.vs:
        v['HC'] = v['ECC'] * v[attr]
    return g.vs['HC']


def voteRank(g, directed=True, f=None, r=1):
    g.vs['s'] = 0
    g.vs['va'] = 1
    g.vs['VR'] = 0
    if f == None:
        f = 1.0 / np.average(g.degree())
    
    spr = []
    spr_set = set()
    if directed:
        mode = 'IN'
    else:
        mode = 'ALL'
    
    times = list()
    while(len(spr) < r):
        max_v = None
        max_s = 0
        
        start_time = time.time()
        #vote phase, each node is voted by its incoming neighbors
        g.vs['s'] = 0 
        for v in g.vs:
            neighs = v.neighbors(mode=mode)
            for w in neighs:
                v['s'] += w['va']
        
            #keep track of most voted node
            if v['s'] > max_s and v not in spr_set:
                max_s = v['s'] 
                max_v = v
        if max_v == None:
            break
        spr.append((max_v,max_v['s']))
        spr_set.add(max_v)
        max_v['va'] = 0
        max_v['s'] = 0
        
        #Weaken the voting ability of neighbors of max_v
        neighs = max_v.neighbors(mode=mode)
        for w in neighs:
            w['va'] = max(w['va'] - f, 0)
    
    del g.vs['s']
    del g.vs['va']
    for v, val in spr:
        v['VR'] = val
    return g.vs['VR']


def IKS(g, directed=True):
    k_N = float(sum(g.vs.degree()))
    g.vs['iks_kcore'] = g.coreness(mode='All')
    g.vs['iks_e'] = 0
    for v in g.vs:
        e_v = 0
        for w in v.neighbors():
            I_w = w.degree() / k_N
            e_v += I_w + log(I_w)
        v['iks_e'] = -e_v


def IKS_Select(g, N):
    g.vs['IKS']=0
    clust = VertexClustering.FromAttribute(g, 'iks_kcore')
    subg_clust = clust.subgraphs()
    subg_clust.sort(key=lambda x:x.vs[0]['iks_kcore'], reverse=True) #sort by k-shell
    subg_clust = [list(sub_g.vs) for sub_g in subg_clust] #change from graph to list of vertexes
    for k_vs in subg_clust:
        k_vs.sort(key=lambda x:x['iks_e']) 
    N = min(N,len(g.vs))
    kcore_size = len(subg_clust)
    k = 0
    spr = []
    while len(spr) < N:
        j = k % kcore_size
        k += 1
        if(len(subg_clust[j])==0):
            continue
        else:
            v = subg_clust[j].pop()
            v = g.vs.find(name=v['name'])
            v['IKS'] = N + 1 - len(spr)
            spr.append( v )        
    return spr


def s_k_w_out(v,w):
    neigh_v = set(v.neighbors())
    neigh_v.add(v)
    neigh_w = set(w.neighbors())
    k_w_out = len(neigh_w-neigh_v)
    return k_w_out


def sc_core(g,alpha):
    g.es['w_score'] = [{}]*len(g.es)
    for e in g.es:
        v = g.vs[e.source]
        w = g.vs[e.target]
        D_vw = g.get_all_simple_paths(v.index,w.index,cutoff=2, mode='ALL')
        D_vw = len(D_vw)-1 #sum( [1 for x in D_vw if len(x)==3] )
        
        k_w_out = s_k_w_out(v,w)
        e['w_score']["%d,%d"%(v.index, w.index)] = 1 + k_w_out * (1 + D_vw/float(2**2))**alpha
        
        w = g.vs[e.source]
        v = g.vs[e.target]
        k_w_out = s_k_w_out(v,w)
        e['w_score']["%d,%d"%(v.index, w.index)] = 1 + k_w_out * (1 + D_vw/float(2**2))**alpha
    
    g.vs['sc_score']=0
    for v in g.vs:
        s = 0
        for w in v.neighbors():
            eid = g.get_eid(v, w, directed=False, error=True)
            s += g.es[eid]['w_score']["%d,%d"%(v.index, w.index)]
        v['sc_score'] = s

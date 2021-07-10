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

import numpy as np
import SprModel.Utilities as ut
import operator
from collections import deque, OrderedDict
import time  


def _checkCommunities(g_comm, clusters, policy_tag):
    g = g_comm
    comm_size = len(clusters)
    for i in range(comm_size):
        _inc = 'A'
        sub = clusters.subgraph(i)
        components = list(sub.components())
        if len(components) == 1:
            continue
        else:
            print ("Community: " + sub.vs[0][policy_tag] + 'broken into ' + str(len(components)) + " communities")
            for comp in components:
                for v_sub in comp:
                    v = g.vs.find(name=sub.vs[v_sub]['name'])
                    v[policy_tag] = v[policy_tag][:-1]+_inc+"]"
                _inc = ut.inc_char(_inc)


def __getNeighbors(v, r=3):
    r_max = r
    r=0
    Q1 = deque() 
    Q2 = deque()
    distance = OrderedDict()
    Q1.append(v)
    distance[v] = 0
    for r in range(r_max):
        while len(Q1) != 0:
            _v = Q1.popleft()
            neis = _v.neighbors()
            for w in neis:
                if not distance.has_key(w):
                    Q2.append(w)
                    distance[w] = distance[_v] + 1 
        tmp = Q1
        Q1 = Q2
        Q2 = tmp        
    distance.pop(v)
    return distance


def __G(g, X):
    gX = [None]*len(g.vs)
    for v in g.vs:
        _sum = 0
        neis = __getNeighbors(v)
        for w in neis.keys():
            _sum += float(v[X] * w[X]) / pow(neis[w],2)
        gX[v.index] = _sum
    return gX


def checkPrepareCommunities(g, c_name=None):
    if c_name == None:
        c_name = 'MLC'
    
    if c_name not in g.vs.attribute_names():
        components = g.components()
        components = components.membership
        
        clusters = g.community_multilevel()
        membership = clusters.membership
        membership = ["["+str(components[i])+"-"+str(membership[i])+"]" for i in range(len(membership)) ]
        g.vs[c_name] = membership # get the membership vector
        _checkCommunities(g, clusters, c_name)
    
    communities = {c:0 for c in set(g.vs[c_name])}
    for v in g.vs:
        communities[ v[c_name] ] += 1
    sorted_comm = sorted(communities.items(), key=operator.itemgetter(1), reverse=False)
    comm, _ = zip(*sorted_comm)
    
    if 'gks_'+c_name not in set(g.vs.attribute_names()):
        g.vs['gks_'+c_name] = None
        for c in comm:
            topics = g.vs.select(lambda x:x[c_name] == c)#['label']
            
            sub_g = g.subgraph(topics.indices)
            ks = sub_g.coreness(mode='ALL')
            sub_g.vs['ks'] = ks
            
            gks = __G(sub_g, 'ks')
            sub_g.vs['gks_'+c_name] = gks
            
            for v in sub_g.vs:
                w = g.vs.select(name=v['name'])
                w['gks_'+c_name] = v['gks_'+c_name]
    return c_name, 'gks_'+c_name


def deflatSpreaders(spr):
    final_local = []
    rows = len(spr)
    i = 0
    while (True):
        if len(spr[i]) > 0:
            tmp = spr[i].popleft()
            final_local.append( tmp )
            i += 1
        else:
            del spr[i]
            rows -= 1
            if rows == 0:
                break
        
        if i % rows == 0:
            i = 0
    return final_local


def detectSpreaders(g, max_spreaders, alg, commSet = None):
    gks_name = 'gks_'+alg
    comms_set = set(g.vs[alg])
    
    if commSet:
        comm = commSet[0]
        communities = commSet[1]
    else:
        communities = {c:0 for c in comms_set }
        for v in g.vs:
            if max_spreaders < len(comms_set):
                communities[ v[alg] ] = max(communities[ v[alg] ], v[gks_name])
            else:
                communities[ v[alg] ] += 1
        
        sorted_comm = sorted(communities.items(), key=operator.itemgetter(1), reverse=False)
        comm, _ = zip(*sorted_comm)
    
    V = len(g.vs)
    spread_done = 0
    nodes_done = 0
    g.vs['PRP'] = 0
    spr = []
    
    if max_spreaders < len(comm):
        comm = list(comm)
        comm.reverse()
        comm = comm[:max_spreaders]
        
        for c in comm: 
            topics = g.vs.select(lambda x:x[alg] == c)
            vertices = list(topics)
            vertices.sort(key = lambda v: v[gks_name], reverse=True)
            spr.append( deque(vertices[:1]) )
    else:
        for c in comm:
            topics = g.vs.select(lambda x:x[alg] == c)
            vertices = list(topics)
            vertices.sort(key = lambda v: v[gks_name], reverse=True)
            k = communities[c] * (max_spreaders-spread_done) / float(V-nodes_done)
            if k < 1:
                k=1
            else:
                k = int(round(k,0))
            spread_done += k
            nodes_done += communities[c]
            vertices = vertices[:k]
            
            spr.append( deque(vertices) )
        spr.reverse()
    spr2 = deflatSpreaders(spr)
    i = max_spreaders
    for v in spr2:
        v['PRP'] = i
        i-=1 
    return spr2

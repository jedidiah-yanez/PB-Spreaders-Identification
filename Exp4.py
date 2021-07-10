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

import os
import operator
import argparse
import matplotlib; 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import SprModel.Utilities as ut
from SprModel.SIR import sim_SIR
from Algorithms.PBSI import __G, checkPrepareCommunities
from Algorithms.SpreadersAlgs import voteRank, HybridRank, IKS, IKS_Select, sc_core


def plot(fig, data, keys=[], names=[], title="", xylabels=['',''], flag=True, vline=None):
    ymax = 0
    for k in keys:
        fig.plot(data[k][0], data[k][1], data[k][2], markersize=5)
        ymax = max(ymax, max(data[k][1]))
    if flag:
        fig.set_title(title, fontsize=12)
        fig.set_xlabel(xylabels[0], fontsize=10)
        fig.set_ylabel(xylabels[1], fontsize=10)
        fig.legend(names, fontsize=10, loc=7)
        fig.grid()
    if vline:
        fig.vlines(vline, 0, ymax, colors='k', linestyles='dashed')


def calculateMetrics(g, metrics, spr_size, weight_attr=None):
    if 'BET' in metrics:
        return g.betweenness(directed=False, weights=weight_attr)
    
    if 'CLO' in metrics:
        return g.closeness(weights=weight_attr)
    
    if 'EIG' in metrics:
        return g.eigenvector_centrality(directed=False, weights=weight_attr)
    
    if 'PGR' in metrics:
        return  g.pagerank(directed=False, weights=weight_attr)
    
    if 'STR' in metrics:
        return g.strength(weights=weight_attr)
    
    if 'DEG' in metrics: 
        return g.vs.degree()
    
    if 'ECC' in metrics:
        return g.eccentricity()
    
    if 'VR' in metrics:
        return voteRank(g, directed=False, r = spr_size)
    
    if 'HC' in metrics:
        return HybridRank(g, directed=False)
    
    if 'PRP' in metrics:
        ks = g.coreness(mode='ALL')
        g.vs['ks'] = ks
        return __G(g, 'ks')
    
    if 'IKS' in metrics:
        IKS(g, directed=False)
        IKS_Select(g, spr_size)
        return g.vs['IKS']
    
    if 'sc_score' in metrics:
        sc_core(g,0.5)
        return g.vs['sc_score']


def computeMetric(g, c_name, metric, spr_size):
    if c_name == None:
        c_name = 'MLC'
    
    communities = {c:0 for c in set(g.vs[c_name])}
    for v in g.vs:
        communities[ v[c_name] ] += 1
    sorted_comm = sorted(communities.items(), key=operator.itemgetter(1), reverse=False)
    comm, _ = zip(*sorted_comm)
    
    if metric+'_'+c_name not in set(g.vs.attribute_names()):
        g.vs[metric + '_' + c_name] = None
        for c in comm:
            topics = g.vs.select(lambda x:x[c_name] == c)
            sub_g = g.subgraph(topics.indices)
            m_values = calculateMetrics(sub_g, [metric], spr_size, weight_attr=None)
            sub_g.vs[metric +'_' + c_name] = m_values
            for v in sub_g.vs:
                w = g.vs.select(name=v['name'])
                w[metric +'_' + c_name] = v[metric +'_' + c_name]
    return metric +'_' + c_name, c_name


def detectSpreaders(g, alg, metric):
    gks_name = metric
    comms_set = set(g.vs[alg])
    spr = []
    for c in comms_set: 
        topics = g.vs.select(lambda x:x[alg] == c)
        vertices = list(topics)
        vertices.sort(key = lambda v: v[gks_name], reverse=True)
        spr.append( vertices[:1][0] )
    return spr


def FSS_Experiment(g, g_path, mu, beta, mc, metrics, spr_size):
    graph_name = 'FSS_'+g_path[g_path.rfind("/")+1:g_path.rfind(".")]+"_spr"+str(spr_size)
    if 'MLC' not in g.vs.attribute_names():
        checkPrepareCommunities(g, 'MLC')
    nodes = list(g.vs)
    res = {metric:[None,[],'.-'] for metric in metrics if metric != "PRP"}
    x = [beta]
    for metric in metrics: 
        if metric != 'PRP':
            print ("Processing Metric: " + metric)
            
            if metric == 'IKS':
                IKS_Select(g, spr_size)
            
            if  metric == 'VR':
                voteRank(g, directed=False, r = spr_size)
            
            nodes.sort(key=lambda x: x[metric], reverse=True)
            seeds = nodes[:spr_size]
            for b in x:
                print ("x: " + str(b))
                spr,_,_ = sim_SIR(g, seeds, b, mu=mu, mc=mc, verbose=False)
                res[metric][1].append(spr)
            res[metric][0] = x
            print (res[metric])
        print ("Processing v2 of Metric: " + metric)
        metric2,_ = computeMetric(g, None, metric, 1)
        seeds = detectSpreaders(g, 'MLC', metric2)
        del g.vs[metric2]
        res[metric2]=[None,[],'.-']
        for b in x:
            print ("x: " + str(b))
            spr,_,_ = sim_SIR(g, seeds, b, mu=mu, mc=mc, verbose=False)
            res[metric2][1].append(spr)
        res[metric2][0] = x
        print (res[metric2])
    
    metrics = sorted(res.keys())
    fig = plt.figure(figsize=(6,6)).add_subplot(111)
    plot(fig, res, keys=metrics, names=metrics, vline = None)
    plt.savefig('plots/'+graph_name+"_mc"+str(mc)+".pdf", dpi=300) 
    ut.save_json('plots/'+graph_name+"_mc"+str(mc)+".json", res)
    print ("number of spreaders = " + str(spr_size))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Graph instance file', dest='graph_path')
    parser.add_argument('-beta', type=float, help='spreading probability', dest='beta', default=0.5)
    parser.add_argument('-mu', type=float, help='spreading probability', dest='mu', default=1.0)
    parser.add_argument('-mc', type=int, help='monte carlo simulations', dest='mc', default=32)
    parser.add_argument('-metrics', nargs='+', type=str, default=['PRP'], dest='metrics')
    # parser.add_argument('-metrics', nargs='+', type=str, default=['PRP', 'BET', 'CLO', 'DEG', 'VR', 'HC'], dest='metrics')
    
    args = parser.parse_args()
    graph_path = args.graph_path
    beta = args.beta
    mu = args.mu
    mc = args.mc
    metrics = args.metrics
    
    if graph_path != None:
        if os.path.exists(graph_path):
            g = ut.openGraph(graph_path)
            if 'name' not in g.vs.attribute_names():
                g.vs['name'] = ['a'+str(v.index) for v in g.vs]
            spr_size = len(set(g.vs['MLC']))
            # metrics = ['BET', 'CLO', 'DEG',  'HC', 'sc_score', 'IKS', 'VR', 'PRP']
            FSS_Experiment(g, graph_path, mu, beta, mc, metrics, spr_size)
        else:
            print ('graph file can not found')
    else: 
        print ('a graph file must be provided')
    print("done...")
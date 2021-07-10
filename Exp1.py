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
import time
import argparse
import numpy as np 
import matplotlib; 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import SprModel.Utilities as ut
from SprModel.SIR import sim_SIR
from Algorithms.PBSI import detectSpreaders, checkPrepareCommunities
from Algorithms.SpreadersAlgs import voteRank, HybridRank, IKS, IKS_Select, sc_core


def plot(fig, data, keys=[], names=[], title="", xylabels=['',''], flag=True, vline=None):
    ymax = 0
    for k in keys:
        fig.plot(data[k][0], data[k][1], data[k][2], markersize=1)
        ymax = max(ymax, max(data[k][1]))
    if flag:
        fig.set_title(title, fontsize=12)
        fig.set_xlabel(xylabels[0], fontsize=10)
        fig.set_ylabel(xylabels[1], fontsize=10)
        fig.legend(names, fontsize=10, loc=7)
        fig.grid()
    if vline:
        fig.vlines(vline, 0, ymax, colors='k', linestyles='dashed')


def calculateMetrics(g,  g_path, metrics, spr_size, weight_attr=None):
    g_attr = set(g.vs.attribute_names())
    times = {metric:0 for metric in metrics}
    if 'BET' in metrics and 'BET' not in g_attr:
        print ('Computing Betweenness, this process could be time consuming...')
        start_time = time.time()
        g.vs['BET'] = g.betweenness(directed=False, weights=weight_attr)
        times['BET'] = time.time() - start_time
        ut.saveGraph(g, g_path)
    
    if 'CLO' in metrics and 'CLO' not in g_attr:
        print ('Computing Closeness...')
        start_time = time.time()
        g.vs['CLO'] = g.closeness(weights=weight_attr)
        times['CLO'] = time.time() - start_time
        ut.saveGraph(g, g_path)
    
    if 'DEG' in metrics and 'DEG' not in g_attr: 
        print ('Computing Degree...')
        start_time = time.time()
        g.vs['DEG'] = g.vs.degree()
        times['DEG'] = time.time() - start_time
        ut.saveGraph(g, g_path)
    
    if 'HC' in metrics and 'HC' not in g_attr:
        print ('Computing HybridRank, this process could be time consuming...')
        start_time = time.time()
        HybridRank(g, directed=False)
        times['HC'] = time.time() - start_time
        ut.saveGraph(g, g_path)
    
    if 'sc_score' in metrics and 'sc_score' not in g_attr:
        print ('Computing SC Score...')
        start_time = time.time()
        sc_core(g,0.5)
        times['sc_score'] = time.time() - start_time
        ut.saveGraph(g, g_path)
    
    if 'IKS' in metrics: 
        print ('Computing IKS...')
        start_time = time.time()
        IKS(g, directed=False)
        times['IKS'] = time.time() - start_time
        ut.saveGraph(g, g_path)
        start_time = time.time()    
        IKS_Select(g, spr_size)
        times['IKS2'] = time.time() - start_time
    
    if 'VR' in metrics: 
        print ('Computing VoteRank, this process could be time consuming...')
        start_time = time.time()
        voteRank(g, directed=False, r = spr_size)
        times['VR'] = time.time() - start_time
    
    if 'PRP' in metrics: 
        print ('Computing PBSI...')
        start_time = time.time()
        c_name,_ = checkPrepareCommunities(g)
        start2_time = time.time()
        ut.saveGraph(g, g_path)
        detectSpreaders(g, spr_size, c_name, commSet = None)
        final_time = time.time()
        times['PRP'] = final_time - start_time
        times['PRP2'] = final_time - start2_time
    return times



def FSS_Experiment(g, g_path, mu, mc, metrics, spr_size):
    graph_name = 'FSS_'+graph_path[graph_path.rfind("/")+1:graph_path.rfind(".")]+"_spr"+str(spr_size)
    times = calculateMetrics(g, g_path, metrics, spr_size, weight_attr=None)
    
    nodes = list(g.vs)
    res = {metric:[None,[],'.-'] for metric in metrics}
    for metric in metrics: 
        print ("Processing Metric: " + metric)
        nodes.sort(key=lambda x: x[metric], reverse=True)
        seeds = nodes[:spr_size]
        x = np.linspace(.01,.15,15)
        for beta in x:
            spr,_,_ = sim_SIR(g, seeds, beta, mu=mu, mc=mc, verbose=False)
            res[metric][1].append(spr)
        res[metric][0] = x
        # print (res[metric])
        
    degrees = np.array(g.degree())
    degrees_2 = np.power(degrees,2)
    epidemic_threshold = np.average(degrees) / np.average(degrees_2)
    
    fig = plt.figure(figsize=(6,6)).add_subplot(111)
    plot(fig, res, keys=metrics, names=metrics, vline = epidemic_threshold)
    plt.savefig('plots/'+graph_name+"_mc"+str(mc)+".pdf", dpi=300) 
    ut.save_json('plots/'+graph_name+"_mc"+str(mc)+".json", res)
    ut.save_json('plots/'+graph_name+"_mc"+str(mc)+"_times.json", times)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Graph instance file', dest='graph_path')
    parser.add_argument('-mu', type=float, help='spreading probability', dest='mu', default=1.0)
    parser.add_argument('-mc', type=int, help='monte carlo simulations', dest='mc', default=32)
    parser.add_argument('-spr', type=int, help='Number of Spreaders', dest='spr_size', default=50)
    parser.add_argument('-metrics', nargs='+', type=str, default=['PRP'], dest='metrics')
    
    args = parser.parse_args()
    graph_path = args.graph_path
    mu = args.mu
    mc = args.mc
    spr_size = args.spr_size
    metrics = args.metrics
    
    if graph_path != None:
        if os.path.exists(graph_path):
            g = ut.openGraph(graph_path)
            if 'name' not in g.vs.attribute_names():
                g.vs['name'] = ['a'+str(v.index) for v in g.vs]
            # metrics = ['BET', 'CLO', 'DEG',  'HC', 'sc_score', 'IKS', 'VR', 'PRP']
            FSS_Experiment(g, graph_path, mu, mc, metrics, spr_size)
        else:
            print ('graph file can not found')
    else: 
        print ('a graph file must be provided')
    print("done...")

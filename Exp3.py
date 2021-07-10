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
import argparse
import matplotlib; 
matplotlib.use('agg')
import numpy as np 
import matplotlib.pyplot as plt
import SprModel.Utilities as ut
from Algorithms.PBSI import detectSpreaders
from Algorithms.SpreadersAlgs import voteRank, HybridRank, IKS_Select 


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


# -i graphs/PGP.graphml -mc 100 -mu 1.0 -betas 0.07 0.10 0.13 -spr 1 100 20

def LS(g, seeds):
    if len(seeds) == 1:
        return 1
    else:
        paths = []
        while(len(seeds) > 0):
            spr = seeds.pop()
            for w in seeds:
                paths.append( g.shortest_paths(source=spr['name'], target=w['name'], weights=None, mode='ALL') )
        paths = np.array(paths + paths)
        return np.average(paths)
    


def avg_sht_pat(g, g_path, metrics, spr_sizes):
    c_name = 'MLC'
    nodes = list(g.vs)
    graph_name = 'JS_'+g_path[g_path.rfind("/")+1:g_path.rfind(".")]
    res = {metric:[None,[],'.-'] for metric in metrics}
    
    x = np.linspace(spr_sizes[0], spr_sizes[1], spr_sizes[2])
    for metric in metrics:
        for spr_size in x:
            spr_size = int(spr_size)
            print ("Processing Metric: " + metric + "    with spr_size: " + str(spr_size))
            
            if metric == 'PRP': 
                detectSpreaders(g, spr_size, c_name, commSet = None)
            
            if  metric == 'VR':
                print ('Computing VoteRank with r = ' + str(spr_size))
                voteRank(g, directed=False, r = spr_size)
            
            if metric == 'IKS':
                print ('Computing IKS with spr = ' + str(spr_size))
                IKS_Select(g, spr_size)
            
            nodes.sort(key=lambda x: x[metric], reverse=True)
            seeds = nodes[:spr_size]
            avg_paths = LS(g, seeds)
            res[metric][1].append(avg_paths)
        res[metric][0] = x
        print (res[metric])
    
    fig = plt.figure(figsize=(6,6)).add_subplot(111)
    plot(fig, res, keys=metrics, names=metrics, vline = None)
    plt.savefig('plots/'+graph_name+"_spr_"+str(spr_sizes[0])+"-"+str(spr_sizes[1])+".pdf", dpi=300) 
    ut.save_json('plots/'+graph_name+"_spr_"+str(spr_sizes[0])+"-"+str(spr_sizes[1])+".json", res)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Graph instance file', dest='graph_path')
    parser.add_argument('-spr', nargs=3, type=int, default=[1,100,20], dest='spr_sizes')
    parser.add_argument('-metrics', nargs='+', type=str, default=['PRP'], dest='metrics')
    
    args = parser.parse_args()
    graph_path = args.graph_path
    spr_sizes = args.spr_sizes
    metrics = args.metrics
    
    if graph_path != None:
        if os.path.exists(graph_path):
            g = ut.openGraph(graph_path)
            if 'name' not in g.vs.attribute_names():
                g.vs['name'] = ['a'+str(v.index) for v in g.vs]
            # metrics = ['BET', 'CLO', 'DEG',  'HC', 'sc_score', 'IKS', 'VR', 'PRP']
            avg_sht_pat(g, graph_path, metrics, spr_sizes)
        else:
            print ('graph file can not found')
    else: 
        print ('a graph file must be provided')
    print("done...")
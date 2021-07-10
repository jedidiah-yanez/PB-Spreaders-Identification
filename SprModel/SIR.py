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

class SIR:
    def __init__(self, graph, beta=0.3, mu=0.08, seed=[]):
        self.g = graph#.copy()   
        self.N = len(graph.vs)
        self.BETA = beta  # infection rate
        self.MU = mu      # recovery rate
        self.keys=['t', 'S', 'I', 'R', 'SI', 'IR']
        self.data_normalized = { k: [] for k in self.keys}
        self.attack_rate=0.0
        self.peak_time=0.0
        self.seed=seed
        
    def run(self, num_steps=1, random_seed=None): 
        self.data = { k: [] for k in self.keys }
        
        if random_seed != None:
            print ('Working with fixed random_seed: ' + str(random_seed))
            np.random.seed(random_seed)
        
        if not len(self.seed):
            self.seed = list(np.random.choice(self.g.vs, size=1, replace=False))
            
        # initialize sets of S/I/R nodes
        I = set(self.seed)
        S = set(self.g.vs).difference(I) 
        R = set()
        t = 0
        
        SI = set(self.seed) #S->I transition
        IR = set()          #I->R transition
        
        while True:
            # generator logic: yield current status every num_steps iterations
            if t % num_steps == 0:
                data=[t, len(S), len(I), len(R), len(SI), len(IR)]
                for i in range(len(self.keys)):
                    k=self.keys[i]
                    self.data[k].append(data[i])
            # stop when there are no infectious nodes left
            
            if not len(I):
                # compute attack rate
                self.attack_rate=self.data["R"][-1]/self.N
                # compute peak time
                ii=np.argmax(self.data["I"])
                self.peak_time=self.data["t"][ii]
                self.data_normalized = {k: [float(v) / self.N for v in self.data[k]] for k in self.data if k not in ["t"]}
                break
            
            SI = set()
            IR = set()
            
            # loop over neighbors of infectious nodes
            tmpI = set(I)
            if random_seed != None: 
                tmpI = sorted(tmpI)
            for i in tmpI: #set(I):
                # TRANSMISSION
                tmpS = S.intersection(i.neighbors())
                if random_seed != None: 
                    tmpS = sorted(tmpS)
                for j in tmpS: #S.intersection(i.neighbors()):
                    # Bernoulli sampling
                    if np.random.uniform() < self.BETA:
                        S.remove(j)
                        I.add(j)
                        SI.add(j)
                        
                # RECOVERY
                # Bernoulli sampling
                if np.random.uniform() < self.MU:
                    I.remove(i)
                    R.add(i)
                    IR.add(i)
            t += 1
        return self.data, self.data_normalized
    
    def getSeed(self):
        return self.seed
    def getData(self):
        return self.data
    def getDataNormalized(self):
        return self.data_normalized
    def getAttackRate(self):
        return self.attack_rate
    def getPeakTime(self):
        return self.peak_time


def sim_SIR(g, seed=[], beta=0.3, mu=0.08, num_steps=1, random_seed=None, mc = 32, verbose=False):
    '''
    perform multiple realization of a SIR simulation
    '''
    spread=[]
    spread_norm=[]
    data_avg = None 
    
    sir = SIR(g, beta, mu, seed)
    for i in range(mc):
        data, data_norm = sir.run(num_steps, random_seed)
        spread.append(data["R"][-1])         
        spread_norm.append(data_norm["R"][-1])
        if i==0:
            data_avg = data
        else:
            for k in data.keys():
                minL = min(len(data[k]), len(data_avg[k]))
                data_avg[k] = np.add(data[k][:minL],data_avg[k][:minL])
        if verbose: print (len(data["R"]), data["R"][-1])
    for k in data.keys():
        data_avg[k] /= mc
    return np.mean(spread), np.mean(spread_norm), data_avg


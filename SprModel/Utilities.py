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

from __future__ import print_function
import scipy.stats as sp
from igraph import Graph
import re 
import sys
import imp
import json
import numpy as np


def saveGraph(g, path, verbose=True):
    if verbose:
        print(g.summary())    
        print("Writing to GRAPHML...")
    graph_file_handler = open(path, mode = "w")        
    g.write_graphml(graph_file_handler)
    graph_file_handler.close()   
        
def openGraph(path):
    graph_file_handler = open(path, mode = 'r')
    g = Graph().Read_GraphML(graph_file_handler)
    del g.vs['id']
    graph_file_handler.close()
    return g


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, 
            np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def save_json(name, data):
    with open(name, "w") as f:
        json.dump(data, f, cls=NumpyEncoder)


def load_json(name):
    with open(name, "r") as f:
        data = json.load(f)
    return data


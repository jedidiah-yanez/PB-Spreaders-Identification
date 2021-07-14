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
import os
import argparse
from igraph import Graph

def error():
    print("Unsupported Input Graph") 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    txt = 'Format of the input graph [EdgeList, NCol, Pajek, mtx]. ' \
        '*EdgeList reads an edge list from a file and creates a graph based on it. \n' \
        '\t Please note that the vertex indices are zero-based. A vertex of zero degree will '\
        '\t be created for every integer that is in range but does not appear in the edgelist. \n' \
        '*NCol. Reads an .ncol file used by LGL. \n' \
        '*Pajek. Reads a Pajek format file and creates a graph based on it. \n\n' \
        'Note: in some cases, first comments lines must be removed prior to read it with this tool.'
    
    parser.add_argument('-type', type=str, help=txt, dest='GraphFormat')
    parser.add_argument('-i', help='input Graph file', dest='graph_path')
    args = parser.parse_args()
    graph_path = args.graph_path
    inputFormat = args.GraphFormat
    
    opc = {'ncol':Graph.Read_Ncol, 'pajek':Graph.Read_Pajek, 'edgelist':Graph.Read_Edgelist, 'mtx':Graph.Read_Ncol}
    
    if graph_path != None:
        if os.path.exists(graph_path):
            graph_name = graph_path[graph_path.rfind("/")+1:graph_path.rfind(".")]+".graphml"
            g = opc.get(inputFormat.lower(), error)(graph_path)
            g.write_graphml(graph_name)
            print ("Graph stored as ", graph_name)
        else:
            print ('graph file can not be found')
    else: 
        print ('a graph file must be provided')
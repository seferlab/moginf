import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import gzip
import cPickle
import myutilities as myutil

foldername="traces_real_static_si_edge"
for tracedata in myutil.listdirectories(foldername):
    path1="{0}".format(tracedata)
    if tracedata.startswith("depo"):
       continue
    for spreadparam in myutil.listdirectories(path1):
        path2="{0}/{1}".format(tracedata,spreadparam)
        for filename in myutil.listfiles(path2):
            print filename
            filepath="{0}/{1}".format(path2,filename)
            file=gzip.open(filepath,"rb")
            spreaddata=cPickle.load(file)
            file.close()
            startnodes=set()
            for node in spreaddata.keys():
                curtime=spreaddata[node]["infect"]
                if curtime==0:
                   startnodes.add(curtime)
            assert len(startnodes)==1  

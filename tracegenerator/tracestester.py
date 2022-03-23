import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import gzip
import cPickle
import myutilities as myutil

#clean those traces that have more than 1 start node

if False:
 filepath="traces_real_static_si_edge/enron-email-clustered-146-1146-2.gml/spreadprob_0.5-s2i_expo_2.0/135_4262905.pkl"
 file=open(filepath,"rb")
 spreaddata=cPickle.load(file)
 file.close()
 startnodes=set()
 for node in spreaddata.keys():
    curtime=spreaddata[node]["infect"]
    if curtime==0:
       startnodes.add(node)
 print startnodes      
 assert len(startnodes)==1
 print startnodes
 print spreaddata.keys()
 print len(spreaddata.keys())
 exit(1)

startfolder="traces_real_static_si_edge"           
for tracedata in myutil.listdirectories(startfolder):
    path1="{0}/{1}".format(startfolder,tracedata)
    if tracedata.startswith("depo"):
       continue
    for spreadparam in myutil.listdirectories(path1):
        print spreadparam
        count=0
        path2="{0}/{1}".format(path1,spreadparam)
        avglen=0
        avgcount=0
        for filename in myutil.listfiles(path2):
            filepath="{0}/{1}".format(path2,filename)
            mytraces="-1"
            #print filepath
            try:
             infile=gzip.open(filepath,"rb")
             spreaddata=cPickle.load(infile) 
             infile.close()
            except IOError:
             #print "normal trace" 
             infile=open(filepath,"rb")
             spreaddata=cPickle.load(infile) 
             infile.close()
            
            startnodes=set()
            alltimes=set()
            for node in spreaddata.keys():
                curtime=spreaddata[node]["infect"]
                alltimes.add(curtime)
                if curtime==0:
                   startnodes.add(node)
            avglen+=len(spreaddata.keys())
            avgcount+=1
            #print startnodes
            if len(startnodes)!=1:
               print startnodes
               print filepath
               print filename
               print len(spreaddata.keys())
               print alltimes
               exit(1)
               count+=1 
            #assert len(startnodes)==1
        avglen/=float(avgcount)       
        print "avg trace count {0}".format(avglen)       
        print "count is {0}".format(count)

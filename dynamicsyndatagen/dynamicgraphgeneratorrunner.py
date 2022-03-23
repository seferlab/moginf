#We will be using graoh6 or sparse6 for storing undirected graphs when node names are numbers(they don't have labels). When there are labels on the nodes and directed graphs, we will use gml file format and then compress it by gzip since networkx can also read compressed graph files with extension .gz
#We will need new random graph growth model implementations for this type of evolution
import networkx as nx
import os
import math
import numpy as np
import itertools
import random
import sys
from optparse import OptionParser
from collections import deque
import myutilities as myutil


def algoparamassigner(nodenum):
    algosparameters={}

    #Random graph models based evolution(node number is always monotonically increasing)
    paramdict1={}
    paramdict1["qcon"]=0.5
    paramdict1["qmod"]=0.6
    paramdict2={}
    paramdict2["qcon"]=0.5
    paramdict2["qmod"]=0.8
    paramdict3={}
    paramdict3["qcon"]=0.3
    paramdict3["qmod"]=0.6
    paramdict4={}
    paramdict4["qcon"]=0.4
    paramdict4["qmod"]=0.8
    algosparameters["dmc"]=[paramdict1,paramdict2,paramdict3,paramdict4]

    paramdict1={}
    paramdict1["k"]=int(nodenum/50)
    paramdict2={}
    paramdict2["k"]=int(nodenum/10)
    algosparameters["lpa"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["p"]=0.05
    paramdict2={}
    paramdict2["p"]=0.1
    algosparameters["rds"]=[paramdict1,paramdict2]
    
    paramdict1={}
    paramdict1["p"]=0.2
    paramdict2={}
    paramdict2["p"]=0.4
    algosparameters["undirected_ff"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["pbackward"]=0.1
    paramdict1["pforward"]=0.3
    paramdict2={}
    paramdict2["pbackward"]=0.2
    paramdict2["pforward"]=0.4
    algosparameters["directed_ff"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["p1"]=0.4
    paramdict1["p2"]=0.8
    paramdict1["p3"]=0.8
    paramdict1["p4"]=0.4
    paramdict2={}
    paramdict2["p1"]=0.5
    paramdict2["p2"]=0.7
    paramdict2["p3"]=0.7
    paramdict2["p4"]=0.4
    algosparameters["kronecker"]=[paramdict1,paramdict2]
    
    #Smoothing based evolution
    #The parameters of smoothing based evolution will be controlled inside the algorithm. They have a fixed distribution assigned.
    paramdict1={}
    paramdict1["intervaldist"]="uniform"
    paramdict1["intervaldistparam"]=0.01 #change probability for uniform case
    paramdict1["startgraphalgo"]="rds"
    paramdict1["startgraphalgoparam"]=0.1
    paramdict2={}
    paramdict2["intervaldist"]="uniform"
    paramdict2["intervaldistparam"]=0.02 
    paramdict2["startgraphalgo"]="rds"
    paramdict2["startgraphalgoparam"]=0.1
    paramdict3={}
    paramdict3["intervaldist"]="uniform"
    paramdict3["intervaldistparam"]=0.05
    paramdict3["startgraphalgo"]="rds"
    paramdict3["startgraphalgoparam"]=0.1
    paramdict4={}
    paramdict4["intervaldist"]="uniform"
    paramdict4["intervaldistparam"]=0.1 #percentage of edges that will change(exist-> nonexist, nonexist,=->exist)
    paramdict4["startgraphalgo"]="rds"
    paramdict4["startgraphalgoparam"]=0.1
    algosparameters["fusedlasso"]=[paramdict1,paramdict2,paramdict3,paramdict4] 
    algosparameters["kernelsmooth"]=[-1]
    algosparameters["exposmooth"]=[-1]
    algosparameters["doubleexposmooth"]=[-1]
    return algosparameters
    
    
def data_generate(startnodenum,stepsize,tcount,samplecount,underfoldername,algo):
    algosparameter=algoparamassigner(startnodenum)
    if algo in directedalgos:
       algotype="directed"
    else:
       algotype="undirected"
        
    params=algosparameter[algo]
    for param in params:
        paramnames=sorted(param.keys())
        strparamlist=[]
        for index in range(0,len(paramnames)):
            paramname=paramnames[index]
            elem="{0}{1}".format(paramname,param[paramname])
            strparamlist.append(elem)
        strsent="_".join(strparamlist)
        strparam=[str(param[paramname]) for paramname in paramnames]
        paramstr=" ".join(strparam)

        code="python dynamicgraphgeneratorpart.py {0} {1} {2} {3} {4} {5} {6} {7} {8}".format(strsent,underfoldername,samplecount,algo,algotype,startnodenum,stepsize,tcount,paramstr)        
        code="runCmd -c \""+code+"\" --nowait"
        os.system(code)
               
        

samplecount=20
graphfolder="synthetic_dynamic_graphs"
#algos=["fusedlasso","dmc","rds","lpa","undirected_ff","directed_ff","kronecker"] #Among the algorithms, smw can not thought as a temporal evolution model.
algos=["dmc","fusedlasso"]
smoothalgos=["fusedlasso","kernelsmooth","exposmooth","doubleexposmooth"]
randomgraphalgos=["dmc","rds","lpa","undirected_ff","directed_ff","kronecker"]
directedalgos=["directed_ff","kronecker"]

#Those variables will be assigned in data_generate depending on algo
allstepsizes="" #nunber of consective steps same graph will be preserved.
tcounts=""
startnodenums="" 
  
#Among randomgraphalgos, directed_ff and kronecker is run by leskovec code so the graphs obtained at consective steps are not progressive. KEEP THIS IN MIND.
#All algorithms except KRONECKER increase one by one whereas kronecker is multiplicative 128,256,512,1024,2048
#Right now we always want to make graph connected. IS this necessary? THINK ABOUT THIS

if __name__ == "__main__":
    for algo in algos:
        if algo in ["kronecker"]: #this is multiplicative so thats why it is different!!
           globals()["allstepsizes"]=[2,3,7] #nunber of consective steps same graph will be preserved.
           globals()["tcounts"]=[10,15]
           globals()["startnodenums"]=[16,32,64] #startnodenums=[126,250,500,1000,2000]
        else:    
           globals()["allstepsizes"]=[1,2,3,5,7,10] #nunber of consective steps same graph will be preserved.
           globals()["tcounts"]=[10,20,50,100]
           globals()["startnodenums"]=[32,64,128,256] #startnodenums=[126,250,500,1000,2000]
        templist=[startnodenums,allstepsizes,tcounts]
        params=list(itertools.product(*templist))
        for startnodenum,stepsize,tcount in params:
            underfoldername="{0}/datan{1}s{2}t{3}".format(graphfolder,startnodenum,stepsize,tcount)
            underfoldername="{0}/{1}".format(underfoldername,algo)
            if not os.path.exists(underfoldername):
               os.makedirs(underfoldername)
            data_generate(startnodenum,stepsize,tcount,samplecount,underfoldername,algo)
           

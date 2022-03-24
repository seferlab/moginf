# 2 STEP CEFER FOR DEFAULT CASE(WITHOUT KERNEL)
import networkx as nx
import numpy as np
import scipy as sp
import scipy.sparse
import random
import os
import math
import operator
import time
from moginfUtil import ceferParams
import moginfUtil
import sys
import moginfKernel
from functools import reduce
sys.path.append("./lib")
import myutilities as myutil
import scipy.optimize
import itertools
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from DiscreteDistribution import DiscreteDistribution as DisDist
from ContinuousDistribution import ContinuousDistribution as ContDist
import moginfcplexrunner

def run2Step(algofields,tracefields,infernodes,allnodes,runfolder):    
    """ CEFER(2 step for partial case) (can be either EM or Kernel)
    Args:
       algofields: fields related to algo
       tracefields: fields related to trace
       infernodes: nodes will be inferred
       allnodes: all nodes
       runfolder : run folder
    """
    traces,evol,smodel,dists,samplerate,iprobmethod = tracefields
    algoinfo,algopar = algofields
    errortype,sparsetype,cover,secondalgo = algoinfo
    assert Trace.IsTraceNoisy(traces)
    count,LOOP = 0, 10
    noise = False
    if Trace.IsTraceNoisy(traces):
       noise = True
    iprobs = {}
    if smodel in ["sir","seir"]:
       dist = DisDist.genPartDist(dists[Trace.I2R][0],dists[Trace.I2R][1],"normal")
       iparam = DisDist.genPdf(dist) 
    else:
       iparam = [] 
    for trin in range(len(traces)):
        sortedtimes = sorted(list(set([time for node in list(traces[trin].keys()) for time in list(traces[trin][node].keys())])))
        iprobs[trin] = {Trace.INFECTED: ceferUtil.getInfectProbs(smodel,iprobmethod,sortedtimes,traces[trin],iparam)}
    while count < LOOP:
        print("iteration {0}".format(count))
        print(iprobs)
        edgecoef,edge2coef,curmax = {} ,{}, 0
        substr = "Subject To\n"
        for trin in range(len(traces)):
            print("tr in {0}".format(trin))
            subalgofields = [algoinfo,algopar]
            print("iprobinfo")
            print(iprobs[trin]) 
            subtracefields = [traces[trin],evol,smodel,dists,samplerate,iprobs[trin]]
            cursubstr,curedgecoef,curedge2coef,curmax = cefercplexrunner.genPartialConst(subalgofields,subtracefields,infernodes,curmax)
            for key in list(curedgecoef.keys()):
                edgecoef.setdefault(key,0.0)
                edgecoef[key] += curedgecoef[key]
            for key in list(curedge2coef.keys()):
                edge2coef.setdefault(key,0.0)
                edge2coef[key] += curedge2coef[key]
            print(len(list(edgecoef.keys())), len(list(edge2coef.keys())))
            substr += cursubstr
        subalgofields = [algoinfo,algopar]
        objstr = cefercplexrunner.genObjConst(subalgofields,edgecoef,edge2coef,curmax,noise)
        boundstr = cefercplexrunner.genBoundConst(edgecoef,curmax)
        outmethod = getattr(InputOutput,"convertCeferOut")
        retvalues = cefercplexrunner.runCode(substr,objstr,boundstr,algopar,evol,runfolder,outmethod)
        print("values at step {0}".format(count))
        print(len(list(retvalues.keys())))
        
        probcoef,probcoef2,curmax = {},{},0
        substr = "Subject To\n"
        if secondalgo == "Default":
           print("inside default part")
           subalgofields = [algoinfo,algopar] 
           for trin in range(len(traces)):
               subtracefields = [traces[trin],trin+1,dists,evol,smodel,samplerate]
               cursubstr,curprobcoef,curprobcoef2 = genDefaultSecondConst(subalgofields,subtracefields,infernodes,retvalues)
               for key in list(curprobcoef.keys()):
                   probcoef.setdefault(key,0.0)
                   probcoef[key] += curprobcoef[key]
               for key in list(curprobcoef2.keys()):
                   probcoef2.setdefault(key,0.0)
                   probcoef2[key] += curprobcoef2[key]    
               substr += cursubstr
               print(len(list(edgecoef.keys())), len(list(edge2coef.keys())))
           objstr = genDefaultSecondObjConst(subalgofields,probcoef,probcoef2)
           boundstr = genDefaultSecondBoundConst(probcoef)
           readoutmethod = getattr(Probs,"returnProbOut")
           tempiprobs = cefercplexrunner.runCode(substr,objstr,boundstr,algopar,evol,runfolder,readoutmethod)
           #for trin in xrange(tempiprobs.keys()):
           #    for  in xrange(tempiprobs[trin].keys()):
               
        elif secondalgo == "Kernel":
           print("inside kernel")
           subalgofields = [algoinfo,algopar] 
           for trin in range(len(traces)):
               trace = traces[trin]
               if len(set([time for node in list(trace.keys()) for time in list(trace[node].keys())])) <= 1:
                  continue 
               subtracefields = [trace,trin+1,dists,evol,smodel,samplerate]
               cursubstr,curprobcoef,curprobcoef2,curmax = ceferKernel.genKernelConst(subalgofields,subtracefields,infernodes,retvalues,curmax)
               for key in list(curprobcoef.keys()):
                   probcoef.setdefault(key,0.0)
                   probcoef[key] += curprobcoef[key]
               for key in list(curprobcoef2.keys()):
                   probcoef2.setdefault(key,0.0)
                   probcoef2[key] += curprobcoef2[key]    
               substr += cursubstr
           objstr = ceferKernel.genKernelObjConst(subalgofields,probcoef,probcoef2,curmax)
           boundstr = ceferKernel.genKernelBoundConst(probcoef)
           readoutmethod = getattr(ceferUtil,"returnKernelProbOut")
           alphavals = cefercplexrunner.runCode(substr,objstr,boundstr,algopar,evol,runfolder,readoutmethod)
           iprobs = ceferKernel.estimateIprob(alphavals,traces)
        print(iprobs)
        print("done iprob")
        exit(1)
        print(runfolder)
        print("iprobs")
        print(iprobs)
        #iprobs = ceferUtil.prepareIprob(iprobs,tracenum)
        count += 1
    print(iprobs)    
    print(edge2val)
    exit(1)
    return edge2val


def genDefaultSecondObjConst(algofields,probcoef,prob2coef):
    """Returns objective function string
    Args:
        algofields:
        probcoef:
        prob2coef:
    Returns:
        objstr: objective function in LP format
    """
    algoinfo,algopar = algofields
    errortype,sparsetype,cover,secondalgo = algoinfo
    assert errortype in ["abse","lse"]
    isPlus = lambda x: "+" if x >= 0 else " "

    objstr = "Minimize\n obj: "
    if errortype in ["abse","lse"]:
       for trin,n,t in list(probcoef.keys()):
           if probcoef[(trin,n,t)] != 0.0:
              objstr += " {0} %.8f i{1}?{2}?{3} ".format(isPlus(probcoef[(trin,n,t)]),trin,n,t) %probcoef[(trin,n,t)] 
    if errortype == "lse":
       objstr += " + [ "
       for trin,n1,t1,n2,t2 in list(prob2coef.keys()):
           if prob2coef[(trin,n1,t1,n2,t2)] != 0.0:
              objstr += " {0} %.8f i{1}?{2}?{3} * i{1}?{4}{5} ".format(isPlus(prob2coef[(trin,n1,t1,n2,t2)]),trin,n1,t1,n2,t2) %(2.0*prob2coef[(trin,n1,t1,n2,t2)])
       objstr += " ] / 2"         
      
    return objstr   
  
def genDefaultSecondBoundConst(probcoef):
    """Returns boundary constraints
    Args:
       probcoef: estimated probcoefs
    Returns:
       boundstr: returns the boundary string
    """
    boundstr = "Bounds\n " + "\n".join(["0 <= i{0}?{1}?{2} <= 1 ".format(index,node,time) for (index,node,time) in list(probcoef.keys()) if round(probcoef[(index,node,time)],3) != 0]) + "\n"
    eqstr = ""
    indices = set([index for index,node,time in list(probcoef.keys())])
    for index in indices:
        allnodes = set([node for tindex,node,time in list(probcoef.keys()) if tindex == index])
        for node in allnodes:
            parts = [" i{0}?{1}?{2} ".format(index,node,time) for tindex,tnode,time in list(probcoef.keys()) if tnode == node and tindex == index]
            if len(parts) != 0: 
               eqstr += " + ".join(parts) + " <= 1.0 \n"
    return boundstr + eqstr + "\n" 


def genDefaultSecondConstPart(node,tracefields,timefields,errortype,retvalues):
    """ generates constraint for each node-time pair
    Args:
       node: node
       tracefields: traceinfo
       timefields: pretime,curtime
       errortype:
       retvalues: edge return values
    Returns:
       substr: constraint string
       probcoef:
       prob2coef:
    """
    pretime,curtime = timefields
    trace,trin,evol,smodel,samplerate,modelparams,sortedtimes = tracefields
    s2icoef = modelparams[0]
    probcoef,prob2coef = {},{}
    if round(trace[node][pretime][Trace.SUSCEPTIBLE],3) == 0:
       return ["",{},{}]
    right = cefercplexrunner.getRight(trace,node,curtime,pretime,smodel)
    for sender in list(trace.keys()):
        if (sender,node) not in retvalues or round(retvalues[(sender,node)],3) == 0.0 or sender == node or round(trace[sender][pretime][Trace.INFECTED],3) == 0:
           continue
        for time in sortedtimes:
            basesum = reduce(operator.mul, [s2icoef[temptime - time] for temptime in range(pretime,curtime) if temptime - time in s2icoef],1.0)
            if round(basesum,3) == 1:
               continue
            if errortype == "abse":
               if math.pow(LOGBASE,right) >= 0.5:  
                  probcoef[(trin,sender,time)] = -1.0 * math.log(basesum,ceferParams.LOGBASE) * retvalues[(sender,node)]
               else:
                  probcoef[(trin,sender,time)] = math.log(basesum,ceferParams.LOGBASE) * retvalues[(sender,node)]
            else:
               probcoef[(trin,sender,time)] = math.log(basesum,ceferParams.LOGBASE) * retvalues[(sender,node)]
    if len(list(probcoef.keys())) == 0:
       return ["",{},{}]
    if errortype == "lse":
       prob2coef = {(trin,n1,t1,n2,t2): probcoef[(trin1,n1,t1)] * probcoef[(trin1,n2,t2)] for (trin1,n1,t1) in list(probcoef.keys()) for (trin1,n2,t2) in list(probcoef.keys())}
       probcoef = {key: -2.0 * right * probcoef[key] for key in list(probcoef.keys())}
    return ("",probcoef,prob2coef)

def genDefaultSecondConst(algofields,tracefields,infernodes,retvalues):
    """ generates infection probability constraint program
    Args:
       algofields:
       tracefields: parameters related to trace
       infernodes: nodes to be run
       retvalues: x values
    Returns:
       substr:
       edgecoef:
       edge2coef:
    """
    algoinfo,algopar = algofields
    trace,trin,dists,evol,smodel,samplerate = tracefields
    substr,probcoef,prob2coef,modelparams = "",{}, {},[]
    errortype,sparsetype,cover,extra = algoinfo
    if smodel == "seir":
       nons2ecoef = DisDist.genPartDist(dists[Trace.S2E][0],dists[Trace.S2E][1],"reverseCdf",dists[Trace.SPROB])
       modelparams.append(nons2ecoef)
    elif smodel in ["si","sir","sis"]:
       nons2icoef = DisDist.genPartDist(dists[Trace.S2I][0],dists[Trace.S2I][1],"reverseCdf",dists[Trace.SPROB])
       modelparams.append(nons2icoef)
    if smodel in ["sir","seir"]:
       i2rcoef = DisDist.genPartDist(dists[Trace.I2R][0],dists[Trace.I2R][1],"normal")
       modelparams.append(i2rcoef)
    if smodel == "seir":
       e2icoef = DisDist.genPartDist(dists[Trace.E2I][0],dists[Trace.E2I][1],"normal")
       modelparams.append(e2icoef)
    sortedtimes = sorted(list(set([time for node in list(trace.keys()) for time in list(trace[node].keys())])))
    for node in set(infernodes).intersection(set(trace.keys())):
        for tindex in range(1,len(sortedtimes)):
            pretime,curtime = sortedtimes[tindex-1:tindex+1]
            subtimefields = [pretime,curtime]
            subtracefields = [trace,trin,evol,smodel,samplerate,modelparams,sortedtimes]
            cursubstr,curprobcoef,curprob2coef = genDefaultSecondConstPart(node,subtracefields,subtimefields,errortype,retvalues)
            for key in list(curprobcoef.keys()):
                probcoef.setdefault(key,0.0)
                probcoef[key] += curprobcoef[key]
            for key in list(curprob2coef.keys()):
                prob2coef.setdefault(key,0.0)
                prob2coef[key] += curprob2coef[key]
    print("inside constraint")            
    print((len(list(probcoef.keys())),len(list(prob2coef.keys()))))            
    if "cover" in algoinfo:
       for node in set(infernodes).intersection(set(trace.keys())):
           for tindex in range(1,len(sortedtimes)):
               curtime = sortedtimes[tindex]
               parts = [" {0} i{1}?{2}?{3} ".format(retvalues[(prenode,node)],trin,prenode,sortedtimes[tindex2]) for prenode in list(trace.keys()) for tindex2 in range(len(sortedtimes)) if (prenode,node) in retvalues and retvalues[(prenode,node)] != 0.0 and tindex2 < tindex] 
               if len(parts) != 0:       
                  tempstr = " + ".join(parts) + " - i{0}?{1}?{2} >= 0\n".format(trin,node,curtime)
                  substr += tempstr
    return (substr,probcoef,prob2coef)


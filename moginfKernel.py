##Cefer Kernel
import networkx as nx
import numpy as np
import scipy as sp
import scipy.sparse
import random
import os
import math
import gzip
import operator
import pickle
import time
import sys
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
from moginfParams import ceferParams

def genKernelObjConst(algofields,probcoef,prob2coef,curmax):
    """ generates kernel obj function
    Args:
        algofields:
        probcoef:
        prob2coef:
        curmax:
    Returns:
        objstr: objective function in LP format
    """
    algoinfo,algopar = algofields
    errortype,sparsetype,cover,secondalgo = algoinfo
    assert errortype in ["abse","lse"]
    isPlus = lambda x: "+" if x >= 0 else " "

    objstr = "Minimize\n obj: "
    if errortype in ["abse","lse"]:
       for trin,alpha in probcoef.keys():
           if probcoef[(trin,alpha)] != 0.0:
              objstr += " {0} %.8f k{1}?{2} ".format(isPlus(probcoef[(trin,alpha)]),trin,alpha) %probcoef[(trin,alpha)] 
    if errortype == "lse":
       objstr += " + [ "
       for t1,a1,t2,a2 in prob2coef.keys():
           if prob2coef[(t1,a1,t2,a2)] != 0.0:
              objstr += " {0} %.8f k{1}?{2} * k{3}?{4} ".format(isPlus(prob2coef[(t1,a1,t2,a2)]),t1,a1,t2,a2) %(2.0*prob2coef[(t1,a1,t2,a2)])
       objstr += " ] / 2"         
    
    return objstr     

def estimateIprob(alphavals,traces):
    """ converts kernel alphas to iprobs
    Args:
       alphavals:
       traces:
    Returns:
       iprobs:
    """
    iprobs = {}
    for trin in xrange(len(traces)):
        iprobs[trin] = {Trace.INFECTED: {}}
        trace = traces[trin]
        alphas = sorted(list(set([alpha for ttrin,alpha in alphavals.keys() if ttrin == trin])))
        sortedtimes = sorted([time for time in trace[trace.keys()[0]] if trace[trace.keys()[0]][time].has_key(Trace.INFECTED)])
        newalphas = []
        for index in xrange(len(sortedtimes)):
            if index in alphas:
               newalphas.append(index)
            else:
               newalphas.append(0.0)
        mat = genKernelMat("gaussian",max(sortedtimes))
        totvec = {node: getKernelCoef(trace,node,sortedtimes,mat) for node in trace.keys()}
        posnodes = set(trace.keys())
        for node in totvec.keys():
            iprobs[trin][Trace.INFECTED][node] = {}
            alphavec = np.array(newalphas)
            print(np.shape(totvec[node]), np.shape(alphavec))
            resvec = np.dot(totvec[node],alphavec)
            for tindex in xrange(len(sortedtimes)):
                time = sortedtimes[tindex]
                iprobs[trin][Trace.INFECTED][node][time] = resvec[tindex]
    return iprobs

def genKernelBoundConst(probcoef):
    """Returns boundary constraints for Kernel case
    Args:
       probcoef: estimated probcoefs
    Returns:
       boundstr: returns the boundary string
    """
    boundstr = "Bounds\n " + "\n".join(["k{0}?{1} >= 0".format(index,alpha) for (index,alpha) in probcoef.keys()]) + "\n"
    return boundstr 

def genKernelConstPart(node,tracefields,timefields,errortype,retvalues,totvec):
    """ generates kernel constraint part
    Args:
       node: node
       tracefields: traceinfo
       timefields: pretime,curtime
       errortype:
       retvalues: edge return values:
       totvec: kernel matrix for node
    Returns:
       substr: constraint string
       probcoef:
       prob2coef:
    """
    pretime,curtime = timefields
    trace,trin,evol,smodel,samplerate,modelparams,sortedtimes = tracefields
    s2icoef = modelparams[0]
    alphas = range(0,np.shape(totvec)[1])
    probcoef = {(trin,alpha): 0.0 for alpha in alphas}
    if round(trace[node][pretime][Trace.SUSCEPTIBLE],3) == 0:
       return ["",{},{}]
    right = cefercplexrunner.getRight(trace,node,curtime,pretime,smodel)
    for sender in trace.keys():
        if not retvalues.has_key((sender,node)) or round(retvalues[(sender,node)],3) == 0.0 or sender == node or round(trace[sender][pretime][Trace.INFECTED],3) == 0:
           continue
        for tindex in xrange(len(sortedtimes)):
            time = sortedtimes[tindex]
            basesum = reduce(operator.mul, [s2icoef[temptime - time] for temptime in xrange(pretime,curtime) if s2icoef.has_key(temptime - time)],1.0)
            if round(basesum,3) == 1:
               continue
            val = math.log(basesum, ceferParams.LOGBASE) * retvalues[(sender,node)]
            if errortype == "abse":
               if math.pow(ceferParams.LOGBASE,right) >= 0.5:
                  val *= -1.0
            for alpha in alphas:      
                probcoef[(trin,alpha)] += (val * totvec[tindex,alpha])
    if len(probcoef.keys()) == 0:
       return ["",{},{}]
    return ("",probcoef,{})


def genKernelMat(type,maxtime):
    """generates kernel matrix
    Args:
       type: type of kernel
       maxtime: kernel leng will be 2*maxtime
    Returns:
       kernmat: kernel matrix
    """
    assert type in ["gaussian"]
    if type == "gaussian":
       rho2 = 1.0 
       kernmat = np.zeros((2*maxtime+1,2*maxtime+1),dtype = np.float64)
       for time1 in xrange(np.shape(kernmat)[0]):
           for time2 in xrange(np.shape(kernmat)[1]):
               kernmat[time1,time2] = math.exp((-0.5 * float(time1 - time2)**2) /rho2)
    return kernmat

def getKernelCoef(trace,node,sortedtimes,mat):
    """ generates coefficients after kernel multiplications
    Args:
       trace: given trace
       node: node
       sortedtimes: sorted seen times
       mat: matrix
    Returns:
       totvec: matrix of size
    """
    timesarr = np.array([trace[node][time][Trace.INFECTED] for time in sortedtimes])
    totvec = np.zeros((len(sortedtimes),len(sortedtimes)),dtype=np.float64) #times x kernelnum
    for index in xrange(len(sortedtimes)):
        meantime = sortedtimes[index]
        addcoef = (-1.0 * meantime) + (np.shape(mat)[0] / 2)
        for index2 in xrange(len(sortedtimes)):
            totvec[index2,index] = sum([mat[sortedtimes[index2] + addcoef,sortedtimes[index3] + addcoef] * timesarr[index3] for index3 in xrange(len(sortedtimes))])
    return totvec


def genKernelConst(algofields,tracefields,infernodes,retvalues,curmax): 
    """ generates kernel constraints
    Args:
       algofields:
       tracefields: parameters related to trace
       infernodes: nodes to be run
       retvalues: x values
       curmax: current maximum
    Returns:
       substr:
       edgecoef:
       edge2coef:
       curmax: current max for prob dist constraint
    """
    isPlus = lambda x: "+" if x > 0 else " "
    algoinfo,algopar = algofields
    trace,trin,dists,evol,smodel,samplerate = tracefields
    substr,probcoef,prob2coef,modelparams = "",{}, {},[]
    errortype,sparsetype,cover,secondalgo = algoinfo
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
    sortedtimes = sorted(list(set([time for node in trace.keys() for time in trace[node].keys()])))
    mat = genKernelMat("gaussian",max(sortedtimes))
    posnodes = set(infernodes).intersection(set(trace.keys()))
    totvec = {node: getKernelCoef(trace,node,sortedtimes,mat) for node in trace.keys()}
    for node in posnodes:
        for tindex in xrange(1,len(sortedtimes)):
            pretime,curtime = sortedtimes[tindex-1:tindex+1]
            subtimefields = [pretime,curtime]
            subtracefields = [trace,trin,evol,smodel,samplerate,modelparams,sortedtimes]
            cursubstr,curprobcoef,curprob2coef = genKernelConstPart(node,subtracefields,subtimefields,errortype,retvalues,totvec[node])
            for key in curprobcoef.keys():
                probcoef.setdefault(key,0.0)
                probcoef[key] += curprobcoef[key]
            for key in curprob2coef.keys():
                prob2coef.setdefault(key,0.0)
                prob2coef[key] += curprob2coef[key]
        alphas = range(0,np.shape(totvec[node])[1])
        parts = [" {0} k{1}?{2} ".format(totvec[node][tindex,a1],int(trin),a1) for a1 in alphas for tindex in xrange(np.shape(totvec[node])[0])]
        if len(parts) != 0:
           eqstr = "probeq{0} : ".format(random.randint(0,100000000)) + " + ".join(parts) +  " <= 1.0 \n"
           substr += eqstr
           #addvar = "add_{0}".format(curmax)
           #curmax += 1
           #addpart = " + 1.0 {0} ".format(addvar)
           #eqstr = "probeq{0} : ".format(random.randint(0,100000000)) + " + ".join(parts) +  " {0} = 1.0 \n".format(addpart)
           #substr += eqstr
        print(len(probcoef.keys()),len(prob2coef.keys()))            
        if "cover" in algoinfo:
            for tindex in xrange(1,len(sortedtimes)):
                curtime = sortedtimes[tindex]
                parts = []
                for alpha in alphas:
                    lcoef = sum([retvalues[(prenode,node)] * totvec[prenode][tindex2,alpha] for prenode in trace.keys() for tindex2 in xrange(len(sortedtimes)) if retvalues.has_key((prenode,node)) and retvalues[(prenode,node)] != 0.0 and totvec[prenode][tindex2,alpha] != 0.0 and tindex2 < tindex])
                    rcoef = totvec[node][tindex,alpha]            
                    parts.append(" {0} {1} k{2}?{3} ".format(isPlus(lcoef-rcoef),lcoef-rcoef,trin,alpha))
                if len(parts) != 0:
                   tempstr = "covereq{0}: ".format(random.randint(0,100000000)) + " ".join(parts) + " >= 0 \n"
                   substr += tempstr
    return (substr,probcoef,prob2coef,curmax)    



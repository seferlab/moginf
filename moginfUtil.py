#Methods related to infection/exposition probabilities
import math
import numpy as np
import scipy as sp
import sys
import os
import random
from InputOutput import InputOutput
from Trace import Trace

class ceferParams():
    """ Cefer parameters
    """
    EPSILON = 0.000001 #epsilon for rhs
    LOGBASE = 2.0 #log base
    INFINITYTIME = 10000000 #infinity time


def randomRoundEdges(edgedict):
    """ round edges randomized way
    Args:
       edgedict: current edge values
    Returns:
       newedges: rounded edge values
    """
    return { key : 1.0 for key in edgedict.keys() if edgedict[key] >= random.random()}

def isRoundingValid(edges,probs):
    """checks whether rounding is valid or not
    Args:
       edges: rounded edges dictionary
       probs: rounded probs dictionary
    Returns:
       bool: boolean value for validness   
    """
    node2time = {node : probs[node].keys()[0] for node in probs.keys() if probs[node].has_key(time)}
    alltimes = sorted(list(set(node2time.values())))
    for node1 in node2time.keys():
        time1 = node2time[node1]
        if time1 == alltimes[0]:
           continue
        flag = False
        for node2 in node2time.keys():
            time2 = node2time[node2]
            if time2 <= time1 and edges.has_key((node1,node2)):
               flag = True
               break
        if not flag:
           return False
    return True

def randomSimulRound(edgedict,allprobs):
    """ random rounds both edges and probabilities until a valid one is found
    Args:
       edgedict:
       allprobs:
    Returns:
       validedges,validprobs: valid rounded probs and edges
    """      
    while True:
       newedges = randomRoundEdges(edgedict)
       newallprobs,flag,index = {}, True, 0
       for probs in allprobs:
           newprobs = randomRoundItimes(probs)
           newallprobs[index] = newprobs
           index += 1 
           if not isRoundingValid(newedges,newprobs):
              flag = False
              break
       if flag:
          break 
    return (newedges,newprobs)                   
       
def randomRoundItimes(probs):
    """ round infection times randomizedly
    Args:
       probs: current probabilistic times
    Returns:
       newprobs: rounded probabilities
    """
    newprobs = {}
    for node in probs.keys():
        newprobs[node] = {}
        for time in probs[node].keys():
            if probs[node][time] >= random.random():
               newprobs[node][time] = 1.0
               break
    return newprobs

def returnKernelProbOut(cplexoutpath,evol):
    """ returns kernel output parameters
    """ 
    retvalues = InputOutput.readCplexOut(cplexoutpath,specific = ["k"])
    alphas = {}
    for varname in retvalues.keys():
        trace,alpha = varname.replace("k","").split("?")
        alphas.setdefault(int(trace),{})
        alphas[int(trace)][int(alpha)] = retvalues[varname] #time might be float
    return alphas

def returnProbOut(cplexoutpath,evol):
    """returns infection probabilities output
        Args:
           cplexoutpath: CPLEX output file
           evol: static/dynamic graph
        Returns:
           retvalues2: edge variables returned
    """
    retvalues = InputOutput.readCplexOut(cplexoutpath,specific = ["i"])
    iprobs = {}
    for varname in retvalues.keys():
        trace,node,time = varname.replace("i","").split("?")
        iprobs.setdefault(int(trace),{})
        iprobs[int(trace)].setdefault(int(node),{})
        iprobs[int(trace)][int(node)][int(time)] = retvalues[varname] #time might be float
    return iprobs

def prepareIprob(prob,tracenum):
    """ prepares Iprob for inference methods
    Args:
       iprob:
       tracenum:
    Returns:
       prob: modified prob
    """
    for index in xrange(tracenum):
        if not prob.has_key(index):
           prob[index] = {ceferParams.INFINITYTIME: 1.0}
        else:   
           val = sum(prob[index].values())   
           if val < 0.9999:
              prob[index][ceferParams.INFINITYTIME] = 1.0 - val
    return prob

def getRecoverLseProb(trace,sortedtimes,probdistcoef,smodel):
    """return infected probabilities
    Args:
       trace: tracedata
       sortedtimes: all times sorted
       probdistcoef: recovering probability distribution as parts
       smodel: spreading model
    Returns:
       returns estimated probabilities
    """
    assert smodel in ["sir","seir"]
    retprobs = {}
    for node in trace.keys():
        A = np.zeros((len(sortedtimes) - 1,len(sortedtimes) - 1),dtype=np.float)
        b = np.zeros((len(sortedtimes) - 1),dtype=np.float)
        for timeindex in xrange(1,len(sortedtimes)):
            pretime,curtime = sortedtimes[timeindex-1:timeindex+1]
            if smodel == "sir":
               rightval = trace[node][curtime][Trace.RECOVERED] - trace[node][pretime][Trace.RECOVERED]
            elif smodel == "seir":
               rightval = trace[node][curtime][Trace.RECOVERED] - trace[node][pretime][Trace.RECOVERED]
               rightval += trace[node][curtime][Trace.INFECTED] - trace[node][pretime][Trace.INFECTED]
            b[timeindex-1] = rightval
            for timeindex2 in xrange(1,timeindex+1):
                if probdistcoef.has_key(timeindex - timeindex2 + 1):
                   A[timeindex - 1,timeindex2 - 1] = probdistcoef[timeindex - timeindex2 + 1]
        optx = scipy.optimize.nnls(A,b)[0] #nonnegative one, better!
        assert len(optx) == len(sortedtimes) - 1
        retprobs[node] = {sortedtimes[timeindex]: optx[timeindex] for timeindex in xrange(len(sortedtimes) - 1)}
        retprobs[node][ceferParams.INFINITYTIME] = max(0.0,1.0-sum(optx))    
        weightedsum = sum(retprobs[node].values())
        retprobs[node] = {time: retprobs[node][time]/float(weightedsum) for time in retprobs[node].keys()}
    return retprobs

def getLseProb(trace,sortedtimes,smodel):
    """returns least square error probabilities
    Args:
       trace: tracedata 
       sortedtimes: alltimes sorted
       smodel: spreading model
    Returns:
       returns the LSE infection probabilities for si model
    """
    assert smodel == "si"
    iprobs = {}
    for node in trace.keys():
        iprobs[node] = {}
        posinfected = [time for time in sortedtimes if trace[node][time][Trace.INFECTED] > 0.01] #possible infection times for node
        if len(posinfected) == 0:
           iprobs[node][ceferParams.INFINITYTIME] = 1.0
           continue
        iprobs[node] = {time: 0.0 for time in sortedtimes if time not in posinfected}
        for itime in posinfected:
            err = sum([trace[node][time][Trace.INFECTED]**2 for time in sortedtimes if time < itime])
            err += sum([(1.0-trace[node][time][Trace.INFECTED])**2 for time in sortedtimes if time >= itime])
            if err == 0.0:
               iprobs[node] = {} 
               iprobs[node][itime] = 1.0
               break 
            iprobs[node][itime] = 1.0 / err
        if err == 0.0:
           continue 
        err = sum([trace[node][time][Trace.INFECTED]**2 for time in sortedtimes])   
        iprobs[node][ceferParams.INFINITYTIME] = 1.0 / err
        totprob = sum(iprobs[node].values())
        iprobs[node] = {key: iprobs[node][key] / float(totprob) for key in iprobs[node].keys()}
    return iprobs

def getInfectProbs(smodel,iprobmethod,sortedtimes,trace,params):
    """ returns estimated i probabilities for noisy data 
    Args:
        smodel: spreading model
        iprobmethod: probability estimation method
        sortedtimes: trace times sorted
        trace: noisy trace info
        params: probability estimation method parameters
    Returns:
        estimated infection/exposed probabilities for uncertain data
    """
    assert iprobmethod in ["lse"]
    if iprobmethod == "lse" and smodel == "si":
       return getLseProb(trace,sortedtimes,smodel)
    elif iprobmethod == "lse" and smodel in ["sir","seir"]:
       return getRecoverLseProb(trace,sortedtimes,params,smodel)
    return   

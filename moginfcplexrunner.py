#  
# Cefer generates graph inference code in LP Format
# This code can easily be modified to run on solvers other than cplex
#
import networkx as nx
import numpy as np
import scipy as sp
import scipy.sparse
import random
import os
import math
import gzip
import pickle
import time
import sys
sys.path.append("./lib")
import myutilities as myutil
import scipy.optimize
import itertools
import moginfEM
import operator
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from DiscreteDistribution import DiscreteDistribution as DisDist
from ContinuousDistribution import ContinuousDistribution as ContDist
import autotracegen
from moginfUtil import ceferParams

class ceferParams():
    """ Cefer parameters
    """
    EPSILON = 0.000001 #epsilon for rhs
    LOGBASE = 2.0 #log base
    INFINITYTIME = 10000000 #infinity time

def data2NonInfectConst(tracefields,node,algofields):
    """writes constaints for noninfected nodes from trace data
    Args:
       tracefields: trace related fields
       node: node 
       algofields: algo related fields
    Returns:
       edgecoef:
       edge2coef:
    """
    trace,dists,statekey = tracefields
    errortype = algofields[0]
    maxitime = max(set([trace[tnode][Trace.INFECTED] for tnode in trace.keys() if Trace.INFECTED in trace[tnode]]))
    diststr = "s2{0}".format(statekey)
    edgecoef,edge2coef = {},{}
    for sender in trace.keys():
        if sender == node or Trace.INFECTED not in trace[sender]:
           continue
        if Trace.RECOVERED in trace[sender]:
           timedif = trace[sender][Trace.RECOVERED] - trace[sender][Trace.INFECTED]
        else:    
           timedif = maxitime - trace[sender][Trace.INFECTED]
        val = 1.0 - ContDist.getContCdf(dists[diststr][0],dists[diststr][1],timedif,dists[Trace.SPROB])
        edgecoef[(sender,node)] = -1.0 * math.log(val,ceferParams.LOGBASE)
    if errortype == "lse":
       edge2coef = {(node1,sender1,sender2): edgecoef[(sender1,node1)] * edgecoef[(sender2,node1)] for (sender1,node1) in edgecoef.keys() for (sender2,node1) in edgecoef.keys()}
       edgecoef = {}
    return edgecoef,edge2coef

def data2InfectConst(tracefields,node,algofields):
    """writes infection constraints from data
    Args:
       tracefields: trace related fields
       node: node
       algofields: algo related fields
    Returns:
       poseffectors: nodes that are infected before so they can affect "node"
       edgecoef: posssible edges from this data
       edge2coef:
    """
    [trace,smodel,dists,samplerate,statekey] = tracefields
    errortype = algofields[0]
    diststr = "s2{0}".format(statekey)
    if samplerate == 0:
       mydif,mydiv = 0.000000001, 0.000000001
    else:
       mydif,mydiv = samplerate, 1.0
    poseffectors,approx,edgecoef,edge2coef = [], 0.001, {}, {}
    for sender in trace.keys():
        if sender == node or Trace.INFECTED not in trace[sender]:
           continue
        if trace[sender][Trace.INFECTED] < trace[node][statekey]:
           timedif = trace[node][statekey] - trace[sender][Trace.INFECTED]  
           if Trace.RECOVERED in trace[sender] and trace[node][statekey] >= trace[sender][Trace.RECOVERED]:
              continue
           poseffectors.append(sender)
           val = 1.0 - ContDist.getContCdf(dists[diststr][0],dists[diststr][1],timedif,dists[Trace.SPROB])
           preval = 1.0 - ContDist.getContCdf(dists[diststr][0],dists[diststr][1],timedif-mydif,dists[Trace.SPROB])
           coef = (math.log(val,ceferParams.LOGBASE) - math.log(preval,ceferParams.LOGBASE)) / mydiv
           assert coef <= 0.00000001
           edgecoef[(sender,node)] = coef
    if errortype == "lse":
       edge2coef = {(node1,sender1,sender2): edgecoef[(sender1,node1)] * edgecoef[(sender2,node1)] for (sender1,node1) in edgecoef.keys() for (sender2,node1) in edgecoef.keys()}
       edgecoef = {key: -2.0 * math.log(ceferParams.EPSILON,ceferParams.LOGBASE) * edgecoef[key] for key in edgecoef.keys()}
    return (poseffectors,edgecoef,edge2coef)

def genPerfectConst(algofields,tracefields,infernodes):       
    """Generates consraints out of given trace data
    Args:
       algofields:
       tracefields:
       infernodes:
    Returns:
       coverstr: constraint string
       edgecoef : edge coefs
       edge2coef: quadratic coef
    """
    algoinfo,algopar = algofields
    trace,evol,smodel,dists,samplerate = tracefields
    if smodel in ["si","sir"]:       
       statekey = Trace.INFECTED
    elif smodel == "seir":
       statekey = Trace.EXPOSED
    errortype = algoinfo[0]    
    coverstr,edgecoef,edge2coef = "",{},{}
    for node in infernodes:
        if node in trace:
           assert(statekey in trace[node])
           subtracefields = [trace,smodel,dists,samplerate,statekey]
           subalgofields = [errortype]
           poseffectors,curedgecoef,curedge2coef = data2InfectConst(subtracefields,node,subalgofields)
           if len(poseffectors) > 0 and "cover" in algoinfo:
              subcoverstr = " + 1.0 ".join(["x{0}?{1}".format(item,node) for item in poseffectors])
              coverstr += "1.0 {0} >= 1.0\n".format(subcoverstr)    
        else:
           subtracefields = [trace,dists,statekey]
           subalgofields = [errortype] 
           curedgecoef,curedge2coef = data2NonInfectConst(subtracefields,node,subalgofields)
        for key in curedgecoef.keys():
            edgecoef[key] = curedgecoef[key]
        for key in curedge2coef.keys():
            edge2coef[key] = curedge2coef[key]
    return (coverstr,edgecoef,edge2coef)

def getRight(trace,node,curtime,pretime,smodel):
    """ estimates right side of constraint
    Args:
       trace: trace 
       node: node
       curtime: current time
       pretime: previous time
       smodel: spreading model
    Returns:
       right: right hand side
    """
    if smodel in ["sir","seir","si"]:
       firstright = float(trace[node][curtime][Trace.SUSCEPTIBLE]) / trace[node][pretime][Trace.SUSCEPTIBLE]
    if smodel in ["sir","seir"]:
       dif = trace[node][curtime][Trace.RECOVERED] - trace[node][pretime][Trace.RECOVERED]
       recovflow = trace[node][pretime][Trace.INFECTED] - dif
       secondright = trace[node][pretime][Trace.SUSCEPTIBLE] - trace[node][curtime][Trace.INFECTED] + recovflow
    elif smodel in ["si"]:
       secondright = 1.0 - (float(trace[node][curtime][Trace.INFECTED])/trace[node][pretime][Trace.SUSCEPTIBLE])
    right = math.log(max(ceferParams.EPSILON,(firstright + secondright)/2.0),ceferParams.LOGBASE)
    return right

def genSeirConst(algofields,tracefields,node,pretime,curtime):   
    """Generates partial constraints for given node for SEIR models
    Args:
       algofields: algo related fields
       tracefields: trace related fields
       node: node
       pretime:
       curtime:
    Returns:
       substr: constraint string
       edgecoef:
       edge2coef:
    """
    errortype = algofields[0]
    [trace,evol,smodel,modelparams,samplerate,iprob] = tracefields
    s2icoef = modelparams[0]
    edgecoef,edge2coef = {},{}
    if round(trace[node][pretime][Trace.SUSCEPTIBLE],3) == 0:
       return ["",{},{}]
    right = getRight(trace,node,curtime,pretime,smodel)
    for sender in trace.keys():
        if sender == node or round(trace[sender][pretime][Trace.INFECTED],3) == 0:
           continue
        totsum = 0.0
        for itime in sorted(iprob[Trace.INFECTED][sender].keys()):
            if itime >= curtime:
               break
            basesum = reduce(operator.mul, [s2icoef[temptime - itime] for temptime in range(pretime,curtime) if s2icoef.has_key(temptime - itime)],1.0)
            if round(basesum,3) == 1: #neighnode has no diffusion affect on node
               continue
            totsum += math.log(basesum,ceferParams.LOGBASE) * iprob[Trace.INFECTED][sender][itime]
        if errortype == "abse":
           if math.pow(ceferParams.LOGBASE,right) >= 0.5:  
              edgecoef[(sender,node)] = -1.0 * totsum
           else:
              edgecoef[(sender,node)] = totsum 
        else:
           edgecoef[(sender,node)] = totsum 
    if len(edgecoef.keys()) == 0:
       return ["",{},{}]
    if errortype == "lse":
       edge2coef = {(tnode,sender1,sender2): edgecoef[(sender1,tnode)] * edgecoef[(sender2,tnode)] for (sender1,tnode) in edgecoef.keys() for (sender2,tnode) in edgecoef.keys()}
       edgecoef = {key: -2.0 * right * edgecoef[key] for key in edgecoef.keys()}
    return ("",edgecoef,edge2coef)


def genPartialConst(algofields,tracefields,infernodes,curmax):
    """Generates consraints out of given partial/undersampled trace data
    Args:
       algofields:
       tracefields:
       infernodes:
       curmax:
    Returns:
       substr: constraint string
       edgecoef:
       edge2coef:
    """
    [algoinfo,algopar] = algofields
    [trace,evol,smodel,dists,samplerate,iprob] = tracefields      
    substr,edgecoef,edge2coef,modelparams= "",{}, {},[],
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
    errortype = algoinfo[0]   
    sortedtimes = sorted(list(set([time for node in trace.keys() for time in trace[node].keys()])))
    for node in set(infernodes).intersection(set(trace.keys())):
        for tindex in range(1,len(sortedtimes)):
            pretime,curtime = sortedtimes[tindex-1:tindex+1]
            subalgofields = [errortype]
            subtracefields = [trace,evol,smodel,modelparams,samplerate,iprob]
            cursubstr,curedgecoef,curedge2coef = genSeirConst(subalgofields,subtracefields,node,pretime,curtime)
            for key in curedgecoef.keys():
                edgecoef.setdefault(key,0.0)
                edgecoef[key] += curedgecoef[key]
            for key in curedge2coef.keys():
                edge2coef.setdefault(key,0.0)
                edge2coef[key] += curedge2coef[key]
    if "cover" in algoinfo:
       print("in cover") 
       for node in set(infernodes).intersection(set(trace.keys())):
           for curtime in sortedtimes[1:]:
               if iprob[Trace.INFECTED][node].has_key(curtime) and round(iprob[Trace.INFECTED][node][curtime],3) != 0.0:
                  parts = ["{0} x{1}?{2}".format(iprob[Trace.INFECTED][prenode][pretime],prenode,node) for prenode in trace.keys() for pretime in sortedtimes[0:tindex] if iprob[Trace.INFECTED][prenode].has_key(pretime) and round(iprob[Trace.INFECTED][prenode][pretime],3) != 0.0 and node != prenode]
                  if len(parts) != 0:
                     leftsum = sum([float(part.split(" ")[0]) for part in parts])
                     if leftsum <= iprob[Trace.INFECTED][node][curtime]:
                        addvar = "add_{0}".format(curmax)
                        curmax += 1
                        addpart = " + 1.0 {0} = ".format(addvar)
                     else:
                        addpart = " >= " 
                     linestr = " + ".join(parts) + " {0} {1} \n".format(addpart,iprob[Trace.INFECTED][node][curtime])
                     substr += linestr 
    return (substr,edgecoef,edge2coef,curmax)
        
def genTraceConst(algofields,tracefields,infernodes,allnodes):
    """Generates Trace constraints from trace data
    Args:
       algofields: algo related fields
       tracefields: trace related fields
       infernodes: nodes to be inferred
       allnodes: allnodes
    Returns:
       substr: generates trace constraints string  
       consnum: number of constraints
       posedges: possible edges from trace data
    """
    algoinfo,algopar = algofields
    traces,evol,smodel,dists,samplerate,iprobmethod = tracefields
    edgecoef,edge2coef,curmax = {} ,{}, 0
    substr = "Subject To\n"
    if not Trace.IsTraceNoisy(traces):    
       for trace in traces:
           subalgofields = [algoinfo,algopar] 
           subtracefields = [trace,evol,smodel,dists,samplerate]
           coverstr,curedgecoef,curedge2coef = genPerfectConst(subalgofields,subtracefields,infernodes)
           for key in curedgecoef.keys():
               edgecoef.setdefault(key,0.0)
               edgecoef[key] += curedgecoef[key]
           for key in curedge2coef.keys():
               edge2coef.setdefault(key,0.0)
               edge2coef[key] += curedge2coef[key]
           substr += coverstr
    else:
       for trace in traces:
           sortedtimes = sorted(list(set([time for node in trace.keys() for time in trace[node].keys()])))
           if smodel in ["sir","seir"]:
              dist = DisDist.genPartDist(dists[Trace.I2R][0],dists[Trace.I2R][1],"normal")
              iparam = DisDist.genPdf(dist) 
           else:
              iparam = []
           iprob = {Trace.INFECTED: ceferUtil.getInfectProbs(smodel,iprobmethod,sortedtimes,trace,iparam)}
           subalgofields = [algoinfo,algopar] 
           subtracefields = [trace,evol,smodel,dists,samplerate,iprob]
           cursubstr,curedgecoef,curedge2coef,curmax = genPartialConst(subalgofields,subtracefields,infernodes,curmax)
           for key in curedgecoef.keys():
               edgecoef.setdefault(key,0.0)
               edgecoef[key] += curedgecoef[key]
           for key in curedge2coef.keys():
               edge2coef.setdefault(key,0.0)
               edge2coef[key] += curedge2coef[key]
           print(len(edgecoef.keys()), len(edge2coef.keys()))
           print("curmax {0}".format(curmax))
           if "cover" in algoinfo:    
              substr += cursubstr    
    return (substr,edgecoef,edge2coef,curmax)

def genObjConst(algofields,edgecoef,edge2coef,curmax,noise):
    """Returns objective function string
    Args:
       algofields:
       edgecoef:
       edge2coef:
       curmax: max of covering constraints
       noise:
    Returns:
        objstr: objective function in LP format
    """
    algoinfo,algopar = algofields
    errortype,sparsetype,cover,secondalgo = algoinfo
    assert errortype in ["abse","lse"]
    isPlus = lambda x: "+" if x >= 0 else " "
    
    objstr = "Minimize\n obj: "
    if errortype in ["abse","lse"]:
       for node1,node2 in edgecoef.keys():
           if node1 != node2 and edgecoef[(node1,node2)] != 0.0:
              objstr += " {0} %.8f x{1}?{2} ".format(isPlus(edgecoef[(node1,node2)]),node1,node2) %edgecoef[(node1,node2)] 
    if errortype == "lse":
       objstr += " + [ "
       for node,sender1,sender2 in edge2coef.keys():
           if edge2coef[(node,sender1,sender2)] != 0.0 and sender1 <= sender2:
              objstr += " {3} %.8f x{1}?{0} * x{2}?{0} ".format(node,sender1,sender2,isPlus(edge2coef[(node,sender1,sender2)])) %(2.0*edge2coef[(node,sender1,sender2)]) 
       objstr += " ] / 2"         
    
    if cover == "cover":
       highval = 1000000.0
       objstr += " ".join([" {0} %.8f x{1}?{2} ".format(isPlus(highval),node1,node2) %highval for node1,node2 in edgecoef.keys() if node1 != node2])
       if noise:
          maxval = 100000000000.0
          objstr += " + ".join([" {0} add_{1} ".format(maxval,index) for index in range(curmax)])
       
    if sparsetype in ["l1","l1l2"]:
       objstr += " + " + " + ".join([" {0} x{1}?{2} ".format(algopar["lambda1"],node1,node2) for node1,node2 in edgecoef.keys() if node1 != node2 and edgecoef[(node1,node2)] != 0.0]) 
    if sparsetype in ["l2","l1l2"]:
       objstr += " + [ " 
       for node1,node2 in edgecoef.keys():
           if node1 != node2:
              objstr += " {0} {1} x{2}?{3} * x{2}?{3} ".format(isPlus(algopar["lambda2"]),algopar["lambda2"],node1,node2)
       objstr += " ] / 2 "        
    return objstr   
  
def genBoundConst(edgecoef,curmax):
    """Returns boundary constraints
    Args:
       edgecoef:
    Returns:
       boundstr: returns the boundary string
    """
    boundstr = "Bounds\n"
    allnodes = set([item for edge in edgecoef.keys() for item in edge ])
    for node1 in allnodes:
        for node2 in allnodes:
            if node1 != node2:  
               boundstr += "0 <= x{0}?{1} <= 1 \n".format(node1,node2)
    for index in range(curmax):
        boundstr += " add_{0} >= 0 \n ".format(index)  
    return boundstr        

def genTempConst(tempsparsetype,algopar,alltimes,posedges):
    """generates temporal constraints for dynamic graphs
    Args:
       tempsparsetype: temporal sparsity type(ex: fused)
       algopar: algorithm parameters
       alltimes: alltimes for dynamic graph case
       posedges: posssible edges 
    Returns:
       tempstr: temporal constraints for dynamic graphs
    """
    tempstr = ""
    if "fused" in algopar.keys():
        alltimes = sorted(list(alltimes))[0:-1]
        for time in alltimes:
            for node1,node2 in posedges:
                varname = "x{0}?{1}?{2}".format(node1,node2,time)
                curvarname = "tempabs{0}".format(varname) #temporal absolute variable for dynamic case
                tempstr += " + {0} {1} ".format(algopar["fused"],curvarname) 
    return tempstr
 
def runCode(consstr,objstr,boundstr,algopar,evol,runfolder,outmethod):
    """Makes CEFER out of given constraints and runs it
    Args:
        consstr: constraint string
        objstr: objective function string
        boundstr: boundary string
        algopar: algorithm parameters
        evol: graph is static/dynamic
        runfolder: run folder
        outmethod = function to be run before returning output
    Returns:
        edge2val : edge value mapping
    """
    PNUM = 1 #processor count
    algoparstr = "init" + "_".join([str(item) for item in set(algopar.keys()).union(algopar.values())])
    outlppath = "{0}/{1}.lp".format(runfolder,algoparstr)
    with open(outlppath,"w") as file:
       file.write(objstr+"\n")
       file.write(consstr)
       file.write(boundstr+"\n")
       file.write("End\n")
    cplexoutpath = "{0}/{1}.lpout".format(runfolder,algoparstr)
    cplexscriptpath = "{0}/{1}.script".format(runfolder,algoparstr)
    with open(cplexscriptpath,"w") as file:
       file.write("read {0}\n".format(outlppath))
       file.write("set threads {0}\n".format(PNUM))
       file.write("optimize\n") 
       file.write("display solution objective\n")
       file.write("display solution variables -\n")  
    t1 = time.time()
    code="cplex < {0} > {1}".format(cplexscriptpath,cplexoutpath)
    os.system(code)
    t2 = time.time()
    print("Graph inferred in {0} seconds".format(t2-t1))
    retvalues = outmethod(cplexoutpath,evol)
    #os.system("rm -rf {0}".format(runfolder))
    return retvalues

def runCefer(traces,vararr,infernodes,allnodes):    
    """runs CEFER and saves inferred graph to resultfilename
    Args:
       traces: traces
       vararr: arguments in array
       infernodes: runs CEFER only on infernodes(important when running CEFER in parallel)
       allnodes: all nodes
    """
    print(vararr)
    
    for var in vararr.keys():
        print(var,vararr[var])
        if type(vararr[var]) == type(""):
           exec('{0}="{1}"'.format(var,vararr[var]), globals())
        else:
           exec('{0}={1}'.format(var,vararr[var]), globals())
    #import sys
    #sys.exit(1)
    print("Running CEFER-{0} {1} {2} {3} in {4} mode".format(errortype,sparsetype,cover,fusedtype,runmode))
    if evol == "static":
       algoinfo = (errortype,sparsetype,cover,secondalgo)
    else:
       algoinfo = (errortype,sparsetype,cover,secondalgo,fusedtype)
    if not os.path.exists(runfolder):
       os.makedirs(runfolder)
    noise = False
    if Trace.IsTraceNoisy(traces):
       noise = True
    if noise and secondalgo in ["Default","Kernel"]:
       tracefields = [traces,evol,smodel,dists,samplerate,iprobmethod]
       algofields = [algoinfo,algopar]
       retvalues = ceferEM.run2Step(algofields,tracefields,infernodes,allnodes,runfolder)
       InputOutput.writeGraph2File(retvalues,resultfilename,evol)
    else:
       if MATRIXMODE:
          consmat,consbmat,covermat,node2loc = trace2Matrix(algoinfo,algopar,traces,infernodes,smodel,dists,iprobmethod)
          constr,objstr,boundstr = matrix2Lp(consmat,consbmat,covermat,node2loc,algoinfo,algopar)
       else:
          tracefields = [traces,evol,smodel,dists,samplerate,iprobmethod]
          algofields = [algoinfo,algopar]
          consstr,edgecoef,edge2coef,curmax = genTraceConst(algofields,tracefields,infernodes,allnodes)
          if evol == "dynamic":
             alltimes = Trace.getAllTimes(traces) 
             consstr += genTempConst(tempsparsetype,algopar,alltimes,posedges)
          subalgofields = [algoinfo,algopar]
          objstr = genObjConst(subalgofields,edgecoef,edge2coef,curmax,noise) 
          boundstr = genBoundConst(edgecoef,curmax)
          outmethod = getattr(InputOutput,"convertCeferOut")
          retvalues = runCode(consstr,objstr,boundstr,algopar,evol,runfolder,outmethod)
          InputOutput.writeGraph2File(retvalues,resultfilename,evol)
                
MATRIXMODE = False

def main():
    """runs CEFER by the parameters provided in config file and outputs the inferred graph
       is called only when code will be run in parallel
    Args:
       configfile: configuration filename
    """
    assert len(sys.argv) == 2
    vararr = readQsubParams(sys.argv[1]) 
    runCefer(vararr)
   
 

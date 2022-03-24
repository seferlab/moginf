import networkx as nx
import random
import os
import sys
sys.path.append("./lib")
import math
import myutilities as myutil
import cPickle
import gzip
import itertools
from Trace import Trace
from InputOutput import InputOutput

def getScore():
    """ gets score for other algos
    """
    return
       

def readNetInfOut(outfile):
    """ reads netinf output
    """
    edges = []
    with open(outfile,"r") as file:
       for line in file:
           line = line.rstrip()
           if line == "":
              continue
           edges.append((int(line.split(",")[0]),int(line.split(",")[1])))
    return edges   

def readMultitreeOut(outfile):
    """ reads multitree out
    """
    edges = []
    with open(outfile,"r") as file:
       for line in file:
           line = line.rstrip()
           if line == "" or line.startswith("src"):
              continue
           splitted = line.split("/")
           edges.append((int(splitted[0]),int(splitted[1])))
    return edges

def runConnie(traces,outrunfolder,traceinfo,lambda1,TOL = 0.1):
    """ runs Connie
    """
    (dist,param) = traceinfo
    assert dist in ["expo","powerlaw","uniform","weibull"]
    if dist == "expo":
       distnum = 2
       assert param[0] == 1.0
    elif dist == "powerlaw":
       distnum = 1
       assert param[0] == -2.0
    elif dist == "uniform":
       distnum = 3
       assert param[0] == 1
    elif dist == "weibull":
       distnum = 4,param[0]
       assert param == (9.749,2.3494)
    allnodes = set([node for trace in traces for node in trace.keys()])
    node2index = {allnodes[index]:index for index in xrange(len(allnodes))}
    difmat = np.zeros((len(traces),len(allnodes)),dtype=np.float).fill(-1)
    for trindex in xrange(len(traces)):
        for node in allnodes:
            if traces[trindex].has_key(node):
               difmat[trindex,node2index[node]] = traces[trindex][node][Trace.INFECTED]
    difmatstr = matrix2MatlabStr(difmat)
    olddir,newdir = os.getcwd(),codedir
    outprefix="-".join(["connie"] + outrunfolder.split("/") + [distnum])
    outfile = "{0}.mat".format(outprefix)
    outfile = "{0}/{1}/{2}".format(olddir,outrunfolder,outfilename)
    #code="opt/stow/matlab-r2012a/bin/matlab -r \" cd('./algos/connie'); connie(0.2,2,[1,2,3;2,3,4;1,2,2],0.2,'../myoutfilename.txt') ; cd('../../'); quit; \"" 
    code="{0} -r \" cd('{1}'); connie({2},{3},{4},{5},'{6}') ; cd('{7}'); quit; \" ".format(MATLABPATH,newdir,lambda1,distnum,difmatstr,TOL,outfile,olddir)
    os.system(code)
    difmat = scipy.io.loadmat(outfile)['A_mle']
    edges = set([(node1,node2) for node1,node2 in list(itertools.product(allnodes,allnodes)) if difmat[node2index[node1],node2index[node2]] != 0.0])
    os.sytem("rm -rf {0}".format(outfile))
    return edges

def runNetrate(traces,outrunfolder,traceinfo,nodenum):
    """ runs Netrate
    """
    (dist,param) = traceinfo
    assert dist in ["expo","powerlaw","rayleigh"]
    if dist == "expo":
       distnum,alpha = 'exp', param[0]
    elif dist == "pl":
       distnum,alpha = 'pl',-1.0-(1.0*param[0])
    elif dist == "rayleigh":
       distnum,alpha = 'rayleigh', 0.5/(param[0]**2)
    olddir,newdir = os.getcwd(),codedir
    outprefix = "-".join(["netrate"] + outrunfolder.split("/") + [distnum,alpha])
    newoutdir = "{0}/{1}".format(olddir,outrunfolder)
    if not os.path.exists(newoutdir):
       os.makedirs(newoutdir) 
    outfile = "{0}/{1}.mat".format(newoutdir,outprefix)
    tracefile = "{0}/{1}.netratetrace".format(codedir,outputprefix)
    node2node = Trace.saveOtherTrace(traces,tracefile,",",maptype="nodemap")
    code = "{0} -r \" cd('{1}'); cd('{2}'); cvx_setup; cd('{3}'); {4}('{5}','{6}',{7},'{8}',{9},'{10}'); cd('{11}'); quit; \" ".format(MATLABPATH,newdir,"cvx","..","netrate","[]",tracefile,WINDOW,distnum,nodenum,outfile,olddir)
    os.system(code)
    difmat = scipy.io.loadmat(outfilepath)['A_hat']
    edges = set([difmat[node1,node2] for node1 in xrange(np.shape(difmat)[0]) for node2 in xrange(np.shape(difmat)[1]) if node1 != node2])
    [os.system("rm -rf {0}".format(file)) for file in [outfile,tracefile]]
    return edges

def saveOtherTrace(traces,tracefile,seperator,maptype=None):
    """ save given traces for other files(leskovec traces)
    """
    allnodes = list(set([node for trace in traces for node in trace.keys()]))
    if maptype == "nodemap":
       sortednodes = sorted(allnodes)
       node2node = {sortednodes[index]: index for index in xrange(len(sortednodes))}
    elif maptype == None:
       node2node = {allnodes[index]: allnodes[index] for index in xrange(len(allnodes))}   
    with open(tracefile,"w") as file:
       file.write( "\n".join(["{0},{0}".format(node2node[node]) for node in allnodes]) + "\n\n")
       for trace in traces:
           time2nodes = {}
           for node in trace.keys():
               time2nodes.setdefault(trace[node][Trace.INFECTED],set()).add(node)
           file.write(seperator.join(["{0}{1}{2}".format(node2node[node],",",time) for time in sorted(time2nodes.keys()) for node in time2nodes[time]]) + "\n")
    return node2node

def runMultitree(traces,outrunfolder,traceinfo,itercount):
    """ runs Multitree
    """
    (dist,param) = traceinfo
    assert dist in ["expo","powerlaw","rayleigh","weibull"]
    if dist == "expo":
       distnum,alpha = 8, 1.0/param[0]
    elif dist == "powerlaw":
       distnum,alpha = 9,-1.0*param[0]
    elif dist == "rayleigh":
       distnum,alpha = 10,1.0/(param[0]**2)
    elif dist == "weibull":
       distnum,alpha = 11,param[0] 
    outprefix = "-".join(["multitree"] + outrunfolder.split("/") + [str(distnum),str(alpha)])
    outfile = "{0}-edge.info".format(outprefix)
    tracefile = "{0}/{1}.multitree".format(codedir,outprefix)
    saveOtherTrace(traces,tracefile,",")
    code = "{0}/network-inference-multitree -i:{1} -o:{2} -e:{3} -a:{4} -d:{5} -nc:{6} -m:{7} -s:{8}".format(codedir,tracefile,outprefix,itercount,alpha,1,-1,distnum,1)
    os.system(code)
    edges = readMultitreeOut(outfile)
    [os.system("rm -rf {0}".format(file)) for file in [outfile,tracefile]]
    return edges
    
def runNetinf(traces,outrunfolder,traceinfo,itercount):
    """ runs Netinf
    """
    (dist,param) = traceinfo
    assert dist in ["expo","powerlaw","rayleigh","weibull"]
    if dist == "expo":
       distnum,alpha = 8, 1.0/param[0]
    elif dist == "powerlaw":
       distnum,alpha = 9,-1.0*param[0]
    elif dist == "rayleigh":
       distnum,alpha = 10,1.0/(param[0]**2)
    elif dist == "weibull":
       distnum,alpha = 11,param[0] 
    outprefix = "-".join(["netinf"] + outrunfolder.split("/") + [str(distnum),str(alpha)])
    outfile = "{0}.txt".format(outprefix)
    tracefile = "{0}/{1}.netinf".format(codedir,outprefix)
    saveOtherTrace(traces,tracefile,";")
    code = "{0}/netinf -i:{1} -o:{2} -e:{3} -a:{4} -m:{5} -s:{6}".format(codedir,tracefile,outprefix,itercount,alpha,distnum,1)
    os.system(code)
    edges = readNetInfOut(outfile)
    [os.system("rm -rf {0}".format("{0}{1}".format(outprefix,extend))) for extend in ["-edge.info",".txt"]]
    os.system("rm -rf {0}".format(tracefile))
    return edges
       
def matrix2MatlabStr(mat):
    """converts matrix to matlab str
    """
    l1,l2 = np.shape(mat)
    matstr = ";".join([" ".join([str(mat[index1,index2]) for index2 in xrange(l2)]) for index1 in xrange(l1)])
    return "[{0}]".format(matstr[0:-1])

def trace2Trace(traces,inmodel,outmodel):
    """converts given trace model to another model
    Args:
       traces:
       outmodel:
    """
    if outmodel == "si":
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED]
    elif outmodel == "sir":
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED]
    elif outmodel == "seir":
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED,Trace.EXPOSED]
    elif outmodel == "sis":
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED]
    print states
    if inmodel == "seir" and outmodel == "si":
       mymap = {Trace.EXPOSED: Trace.INFECTED, Trace.SUSCEPTIBLE: Trace.SUSCEPTIBLE}
    else:
       mymap = {state: state for state in states}
    newtraces = [{node: {mymap[state]: trace[node][state] for state in trace[node].keys() if mymap.has_key(state) and mymap[state] in states} for node in trace.keys()} for trace in traces]
    for trace in newtraces:
        for node in trace.keys():
            if len(trace[node].keys()) == 0:
               del trace[node]
    return newtraces        

def runOther(algo,algoparams,tracefolder,resultfilename,tracecount,spreadparam,iprobmethod,runfolder,noise,graphfile,smodel):
    """ run other algos
    """
    G = InputOutput.readGraphAndParams(graphfile)
    (traces,allnodes) = InputOutput.readTraces(tracefolder,tracecount,smodel)
    if noise != 0.0:
       traces = Trace.roundTrace(traces,"random",{})
    assert not Trace.IsTraceNoisy(traces)
    traces = trace2Trace(traces,smodel,"si")
    methodname = "run{0}".format(algo[0])
    method = globals()[methodname]
    if smodel in ["sir","si"]:
       traceinfo = spreadparam[Trace.S2I] 
    elif smodel == "seir":
       traceinfo = spreadparam[Trace.S2E]
    methodvars = [traces,runfolder,traceinfo]
    if algo[0] in ["Netinf","Multitree"]:
       methodvars.append(int(algoparams["degree"]*G.number_of_edges()))
    elif algo[0] == "Connie":
       methodvars.append(algoparams["lambda1"])   
    elif algo[0] in "Netrate":
       methodvars.append(algoparams["nodenum"])
    edges = list(method(*methodvars))
    with open(resultfilename,"w") as file:
       file.write("\n".join(["{0} {1} 1.0".format(edge[0],edge[1]) for edge in edges]))
    

MATLABPATH = "/usr/bin/matlab"
codedir = "otheralgos"
WINDOW = 10

def main():
    with gzip.open(sys.argv[1],"rb") as file:
        algo = cPickle.load(file)
        algoparams = cPickle.load(file)
        tracefolder = cPickle.load(file)
        resultfilename = cPickle.load(file)
        tracecount = cPickle.load(file)
        spreadparam = cPickle.load(file)
        iprobmethod = cPickle.load(file)
        runfolder = cPickle.load(file)
        noise = cPickle.load(file)
        path = cPickle.load(file)
        smodel = cPickle.load(file)
    runOther(algo,algoparams,tracefolder,resultfilename,tracecount,spreadparam,iprobmethod,runfolder,noise,path,smodel)

if __name__ == "__main__":
    main()

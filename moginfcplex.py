#  
# Cefer generates graph inference code in LP Format
# This code can easily be modified to run on solvers other than cplex
#
import os
import gzip
import pickle
import sys
import networkx as nx
sys.path.append("./lib")
import myutilities as myutil
import math
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import deepcopy
import moginfcplexrunner
from Trace import Trace
from InputOutput import InputOutput
from Score import Score

def assignDefaultCeferParameters():
    """assigns default parameter values
    Returns:
       returns default value array
    """
    vararr = {}
    vararr["tracefolder"] = "traces"
    vararr["resultfilename"] = "ceferresults.edg"
    vararr["runfolder"] = "temprun"
    vararr["model"] = "si"
    vararr["evol"] = "static"
    vararr["errortype"] = "abse"
    vararr["sparsetype"] = "None"
    vararr["cover"] = "cover"
    vararr["secondalgo"] = "Kernel"
    vararr["fusedtype"] = "fused"
    vararr["iprobmethod"] = "lse"
    vararr["algopar"] = {}
    vararr["dists"] = {}
    vararr["dist"] = "cont"
    vararr["runmode"] = "serial"
    vararr["parallelcount"] = 1
    vararr["graphfilename"] = "graph.gml" 
    vararr["printscore"] = "f1"
    vararr["parstartin"] = 0
    return vararr

def readCeferConfigFile(filename):
    """reads CEFER configuration file and returns parameters
    Args:
       filename: configuration filename
    Returns:
       returns parameters from configuration file
    """
    vararr = assignDefaultCeferParameters()
    with open(filename,"r") as file:
        for line in file:
            varname,value = line.rstrip().split(":")
            value = value.lstrip()
            if varname in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
               dist = value.split(" ")[0]
               param = [float(item) for item in value.split(" ")[1:]]
               vararr["dists"][varname] = (dist,tuple(param))
            elif varname in [Trace.SPROB]:
               vararr["dists"][varname] = float(value)   
            elif varname in ["lambda1","lambda2","fusedlambda"]:
               vararr["algopar"][varname] = float(value)
            elif varname == "noise":
               vararr[varname] = float(value)
            elif varname in ["parallelcount","tracecount","parstartin","samplerate"]:
               vararr[varname] = int(value) 
            else:
               vararr[varname] = value
    assert vararr["runmode"] in ["serial","parallel"] and vararr["smodel"] in ["si","sir","seir","sis"] 
    return vararr

def plotgraph(inferdict,evol,resultfile):
    """plots inferred graph
    Args:
       inferdict: inferred dict
       evol: static/dynamic
       resultfile: resultfile
    """
    if evol == "static":
       G = nx.DiGraph()
       for node1,node2 in inferdict[0]:
           G.add_edge(node1,node2)
       plt.clf()
       nx.draw(G,pos=nx.spring_layout(G))
       plotfile = resultfile.replace(".edg",".png")
       plt.savefig(plotfile)
         
def getScore(graphfilename,resultfile,scorename,evol):
    """estimates given score type
    Args:
       graphfilename: graphfilename
       resultfile: resultfilename
       scorename: score type
       evol: evolution
    Returns:
       score: returns score of type scorename
    """
    G = InputOutput.readGraphAndParams(graphfilename)
    truthdict = {0: set(G.edges())}
    roundmethods = [("epsilon",value) for value in [0.2,0.999999999999,0.00001,0.05,0.1,0.15,0.25,0.35,0.3,0.4]]
    allscores,inferdicts = [],[]
    for method,params in roundmethods:
        edge2val = InputOutput.readCeferOut(resultfile,evol)
        inferdict = Score.roundScore2Set(edge2val,method,params,evol)
        allscores.append(Score.getBinaryScores(truthdict,inferdict,[scorename])[scorename])
        inferdicts.append(inferdict)
    maxscore = max(allscores)
    plotgraph(inferdicts[allscores.index(maxscore)],evol,resultfile)
    return maxscore

def makeSingleResult(resultfile,parallelcount):
    """ makes single resultfile out of parallel results
    Args:
       resultfile: resultfile prefix
       parallelcount: number of parallel files
    """
    with open(resultfile,"w") as outfile:
       for index in xrange(parallelcount):
           resultpath = "{0}_par{1}".format(resultfile,index+1)
           with open(resultpath,"r") as infile:
              for line in infile:
                  outfile.write(line)
    [os.system("rm -rf {0}_par{1}".format(resultfile,index+1)) for index in xrange(parallelcount)]

def main():
   """runs CEFER by the parameters provided in config file and outputs the inferred graph
   Args:
      configfile: configuration filename
   """
   if len(sys.argv) != 2:
      print("Configuration file must be given\n Look at example.config for format details!!")
      exit(1)
   code = "which cplex"
   checkval = os.system(code)
   if checkval != 0:
      print("Cplex is not installed in this system!!")
      print("Terminating Program!!")
      exit(1)  
   vararr = readCeferConfigFile(sys.argv[1])
   (traces,allnodes) = InputOutput.readTraces(vararr["tracefolder"],vararr["tracecount"],vararr["smodel"])
   print("{0} Traces read under folder {1}".format(len(traces),vararr["tracefolder"]))
   if vararr["runmode"] == "serial":
      t1 = time.time()
      moginfcplexrunner.runCefer(traces,vararr,allnodes,allnodes)
      t2 = time.time()
      print("Total inference done in {0} seconds".format(t2 - t1))
   elif vararr["runmode"] == "parallel":
      inferblocks,allnodes = [], sorted(list(allnodes))
      for index in xrange(vararr["parallelcount"]):
          frac = int(math.floor(len(allnodes)/float(vararr["parallelcount"])))
          start, end = index * frac, (index + 1) * frac
          if index == vararr["parallelcount"] - 1:
             end = len(allnodes)
          inferblocks.append(allnodes[start:end])
      print("{0} parallel jobs".format(vararr["parallelcount"]))    
      for index in range(vararr["parstartin"],len(inferblocks)):
          newvararr = deepcopy(vararr)
          newvararr["resultfilename"] = vararr["resultfilename"] + "_par{0}".format(index+1) 
          moginfcplexrunner.runCefer(traces,newvararr,inferblocks[index],allnodes)
      makeSingleResult(vararr["resultfilename"],vararr["parallelcount"])    
   if vararr["printscore"] != "None":
      score = getScore(vararr["graphfilename"],vararr["resultfilename"],vararr["printscore"],vararr["evol"])
      print("Comparing Cefer graph with truth {0} :".format(vararr["graphfilename"]))
      print("{0} score: {1} given {2} traces".format(vararr["printscore"],score,vararr["tracecount"]))       
        
if  __name__ == '__main__':
    main()

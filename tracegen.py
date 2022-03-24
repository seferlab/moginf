#  
# This code generates traces for SEIR models
#
import networkx as nx
import random
import os
import sys
sys.path.append("./lib")
sys.path.append("../lib")
import math
import myutilities as myutil
import pickle
import gzip
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace


def assignDefaultTraceParameters():
    """assigns default parameter values for tracegen
    Returns:
       returns default value array
    """
    vararr={}
    vararr["traceoutfolder"] = "traces"
    vararr["graphfilename"] = "graph.gml"
    vararr["tracecount"] = 100
    vararr["samplerate"] = 1
    vararr["startnodecount"] = 1
    vararr["noise"] = 0.5
    vararr["noiseshape"] = "uniform"
    vararr["smodel"] = "si"
    vararr["ongraph"] = False
    vararr["evol"] = "static"
    vararr["dists"] = {}
    vararr["dist"] = "cont"
    vararr["modelparams"] = {}
    vararr["startnode"] = None
    return vararr

def readTraceConfigFile(filename):
    """reads trace configuration file
    Args:
       filename: configuration filename
    Returns:
       returns parameters from configuration file
    """
    vararr = assignDefaultTraceParameters()
    with open(filename,"r") as file:
        for line in file:
            varname,value = line.rstrip().split(":")
            value = value.lstrip()
            if varname in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
               dist = value.split(" ")[0]
               param = [float(item) for item in value.split(" ")[1:]]
               vararr["dists"][varname] = (dist,tuple(param))
            elif varname == Trace.SPROB:
               vararr["dists"][varname] = float(value)
            elif varname == "ongraph":
               vararr[varname] = value == "True"   
            elif varname == "noise":
               vararr[varname] = float(value)
            elif varname in ["samplerate","tracecount","startnodecount","startnode"]:
               vararr[varname] = int(value)
            elif varname == "sismaxtime":
               vararr["modelparams"][varname] = int(value)
            else:
               vararr[varname] = value
    assert vararr["smodel"] in ["si","sir","seir","sis"] 
    return vararr

def main():
   """generates traces by the parameters provided in config file
   Args:
      configfile: configuration file
   """
   if len(sys.argv)!=2:
      print("Trace Configuration file must be given\n Look at example.config for format details!!")
      exit(1)   
   vararr = readTraceConfigFile(sys.argv[1])
   for var in vararr.keys():
       if type(vararr[var])==type(""):
          exec('{0}="{1}"'.format(var,vararr[var]), globals())
       else:
          exec('{0}={1}'.format(var,vararr[var]), globals())      
   assert dist in ["dis","cont"]
   if graphfilename.endswith(".gml"):
      G = InputOutput.readGraphAndParams(graphfilename,"gml")
   else:
      G = InputOutput.readGraphAndParams(graphfilename)
   if not ongraph:
      InputOutput.addAttribute2Graph(G,dists)
   if dist == "dis":       
      traces = DisTrace.genTraces(G,startnodecount,evol,tracecount,smodel,modelparams,startnode)
   elif dist == "cont":
      traces = ContTrace.genTraces(G,startnodecount,evol,tracecount,smodel,modelparams,startnode)  
   if samplerate != 0 or noise != 0.0:
      noiseinfo = (noiseshape,noise)
      traces = Trace.modifyTraces(traces,samplerate,noiseinfo,smodel)
   if noise == 0.0:
      noisesamplestr = "{0}-{1}".format(noise,samplerate)
   else:
      noisesamplestr = "{0}-{1}-{2}".format(noiseshape,noise,samplerate)  
   if not ongraph:
      sfolder = Trace.getSpreadFolder(smodel,dists,dist)
      tracefolder = "{0}/{1}/{2}".format(traceoutfolder,sfolder,noisesamplestr)
   else: 
      tracefolder = "{0}/{1}".format(traceoutfolder,noisesamplestr)
   if not os.path.exists(tracefolder):
      os.makedirs(tracefolder)
   maxindex = max([int(filename.replace(".plain","")) for filename in  myutil.listfiles(tracefolder)] + [-1])
   for index in range(len(traces)):
       tracefile = "{0}.plain".format(maxindex + index + 1)
       tracefilepath = "{0}/{1}".format(tracefolder,tracefile)
       InputOutput.writePlainTrace(traces[index],tracefilepath)
       
if __name__ == "__main__":
    main()

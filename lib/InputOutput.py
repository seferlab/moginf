import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import sys
import os
import myutilities as myutil
import gzip
import string
import pickle
from Trace import Trace

class InputOutput():     
    """ Methods related to Input/Output(graph/trace)
    """
    def __init__(self):
       return
    
    @staticmethod 
    def writeConfigFile(vararr,filename): 
       """write given varar to configuration file
       Args:
          vararr: config file parameters
          filename: configuration filename
       """
       assert vararr["smodel"] in ["si","sir","seir","sis"] 
       with open(filename,"w") as outfile:
          outfile.write("\n".join(["{0}: {1}".format(varname,value) for varname,value in vararr.items()]) + "\n")

    @staticmethod
    def writeGraph2File(retvalues,resultfile,evol):
        """stores retvalues as edgelist to file
        Args:
           retvalues: edge values in dictionary
           resultfilename: graph resultfilename
           evol: static/dynamic graph
        """
        if evol == "static":
           with open(resultfile,"w") as file:
              file.write("".join(["{0} {1} {2}\n".format(node1,node2,retvalues[(node1,node2)]) for node1,node2 in retvalues.keys()]))
        elif evol == "dynamic":
           alltimes = {}
           for node1,node2,time in retvalues.keys():
               alltimes.setdefault(time,set())
               alltimes[time].add((node1,node2))   
           with open(resultfile,"w") as file:
               for time in alltimes.keys():
                   file.write("Time: {0}\n".format(time))
                   file.write("".join(["{0} {1} {2}\n".format(node1,node2,retvalues[(node1,node2)]) for node1,node2 in retvalues[time].keys()]))
               
    @classmethod
    def convertCeferOut(self,cplexoutpath,evol):
        """reads output solution and returns the edges with values
        Args:
           cplexoutpath: CPLEX output file
           evol: static/dynamic graph
        Returns:
           retvalues2: edge variables returned
        """
        retvalues = self.readCplexOut(cplexoutpath,specific = ["x"])
        retvalues2 = {}
        if evol == "dynamic":
           for varname in retvalues.keys():
               node1,node2,time = [int(item) for item in varname.replace("x","").split("?")]
               retvalues2[(node1,node2,time)] = retvalues[varname]
        elif evol == "static":
           for varname in retvalues.keys():
              node1,node2 = [int(item) for item in varname.replace("x","").split("?")]
              retvalues2[(node1,node2)] = retvalues[varname]
        return retvalues2

    @classmethod
    def readCplexOut(self,outfile,specific=[]):
        """reads CPLEX output file and returns only SPECIFIC variable values as dictionary
        Args:
           outfile: CPLEX output file
           specific: specific variable prefixes such as x
        Returns:
           retvalues: variable-value dictionary
        """
        retvalues = {}
        varflag = False
        with open(outfile,"r") as file:
           for line in file:
               line = line.rstrip()
               if not varflag and line.find("CPLEX> Variable Name")!=-1 and line.find("Solution Value")!=-1:
                  varflag=True
                  continue
               if varflag:
                  for varname in specific: 
                      if line.startswith(varname):
                         key,value = line.split()
                         retvalues[key] = float(value)
                         break
        return retvalues

    @classmethod
    def addAttribute2Graph(self,G,dists):
        """ appends dist info to G as attributes(if there is interval, assigns one)
        Args:
          G: graph
          dists: distributions
        """
        for key in dists.keys():
            if key in [Trace.I2R,Trace.E2I,Trace.I2S]:
               if type(dists[key]) == dict:
                  for v in G.nodes():
                      G.nodes[v][key] = (dists[key]["dist"],tuple([random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]))
               else:
                  for v in G.nodes():
                      G.nodes[v][key] = dists[key]
            elif key in [Trace.S2I,Trace.S2E]:
               if type(dists[key]) == dict:
                  for u,v in G.edges():
                      G[u][v][key] = (dists[key]["dist"],tuple([random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]))
               else:
                  for u,v in G.edges():
                      G[u][v][key] = dists[key]
            elif key == Trace.SPROB:
               if type(dists[key]) == dict:
                  for u,v in G.edges():
                      G[u][v][key] = random.uniform(dists[key]["start"],dists[key]["end"]) 
               else:
                  for u,v in G.edges():
                      G[u][v][key] = dists[key]
                                           
    @staticmethod
    def roundCeferOut(edge2val):
        """Reads Graph combined with parameters
        Args:
          edge2val: edges to values
        Returns:
          inferedges: set of inferred edges
        """
        return set([edge for edge in edge2vale.keys() if edge2val[edge] >= 0.000001])
    
    @staticmethod
    def readCeferOut(outfile,evol):
        """Reads CEFER out
        Args:
          outfile: CEFER outfile
          evol: evol
        Returns:
          resdict: result hash
        """
        resdict = {}
        with open(outfile,"r") as infile:
           for line in infile:
               node1,node2,val = line.rstrip().split(" ")
               resdict[(int(node1), int(node2))] = float(val)
        return resdict

    @staticmethod
    def writeCeferOut(outfile,resdict):
        with open(outfile,"w") as outfile:
             outfile.write("\n".join(["{0} {1} {2}".format(node1,node2,val) for (node1,node2),val in resdict.items()]) + "\n")

    @staticmethod
    def writeHistoryResult(history,infermode,resultfile,wformat="pkl"):
        """writes discrete greedy solution to file 
          Args:
            history:
            infermode:
            resultfile:
            format:
        """
        print(resultfile)
        exit(1)
        if wformat == "pkl":
           with gzip.open(resultfile,"wb") as outfile:
                cPickle.dump(history)
        elif wformat == "plain":
           pass 

    @staticmethod
    def readHistoryResult(resultfile,infermode,rformat):
        """reads discrete greedy output from file 
          Args:
            resultfile: resultfile to be read
            infermode:
            rformat: read format
          Returns:
            history: 
        """
        if rformat == "pkl":
           with open(resultfile,"rb") as infile:
                history = cPickle.load(infile)
        return history

    @staticmethod
    def readGraphAndParams(graphpath,outformat="pickle"):
       """Reads Graph combined with parameters
       Args:
          graphpath: graph file
       Returns:
          G: directed graph with attributes
       """
       if graphpath.endswith("pickle"):
          outformat = "pickle"
       elif graphpath.endswith(".gml"):
          outformat = "gml" 
       if outformat == "pickle": 
          return nx.DiGraph(nx.read_gpickle(graphpath))
       elif outformat == "gml":
          return nx.DiGraph(nx.read_gml(graphpath))

    @classmethod  
    def WriteGraphAndParams(self,G,dists,outfile,smodel,format="pickle"):
       """Writes Graph combined with parameters as weighted graph to outfile
       Args:
         G: graph
         dists: model distributions
         outfile: graph filename
         smodel: spreading model
       """
       self.addAttribute2Graph(G,dists)
       if format == "pickle": 
          nx.write_gpickle(G,outfile)
       elif format == "gml":
          nx.write_gml(G,outfile)

    @staticmethod
    def readHistoryOutput(outfile):
        with gzip.open(outfile,"rb") as outfile:
             return cPickle.load(outfile)

    @staticmethod    
    def writeHistoryOutput(outfile,scores):
        with gzip.open(outfile,"wb") as outfile:
             cPickle.dump(outfile,scores)

    @staticmethod
    def writeTrace(trace,tracefilename,outformat="plain"):
        """writes traces to outfilename
        Args:
           trace: trace
           outfilename: trace outfilename
           smodel:
        """
        methodname = "write{0}Trace".format(outformat.capitalize())
        getattr(InputOutput,methodname)(trace,outfilename)
    
    @staticmethod
    def writePklTrace(trace,outfilename):
       """writes traces to outfilename
       Args:
          trace: trace
          outfilename: trace outfilename
          smodel:
       """
       with gzip.open(outfilename,"wb") as outfile:
           cPickle.dump(trace,outfile)

    @staticmethod
    def writePlainTrace(trace,outfilename):
       """Writes traces to readable text file as plain
       Args:
          trace: trace
          outfilename:
          smodel:
       """
       with open(outfilename,"w") as file:
          if Trace.IsTraceNoisy([trace]):
             sortedtimes = sorted(trace[trace.keys()[0]].keys())
             tracestr = "\n".join(["{0} {1} {2} {3}".format(node,time,state,trace[node][time][state]) for node in trace.keys() for time in sortedtimes for state in trace[node][time].keys()])
          else:
             unsortednodes = []
             for node in trace.keys():
                 if Trace.INFECTED in trace[node]:
                    unsortednodes.append((node,trace[node][Trace.INFECTED]))
                 elif Trace.EXPOSED in trace[node]:
                    unsortednodes.append((node,trace[node][Trace.EXPOSED]))
             sortednodes = [node for node,time in sorted(unsortednodes, key = lambda element : element[1])]
             assert len(set(sortednodes).intersection(set(trace.keys()))) == len(sortednodes)
             tracestr = "\n".join(["{0} {1} {2}".format(node,state,trace[node][state]) for node in sortednodes for state in trace[node].keys()])  
          file.write("{0}\n".format(tracestr))

    @staticmethod
    def readPlainTrace(tracefile):
        """reads given tracefile
        Args:
          tracefile: trace file
        """
        trace = {}
        with open(tracefile,"r") as file:
          for line in file:
              if len(line.rstrip().split("\n")[0].split(" ")) == 3:
                 noisy = False 
              elif len(line.rstrip().split("\n")[0].split(" ")) == 4:
                 noisy = True
        with open(tracefile,"r") as file:
          for line in file:
              for item in line.rstrip().split("\n"):
                  splitted = item.split(" ")
                  if noisy:
                     trace.setdefault(int(splitted[0]),{})
                     trace[int(splitted[0])].setdefault(int(splitted[1]),{})
                     trace[int(splitted[0])][int(splitted[1])][splitted[2]] = float(splitted[3])
                  else:
                     trace.setdefault(int(splitted[0]),{})
                     trace[int(splitted[0])][splitted[1]] = float(splitted[2])
        return trace

    @staticmethod
    def readPklTrace(tracefile):
        """ reads given tracefile
        """
        with gzip.open(tracefile,"rb") as file:
            return cPickle.load(file)

    @staticmethod
    def readTraces(tracefolder,tracecount,smodel,mininfected=-1,informat="plain",rettype="list"):
        """reads tracecount many traces
        Args:
          tracefolder: tracefolder
          tracecount: number of traces wanted
          smodel: spreading model
          mininfected: minimum number of infected nodes in each trace
          informat: trace format
          rettype = return type
        Returns:
          traces: traces as dict
        """
        assert rettype in ["list","dict"] 
        if smodel in ["si","sir","sis"]:
           statekey = Trace.INFECTED
        elif smodel == "seir":
           statekey = Trace.EXPOSED
        allnodes = set()
        if rettype == "list":
           traces = []
        elif rettype == "dict":
           traces = {} 
        filenames = myutil.listfiles(tracefolder)
        random.shuffle(filenames)  
        for filename in filenames:
            tracefilepath = "{0}/{1}".format(tracefolder,filename)
            method = getattr(InputOutput,"read{0}Trace".format(informat.capitalize()))
            trace = method(tracefilepath)
            if len([node for node in trace.keys() if statekey in trace[node]]) >= mininfected:
               if rettype == "list": 
                  traces.append(trace)
               elif rettype == "dict":
                  traces[filename] = trace 
               allnodes |= set(trace.keys())
            if len(traces) == tracecount:
               return (traces, allnodes) 
        print("Error!! NOT ENOUGH TRACES under {0}".format(tracefolder))

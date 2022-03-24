#  
# This code generates traces for SEIR models automatically
#
import networkx as nx
import random
import os
import sys
sys.path.append("./lib")
import math
import myutilities as myutil
import pickle
import gzip
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from Pbs import Pbs

def genTraceConfigFile(configfile,vararr):
    """generates trace config file
    Args:
       configfile:
       vararr:
    """
    with open(configfile,"w") as file:
        for var in vararr.keys():
            if var == "dists":
               for key in vararr[var].keys():
                   if key in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
                      paramstr = " ".join([str(item) for item in vararr[var][key][1]])
                      file.write("{0}: {1} {2}\n".format(key,vararr[var][key][0],paramstr))
                   elif key == Trace.SPROB:
                      file.write("{0}: {1}\n".format(key,vararr[var][key])) 
            else:
               file.write("{0}: {1}\n".format(var,vararr[var]))

def getDistParams(smodel,prob,realdata,evol):
    spreadparams = []
    if smodel == "si":
       if realdata == "real": 
          sparam = {Trace.S2I: ("weibull", (9.5, 2.3)) , Trace.SPROB: 0.08}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull", (9.5, 2.3)) , Trace.SPROB: 0.1}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull", (9.5, 2.3)) , Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull",(9.5,2.3)) , Trace.SPROB: 0.3}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull",(9.5,2.3)) , Trace.SPROB: 0.5}
          spreadparams.append(sparam)
       elif realdata == "syn":
          sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (1.0,)) , Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (0.5,)) , Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (0.5,)) , Trace.SPROB: 0.3}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.SPROB: 0.3}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (1.0,)) , Trace.SPROB: 0.3}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (0.5,)) , Trace.SPROB: 0.4}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.SPROB: 0.4}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (1.0,)) , Trace.SPROB: 0.4}
          spreadparams.append(sparam)
    elif smodel == "sir":
       if realdata == "real":
          sparam = {Trace.S2I: ("weibull", (9.5, 2.3)) , Trace.I2R: ("expo", (0.1,)), Trace.SPROB: 0.1}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull", (9.5, 2.3)) , Trace.I2R: ("expo", (0.1,)), Trace.SPROB: 0.2}
          spreadparams.append(sparam)
       elif realdata == "syn":
          sparam = {Trace.S2I: ("expo", (1.0,)) , Trace.I2R: ("expo", (1.0,)) , Trace.SPROB: 0.08}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.I2R: ("expo", (1.5,)) , Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.I2R: ("expo", (1.0,)) , Trace.SPROB: 0.1}
          spreadparams.append(sparam)
    elif smodel == "seir":
       if realdata == "real": 
          sparam = {Trace.S2E: ("weibull", (9.5, 2.3)) , Trace.E2I: ("weibull", (1.1, 2.1)) , Trace.I2R: ("expo", (0.05,)), Trace.SPROB: 0.1}
          spreadparams.append(sparam)
          sparam = {Trace.S2E: ("weibull", (9.5, 2.3)) , Trace.E2I: ("weibull", (1.1, 2.1)) , Trace.I2R: ("expo", (0.05,)), Trace.SPROB: 0.2}
          spreadparams.append(sparam)
          sparam = {Trace.S2E: ("weibull", (9.5, 2.3)) , Trace.E2I: ("weibull", (1.1, 2.1)) , Trace.I2R: ("expo", (0.05,)), Trace.SPROB: 0.08}
          spreadparams.append(sparam)
       elif realdata == "syn":
          sparam = {Trace.S2E:("expo", (1.0,)) , Trace.I2R: ("expo", (1.0,)) , Trace.E2I: ("expo", (2.0,)), Trace.SPROB: 0.1}
          spreadparams.append(sparam)
          sparam = {Trace.S2E: ("expo", (1.0,)) , Trace.I2R: ("expo", (1.0,)) , Trace.E2I: ("expo", (2.0,)), Trace.SPROB: 0.08}
          spreadparams.append(sparam)
          sparam = {Trace.S2E: ("expo", (1.0,)) , Trace.I2R: ("expo", (1.0,)) , Trace.E2I: ("expo", (2.0,)), Trace.SPROB: 0.2}
          spreadparams.append(sparam)
    elif smodel == "sis":
       sparam = {Trace.S2I: ("expo", (2.0,)) , Trace.I2S: ("expo", (1.5,)) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
       sparam = {Trace.S2I: ("expo", (1.0,)) , Trace.I2S: ("expo", (1.5,)) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
       sparam = {Trace.S2I: ("expo", (0.5,)) , Trace.I2S: ("expo", (1.5,)) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
    return spreadparams

def genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,configfolderpref,configpref,realdata,smodel,evol,prob,pbsfolder):
    graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,graphpref)
    tracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(tracefolderpref,tracepref,realdata,evol,smodel,"edge",prob)
    configfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(configfolderpref,configpref,realdata,evol,smodel,"edge",prob)
    [os.makedirs(folder) for folder in [pbsfolder,configfolder,tracefolder] if not os.path.exists(folder)]
    return graphfolder,tracefolder,configfolder

def returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount):    
    path2info = {}
    if realdata == "real" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
           filepath = "{0}/{1}".format(graphfolder,filename)
           filetracefolder = "{0}/{1}".format(tracefolder,filename)
           path2info[filepath] = (filename,filetracefolder)       
    elif realdata == "syn" and evol == "static":
       for data in myutil.listdirectories(graphfolder):
           datadirname = "{0}/{1}".format(graphfolder,data)
           for graphalgo in myutil.listdirectories(datadirname):
               graphdirname = "{0}/{1}".format(datadirname,graphalgo)
               filenames = myutil.listfiles(graphdirname) 
               assert len(filenames) == filesamplecount
               for filename in filenames:
                   filepath = "{0}/{1}".format(graphdirname,filename)
                   filetracefolder = "{0}/{1}/{2}/{3}".format(tracefolder,data,graphalgo,filename)
                   path2info[filepath] = (filename,filetracefolder)
    elif realdata == "real" and evol == "dynamic":
       for dirname in myutil.listdirectories(graphfolder):
           dirpath = "{0}/{1}".format(graphfolder,dirname)
           dirtracefolder = "{0}/{1}".format(tracefolder,dirname)
           path2info[dirpath] = (dirname,dirtracefolder)
    elif realdata == "syn" and evol == "dynamic":
       for data in myutil.listdirectories(graphfolder):
           dirname = "{0}/{1}".format(graphfolder,data)
           for algo in myutil.listdirectories(dirname):
               algodirname = "{0}/{1}".format(dirname,algo)
               for param in myutil.listdirectories(algodirname):
                   paramdirname = "{0}/{1}".format(algodirname,param)
                   for sampledirname in myutil.listdirectories(paramdirname):
                       samplepath = "{0}/{1}".format(paramdirname,sampledirname)
                       dirtracefolder = "{0}/{1}".format(tracefolder,sampledirname)
                       path2info[samplepath] = (sampledirname,dirtracefolder)
    return path2info


def main():    
    """automatically generates traces by calling tracegen.py(assumes trace info is given in graphfile)
    """
    graphpref = "graphs"
    tracepref = "traces"
    configpref = "config"
    realdata = "real"
    evol = "static"
    noise = 0.0
    noiseshape = "uniform" #"add"
    samplerate = 0
    tracefolderpref = "."
    graphfolderpref = "."
    configfolderpref = "."
    pbsfolder = "pbsfolder"
    smodel = "si" #"seir","sir","sis","model2"
    prob = "cont" #"dis"
    filesamplecount = 10 #for syn data
    startcount = 1
    tracecount = 300
    sismaxtime = 100
    
    if prob == "cont" and samplerate == 0:
       assert noise == 0.0
    graphfolder,tracefolder,configfolder = genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,configfolderpref,configpref,realdata,smodel,evol,prob,pbsfolder)
    path2info = returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount)
    spreadparams = getDistParams(smodel,prob,realdata,evol)
    for path in path2info.keys():
        Gname,filetracefolder = path2info[path]
        extensionstr = "-".join(path.replace("../","").split("/")[1:])
        for spreadparam in spreadparams:
            sparamstr = Trace.getSpreadFolder(smodel,spreadparam,prob)
            vararr = {}
            vararr["ongraph"] = False
            vararr["traceoutfolder"] = filetracefolder
            vararr["graphfilename"] = path
            vararr["tracecount"] = tracecount
            vararr["samplerate"] = samplerate
            vararr["startnodecount"] = startcount
            vararr["noise"] = noise
            vararr["noiseshape"] = noiseshape
            vararr["smodel"] = smodel
            vararr["evol"] = evol
            vararr["dists"] = spreadparam
            vararr["dist"] = prob
            if smodel == "sis":
               vararr["sismaxtime"] = sismaxtime
            configpath = "{0}/{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}.config".format(configfolder,extensionstr,realdata,evol,smodel,prob,tracecount,sparamstr,samplerate,noise)
            genTraceConfigFile(configpath,vararr)
            code = "python tracegen.py {0}".format(configpath)
            #os.system(code)
            #continue
            pbsfilename = "{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}.pbs".format(realdata,evol,smodel,prob,extensionstr,tracecount,sparamstr,samplerate,noise)
            Pbs.submitPbs(code,pbsfolder,pbsfilename,"pool2")
   
if __name__ == "__main__":
    main()

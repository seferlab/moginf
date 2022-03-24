#  
# This code runs cefer automatically
#
import networkx as nx
import random
import os
import sys
sys.path.append("./lib")
import math
import numpy as np
import myutilities as myutil
import cPickle
import gzip
import itertools
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
import OtherAlgos
from Pbs import Pbs

def genOtherConfigFile(configpath,algo,algoparams,tracefolder,resultfilename,tracecount,spreadparam,iprobmethod,runfolder,noise,path,smodel):
    """generates config file for other algo
    """
    with gzip.open(configpath,"wb") as file:
        cPickle.dump(algo,file)
        cPickle.dump(algoparams,file)
        cPickle.dump(tracefolder,file)
        cPickle.dump(resultfilename,file)
        cPickle.dump(tracecount,file)
        cPickle.dump(spreadparam,file)
        cPickle.dump(iprobmethod,file)
        cPickle.dump(runfolder,file)
        cPickle.dump(noise,file)
        cPickle.dump(path,file)
        cPickle.dump(smodel,file)

def genCeferConfigFile(configfile,vararr):
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

def returnalgos(smodel,evol,noise,samplerate):
    """returns all algos to be run
    """
    lambda1s = [0.01,0.05,0.1,0.15,0.2,0.3,0.5,0.75,1.0]
    lambda2s = [100000.0]
    fusedlambdas = [0.1]
    edgeratios = [1.0]
    algos = []
    if evol == "static":
       if noise == 0.0: 
           edgestaticoptmyalgos = [("abse",None,"cover",None)]
       else:
          edgestaticoptmyalgos = [("abse",None,"cover","Kernel")]
       for algo in edgestaticoptmyalgos:
           if algo[1] == "l1":
              algos.extend([(algo,{"lambda1": lambda1}) for lambda1 in lambda1s])  
           elif algo[1] == "l2":
              algos.extend([(algo,{"lambda2": lambda2}) for lambda2 in lambda2s])
           elif algo[1] == "l1l2":
              algos.extend([(algo,{"lambda1": lambda1, "lambda2": lambda2}) for lambda1,lambda2 in list(itertools.product(lambda1s,lambda2s))]) 
	   elif algo[1] == None:
              algos.append((algo,{}))
       edgestaticcovermyalgos = []
       edgestaticotheralgos = ["Multitree","Netinf"]
       for algo in edgestaticotheralgos:
           if algo in ["Multitree","Netinf"]:
              algos.extend([((algo,),{"degree": param}) for param in edgeratios]) 
           elif algo == "Netrate":
              algos.extend([((algo,),{})])
    return algos

def getRunParams(smodel,prob,realdata,evol):
    spreadparams = []
    if smodel == "si":
       if realdata == "real":
	  sparam = {Trace.S2I: ("weibull",(9.5,2.3)) , Trace.SPROB: 0.08}
          spreadparams.append(sparam) 
          sparam = {Trace.S2I: ("weibull",(9.5,2.3)) , Trace.SPROB: 0.1}
          spreadparams.append(sparam)
          sparam = {Trace.S2I: ("weibull",(9.5,2.3)) , Trace.SPROB: 0.2}
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
       sparam = {Trace.S2I: ("expo", 2.0) , Trace.I2S: ("expo", 1.5) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
       sparam = {Trace.S2I: ("expo", 1.0) , Trace.I2S: ("expo", 1.5) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
       sparam = {Trace.S2I: ("expo", 0.5) , Trace.I2S: ("expo", 1.5) , Trace.SPROB: 0.1}
       spreadparams.append(sparam)
    return spreadparams

def getAlgostr(algo,algoparams):
    """ get algo str
    """
    return "_".join([str(item) for item in list(algo)] + [str(item) for items in algoparams.items() for item in items])


def runCefer(tracefields,algofields,graphfields,iofields,extrafields):        
    """ infers by CEFER 
    """
    frac,spreadparam,smodel,evol,prob,samplerate,noise,iprobmethod = tracefields
    algo,algoparams = algofields
    (nodecount,) = graphfields
    extensionstr,filetracefolder,resultfolder,runfolder,configfolder,path,pbsfolder = iofields
    complete,completeupto,runmode,printscore,parallelcount = extrafields
    tracecount = int(math.floor(nodecount*frac))
    sparamstr = Trace.getSpreadFolder(smodel,spreadparam,prob)
    noisestr = "{0}-{1}".format(noise,samplerate)
    algostr = getAlgostr(algo,algoparams)
    tresultfolder = "{0}/{1}/{2}/{3}/{4}/{5}".format(resultfolder,extensionstr,sparamstr,algostr,frac,noisestr)
    if not os.path.exists(tresultfolder):
       os.makedirs(tresultfolder)
    indices = set([-1] + [int(myfile.replace(".edg","")) for myfile in myutil.listfiles(tresultfolder) if myfile.find("edg_par")==-1])
    if complete and len(indices) == completeupto + 1 :
       return
    uncomp = [(int(myfile.split("_outof_")[0].split("edg_par")[1]),int(myfile.split("_outof_")[1]),myfile.split("_par")[0],int(myfile.split(".edg_")[0])) for myfile in myutil.listfiles(tresultfolder) if myfile.find("edg_par")!=-1]
    if len(uncomp) != 0:
       assert runmode == "parallel"
       maxpart = max([index for index,count,fname,code in uncomp])
       parcount = list(set([count for index,count,fname,code in uncomp]))[0]
       resultfilename = "{0}/{1}".format(tresultfolder,list(set([fname for index,count,fname,code in uncomp]))[0])
       curindex = list(set([code for index,count,fname,code in uncomp]))[0]
    else:
       resultfilename = "{0}/{1}.edg".format(tresultfolder,max(indices)+1)
       curindex = max(indices)+1
    runfolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(runfolder,extensionstr,sparamstr,algostr,frac,noisestr,curindex)
    if not os.path.exists(runfolder):
       os.makedirs(runfolder) 
    tracefolder = "{0}/{1}/{2}".format(filetracefolder,sparamstr,noisestr)
    vararr = {"graphfilename": path, "dist": prob, "dists": spreadparam}
    for item in ["tracefolder","tracecount","samplerate","noise","smodel","runfolder","evol","iprobmethod","printscore","runmode","resultfilename","parallelcount"]:
        exec('vararr["{0}"] = {0}'.format(item)) in locals(),globals()
    if evol == "static":    
       vararr["errortype"],vararr["sparsetype"],vararr["cover"],vararr["secondalgo"] = algo
    elif evol == "dynamic":
       vararr["errortype"],vararr["sparsetype"],vararr["cover"],vararr["secondalgo"],vararr["fusedtype"] = algo  
    for param in algoparams.keys():
        vararr[param] = algoparams[param]
    if len(uncomp) != 0:
       vararr["parallelcount"] = parcount
       vararr["parstartin"] = maxpart
    configpath = "{0}/{1}.config".format(configfolder,"_".join([str(key) for key in [extensionstr,sparamstr,algostr,noisestr,tracecount,curindex,iprobmethod]]))
    genCeferConfigFile(configpath,vararr)
    code = "python cefercplex.py {0}".format(configpath)
    #os.system(code)
    #return            
    pbsfilename = "{0}".format("_".join([str(key) for key in [extensionstr,sparamstr,algostr,noisestr,tracecount,curindex,iprobmethod]]))
    Pbs.submitPbs(code,pbsfolder,pbsfilename,"pool1")

def runOther(tracefields,algofields,graphfields,iofields,extrafields):    
    """ runs other algorithms
    """
    frac,spreadparam,smodel,evol,prob,samplerate,noise,iprobmethod = tracefields
    algo,algoparams = algofields
    (nodecount,) = graphfields
    extensionstr,filetracefolder,resultfolder,runfolder,configfolder,path,pbsfolder = iofields
    complete,completeupto = extrafields
               
    tracecount = int(math.floor(nodecount*frac))
    sparamstr = Trace.getSpreadFolder(smodel,spreadparam,prob)
    noisestr = "{0}-{1}".format(noise,samplerate)
    algostr = getAlgostr(algo,algoparams)
    tresultfolder = "{0}/{1}/{2}/{3}/{4}/{5}".format(resultfolder,extensionstr,sparamstr,algostr,frac,noisestr)
    if not os.path.exists(tresultfolder):
       os.makedirs(tresultfolder)
    indices = set([-1] + [int(myfile.replace(".edg","")) for myfile in myutil.listfiles(tresultfolder)])
    if complete and len(indices) == completeupto + 1:
       return
    runfolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(runfolder,extensionstr,sparamstr,algostr,frac,noisestr,max(indices)+1)
    if not os.path.exists(runfolder):
       os.makedirs(runfolder) 
    resultfilename = "{0}/{1}.edg".format(tresultfolder,max(indices)+1)
    tracefolder = "{0}/{1}/{2}".format(filetracefolder,sparamstr,noisestr)
    configpath = "{0}/{1}.config".format(configfolder,"_".join([str(key) for key in [extensionstr,sparamstr,algostr,noisestr,tracecount,max(indices)+1,iprobmethod]]))
    genOtherConfigFile(configpath,algo,algoparams,tracefolder,resultfilename,tracecount,spreadparam,iprobmethod,runfolder,noise,path,smodel)
    code = "python OtherAlgos.py {0}".format(configpath)
    #os.system(code)
    #return
    pbsfilename = "{0}".format("_".join([str(key) for key in [extensionstr,sparamstr,algostr,noisestr,tracecount,max(indices)+1,iprobmethod]]))
    Pbs.submitPbs(code,pbsfolder,pbsfilename,"pool1")
                
def genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,configfolderpref,configpref,resultfolderpref,resultpref,runfolderpref,runpref,realdata,evol,smodel,infertype,prob,pbsfolder):
    graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,graphpref)
    tracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(tracefolderpref,tracepref,realdata,evol,smodel,infertype,prob)
    configfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(configfolderpref,configpref,realdata,evol,smodel,infertype,prob)
    resultfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(resultfolderpref,resultpref,realdata,evol,smodel,infertype,prob)
    runfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(runfolderpref,runpref,realdata,evol,smodel,infertype,prob)
    [os.makedirs(folder) for folder in [pbsfolder,configfolder,tracefolder,runfolder,resultfolder] if not os.path.exists(folder)]
    return graphfolder,configfolder,tracefolder,runfolder,resultfolder

def returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount):      
    path2info = {}
    if realdata == "real" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
	   #if filename != "trial.gml":
           #   continue		   
           filepath = "{0}/{1}".format(graphfolder,filename)
           filetracefolder = "{0}/{1}".format(tracefolder,filename)
           G = InputOutput.readGraphAndParams(filepath)
           path2info[filepath] = (filename,filetracefolder,G.number_of_nodes())       
    elif realdata == "syn" and evol == "static":
       for data in myutil.listdirectories(graphfolder):
           datadirname = "{0}/{1}".format(graphfolder,data)
           for graphalgo in myutil.listdirectories(datadirname):
               graphdirname = "{0}/{1}".format(datadirname,graphalgo)
               filenames = myutil.listfiles(graphdirname) 
               assert len(filenames) == filesamplecount
               for filename in filenames:
                   myid = int(filename.replace(".gml",""))
                   if myid >= 5:
                      continue 
                   filepath = "{0}/{1}".format(graphdirname,filename)
                   filetracefolder = "{0}/{1}/{2}/{3}".format(tracefolder,data,graphalgo,filename)
                   G = InputOutput.readGraphAndParams(filepath)
                   path2info[filepath] = (filename,filetracefolder,G.number_of_nodes())
    elif realdata == "real" and evol == "dynamic":
       for dirname in myutil.listdirectories(graphfolder):
           dirpath = "{0}/{1}".format(graphfolder,dirname)
           dirtracefolder = "{0}/{1}".format(tracefolder,dirname)
           G = {int(filename.split("_")[1]): InputOutput.readGraphAndParams("{0}/{1}".format(dirtracefolder,filename)) for filename in myutil.listfiles(dirtracefolder)}
           nodecount = G[sorted(G.keys())[0]].number_of_nodes()    
           path2info[dirpath] = (G,dirname,dirtracefolder,nodecount)
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
                       G = {int(filename.split("_")[1]): InputOutput.readGraphAndParams("{0}/{1}".format(dirtracefolder,filename)) for filename in myutil.listfiles(dirtracefolder)}
                       nodecount = G[sorted(G.keys())[0]].number_of_nodes()    
                       path2info[samplepath] = (G,sampledirname,dirtracefolder,nodecount)
    return path2info

def assignProbMethod(noise,smodel):
    iprobmethods = []
    if noise == 0.0:
       iprobmethods.append(None)  
    elif smodel == "sir":
       iprobmethods.append("lse")
    elif smodel == "si":
       iprobmethods.append("lse") 
    elif smodel == "seir":
       iprobmethods.append("lse") 
    elif smodel == "sis":
       iprobmethods.append(None) 
    return iprobmethods

def checkVars(prob,samplerate,noise):
    if prob == "cont" and samplerate == 0:
       assert noise == 0.0

def main():   
    """automatically generates traces by calling tracegen.py or other algos
    """
    graphpref = "graphs"
    tracepref = "traces"
    configpref = "runconfig"
    runpref = "run"
    resultpref = "result"
    realdata = "real"
    evol = "static"
    noise = 0.0
    samplerate = 0
    resultfolderpref = "."
    runfolderpref = "."
    tracefolderpref = "."
    graphfolderpref = "."
    configfolderpref = "."
    pbsfolder = "pbsfolder"
    smodel = "si" #"seir","sir","sis","model2"
    prob = "cont" #"dis"
    filesamplecount = 10 #for syn data
    infertype = "edge"
    tracefractions = list(np.arange(0.005,0.05,0.005)) + list(np.arange(0.05,0.1,0.01)) + list(np.arange(0.1,0.5,0.1)) + list(np.arange(0.5,1.6,0.2))  
    #tracefractions = [0.005,0.01,0.015,0.02,0.025,0.03,0.35,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.5,0.6,0.75]
    #tracefractions = [0.1,0.2,0.4,0.5]
    complete = True
    completeupto = 5
    runmode = "serial"
    parallelcount = 5
    printscore = None

    checkVars(prob,samplerate,noise)   
    iprobmethods = assignProbMethod(noise,smodel)
    (graphfolder,configfolder,tracefolder,runfolder,resultfolder) = genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,configfolderpref,configpref,resultfolderpref,resultpref,runfolderpref,runpref,realdata,evol,smodel,infertype,prob,pbsfolder)
    path2info = returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount)
    spreadparams = getRunParams(smodel,prob,realdata,evol)
    algos = returnalgos(smodel,evol,noise,samplerate)
    paramcombs = list(itertools.product(*[algos,spreadparams,tracefractions,iprobmethods]))
    for path in path2info.keys():
        Gname,filetracefolder,nodecount = path2info[path]
        extensionstr = "-".join(path.replace("./","").split("/")[1:])
        for (algo,algoparams),spreadparam,frac,iprobmethod in paramcombs:
            if algo[0] in ["Multitree","Netinf","Netrate","Connie"]:
               tracefields = [frac,spreadparam,smodel,evol,prob,samplerate,noise,iprobmethod]
               algofields = [algo,algoparams]
               graphfields = [nodecount]
               iofields = [extensionstr,filetracefolder,resultfolder,runfolder,configfolder,path,pbsfolder]
               extrafields = [complete,completeupto]
               runOther(tracefields,algofields,graphfields,iofields,extrafields)
            else:
               tracefields = [frac,spreadparam,smodel,evol,prob,samplerate,noise,iprobmethod]
               algofields = [algo,algoparams]
               graphfields = [nodecount]
               iofields = [extensionstr,filetracefolder,resultfolder,runfolder,configfolder,path,pbsfolder]
               extrafields = [complete,completeupto,runmode,printscore,parallelcount]
               runCefer(tracefields,algofields,graphfields,iofields,extrafields)
              

if __name__ == "__main__":
    main()




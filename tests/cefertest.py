import random
import unittest
import sys
sys.path.append("../lib")
sys.path.append("..")
import networkx as nx
import myutilities as myutil
import os
import math
import itertools
import operator
import autotracegen
import autoceferrun
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from DiscreteDistribution import DiscreteDistribution as DisDist
from ContinuousDistribution import ContinuousDistribution as ContDist
from Score import Score

class CeferTest(unittest.TestCase):
    """ Cefer test class
    """
    def __init__(self,testname,traces=None,smodel=None,noise=None,samplerate=None,G=None,spreadparam=None,algo=None,algoparams=None,iprobmethod=None,resultfolder=None,evol=None): 
       """constructor
         Args:
            testname: test to run
            traces: trace datas dictionary
            smodel: spreading model
            distinfo: distribution information
            G: graph
            noise:
            samplerate:
            prob:
       """
       super(TraceTest, self).__init__(testname)
       self.traces = traces
       self.smodel = smodel
       self.spreadparam = spreadparam
       self.noise = noise
       self.samplerate = samplerate
       self.prob = prob
       self.G = G
       self.algo = algo
       self.evol = evol
       self.algoparams = algoparams
       self.iprobmethod = iprobmethod
       self.normalresultfolder = normalresultfolder
       self.matrixresultfolder = matrixresultfolder

    def testParallelComp(self):
       """ tests parallel execution generates same results with serial one
       """
       return

    def testCeferValidity(self):
       """ checks validity of results
       """
       #read scores from normal output 
       normalresultfiles = [myfile for myfile in myutil.listfiles(self.normalresultfolder)]
       self.assertEqual(len(normalresultfiles),1)
       normalresultfile = "{0}/{1}".format(self.normalresultfolder,normalresultfiles[0])
       matrixresultfiles = [myfile for myfile in myutil.listfiles(self.matrixresultfolder)]
       self.assertEqual(len(matrixresultfiles),1)
       matrixfile = "{0}/{1}".format(self.matrixresultfolder,matrixresultfiles[0])
       truthdict = {0: set(self.G.edges())}
       normedge2val = InputOutput.readCeferOut(normalresultfile,self.evol)
       matedge2val = InputOutput.readCeferOut(matrixresultfile,self.evol)
       roundmethods = [("epsilon",value) for value in [0.2,0.999999999999,0.00001,0.05,0.1,0.15,0.25,0.35,0.3,0.4]]
       normallscores,inferdicts = {},[]
       for method,params in roundmethods:
           norminferdict = Score.roundScore2Set(normedge2val,method,params,self.evol)
           matinferdict = Score.roundScore2Set(matedge2val,method,params,self.evol)
           self.assertEqual( len(set(norminferdict[0]).difference(set(matinferdict[0]))), 0)
           self.assertEqual( len(set(matinferdict[0]).difference(set(norminferdict[0]))), 0)
           normcurscores = Score.getBinaryScores(truthdict,norminferdict,"all")
           for score in normcurscores.keys():
               normallscores.setdefault(score,set()).add(normcurscore[key])
       inferdicts.append(norminferdict) #equality
       normmaxscores = {score: max(list(normallscores[score])) for score in normallscores.keys()}    
       matmaxscores = {score: max(list(matallscores[score])) for score in matallscores.keys()}
      
    def testParallelReturn(self):
       """ checks validity of results
       """
       return
       
    def testOtherValidity(self):
       """ checks validity of results
       """
       return

    def testOtherValidity(self):
       return

    def testScoreEstimation(self):
       return
   
    def testCompleteUpto(self):
       return

def getManyDist(smodel,prob):
    if prob == "cont": 
       dists = [("expo",(1.0,)),("expo", (2.0,)),("expo", (3.0,)),("expo", (4.0,)),("expo", (5.0,)),("expo",(0.5,)),("weibull", (2.5, 9.5)) ,("weibull", (1.5, 4.5)) ,("weibull", (2.5, 4.5)),("weibull", (3.5, 5.5)), ("rayleigh",(0.5,)),("rayleigh",(1.0,)),("rayleigh",(1.5,)),("rayleigh",(2.0,))]
    elif prob == "dis":
       dists = [("expo",(1.0,)),("expo", (2.0,)),("expo", (3.0,)),("expo", (4.0,)),("expo", (5.0,)),("expo",(0.5,)),("weibull", (2.5, 9.5)) ,("weibull", (1.5, 4.5)) ,("weibull", (2.5, 4.5)),("weibull", (3.5, 5.5)), ("rayleigh",(0.5,)),("rayleigh",(1.0,)),("rayleigh",(1.5,)),("rayleigh",(2.0,))]
    sprobs = [0.05,0.08,0.1,0.2,0.25,0.5,0.75]
    random.shuffle(dists)
    random.shuffle(sprobs)
    spreadparams = []
    if smodel == "si":
       for s2i,sprob in itertools.product(dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "sir":
       for s2i,i2r,sprob in itertools.product(dists,dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.I2R: i2r, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "seir":
       for s2e,e2i,i2r,sprob in itertools.product(dists,dists,dists,sprobs):
           sparam = {Trace.S2E: s2e, Trace.E2I: e2i, Trace.I2R: i2r, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "sis":
       for s2i,i2s,sprob in itertools.product(dists,dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.I2S: i2s,Trace.SPROB: sprob}
           spreadparams.append(sparam)
    return spreadparams

  
def getTests(traces,smodel,noise,samplerate,prob,G,spreadparam,algo,algoparams,iprobmethod,resultfolder):
    """ generates tests for CEFER
    """
    suiteFew = unittest.TestSuite()
    suiteFew.addTest(CeferTest("testParallelComp",traces,smodel,noise,samplerate,prob,G,spreadparam,algo,algoparams,iprobmethod,resultfolder))
    suiteFew.addTest(CeferTest("testParallelReturn",traces,smodel,distinfo,G,noise,samplerate,prob))
    suiteFew.addTest(TraceTest("testCeferValidity",traces,smodel,distinfo,G,noise,samplerate,prob)) #via matrix
    suiteFew.addTest(TraceTest("testOtherValidity",traces,smodel,distinfo,G,noise,samplerate,prob))
    suiteFew.addTest(TraceTest("testScoreEstimation",traces,smodel,distinfo,G,noise,samplerate,prob))
    suiteFew.addTest(TraceTest("testCompleteUpto",traces,smodel,distinfo,G,noise,samplerate,prob))
    return suiteFew
    

def execTrace(path,traceconfigfolder,sparamstr,realdata,smodel,evol,prob,pbsfolder,frac,samplerate,noise,sismaxtime,startcount,spreadparam,filetracefolder):       
    if path.endswith(".gml"):
       G = InputOutput.readGraphAndParams(path,"gml")
    else:
       G = InputOutput.readGraphAndParams(path)
    tracecount = int(round(G.number_of_nodes() * frac))  
    extensionstr = "-".join(path.replace("../","").split("/")[1:])
    vararr = {}
    vararr["ongraph"] = False
    vararr["traceoutfolder"] = filetracefolder
    vararr["graphfilename"] = path
    vararr["tracecount"] = tracecount
    vararr["samplerate"] = samplerate
    vararr["startnodecount"] = startcount
    vararr["noise"] = noise
    vararr["smodel"] = smodel
    vararr["evol"] = evol
    vararr["dists"] = spreadparam
    vararr["dist"] = prob
    if smodel == "sis":
       vararr["sismaxtime"] = sismaxtime
    configpath = "{0}/{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}.config".format(traceconfigfolder,extensionstr,realdata,evol,smodel,prob,tracecount,sparamstr,samplerate,noise)
    autotracegen.genTraceConfigFile(configpath,vararr)
    code = "python ../tracegen.py {0}".format(configpath)
    os.system(code)
    [os.system("rm -rf {0}".format(folder)) for folder in [configpath]]
     

def returnalgos(smodel,evol):
    """returns all algos and parameters to be tested
    """
    lambda1s = [0.01,0.05,0.1,0.15,0.2,0.3,0.5,0.75,1.0]
    lambda2s = [100000.0]
    fusedlambdas = [0.1]
    edgeratios = [0.05,0.2,0.5,1.0,2.5]
    algos = []
    if evol == "static":
       edgestaticoptmyalgos = [("abse","l1","cover"),("abse","l1",None),("abse",None,"cover")]
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

def main():
    """ main function for testing CEFER
    """
    samplerates = [0,1,2,4,8]
    noises = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,0.75]
    realdatas = ["real"]
    evols = ["static"]
    probs = ["cont","dis"]
    smodels = ["si"]
    graphpref = "graphs"
    tracepref = "testcefertraces"
    traceconfigpref = "testceferconfig"
    ceferconfigpref = "testceferconfig"
    runpref = "testrun"
    resultpref = "testresult"
    tracefolderpref = "."
    graphfolderpref = ".."
    traceconfigfolderpref = "."
    ceferconfigfolderpref = "."
    resultfolderpref = "."
    runfolderpref = "."
    pbsfolder = "pbsfolder"
    startcount = 1
    sismaxtime = 100
    filesamplecount = 20 #for syn data
    infertype = "edge"
    fracs = [0.005,0.01,0.015,0.02,0.025,0.05,0.1,0.2,0.25,0.5] 
    completes = [True,False]
    completeuptos = [1,2,5,10]
    runmodes = ["serial","parallel"]
    parallelcounts = [1,5,10,50,100]
    printscore = None
    
    runblocks = [("serial",None)] + [("parallel",par) for par in parallelcounts]
    parts = list(itertools.product(samplerates,noises,realdatas,evols,fracs,probs,smodels,completes,completeuptos,runblocks))
    random.shuffle(parts)
    for (samplerate,noise,realdata,evol,frac,prob,smodel,complete,completeupto,runblock) in parts:
        tracecount = 750 * frac
        if samplerate == 0 and noise != 0.0:
           continue
        for spreadparam in getManyDist(smodel,prob):
            sparamstr = Trace.getSpreadFolder(smodel,spreadparam,prob)
            graphfolder,tracefolder,traceconfigfolder = autotracegen.genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,traceconfigfolderpref,traceconfigpref,realdata,smodel,evol,prob,pbsfolder)
            path2info = autotracegen.returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount)
            for path in path2info.keys():
                Gname,filetracefolder = path2info[path]
                execTrace(path,traceconfigfolder,sparamstr,realdata,smodel,evol,prob,pbsfolder,frac,samplerate,noise,sismaxtime,startcount,spreadparam,filetracefolder)
                runmode,parallelcount = runblock
                #checkVars(prob,samplerate,noise)   
                iprobmethods = autoceferrun.assignProbMethod(noise,smodel)
                (graphfolder,ceferconfigfolder,tracefolder,runfolder,resultfolder) = autoceferrun.genMainFolders(graphfolderpref,graphpref,tracefolderpref,tracepref,ceferconfigfolderpref,ceferconfigpref,resultfolderpref,resultpref,runfolderpref,runpref,realdata,evol,smodel,infertype,prob,pbsfolder)
                path2info = autoceferrun.returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount)
                Gname,filetracefolder,nodecount = path2info[path]
                algos = returnalgos(smodel,evol)
                paramcombs = list(itertools.product(*[algos,iprobmethods]))
                extensionstr = "-".join(path.replace("./","").split("/")[1:])
                for (algo,algoparams),iprobmethod in paramcombs:
                    if algo[0] in ["Multitree","Netinf","Netrate","Connie"]:
                       autoceferrun.runOther(frac,spreadparam,nodecount,algo,algoparams,extensionstr,filetracefolder,iprobmethod,resultfolder,runfolder,ceferconfigfolder,path,smodel,evol,prob,complete,completeupto,samplerate,noise)
                    else:
                       autoceferrun.runCefer(frac,spreadparam,nodecount,algo,algoparams,extensionstr,filetracefolder,iprobmethod,resultfolder,runfolder,ceferconfigfolder,path,smodel,evol,prob,complete,completeupto,runmode,printscore,parallelcount,samplerate,noise)
                    #[os.system("rm -rf {0}".format(folder)) for folder in [configpath]]
                    noisestr = "{0}-{1}".format(noise,samplerate)
                    tresultfolder = "{0}/{1}/{2}/{3}/{4}/{5}".format(resultfolder,extensionstr,sparamstr,algostr,frac,noisestr)
                    suiteFew = getTests(traces,smodel,noise,samplerate,prob,G,spreadparam,algo,algoparams,iprobmethod,tresultfolder)
                    unittest.TextTestRunner(verbosity = 2).run(suiteFew)
                
    
if __name__ == '__main__':
    main()

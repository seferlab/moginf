#Our difference from other papers is that we give the formulation from system dynamics point of view which can be then be easily be solved by compressed sensing algorithms.
#Effectiveness of compressed sensing algorithms are proven by the quality of solution with less sample. 
#Compressed sensing's effectiveness can also be explained by its affectiveness over random matrices(in offline case, we might assume the traces will produce a random matrix).
#Better solutions(with less sample) can be obtained by constructing the compressed sensing matrix online via RIP or some other property.
#We will also discuss compressed sensing's affectiveness on temporally evolving networks which we will construct with less samples. Importance of less sample construction will be more important especially when the network is changing over time since without compressed sengsing we will theoretically need a lot of samples.

# cost of sampling: uniform cost
# nonuniform cost: 1- directly related to expected time of spreading(t). After t time steps, cascade will die out.
# 2- Directly proportional to expected to number of nodes affected. Total of number of nodes infected sometime during the cascade.
# Those will be also be estimated by cascade's parameters.

# !! We will treat the graph as a signal to recover by compressed sensing.

import networkx as nx
import numpy as np
import scipy as sp
import random
import os
import sys
import math
import myutilities as myutil
import operator
import itertools
import community
import gzip
import cPickle
sys.path.append("../tracegenerator/")
import tracegeneratorpart

#Offline plots
#1- Performance vs number of samples(either selected uniformly or by other means) for a given graph
#We can also discuss effectiveness of graph properties based rounding method in this plot.
#2- Recoverability of various growth models(dmc,lpa): This will again be performance vs number of samples plot(similar to figure 1 above)  
#3- Recoverability in terms of trace sample parameters. Recovery performance plot from just specific parameter samples.
  #For SIR, this will be heatmap showing where axises show different parameter combinations and value is the recovery performance
  #For SEIR, this will be 3-d heatmap(how to do it??) 
  #In this plot, our generalized aim is to understand which parameters are the best in terms of recovery(high recovery, high spread etc)
#Another interesting case of the problem: Assume we can see the traces but don't know their associated parameters. In this case, we should also estimate the parameters from sample trace. THIS CASE IS QUITE TOUGH so we won't focus on this right now.
#4- The case where the data is not fully observable. Performance of algorithms vs percentage 


def pbsrunner(code2run,pbsfolder,pbsfilename):
    pbsfilepath="{0}/{1}.pbs".format(pbsfolder,pbsfilename)
    errorpath="{0}/{1}.err".format(pbsfolder,pbsfilename)
    outpath="{0}/{1}.out".format(pbsfolder,pbsfilename)
    print pbsfilepath
    print errorpath
    print outpath
    file=open(pbsfilepath,"w")
    file.write("#!/bin/sh\n")
    num=random.randrange(0,2)
    num=2
    if num==0:
       queuename="workstation"
    elif num==1:
       queuename="throughput"
    elif num==2:
       queuename="high_throughput"
    else:
       print "error"
       exit(1)
    file.write("#PBS -q {0}\n".format(queuename))
    file.write("#PBS -r n\n")
    file.write("#PBS -V\n")
    #file.write("#PBS -l walltime=96:00:00\n")
    file.write("#PBS -o {0}\n".format(outpath))
    file.write("#PBS -e {0}\n".format(errorpath))
    file.write("cd $PBS_O_WORKDIR\n")
    file.write(code2run+"\n")
    file.close()
    code="qsub {0}".format(pbsfilepath)
    os.system(code)


def returnspecificfolder(specificparam):   
    if len(specificparam)==1 and specificparam[0][0]=="all":
       specificfolder="{0}".format(specificparam[0][0])
    else:
       specificfolder=""
       allstrs=[]
       for mytype,distinfo in specificparam:
           params=[mytype]
           if type(distinfo)!=tuple:
              params.append(distinfo)
           else:   
              for elem in distinfo:
                  if type(elem)==tuple:
                     params.extend(elem)
                  params.append(elem)   
           tempparamlist=[]
           for param in params:
               tempparamlist.append(str(param))
           allstrs.append("_".join(tempparamlist))
       specificfolder="-".join(allstrs)
    return specificfolder


#Global variables
pbsfolder="pbsfolder"
spreadmodel="si" #"seir", "si", "sis"
infertype="edge" #difprobpartial","difprobunknown","spreadprobunknown
datatype=("synthetic","graphwise") #datatype=[("real","dynamic"),("real","static"),("synthetic","graphwise"),("synthetic,""paramwise",),("synthetic","dynamic")]
#setcover can also be run on SEIR model!! However, each constraint will involve many variables.
fixededgestaticoptmyalgos=[("lse","","sum"),("lse","l1","sum"),("lse","l2","sum"),("lse","l1l2","sum"),("abse","","sum"),("abse","l1","sum"),("abse","l2","sum"),("abse","l1l2","sum"),("lse","","mul"),("lse","l1","mul"),("lse","l2","mul"),("lse","l1l2","mul"),("abse","","mul"),("abse","l1","mul"),("abse","l2","mul"),("abse","l1l2","mul")]
fixededgestaticcovermyalgos=["setcover"]
fixededgestaticotheralgos=["connie","netrate","netinf","multitree"]
fixededgestaticcsalgos=["elasticnet","gpsr","nuirls","l1ls","nonl1ls"] #nuirls won't be used that much
#fixededgedynamicoptmyalgos=[("lse",""),("lse","l1"),("lse","l2"),("lse","l1l2"),("abse",""),("abse","l1"),("abse","l2"),("abse","l1l2")]
fixededgedynamicoptmyalgos=[("lse","","sum","fused"),("lse","l1","sum","fused"),("lse","l2","sum","fused"),("lse","l1l2","sum","fused"),("abse","","sum","fused"),("abse","l1","sum","fused"),("abse","l2","sum","fused"),("abse","l1l2","sum","fused"),("lse","","mul","fused"),("lse","l1","mul","fused"),("lse","l2","mul","fused"),("lse","l1l2","mul","fused"),("abse","","mul","fused"),("abse","l1","mul","fused"),("abse","l2","mul","fused"),("abse","l1l2","mul","fused")]
fixededgedynamiccovermyalgos=["setcover"]
fixededgedynamicotheralgos=[]
fixededgedynamiccsalgos=["elasticnet","gpsr","nuirls","l1ls","nonl1ls"] #nuirls won't be used that much
#difprobunknown is when nothing is unknown about prob distribution(neither dist nor parameters). difprobpartial is the case where prob dist is known but parameters is unknown.
fixeddifprobpartialoptmyalgos=[("lse",""),("lse","l1"),("lse","l2"),("lse","l1l2"),("abse",""),("abse","l1"),("abse","l2"),("abse","l1l2")]
fixeddifprobpartialcovermyalgos=[]
fixeddifprobpartialotheralgos=["netrate"]
fixeddifprobpartialcsalgos=[]
fixeddifprobunknownoptmyalgos=[("lse",""),("lse","l1"),("lse","l2"),("lse","l1l2"),("abse",""),("abse","l1"),("abse","l2"),("abse","l1l2")]
fixeddifprobunknowncovermyalgos=[]
fixeddifprobunknownotheralgos=[]
fixeddifprobunknowncsalgos=[]
fixedspreadprobunknownoptmyalgos=[("lse",""),("lse","l1"),("lse","l2"),("lse","l1l2"),("abse",""),("abse","l1"),("abse","l2"),("abse","l1l2")]
fixedspreadprobunknowncovermyalgos=[]
fixedspreadprobunknownotheralgos=["connie"]
fixedspreadprobunknowncsalgos=[]


edgestaticoptmyalgos=[("lse","l1","sum"),("abse","l1","sum")]
edgestaticcovermyalgos=[]
edgestaticotheralgos=["multitree","netinf","netrate"]
edgestaticcsalgos=[]
edgedynamicoptmyalgos=[("lse","l1","sum","fused"),("abse","l1","sum","fused")] #("abse","l1","mul","fused")
edgedynamiccovermyalgos=[]
edgedynamicotheralgos=[]
edgedynamiccsalgos=[]
difprobpartialoptmyalgos=[("lse","l1"),("lse","l1l2")]
difprobpartialcovermyalgos=[]
difprobpartialotheralgos=["netrate","connie"]
difprobpartialcsalgos=[]
difprobunknownoptmyalgos=[("lse","l1"),("lse","l1l2")]
difprobunknowncovermyalgos=[]
difprobunknownotheralgos=[]
difprobunknowncsalgos=[]
spreadprobunknownoptmyalgos=[("lse","l1"),("lse","l1l2")]
spreadprobunknowncovermyalgos=[]
spreadprobunknownotheralgos=["connie"]
spreadprobunknowncsalgos=[]

#everyhing includding netrate until 0.2
#dmc is run exept 0.2 and 0.9

#tracefractions=[0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.7,0.1,0.15,0.2] #this fraction is percentage of node fractions that are available fraction count
#tracefractions=[0.03,0.04,0.05,0.1,0.12,0.15,0.2,0.25,0.3,0.5,0.75,1.0] #025,0.3
#tracefractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.4,0.45,0.5,0.6,0.8,1.0,1.2] #025,0.3

tracefractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,1.0] #for dynamic #,1.5,2.0,2.5
#tracefractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2] #for static
#tracefractions=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5] #for static

#tracefractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.4,0.45,0.5]
#tracefractions=[0.1,0.2,0.5,0.8] #025,0.3
sampleintervals=[1,2,4] #,2,4,8,12 length of sampling interval over a given diffusion sample

#length of sampling interval will be 1 for dynamic graphs
sampleselections=["uniform"] #sampleselections=["uniform","positivedegreecorrelation","negativedegreecorrelation"] #offline samples can be obtained in many ways
#One way is uniform sampling over each nodes.
#positivedegreecorrelation: high degree nodes have more samples starting from them
#negativedegreecorrelation: low degree nodes have more samples starting from them 
logbase=2.0 #will be necessary since there might be different bases for nodes affection
lognormalmean=0.0
lognormalstddev=1.0
lognormalexpocoef=2.0 #??
#rightmethods=[("rightnoise",(lognormalmean,lognormalstddev,lognormalexpocoef)),("epsilon",0.01)]
rightmethods=[("epsilon",0.01)]
graphpropertieserror1={} #graphproperties error: we will use clustering coef,scale-free exponent and modularity value
graphpropertieserror1["modularity"]=0.05
graphpropertieserror1["clusteringcoef"]=0.1
graphpropertieserror1["degreedist"]=0.1
#roundingmethods=[("epsilon",0.2),("epsilon",0.999999999999),("graphproperties",graphpropertieserror1)] #epsilonrounding rounds values below epsilon to 0 and values 1-epsilon to 1
roundingmethods=[("epsilon",0.2),("epsilon",0.999999999999),("epsilon",0.00001),("epsilon",0.05),("epsilon",0.1),("epsilon",0.15),("epsilon",0.25),("epsilon",0.35),("epsilon",0.3),("epsilon",0.4),("epsilon",0.5),("epsilon",0.45),("epsilon",0.55),("epsilon",0.6),("epsilon",0.65),("epsilon",0.7),("epsilon",0.75),("epsilon",0.8),("epsilon",0.85),("epsilon",0.9),("epsilon",0.95)]
#roundingmethods=[("epsilon",0.999999999999)]
#lambda1s=[0.1,1.0,5.0]
#lambda2s=[0.1,1.0,5.0]
#fusedlambdas=[0.1,1.0,5.0]

#lambda1s=[0.001,0.005,0.01,0.05,0.1,0.3,0.5,1.0,1.5,15.0,30.0]
lambda1s=[0.0001,0.001,0.005,0.01,0.05,0.1,0.5,1.0,1.5,5.0,10.0,15.0,30.0,50.0,100.0]
#lambda1s=[0.1]
lambda2s=[0.1]
fusedlambdas=[0.1]

#we might have another lambda3 for fused lasso penalty for dynamic graph inference but for now ignore it
if spreadmodel in ["sir","seir"]:
   infectedprobabilityassignermethods=["recoverlse"] # "flowdifference", this will be in our algoritm as a method of assignng probability to infected time.
   #recoverlse is nice general noise-independent probability assignment method(this is based on i-r flow)
elif spreadmodel in ["si","sis"]:
   infectedprobabilityassignermethods=["flowdifference"] # "flowdifference" is based to s-i flow
   
setcoverweightfunction="expodecay" #expodecay,difference,logsum
traceroundingmethods=[("categorical","")]
dumpfolderstatic="offlineinferencedumps"
resultfolderstatic="offlineresults"
runfolderstatic="offlineruns"
tracefolderstatic="traces"
filesamplecount=20 #for synthetic data 
algoruncount=3 #number of times inference algorithm will be run. Scores posted on the plot will be average of algorithm over this many counts. In each time, samples will be selected randomly and algorithm will be run.
tracestartnodenums=["all"] #only trace samples which start from that many nodes will be considered. "all" means everything will be considered.

#we will mode observability here under those variables
#partially observable degree, related to signal smoothness, right now lets also use this to model noisestatus
#if this is low, it means we are tracing(can trace, can find such as word) a word that is nice representative of infected state(get a sharper picture). If high, we can't find such thing so what we observe is quite smoothed picture of reality.
partialdegree=0.8 #1.0 is perfect case and 0.0 is the worst case(random noise)
partialcorrelationfunction="randomnoise" 
missingcompletemode=True
missingcompleteuptocount=3 #there must be 3 repetitions of each
#netinfdegreeratios=[5,8,10,15,20,25,30,50] #degree of graph for edgepercentages for netinf(iteration count)
#multitreedegreeratios=[5,8,10,15,20,25,30,50] #same as for netinf, degree of graph for edgepercentages for netinf(iteration count)

#netinfdegreeratios=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,4.5,10.0,20.0,40.0]
#multitreedegreeratios=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,4.5,10.0,20.0,40.0]
netinfdegreeratios=[0.2,0.3,0.6,0.8,1.2,3.0,5.0,10.0,20.0,30.0,40.0,50.0,100.0]
multitreedegreeratios=[0.2,0.3,0.6,0.8,1.2,3.0,5.0,10.0,20.0,30.0,40.0,50.0,100.0]

def specificparamreturner(spreadmodel,infertype,realdata,graphevolution):
    specificparamlist=[]
    #specificparam is list of params we are interested in. It may have any number of elements, specificparam gives us constraints that must be satisfied!
    if infertype=="edge":
       if spreadmodel=="si":
          #specificparam1=[("spreadprob",0.5),("s2i",("expo",4.0))]
          #specificparamlist.append(specificparam1)
          #specificparam2=[("all","all")]
          #specificparamlist.append(specificparam2)

          #specificparam5=[("spreadprob",0.1),("s2i",("powerlaw",-1.5))]
          #specificparamlist.append(specificparam5)
          #specificparam6=[("spreadprob",0.1),("s2i",("expo",4.0))]
          #specificparamlist.append(specificparam6)
          
          #specificparam6=[("spreadprob",0.05),("s2i",("expo",2.0))]
          #specificparamlist.append(specificparam6)

          #weibull here
          #specificparam7=[("spreadprob",0.5),("s2i",("weibull",3.0,2.5))]
          #specificparamlist.append(specificparam7)

          #specificparam7=[("spreadprob",0.1),("s2i",("weibull",3.0,2.5))]
          #specificparamlist.append(specificparam7)
          
          #here!! 
          specificparam7=[("spreadprob",0.5),("s2i",("expo",2.0))]
          specificparamlist.append(specificparam7)

          specificparam7=[("spreadprob",0.1),("s2i",("expo",2.0))]
          specificparamlist.append(specificparam7)

          #specificparam7=[("spreadprob",0.1),("s2i",("expo",3.0))]
          #specificparamlist.append(specificparam7)

          #specificparam7=[("spreadprob",0.1),("s2i",("expo",4.0))]
          #specificparamlist.append(specificparam7)

          #specificparam7=[("spreadprob",0.1),("s2i",("expo",6.0))]
          #specificparamlist.append(specificparam7)

          #specificparam7=[("spreadprob",0.1),("s2i",("expo",8.0))]
          #specificparamlist.append(specificparam7)
          
          #specificparam5=[("spreadprob",0.1),("s2i",("powerlaw",-1.0))]
          #specificparamlist.append(specificparam5)
          
          #specificparam8=[("spreadprob",0.2),("s2i",("expo",2.0))]
          #specificparamlist.append(specificparam8)

          #specificparam9=[("spreadprob",0.1),("s2i",("expo",4.0))]
          #specificparamlist.append(specificparam9)

          #specificparam10=[("spreadprob",0.2),("s2i",("expo",2.0))]
          #specificparamlist.append(specificparam10)
          
          #specificparam9=[("spreadprob",0.2),("s2i",("expo",2.0))]
          #specificparamlist.append(specificparam9)
          #specificparam12=[("spreadprob",0.05),("s2i",("expo",2.0))]
          #specificparamlist.append(specificparam12)
       elif spreadmodel=="sir":
          #specificparam1=[("spreadprob",0.2),("s2i",("powerlaw",-1.5)),("i2r",("expo",4))]
          #specificparamlist.append(specificparam1)
          #specificparam2=[("spreadprob",0.3),("s2i",("powerlaw",-2.0)),("i2r",("expo",2))]
          #specificparamlist.append(specificparam2)
           
          #specificparam2=[("spreadprob",0.1),("s2i",("expo",2)),("i2r",("m",(("expo",2.0),("expo",4.0)),(0.5,0.5)))]
          #specificparamlist.append(specificparam2)

          specificparam3=[("spreadprob",0.1),("s2i",("expo",2.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.1),("s2i",("expo",3.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.1),("s2i",("expo",4.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.1),("s2i",("expo",6.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.1),("s2i",("expo",8.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.5),("s2i",("expo",2.0)),("i2r",("expo",2.0))]
          specificparamlist.append(specificparam3)
          specificparam3=[("spreadprob",0.5),("s2i",("expo",2.0)),("i2r",("expo",4.0))]
          specificparamlist.append(specificparam3)
       elif spreadmodel=="seir":
          specificparam1=[("spreadprob",0.1),("s2e",("expo",2)),("i2r",("expo",2)),("e2i",("expo",2))]
          specificparamlist.append(specificparam1)
       elif spreadmodel=="sis":
          specificparam1=[("spreadprob",0.1),("s2i",("expo",2)),("i2s",("powerlaw",-1.5))]
          specificparamlist.append(specificparam1)   
       else:
          print "ERROR: No specificparamlist for spreadmodel {0}".format(spreadmodel)
          exit(1)
    else:
       print "ERROR: No specificparamlist for infertype {0}".format(infertype)
       exit(1) 
    return specificparamlist


if __name__ == "__main__":
    edgestaticallalgos=[]
    edgestaticallalgos.extend(edgestaticoptmyalgos)
    #edgestaticallalgos.extend(edgestaticcovermyalgos)
    edgestaticallalgos.extend(edgestaticotheralgos)

    edgedynamicallalgos=[]
    edgedynamicallalgos.extend(edgedynamicoptmyalgos)
    #edgedynamicallalgos.extend(edgedynamiccovermyalgos)
    #edgedynamicallalgos.extend(edgedynamicotheralgos)

    difprobpartialallalgos=[]
    difprobpartialallalgos.extend(difprobpartialoptmyalgos)
    difprobpartialallalgos.extend(difprobpartialcovermyalgos)
    difprobpartialallalgos.extend(difprobpartialotheralgos)

    difprobunknownallalgos=[]
    difprobunknownallalgos.extend(difprobunknownoptmyalgos)

    spreadprobunknownallalgos=[]
    spreadprobunknownallalgos.extend(spreadprobunknownoptmyalgos)   
    spreadprobunknownallalgos.extend(spreadprobunknownotheralgos)
    
    realdata=datatype[0]
    graphevolution=datatype[1]
    if realdata=="real" and graphevolution=="static":
       graphfolderprefix="../staticrealdata"
       tracefolderprefix="../tracegenerator"
    elif realdata=="real" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicrealdata"
       tracefolderprefix="../tracegenerator"
    elif realdata=="synthetic" and graphevolution=="graphwise":
       graphfolderprefix="../syndatagen"
       tracefolderprefix="../tracegenerator"
    elif realdata=="synthetic" and graphevolution=="paramwise":
       graphfolderprefix="../syndatagen"
       tracefolderprefix="../tracegenerator"
    elif realdata=="synthetic" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicsyndatagen"
       tracefolderprefix="../tracegenerator"

    if realdata=="synthetic":
       graphfolder="{0}_{1}_{2}".format(realdata,graphevolution,"graphs")
       graphfolder="{0}/{1}".format(graphfolderprefix,graphfolder)
    elif realdata=="real":
       graphfolder="{0}".format(graphfolderprefix)
       
    tracefolder="{0}_{1}_{2}_{3}_{4}".format(tracefolderstatic,realdata,graphevolution,spreadmodel,infertype)
    tracefolder="{0}/{1}".format(tracefolderprefix,tracefolder)
    dumpfolder="{0}_{1}_{2}_{3}_{4}".format(dumpfolderstatic,realdata,graphevolution,spreadmodel,infertype)
    resultfolder="{0}_{1}_{2}_{3}_{4}".format(resultfolderstatic,realdata,graphevolution,spreadmodel,infertype)
    runfolder="{0}_{1}_{2}_{3}_{4}".format(runfolderstatic,realdata,graphevolution,spreadmodel,infertype)
    
    if not os.path.exists(pbsfolder):
       os.makedirs(pbsfolder)
    if not os.path.exists(dumpfolder):
       os.makedirs(dumpfolder)
    if not os.path.exists(resultfolder):
       os.makedirs(resultfolder)
    if not os.path.exists(runfolder):
       os.makedirs(runfolder)

    #IMPORTANT: Dynamic graph inference can not be run sampleinterval other than 1!!!
    if graphevolution in ["dynamic"]:
       assert len(sampleintervals)==1 and sampleintervals[0]==1
       
    specificparamlist=specificparamreturner(spreadmodel,infertype,realdata,graphevolution)
    
    path2info={}
    if realdata=="real" and graphevolution=="static":
       filenames=myutil.listfiles(graphfolder)
       for filename in filenames:
           #if filename.find("146")==-1 and filename.find("168")==-1 and filename.find("212")==-1 and filename.find("248")==-1:
           #   continue
           if filename.find("stanford")!=-1 and filename.find("207")!=-1:
              pass
           elif filename.find("enron")!=-1 and filename.find("261")!=-1:
              continue
           else:
              continue 
           filepath="{0}/{1}".format(graphfolder,filename)
           print filepath
           if filename.find(".gml")!=-1:
              G=nx.read_gml(filepath,relabel=True)
           elif filename.find(".sparse6")!=-1:
              G=nx.Graph(nx.read_sparse6(filepath))
           elif filename.find(".graph6")!=-1:
              G=nx.Graph(nx.read_graph6(filepath))
           else:
              print "{0} extension is unknown!!".format(filename)
              exit(1)
           
           if type(G)==nx.Graph:
              graphtype="undirected"
           elif type(G)==nx.DiGraph:
              graphtype="directed"
           Gname=filename

           senttracefolder="{0}/{1}".format(tracefolder,Gname)
           assert os.path.exists(senttracefolder)
           path2info[filepath]=(G,Gname,graphtype,senttracefolder)
    elif realdata=="real" and graphevolution=="dynamic":
       graphdirnames=myutil.listdirectories(graphfolder)
       for graphdirname in graphdirnames:
           graphdirpath="{0}/{1}".format(graphfolder,graphdirname)
           Gall={}
           Gname=graphdirname
           for filename in myutil.listfiles(graphdirpath):
               filepath="{0}/{1}".format(graphdirpath,filename)
               if filename.find(".gml.gz")!=-1:
                  time=int(filename.replace(".gml.gz",""))
                  G=nx.read_gml(filepath,relabel=True) 
               elif filename.find(".gml")!=-1:
                  time=int(filename.replace(".gml",""))
                  G=nx.read_gml(filepath,relabel=True)
               elif filename.find(".sparse6")!=-1:
                  time=int(filename.replace(".sparse6",""))
                  G=nx.Graph(nx.read_sparse6(filepath))
               elif filename.find(".graph6")!=-1:
                  time=int(filename.replace(".graph6",""))
                  G=nx.Graph(nx.read_graph6(filepath))
               else:
                  print "{0} extension is unknown!!".format(filename)
                  exit(1)

               if type(G)==nx.Graph:
                  graphtype="undirected"
                  Gall[time]=nx.Graph(G)
               elif type(G)==nx.DiGraph:
                  graphtype="directed"
                  Gall[time]=nx.DiGraph(G)   
                
           senttracefolder="{0}/{1}".format(tracefolder,Gname)
           assert os.path.exists(senttracefolder)    
           path2info[graphdirpath]=(Gall,Gname,graphtype,senttracefolder)
    elif realdata=="synthetic" and graphevolution=="graphwise":
       datas=myutil.listdirectories(graphfolder)
       for data in datas:
           print data
           if data.find("250deg10")==-1:
              continue 
           datadirname="{0}/{1}".format(graphfolder,data)
           graphalgos=myutil.listdirectories(datadirname)
           for graphalgo in graphalgos:
               print graphalgo
               if graphalgo in ["kronecker","directed_ff"]:
                  graphtype="directed"
               else:
                  graphtype="undirected"
               graphdirname="{0}/{1}".format(datadirname,graphalgo)
               filenames=myutil.listfiles(graphdirname) 
               assert len(filenames)==filesamplecount
               for filename in filenames:
                   if (graphalgo in ["dmc"]) and filename.startswith("13"):
                      continue
                   elif (graphalgo in ["rds"]) and filename.startswith("13"):
                      pass
                   elif (graphalgo in ["lpa"]) and filename.startswith("13"):
                      pass
                   elif (graphalgo in ["undirected_ff"]) and filename.startswith("13"):
                      pass 
                   elif (graphalgo in ["smw"]) and filename.startswith("13"):
                      continue
                   else:
                      continue
                   filepath="{0}/{1}".format(graphdirname,filename)
                   if filename.find(".gml")!=-1:
                      G=nx.read_gml(filepath,relabel=True)
                   elif filename.find(".sparse6")!=-1:
                      G=nx.Graph(nx.read_sparse6(filepath))
                   elif filename.find(".graph6")!=-1:
                      G=nx.Graph(nx.read_graph6(filepath))
                   else:
                      print "{0} extension is unknown!!".format(filename)
                      exit(1)
                   Gname=filename

                   tracefolder2="{0}/{1}/{2}".format(tracefolder,data,graphalgo)
                   senttracefolder="{0}/{1}".format(tracefolder2,Gname)
                   assert os.path.exists(senttracefolder)
                   path2info[filepath]=(G,Gname,graphtype,senttracefolder)
    elif realdata=="synthetic" and graphevolution=="dynamic":
       datas=myutil.listdirectories(graphfolder)
       for data in datas:
           datadirname="{0}/{1}".format(graphfolder,data)
           nodenum,temp=data.replace("datan","").split("s")
           stepsize,tcount=temp.split("t")
           nodenum=int(nodenum)
           stepsize=int(stepsize)
           tcount=int(tcount)
           if nodenum==32 and stepsize==2 and tcount==20:
              pass
           else:
              continue
           datadirname="{0}/{1}".format(graphfolder,data)
           graphalgos=myutil.listdirectories(datadirname)
           for graphalgo in graphalgos:
               if graphalgo in ["kronecker","directed_ff"]:
                  graphtype="directed"
               else: #fused lasso is also undirected
                  graphtype="undirected"
               if graphalgo not in ["fusedlasso"]:
                  continue
               algodirname="{0}/{1}".format(datadirname,graphalgo)
               params=myutil.listdirectories(algodirname)
               for param in params:
                   if param.find("intervaldistuniform_intervaldistparam0.1_startgraphalgords_startgraphalgoparam0.1")!=-1:
                      pass
                   elif param.find("intervaldistuniform_intervaldistparam0.02_startgraphalgords_startgraphalgoparam0.1")!=-1:
                      pass
                   else:
                      continue 
                   paramdirname="{0}/{1}".format(algodirname,param)
                   sampledirnames=myutil.listdirectories(paramdirname)
                   print "info start"
                   print param
                   print len(sampledirnames)
                   assert len(sampledirnames)==filesamplecount #??
                   for sampledirname in sampledirnames:
                       if sampledirname.find("20")==-1:
                          continue 
                       samplepath="{0}/{1}".format(paramdirname,sampledirname)
                       Gall={}
                       filenames=myutil.listfiles(samplepath)
                       for filename in filenames:
                           filepath="{0}/{1}".format(samplepath,filename)
                           if filename.find(".gml.gz")!=-1:
                              time=int(filename.replace(".gml.gz",""))
                              G=nx.read_gml(filepath,relabel=True) 
                           elif filename.find(".gml")!=-1:
                              time=int(filename.replace(".gml",""))
                              G=nx.read_gml(filepath,relabel=True)
                           elif filename.find(".sparse6")!=-1:
                              time=int(filename.replace(".sparse6",""))
                              G=nx.Graph(nx.read_sparse6(filepath))
                           elif filename.find(".graph6")!=-1:
                              time=int(filename.replace(".graph6",""))
                              G=nx.Graph(nx.read_graph6(filepath))
                           else:
                              print "{0} extension is unknown!!".format(filename)
                              exit(1)
                            
                           if graphtype=="undirected":
                              assert type(G)==nx.Graph
                              Gall[time]=nx.Graph(G)
                           elif graphtype=="directed":
                              assert type(G)==nx.DiGraph
                              Gall[time]=nx.DiGraph(G)

                       Gname=sampledirname
                       tracefolder2="{0}/{1}/{2}/{3}".format(tracefolder,data,graphalgo,param)
                       senttracefolder="{0}/{1}".format(tracefolder2,Gname)
                       assert os.path.exists(senttracefolder)
                       path2info[samplepath]=(Gall,Gname,graphtype,senttracefolder)
                                

    algoplusroundmethodlist=[]
    sentroundingmethods={}
    if infertype=="difprobpartial":
       sentroundingmethods=["-1"]
       tempparam=["-1"]
       sentalgos=difprobpartialallalgos
       algoplusroundmethodlist=list(itertools.product(sentalgos,tempparam))
    elif infertype=="difprobunknown":
       sentroundingmethods=["-1"]
       tempparam=["-1"]
       sentalgos=difprobunknownallalgos
       algoplusroundmethodlist=list(itertools.product(sentalgos,tempparam))
    elif infertype=="spreadprobunknown":
       sentroundingmethods=["-1"]
       tempparam=["-1"]
       sentalgos=spreadprobunknownallalgos
       algoplusroundmethodlist=list(itertools.product(sentalgos,tempparam))
    elif infertype=="edge":
       if graphevolution in ["static","graphwise","paramwise"]:
          for algo in edgestaticallalgos:
              if algo in edgestaticoptmyalgos:
                 sentroundingmethods[algo]=roundingmethods
                 tempparam=["-1"]
                 algoplusroundmethodlist+=list(itertools.product([algo],tempparam))
              elif algo in edgestaticcovermyalgos:
                 sentroundingmethods[algo]=["-1"]
                 tempparam=["-1"]
                 algoplusroundmethodlist+=list(itertools.product([algo],tempparam))
              elif algo in edgestaticotheralgos:
                 if algo in ["netinf","multitree"]:
                    sentroundingmethods[algo]=["-1"]
                    tempparam=["-1"]
                    algoplusroundmethodlist+=list(itertools.product([algo],tempparam)) 
                 elif algo in ["netrate","connie"]:
                    sentroundingmethods[algo]=[("epsilon",0.999999999999)]
                    tempparam=["-1"]
                    algoplusroundmethodlist+=list(itertools.product([algo],tempparam))
              else:
                  print "algo {0} is unknown for dynamic case!!".format(algo)
                  exit(1)       
       elif graphevolution in ["dynamic"]:
          for algo in edgedynamicallalgos:
              if algo in edgedynamicoptmyalgos:
                  sentroundingmethods[algo]=roundingmethods
                  tempparam=["-1"]
                  algoplusroundmethodlist+=list(itertools.product([algo],tempparam))
              elif algo in edgedynamiccovermyalgos:
                  sentroundingmethods[algo]=["-1"]
                  tempparam=["-1"]
                  algoplusroundmethodlist+=list(itertools.product([algo],tempparam))
              else:
                  print "algo {0} is unknown for dynamic case!!".format(algo)
                  exit(1)
    paramlist=list(itertools.product(algoplusroundmethodlist,sampleselections,specificparamlist))
   
    for graphdirorfilepath in path2info.keys():
        G,Gname,graphtype,senttracefolder=path2info[graphdirorfilepath]
        for (algo,tempparam),sampleselection,specificparam in paramlist:
            partialstr="{0}-{1}".format(partialdegree,partialcorrelationfunction)
            specificfolder=tracegeneratorpart.returnspecificfolder(specificparam)
            
            #we could have shrinked the parametrization below but we didnot since in the future we might need different parametrization
            if algo in fixededgestaticoptmyalgos:
               typeerror=algo[0]
               typesparse=algo[1]
               if typesparse=="":
                  mylambda1s=[""]
                  mylambda2s=[""]
               elif typesparse=="l1":
                  mylambda1s=list(lambda1s)
                  mylambda2s=[""]
               elif typesparse=="l2":
                  mylambda1s=[""]
                  mylambda2s=list(lambda2s) 
               elif typesparse=="l1l2":
                  mylambda1s=list(lambda1s)
                  mylambda2s=list(lambda2s)
               parameterlist=list(itertools.product(mylambda1s,mylambda2s,[logbase],rightmethods,infectedprobabilityassignermethods))
               algostr="".join(algo)
            elif algo in fixededgestaticcsalgos:
               parameterlist=list(itertools.product(lambda1s,[logbase],rightmethods))
               algostr=algo
	    elif algo in fixededgestaticcovermyalgos:
               parameterlist=list(itertools.product(traceroundingmethods,[weightfunction]))
               algostr=algo
            elif algo in fixededgestaticotheralgos:
               if algo in ["netinf"]: 
                  parameterlist=list(itertools.product(traceroundingmethods,netinfdegreeratios))
               elif algo in ["multitree"]: 
                  parameterlist=list(itertools.product(traceroundingmethods,multitreedegreeratios))
               elif algo in ["netrate","connie"]:
                  parameterlist=list(itertools.product(traceroundingmethods,["-1"]))
               else:
                  print "other algo {0} is unknown for parameters".format(algo)
                  exit(1)
               algostr=algo
            elif algo in fixededgedynamicoptmyalgos:
               typeerror=algo[0]
               typesparse=algo[1]
               if typesparse=="":
                  mylambda1s=[""]
                  mylambda2s=[""]
                  myfusedlambdas=list(fusedlambdas)
               elif typesparse=="l1":
                  mylambda1s=list(lambda1s)
                  mylambda2s=[""]
                  myfusedlambdas=list(fusedlambdas)
               elif typesparse=="l2":
                  mylambda1s=[""]
                  mylambda2s=list(lambda2s)
                  myfusedlambdas=list(fusedlambdas)
               elif typesparse=="l1l2":
                  mylambda1s=list(lambda1s)
                  mylambda2s=list(lambda2s)
                  myfusedlambdas=list(fusedlambdas)
               parameterlist=list(itertools.product(mylambda1s,mylambda2s,myfusedlambdas,[logbase],rightmethods,infectedprobabilityassignermethods))
               algostr="".join(algo)
            elif algo in fixededgedynamiccsalgos:
               parameterlist=list(itertools.product(lambda1s,[logbase],rightmethods))
               algostr=algo
	    elif algo in fixededgedynamiccovermyalgos:
               parameterlist=list(itertools.product(traceroundingmethods,[weightfunction]))
               algostr=algo
            elif algo in fixededgedynamicotheralgos:
               parameterlist=list(itertools.product(traceroundingmethods))
               algostr=algo
            elif algo in fixeddifprobpartialoptmyalgos:
               parameterlist=list(itertools.product(lambda1s,lambda2s,[logbase],rightmethods,infectedprobabilityassignermethods))
               algostr="{0}{1}".format(algo[0],algo[1])
            elif algo in fixeddifprobpartialcsalgos:
               parameterlist=list(itertools.product(lambda1s,[logbase],rightmethods))
               algostr=algo
	    elif algo in fixeddifprobpartialcovermyalgos:
               parameterlist=list(itertools.product(traceroundingmethods,[weightfunction]))
               algostr=algo
            elif algo in fixeddifprobpartialotheralgos:
               parameterlist=list(itertools.product(traceroundingmethods))
               algostr=algo 
            elif algo in fixeddifprobunknownoptmyalgos:
               parameterlist=list(itertools.product(lambda1s,lambda2s,[logbase],rightmethods,infectedprobabilityassignermethods))
               algostr="{0}{1}".format(algo[0],algo[1])
            elif algo in fixeddifprobunknowncsalgos:
               parameterlist=list(itertools.product(lambda1s,[logbase],rightmethods))
               algostr=algo
	    elif algo in fixeddifprobunknowncovermyalgos:
               parameterlist=list(itertools.product(traceroundingmethods,[weightfunction]))
               algostr=algo
            elif algo in fixeddifprobunknownotheralgos:
               parameterlist=list(itertools.product(traceroundingmethods))
               algostr=algo      
            elif algo in fixedspreadprobunknownoptmyalgos:
               parameterlist=list(itertools.product(lambda1s,lambda2s,[logbase],rightmethods,infectedprobabilityassignermethods))
               algostr="{0}{1}".format(algo[0],algo[1])
            elif algo in fixedspreadprobunknowncsalgos:
               parameterlist=list(itertools.product(lambda1s,[logbase],rightmethods))
               algostr=algo
	    elif algo in fixedspreadprobunknowncovermyalgos:
               parameterlist=list(itertools.product(traceroundingmethods,[weightfunction]))
               algostr=algo
            elif algo in fixedspreadprobunknownotheralgos:
               parameterlist=list(itertools.product(traceroundingmethods))
               algostr=algo                      
            else:
               print "here error!! algo is unknown {0}".format(algo)
               exit(1)

            print "here"
            extensionparts=graphdirorfilepath.replace(graphfolder+"/","").split("/")[0:-1]
            extensionstr=""
            if len(extensionparts)!=0:
               extensionstr="-".join(extensionparts)
            print extensionstr
            extensionstr2="/".join(extensionstr.split("-"))
            
            for parameter in parameterlist:
                print "inside"
                parameterstrlist=[]
                for elem in parameter:
                    parameterstrlist.append(str(elem).replace(" ","").replace("'","").replace("(","").replace(")","").replace(",",""))
                parameterstr="?".join(parameterstrlist)
                print parameterstr
                print algo
                
                if algo in fixededgestaticoptmyalgos:
                   lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod=parameter
                   algoparams=[lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod]
                elif algo in fixededgestaticcsalgos:
                   lambda1,logbase,rightmethod=parameter
                   algoparams=[lambda1,logbase,rightmethod]
                elif algo in fixededgestaticcovermyalgos:
                   traceroundingmethod,weightfunction=parameter
                   algoparams=[traceroundingmethod,weightfunction]
                elif algo in fixededgestaticotheralgos:
                  if algo in ["netinf"]:
                     traceroundingmethod,netinfdegreeratio=parameter
                     algoparams=[traceroundingmethod,netinfdegreeratio]
                  elif algo in ["multitree"]:
                     traceroundingmethod,multitreedegreeratio=parameter
                     algoparams=[traceroundingmethod,multitreedegreeratio]
                  elif algo in ["netrate","connie"]:
                     traceroundingmethod,tempparam=parameter
                     assert tempparam=="-1"
                     algoparams=[traceroundingmethod]
                  else:
                     print "other algo {0} is unknown for parameters".format(algo)
                     exit(1)
                elif algo in fixededgedynamicoptmyalgos:
                   lambda1,lambda2,fusedlambda,logbase,rightmethod,infectedprobabilityassignermethod=parameter
                   algoparams=[lambda1,lambda2,fusedlambda,logbase,rightmethod,infectedprobabilityassignermethod]
                elif algo in fixededgedynamiccsalgos:
                   lambda1,logbase,rightmethod=parameter
                   algoparams=[lambda1,logbase,rightmethod]
                elif algo in fixededgedynamiccovermyalgos:
                   traceroundingmethod,weightfunction=parameter
                   algoparams=[traceroundingmethod,weightfunction]
                elif algo in fixededgedynamicotheralgos:
                   traceroundingmethod=parameter
                   algoparams=[traceroundingmethod]
                elif algo in fixeddifprobpartialoptmyalgos:
                   lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod=parameter
                   algoparams=[lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod]
                elif algo in fixeddifprobpartialcsalgos:
                   lambda1,logbase,rightmethod=parameter
                   algoparams=[lambda1,logbase,rightmethod]
                elif algo in fixeddifprobpartialcovermyalgos:
                   traceroundingmethod,weightfunction=parameter
                   algoparams=[traceroundingmethod,weightfunction]
                elif algo in fixeddifprobpartialotheralgos:
                   traceroundingmethod=parameter
                   algoparams=[traceroundingmethod]
                elif algo in fixeddifprobunknownoptmyalgos:
                   lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod=parameter
                   algoparams=[lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod]
                elif algo in fixeddifprobunknowncsalgos:
                   lambda1,logbase,rightmethod=parameter
                   algoparams=[lambda1,logbase,rightmethod]
                elif algo in fixeddifprobunknowncovermyalgos:
                   traceroundingmethod,weightfunction=parameter
                   algoparams=[traceroundingmethod,weightfunction]
                elif algo in fixeddifprobunknownotheralgos:
                   traceroundingmethod=parameter
                   algoparams=[traceroundingmethod]
                elif algo in fixedspreadprobunknownoptmyalgos:
                   lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod=parameter
                   algoparams=[lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod]
                elif algo in fixedspreadprobunknowncsalgos:
                   lambda1,logbase,rightmethod=parameter
                   algoparams=[lambda1,logbase,rightmethod]
                elif algo in fixedspreadprobunknowncovermyalgos:
                   traceroundingmethod,weightfunction=parameter
                   algoparams=[traceroundingmethod,weightfunction]
                elif algo in fixedspreadprobunknownotheralgos:
                   traceroundingmethod=parameter
                   algoparams=[traceroundingmethod]   
                else:
                   print "error!! algo is unknown {0}".format(algo)
                   exit(1)

                print "resultfolder is {0}".format(resultfolder)
                print extensionstr2
                if extensionstr2!="":
                   sentresultfolder="{0}/{1}/{2}/{3}/{4}".format(resultfolder,extensionstr2,Gname,specificfolder,algostr)
                else:
                   sentresultfolder="{0}/{1}/{2}/{3}".format(resultfolder,Gname,specificfolder,algostr) 
                if not os.path.exists(sentresultfolder):
                   os.makedirs(sentresultfolder)
                if extensionstr2!="":
                   sentrunfolder="{0}/{1}/{2}/{3}/{4}".format(runfolder,extensionstr2,Gname,specificfolder,algostr)
                else:
                   sentrunfolder="{0}/{1}/{2}/{3}".format(runfolder,Gname,specificfolder,algostr) 
                if not os.path.exists(sentrunfolder):
                   os.makedirs(sentrunfolder)
                print "useful resultfolder is {0}".format(sentresultfolder)
                print "useful runfolder is {0}".format(sentrunfolder)
                print extensionstr
                print extensionstr2

                print "running for {0} at location {1}".format(Gname,graphdirorfilepath)
                print specificfolder
                dumpfile="offlinedump_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pkl".format(extensionstr,Gname,algostr,specificfolder,sampleselection,algoruncount,partialstr,parameterstr) #we used spreadmodel,syn and other parameters on the name of dumpfolder
                dumppath="{0}/{1}".format(dumpfolder,dumpfile)
                outfile=gzip.open(dumppath,"wb")
                cPickle.dump(G,outfile) #this will be hash for dynamic graphs  
                cPickle.dump(Gname,outfile)
                cPickle.dump(graphtype,outfile)
                
                cPickle.dump(infertype,outfile)       
                cPickle.dump(algo,outfile)
                cPickle.dump(algoparams,outfile)
                cPickle.dump(parameterstr,outfile)
                cPickle.dump(sentroundingmethods[algo],outfile)
                
                cPickle.dump(sampleselection,outfile)
                cPickle.dump(tracefractions,outfile)
                cPickle.dump(algoruncount,outfile)
                cPickle.dump(spreadmodel,outfile)
                cPickle.dump(realdata,outfile)
                cPickle.dump(graphevolution,outfile)
                         
                cPickle.dump(partialdegree,outfile)
                cPickle.dump(partialcorrelationfunction,outfile)
                cPickle.dump(partialstr,outfile)
                cPickle.dump(sentresultfolder,outfile)
                cPickle.dump(sentrunfolder,outfile)
                cPickle.dump(senttracefolder,outfile)

                cPickle.dump(fixededgestaticoptmyalgos,outfile)
                cPickle.dump(fixededgestaticcovermyalgos,outfile)
                cPickle.dump(fixededgestaticotheralgos,outfile)
                cPickle.dump(fixededgestaticcsalgos,outfile)
                cPickle.dump(fixededgedynamicoptmyalgos,outfile)
                cPickle.dump(fixededgedynamiccovermyalgos,outfile)
                cPickle.dump(fixededgedynamicotheralgos,outfile)
                cPickle.dump(fixededgedynamiccsalgos,outfile)

                cPickle.dump(fixeddifprobpartialoptmyalgos,outfile)
                cPickle.dump(fixeddifprobpartialcovermyalgos,outfile)
                cPickle.dump(fixeddifprobpartialotheralgos,outfile)
                cPickle.dump(fixeddifprobpartialcsalgos,outfile)
                cPickle.dump(fixeddifprobunknownoptmyalgos,outfile)
                cPickle.dump(fixeddifprobunknowncovermyalgos,outfile)
                cPickle.dump(fixeddifprobunknownotheralgos,outfile)
                cPickle.dump(fixeddifprobunknowncsalgos,outfile)
                cPickle.dump(fixedspreadprobunknownoptmyalgos,outfile)
                cPickle.dump(fixedspreadprobunknowncovermyalgos,outfile)
                cPickle.dump(fixedspreadprobunknownotheralgos,outfile)
                cPickle.dump(fixedspreadprobunknowncsalgos,outfile)

                cPickle.dump(specificparam,outfile)
                cPickle.dump(specificfolder,outfile)
                cPickle.dump(sampleintervals,outfile)
                cPickle.dump(missingcompletemode,outfile)
                cPickle.dump(missingcompleteuptocount,outfile)
                outfile.close()

                #print "submitting {0}".format(dumppath)    
                code="python offlinegraphinferencepart.py {0}".format(dumppath)
                import sys
                sys.exit(1)
                
                #code="runCmd --nowait -- " + code
                #os.system(code)
                #continue
                
                pbsfilename="{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}-{10}-{11}".format(realdata,graphevolution,infertype,spreadmodel,extensionstr,Gname,algostr,specificfolder,sampleselection,algoruncount,partialstr,parameterstr)
                pbsrunner(code,pbsfolder,pbsfilename) 

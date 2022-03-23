#We will need another trace generator when we want to see the performance of the algorithms for given fixed spreading parameters!!!
#Questioms: Which distributions are nice to model exposed duration? Incubation period for diseases is in general right skewed distribution.
#1-Alternative methods of estimating an incubation distribution: examples from severe acute respiratory syndrome
#2-The distribution of incubation periods of infectious disease.
#3-Estimating incubation period distributions with coarse data #no bibtex found
#4-Incubation periods of acute respiratory viral infections: a systematic review
#5-A comparative epidemiologic analysis of SARS in Hong Kong, Beijing and Taiwan
#6-Generation interval contraction and epidemic data analysis
#7-The estimation of SARS incubation distribution from serial interval data using a convolution likelihood.
#8-Modeling the Incubation Period of Inhalational Anthrax #no bibtex found
#All of the above papers discuss log normal behaviour of incubation period
#Other than that, we can use various other distributions to model other type of spreading(such as idea spreading, etc).

#!!Graph recovery under fixed spreading parameter will also be run with the traces created here. 
#This generates sample traces for a given graph. This calls runner for each graph type with different parameter
import networkx as nx
import numpy as np
import scipy as sp
import random
import os
import sys
import math
import myutilities as myutil
import operator
import cPickle
import gzip
import itertools

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
    num=0
    if num==0:
       queuename="workstation"
    elif num==1:
       queuename="throughput"
    elif num==2:
       queuename="workstation"
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


pbsfolder="pbsfolder"
tracefolderstatic="traces"
dumpfolderstatic="tracedumps"
tracestartnodenum=1  #number of starting nodes at each trace. Multiple starting nodes don't always mean that trace have started from multiple nodes at the same time. It might also mean we're only observing some proportion of trace after it has started(We won't be able to start following the trace after it has started) 
#spreadprobunknown and difprobpartial will use static coefs below
staticsitracecoef=200.0
staticsistracecoef=100.0 #we also need observation interval since otherwise it might loop forever(will be modeled by maxspreadtime)
staticsirtracecoef=100.0
staticseirtracecoef=100.0
dynamicsitracecoef=200.0
dynamicsistracecoef=200.0
dynamicsirtracecoef=200.0
dynamicseirtracecoef=200.0

#spreading time will be modeled by spread2infectedparams and whether an infected node eventually will infect some other susceptible node is modeled by spreadprob.
spreadproblist=[0.05,0.1,0.2,0.3,0.4,0.5]
#("lognormal",1,0.25),("lognormal",2,0.25)
#weibull distribution is (scale,shape) and gauss is (mean,variance). exponential=geometric distribution, lognormal(mu,sigma), rayleigh()
#We expect average mean time to be somewhat around 2 and 4.
#weibull can be used to model various shapes
#("m",(("expo",2.0),("rayleigh",3.0)),(0.5,0.5)),("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))

#s2iparams=[("expo",1.0),("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("powerlaw",-2.5),("powerlaw",-3.0),("powerlaw",-4.0),("powerlaw",-5.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
s2iparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("weibull",3.0,2.5),("weibull",1.5,2.0),("weibull",1.5,3.0),("weibull",3.0,2.0),("weibull",5.0,1.5)]
#s2iparams=[("lognormal",2.0,),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0)]
#i2rparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("powerlaw",-3.0),("powerlaw",-4.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
i2rparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0)]
#[("expo",2.0),("expo",4.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
s2eparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("powerlaw",-3.0),("powerlaw",-4.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
#e2iparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("powerlaw",-3.0),("powerlaw",-4.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
e2iparams=[("expo",2.0),("powerlaw",-1.5)]
i2sparams=[("expo",2.0),("expo",3.0),("expo",4.0)] #sis only works for expo(geometric)!!
r2sparams=[("expo",2.0),("expo",3.0),("expo",4.0),("expo",6.0),("expo",8.0),("rayleigh",3.0),("powerlaw",-1.5),("powerlaw",-2.0),("powerlaw",-3.0),("powerlaw",-4.0),("weibull",3.0,2.5),("weibull",5.0,1.5)]
#s2imultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
i2rmultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
s2emultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
e2imultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
i2smultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
r2smultimodaldist=[("m",(("expo",2.0),("expo",4.0)),(0.5,0.5))]
#s2iparams.extend(s2imultimodaldist)
i2rparams.extend(i2rmultimodaldist)
s2eparams.extend(s2emultimodaldist)
e2iparams.extend(e2imultimodaldist)
i2sparams.extend(i2smultimodaldist)
r2sparams.extend(r2smultimodaldist)
#All of those distributions except the ones in multimodaldist are unimodal.
#Those distributions define variety of different characteristics such as modality,skewness, etc. Multimodal distributions correspond to more than one clusters of distributions which seems realistic. Different people groups might have different spreading characteritics(clusters) which may not be explained by a unimodal distribution. For example, some people are good responders(exposed duration is short) and some other people are bad responders(wait a lot until spreading to others). In this case, using a unimodal distribution and explanining the good and bad responders as left and right tails of a unimodal distribution with a single mean may not be the best way.


#This will effect only spread to infection transition probabilities not the rest!! 
#difprobpartial(netrateish case)
difprobpartialparams=[] #(dist,range)
unimodaldists=[("expo",(2.0,4.0)),("powerlaw",(-1.5,-2.5))]
#im not sure whether ("rayleigh",(3.0,5.0)) and ("weibull",(3.0,5.0),(1.5,2.5)) can be solved with this method
difprobpartialparams.extend(unimodaldists)
multimodaldists=[] #not using right now, it will be ok if we dont use(really tough)
difprobpartialparams.extend(multimodaldists)
difprobpartialselectionrule=["uniform"]

#difprobunknown
#random distribution type same for all edges
#this can be done by using edge inference information(no additional trace generation is required)

#spreadprobunknown, must guarantee random variable will be between 0 and 1. right now just use uniform distribution
spreadprobdists=[]
unimodaldists=[("uniform",1.0)]
spreadprobdists.extend(unimodaldists)
multimodaldists=[]
spreadprobdists.extend(multimodaldists)
tracetype="edge" #"edge","difprobpartial","difprobunknown","spreadprobunknown"  #difprobpartial,difprobunknown,spreadprobunknown will only run for static graphs
#do not thing to run tracegenerator for difprobunknown case!!!!!
#additional case will be :difprobunknown where each edge has different unknown distributions(very difficult)!! NOT IMPLEMENTED

maxspreadtime=1000000 #maximum number of temporal observation steps
samplecount=20 #samplecount for synthetic data
sisboundedspreadcoef=0.4
model1boundedspreadcoef=0.5
model2boundedspreadcoef=0.5
model5boundedspreadcoef=0.5
datatype=("real","dynamic")
#datatype=[("real","dynamic"),("real","static"),("synthetic","graphwise"),("synthetic,""paramwise",),("synthetic","dynamic")]
spreadmodel="si" #"seir", "si", "sir", "sis", "sirs", "model1", "model2", "model3", "model4", "model5"
missingcompletemode=True #just completes missing file traces if same param exists
#There will be 4 extra models
#In all of those cases, complicatednesss from graph inference perspective comes from systems edge dependence on more than one state.
#1- SIS-new(loopy): 
#2- SIR-new(not loopy): Infected nodes recover based on other recoveries (testing) (state changes may be loopy rather than linear like SIR or SEIR)
#3- SIRS-new: condition 1 + recover nodes can again becomes susceptible depending on neighbours(state changes may be loopy rather than linear like SIR or SEIR)
#4- SIRI-new(loopy): Recovered nodes has tendency to get infection again without being susceptible again
#5- SIS-new2(loopy): The model is SIS and nodes each time becomes infected with higher and higher probability(more they stayed infected, they might also decrease chance of being infected next time(bagisiklik))
# In all of those cases, probability distributions can be arbitrary and they may change over time. We can easily model all of those
#1,2,5 are qute IMPORTANT
#asyn can also be modeled in the probability distribution(first few parts will be 0).
#If we know about a delay duration in epidemic reaching other node, we can make sure first few parts of prob distribution is by assigning variables to 0
#spreading prob and spreading time dist is different
#graph is unweighted, since weight does not mean anything. Strongness of interaction or effect can be better modeled the diffusion parameter at that time step

def tracerunner(tracetype):
    if tracetype=="difprobunknown":
       print "this is difprobunknown, IM NOT RUNNING ANYTHING"
       print "This is wrong no need to run anything"
       exit(1)

    if not os.path.exists(pbsfolder):
       os.makedirs(pbsfolder)
     
    realdata=datatype[0]
    graphevolution=datatype[1]
    if tracetype!="edge" and graphevolution in ["dynamic"]:
       print "Tracetypes other than edge can not be run for dynamic graphs!!"
       print "Error"
       print "IMPORTANT WARNING: In the future, we might implement dynamic spreadprobunknown case where each edge spreadprob changes over time"
       print "This can either due to some given probability distribution or random"
       exit(1)
       
    if realdata=="real" and graphevolution=="static":
       graphfolderprefix="../staticrealdata"
    elif realdata=="real" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicrealdata"
    elif realdata=="synthetic" and graphevolution=="graphwise":
       graphfolderprefix="../syndatagen"
    elif realdata=="synthetic" and graphevolution=="paramwise":
       graphfolderprefix="../syndatagen"
    elif realdata=="synthetic" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicsyndatagen"   
    
    if realdata=="synthetic":
       graphfolder="{0}_{1}_{2}".format(realdata,graphevolution,"graphs")
       graphfolder="{0}/{1}".format(graphfolderprefix,graphfolder)
    elif realdata=="real":
       graphfolder="{0}".format(graphfolderprefix)
    tracefolder="{0}_{1}_{2}_{3}_{4}".format(tracefolderstatic,realdata,graphevolution,spreadmodel,tracetype)
    dumpfolder="{0}_{1}_{2}_{3}_{4}".format(dumpfolderstatic,realdata,graphevolution,spreadmodel,tracetype)

    path2info={}
    if realdata=="real" and graphevolution=="static":
       if not os.path.exists(dumpfolder):
          os.makedirs(dumpfolder)
       if not os.path.exists(tracefolder):
          os.makedirs(tracefolder)   
       filenames=myutil.listfiles(graphfolder)
       for filename in filenames:
           #if filename.find("enron")!=-1 and filename.find("261")!=-1:
           #   pass
           #else:
           #   continue 
           print filename 
           #nodenum=int(filename.split("-")[3])
           #if nodenum<250:
           #   continue 
           #if filename.find("stanford")==-1:
           #   continue 
           filepath="{0}/{1}".format(graphfolder,filename)
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
           nodenum=G.number_of_nodes()
           Gname=filename
           path2info[filepath]=(G,Gname,graphtype,nodenum)
    elif realdata=="real" and graphevolution=="dynamic":
       if not os.path.exists(dumpfolder):
          os.makedirs(dumpfolder)
       if not os.path.exists(tracefolder):
          os.makedirs(tracefolder)   
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
           firsttime=sorted(Gall.keys())[0]
           nodenum=Gall[firsttime].number_of_nodes()       
           path2info[graphdirpath]=(Gall,Gname,graphtype,nodenum)
    elif realdata=="synthetic" and graphevolution=="graphwise":
       if not os.path.exists(dumpfolder):
          os.makedirs(dumpfolder)
       datas=myutil.listdirectories(graphfolder)
       for data in datas:
           #if data!="datan125deg26" or data!="datan125deg10":
           #   continue 
           dirname="{0}/{1}".format(graphfolder,data)
           nodenum,degree=data.replace("datan","").split("deg")
           nodenum=int(nodenum)
           degree=int(degree)
           if nodenum==250 and degree==10:
              pass
           else: 
              continue
           algos=myutil.listdirectories(dirname)
           for algo in algos:
               if algo in ["kronecker","directed_ff"]:
                  graphtype="directed"
               else:
                  graphtype="undirected"
               dirname2="{0}/{1}".format(dirname,algo)
               filenames=myutil.listfiles(dirname2)
               assert len(filenames)==samplecount
               for filename in filenames:
                   if filename.find("13")!=-1:
                      print "running for {0}".format(filename)
                      pass
                   else:
                      continue
                   #if filename.find("13")!=-1 or filename.find("15")!=-1 or filename.find("12")!=-1:
                   #   print "running for {0}".format(filename)
                   #  pass
                   #else:
                   #   continue
                   filepath="{0}/{1}".format(dirname2,filename)
                   if filename.find(".gml")!=-1:
                      G=nx.read_gml(filepath,relabel=True)
                   elif filename.find(".sparse6")!=-1:
                      G=nx.Graph(nx.read_sparse6(filepath))
                   elif filename.find(".graph6")!=-1:
                      G=nx.Graph(nx.read_graph6(filepath))
                   else:
                      print "{0} extension is unknown!!".format(filename)
                      exit(1)
                   
                   if graphtype=="undirected":
                      assert type(G)==nx.Graph
                   elif graphtype=="directed":
                      assert type(G)==nx.DiGraph
                   
                   Gname=filename
                   nodenum=G.number_of_nodes()
                   path2info[filepath]=(G,Gname,graphtype,nodenum)
    elif realdata=="synthetic" and graphevolution=="dynamic":
       if not os.path.exists(dumpfolder):
          os.makedirs(dumpfolder)
       datas=myutil.listdirectories(graphfolder)
       for data in datas:
           dirname="{0}/{1}".format(graphfolder,data)
           nodenum,temp=data.replace("datan","").split("s")
           stepsize,tcount=temp.split("t")
           nodenum=int(nodenum)
           stepsize=int(stepsize)
           tcount=int(tcount)
           #if nodenum==128 and stepsize==2 and tcount==20:
           #   pass
           if nodenum==32 and stepsize==2 and tcount==20:
              pass
           elif nodenum==64 and stepsize==2 and tcount==20:
              pass
           else:
              continue 
           algos=myutil.listdirectories(dirname)
           for algo in algos:
               if algo in ["kronecker","directed_ff"]:
                  graphtype="directed"
               else: #this also involves fusedlasso!!!
                  graphtype="undirected"
               if algo not in ["fusedlasso"]:
                  continue 
               algodirname="{0}/{1}".format(dirname,algo)
               params=myutil.listdirectories(algodirname)
               for param in params:
                   #if param.find("intervaldistuniform_intervaldistparam0.02")==-1:
                   #   continue
                   paramdirname="{0}/{1}".format(algodirname,param)
                   sampledirnames=myutil.listdirectories(paramdirname)
                   print param
                   print len(sampledirnames)
                   assert len(sampledirnames)==samplecount
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
                       firsttime=sorted(Gall.keys())[0]
                       nodenum=Gall[firsttime].number_of_nodes()
                       path2info[samplepath]=(Gall,Gname,graphtype,nodenum)
    elif realdata=="synthetic" and graphevolution=="paramwise":
       if not os.path.exists(dumpfolder):
          os.makedirs(dumpfolder)
       datas=myutil.listdirectories(graphfolder)
       for data in datas:
           dirname="{0}/{1}".format(graphfolder,data)
           nodenum=int(data.replace("datan",""))
           if nodenum!=500:
              continue 
           algos=myutil.listdirectories(dirname)
           for algo in algos:
               if algo in ["kronecker","directed_ff"]:
                  graphtype="directed"
               else:
                  graphtype="undirected"
               algodirname="{0}/{1}".format(dirname,algo)
               params=myutil.listdirectories(algodirname)
               for param in params:
                   paramdirname="{0}/{1}".format(algodirname,param)
                   filenames=myutil.listfiles(paramdirname)
                   assert len(filenames)==samplecount
                   for filename in filenames:
                       filepath="{0}/{1}".format(paramdirname,filename)
                       if filename.find(".gml")!=-1:
                          G=nx.read_gml(filepath,relabel=True)
                       elif filename.find(".sparse6")!=-1:
                          G=nx.Graph(nx.read_sparse6(filepath))
                       elif filename.find(".graph6")!=-1:
                          G=nx.Graph(nx.read_graph6(filepath))
                       else:
                          print "{0} extension is unknown!!".format(filename)
                          exit(1)
                       
                       if graphtype=="undirected":
                          assert type(G)==nx.Graph
                       elif graphtype=="directed":
                          assert type(G)==nx.DiGraph
                   
                       Gname=filename
                       nodenum=G.number_of_nodes()
                       path2info[filepath]=(G,Gname,graphtype,nodenum)
    else: 
       print "error for UNKNOWN!!"
       print realdata
       print graphevolution
       exit(1)
       
    if tracetype=="edge":
       spreadprobparams=spreadproblist
       sents2iparams=s2iparams
    elif tracetype=="difprobpartial":
       spreadprobparams=spreadproblist 
       sents2iparams=list(itertools.product(difprobpartialparams,difprobpartialselectionrule))
    elif tracetype=="spreadprobunknown":
       spreadprobparams=spreadprobdists
       sents2iparams=s2iparams
    else:
       print "tracetype {0} is unknonn!!".format(tracetype)
       exit(1)
    
    if spreadmodel=="si":
       spreadmodelparams=[spreadprobparams,sents2iparams]
    elif spreadmodel=="sis":
       spreadmodelparams=[spreadprobparams,sents2iparams,i2sparams]
    elif spreadmodel=="sir":
       spreadmodelparams=[spreadprobparams,sents2iparams,i2rparams]
    elif spreadmodel=="seir":
       spreadmodelparams=[spreadprobparams,sents2iparams,e2iparams,i2rparams]
 
    print spreadprobparams
    print sents2iparams
    print "s"
    print "starting trace generation"
    for graphdirorfilepath in path2info.keys():
        print graphfolder
        print graphdirorfilepath
        G,Gname,graphtype,nodenum=path2info[graphdirorfilepath] #G is hash for dynamic graphs
        maxtracecountperparam=nodenum*20 #max trace count per param!!
        if graphevolution in ["static","graphwise","paramwise"]:
           if spreadmodel=="si":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*staticsitracecoef)
           elif spreadmodel=="sis":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*staticsistracecoef)
              globals()["maxspreadtime"]=int(nodenum*sisboundedspreadcoef)
           elif spreadmodel=="sir":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*staticsirtracecoef)
           elif spreadmodel=="seir":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*staticseirtracecoef)
        elif graphevolution=="dynamic":
           if spreadmodel=="si":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*dynamicsitracecoef)
           elif spreadmodel=="sis":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*dynamicsistracecoef)
              globals()["maxspreadtime"]=int(nodenum*sisboundedspreadcoef)
           elif spreadmodel=="sir":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*dynamicsirtracecoef)
           elif spreadmodel=="seir":
              tracesamplenum=int((float(nodenum)/tracestartnodenum)*dynamicseirtracecoef)
        
        extensionparts=graphdirorfilepath.replace(graphfolder+"/","").split("/")[0:-1]
        print extensionparts
        extensionstr=""
	if len(extensionparts)!=0:
           extensionstr="/".join(extensionparts)
           tracefolder2="{0}/{1}".format(tracefolder,extensionstr)
           if not os.path.exists(tracefolder2):
              os.makedirs(tracefolder2)
        else:
           tracefolder2=tracefolder 
        print tracefolder2      
        print extensionstr
        extensionstr2="_".join(extensionstr.split("/"))
        print extensionstr2
        senttracefolder2="{0}/{1}".format(tracefolder2,Gname)
        dumpfile="dump_{0}_{1}_{2}_{3}.pkl".format(extensionstr2,Gname,tracestartnodenum,tracesamplenum)
        dumppath="{0}/{1}".format(dumpfolder,dumpfile)
        outfile=gzip.open(dumppath,"wb")
        cPickle.dump(G,outfile)
        cPickle.dump(senttracefolder2,outfile)
        cPickle.dump(Gname,outfile)
        cPickle.dump(graphtype,outfile)
        cPickle.dump(tracesamplenum,outfile)
        cPickle.dump(tracestartnodenum,outfile)
        cPickle.dump(globals()["maxspreadtime"],outfile)
        cPickle.dump(tracetype,outfile)
        cPickle.dump(spreadmodel,outfile)
        cPickle.dump(spreadmodelparams,outfile)
        cPickle.dump(missingcompletemode,outfile)
        cPickle.dump(maxtracecountperparam,outfile)
        outfile.close()
         
        code="python tracegeneratorpart.py {0}".format(dumppath)
        ##code="runCmd -c \""+code+"\" --nowait"
        #code="runCmd --nowait -- " + code
        os.system(code)
        continue

        print "running pbs"
        pbsfilename="{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}".format(realdata,graphevolution,Gname,spreadmodel,extensionstr2,tracestartnodenum,tracesamplenum,tracetype)
        print code
        print pbsfolder
        print pbsfilename
        pbsrunner(code,pbsfolder,pbsfilename)
                                      
    
if __name__ == "__main__":
    tracerunner(tracetype)
    

#This generates sample traces for a given graph. This calls runner for each graph type with different parameter
import networkx as nx
import numpy as np
import scipy as sp
import scipy.stats
import random
import os
import sys
import math
import myutilities as myutil
import operator
import cPickle
import gzip
import itertools


def sis_statebased_spreader(startinfectednodes,allnodes,graphevolution,nonspreadpertimeedgehash,notspread2ipertime):
    print "running sis"
    nodestate={}
    nextnodestate={}
    infecttimes={}
    sustimes={}
    for node in allnodes:
        nodestate[node]=("s",1.0)
        nextnodestate[node]=""
        sustimes[node]=[0]
    for node in startinfectednodes:
        nodestate[node]=("i",1.0)
        infecttimes[node]=[0]
        del sustimes[node]
    infectednodes=set(startinfectednodes)
    #dynamic graphs are bounded by their temporal duration!!
    if graphevolution in ["static","paramwise","graphwise"]:
       maxtemporaltime=globals()["maxspreadtime"]
    elif graphevolution in ["dynamic"]:
       maxtemporaltime=max(Gall.keys())+1
    for time in range(0,maxtemporaltime): #allnodes are not guaranted to be infected, since we also have spreadprob.
        if len(infectednodes)==0:
           print "pre breaking at time {0}".format(time) 
           break
        for node in allnodes:
            if nodestate[node][0]=="s":
               psus=nodestate[node][1] #probability of staying at sus states 
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
                      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:    
                      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="i":
                      timedif=time-infecttimes[spreadernode][-1]+1
                      if timedif<=10: #for values greater than this, we assume psus have already converged to 1-spreaprob
                         psus *= nonspreadpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-psus):
                  assert node not in infectednodes
                  infectednodes.add(node)
                  nextnodestate[node]=("i",1.0)
                  if not infecttimes.has_key(node):
                     infecttimes[node]=[]
                  infecttimes[node].append(time+1)
               else:
                  nextnodestate[node]=("s",psus)
            elif nodestate[node][0]=="i":
               pinfect=nodestate[node][1] #probability of staying at infected states 
               timedif=time-infecttimes[node][-1]+1
               if timedif<=10: #for values greater than this,
                  pinfect *= notspread2ipertime[timedif]
               if random.random() <= (1.0-pinfect):
                  assert node in infectednodes
                  infectednodes.remove(node)
                  nextnodestate[node]=("s",1.0)
                  if not sustimes.has_key(node):
                     sustimes[node]=[]
                  sustimes[node].append(time+1)
               else:
                  nextnodestate[node]=("i",pinfect)
        nodestate=dict(nextnodestate)
    
    #error checking part
    print "error checking part"
    print "simulation length {0}".format(maxtemporaltime)
    avginfcount=0.0
    avgsuscount=0.0
    avginfectduration=0.0
    avgsusduration=0.0
    tempallnodes=set()
    tempallnodes|=set(infecttimes.keys())
    tempallnodes|=set(sustimes.keys())
    for node in tempallnodes:
        if infecttimes.has_key(node):
           avginfcount+=len(infecttimes[node])
        if sustimes.has_key(node):
           avgsuscount+=len(sustimes[node])
    avginfcount/=float(len(infecttimes.keys()))
    avgsuscount/=float(len(sustimes.keys()))
    print "avg number of times becoming infected is {0}".format(avginfcount)
    print "avg number of times becoming sus is {0}".format(avgsuscount)
    print "per node average count:"
    for node in infecttimes.keys():
        if not sustimes.has_key(node):
           assert len(infecttimes[node])==1 
           continue
        nodeinfecttimes=set(infecttimes[node])
        nodesustimes=set(sustimes[node])
        alltimes=set(nodeinfecttimes)
        alltimes |= set(nodesustimes)
        sortedalltimes=sorted(list(alltimes))
        lastflag=-1
        for time in sortedalltimes:
            if time in nodeinfecttimes:
               assert lastflag!=1 
               lastflag=1
            elif time in nodesustimes:
               assert lastflag!=2 
               lastflag=2
    
    #convert to nodestate vector,
    print "converting part!!"           
    nodestate={}
    alltimes=set()
    for node in infecttimes.keys():
        alltimes |= set(infecttimes[node])
    for node in sustimes.keys():
        alltimes |= set(sustimes[node])
    globalmaxtime=sorted(alltimes)[-1]    
    globalallnodes=set()
    globalallnodes|=set(infecttimes.keys())
    globalallnodes|=set(sustimes.keys())
    for node in globalallnodes:                             
        nodestate[node]={}
        for time in alltimes:
            nodestate[node][time]={}
            nodestate[node][time]["s"]=-1
            nodestate[node][time]["i"]=-1
        if not sustimes.has_key(node):
           assert len(infecttimes[node])==1
           for time in alltimes:
               nodestate[node][time]["i"]=1
               nodestate[node][time]["s"]=0
           continue
        if not infecttimes.has_key(node):
           assert len(sustimes[node])==1
           for time in alltimes:
               nodestate[node][time]["i"]=0
               nodestate[node][time]["s"]=1
           continue       
        nodeinfecttimes=set(infecttimes[node])
        nodesustimes=set(sustimes[node])
        nodealltimes=set(nodeinfecttimes)
        nodealltimes |= set(nodesustimes)
        sortednodealltimes=sorted(list(nodealltimes))
        for timeindex in range(1,len(sortednodealltimes)):
            curtime=sortednodealltimes[timeindex]
            prevtime=sortednodealltimes[timeindex-1]
            if curtime in nodeinfecttimes:
               for time in alltimes:
                   if time<curtime and time>=prevtime:
                      nodestate[node][time]["i"]=0
                      nodestate[node][time]["s"]=1
            elif curtime in nodesustimes:
               for time in alltimes:
                   if time<curtime and time>=prevtime:
                      nodestate[node][time]["i"]=1
                      nodestate[node][time]["s"]=0
        #do not forget to assign the last state
        nodelasttime=sortednodealltimes[-1]
        if nodelasttime in nodeinfecttimes:
           for time in alltimes:
               if time>=nodelasttime and time<=globalmaxtime:
                  nodestate[node][time]["i"]=1
                  nodestate[node][time]["s"]=0
        elif nodelasttime in nodesustimes:
           for time in alltimes:
               if time>=nodelasttime and time<=globalmaxtime:
                  nodestate[node][time]["i"]=0
                  nodestate[node][time]["s"]=1          
    #nodestates assignment error checking part!!       
    for node in nodestate.keys():
        for time in nodestate[node].keys():
            assert nodestate[node][time]["i"]!=-1 and nodestate[node][time]["s"]!=-1
    globals()["infectiontime"]={} #THAT IS VERY IMPORTANT!!!!!
    globals()["sustime"]={} #THAT IS VERY IMPORTANT!!!!! 
    for node in infecttimes.keys():
        globals()["infectiontime"][node]=list(infecttimes[node])
    for node in sustimes.keys():
        globals()["sustime"][node]=list(sustimes[node])     
    return


#nonspreadpertime is hash and starts from 1(not 0)!.
#nonspreadpertime is estimated from 1-f1s-f2s,1-f3s etc
def si_statebased_spreader(startinfectednodes,allnodes,graphevolution,nonspreadpertimeedgehash):
    nodestate={}
    nextnodestate={}
    infecttimes={}
    susnodes=set(allnodes)
    for node in allnodes:
        nodestate[node]=("s",1.0)
        nextnodestate[node]=""
    for node in startinfectednodes:
        nodestate[node]=("i",1.0)
        infecttimes[node]=0
        susnodes.remove(node)
    fixedcount=0
    lastcount=len(susnodes)+1
    #dynamic graphs are bounded by their temporal duration!!
    if graphevolution in ["static","paramwise","graphwise"]:
       maxtemporaltime=globals()["maxspreadtime"]
    elif graphevolution in ["dynamic"]:
       maxtemporaltime=max(Gall.keys())+1
    for time in range(0,maxtemporaltime): #allnodes are not guaranted to be infected, since we also have spreadprob.
        if len(susnodes)==lastcount:
           fixedcount+=1
        else:
           fixedcount=0 
        if fixedcount==15:
           #print "fixedcount reached!!"
           #print len(susnodes)
           #print time
           break 
        lastcount=len(susnodes) 
        for node in allnodes:
            if nodestate[node][0]=="s":
               psus=nodestate[node][1] 
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
		      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:
                      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="i":
                      timedif=time-infecttimes[spreadernode]+1
                      if timedif<=10: #for values greater than this, we assume psus have already converged to 1-spreaprob
                         psus *= nonspreadpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-psus):
                  susnodes.remove(node)
                  nextnodestate[node]=("i",1.0)
                  assert node not in infecttimes.keys()
                  infecttimes[node]=time+1
               else:
                  nextnodestate[node]=("s",psus)   
            elif nodestate[node][0]=="i":
               nextnodestate[node]=("i",1.0)  
        nodestate=dict(nextnodestate)
    #There can only be single startnode at time=0    
    startnodes=set()    
    for node in infecttimes.keys():
        if infecttimes[node]==0:
           startnodes.add(node)
    print startnodes
    assert len(startnodes)==len(startinfectednodes)
    globals()["infectiontime"]={} #THAT IS VERY IMPORTANT!!!!!
    for node in infecttimes.keys():
        globals()["infectiontime"][node]=infecttimes[node]       
    return

#nonspreadpertime and recoverpertime are hash and they start from 1(not 0)!.
#nonspreadpertime is estimated from 1-f1s-f2s,1-f1s etc
def sir_statebased_spreader(startinfectednodes,allnodes,graphevolution,nonspreadpertimeedgehash,recoverpertime):
    print "running sir"
    nodestate={}
    nextnodestate={}
    infecttimes={}
    recovertimes={}
    for node in allnodes:
        nodestate[node]=("s",1.0)
        nextnodestate[node]=""
    for node in startinfectednodes:
        nodestate[node]=("i",1.0)
        infecttimes[node]=0
    #dynamic graphs are bounded by their temporal duration!!
    if graphevolution in ["static","paramwise","graphwise"]:
       maxtemporaltime=globals()["maxspreadtime"]
    elif graphevolution in ["dynamic"]:
       maxtemporaltime=max(Gall.keys())+1
    for time in range(0,maxtemporaltime): #allnodes are not guaranted to be infected, since we also have spreadprob.    
        flag=False
        for node in allnodes:
            if nodestate[node][0]=="i":
               flag=True
               break
        if not flag:
           break 
        for node in allnodes:
            if nodestate[node][0]=="s":
               psus=nodestate[node][1] 
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
                      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:    
                      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="i":
                      timedif=time-infecttimes[spreadernode]+1
                      if timedif<=10: #for values greater than this, we assume psus have already converged to 1-spreaprob
                         psus *= nonspreadpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-psus):
                  nextnodestate[node]=("i",1.0)
                  assert node not in infecttimes.keys()
                  infecttimes[node]=time+1
               else:
                  nextnodestate[node]=("s",psus)   
            elif nodestate[node][0]=="i":
               pinfect=nodestate[node][1] 
               timedif=time-infecttimes[node]+1
               if timedif<=10: #for values greater than this, it wwill be converging to almost 0
                  pinfect *= recoverpertime[timedif]
               if random.random() <= (1.0-pinfect):
                  nextnodestate[node]=("r",1.0)
                  #notrecoveredset.remove(node)
                  assert node not in recovertimes.keys()
                  recovertimes[node]=time+1
               else:
                  nextnodestate[node]=("i",pinfect) 
            elif nodestate[node][0]=="r":    
               nextnodestate[node]=("r",1.0)
        nodestate=dict(nextnodestate)
    globals()["infectiontime"]={} #THAT IS VERY IMPORTANT!!!!!
    globals()["recovertime"]={} #THAT IS VERY IMPORTANT!!!!!
    for node in infecttimes.keys():
        globals()["infectiontime"][node]=infecttimes[node]
    for node in recovertimes.keys():
        globals()["recovertime"][node]=recovertimes[node]          
    return

#nonspreadpertime, exposedpertime and recoverpertime starts from 1(not 0)!
#nonspreadpertime is estimated from 1-f1s-f2s,1-f1s etc
def seir_statebased_spreader(startinfectednodes,allnodes,graphevolution,nonspreadpertimeedgehash,exposedpertime,recoverpertime):
    nodestate={}
    nextnodestate={}
    exposedtimes={}
    infecttimes={}
    recovertimes={}
    for node in allnodes:
        nodestate[node]=("s",1.0)
        nextnodestate[node]=""
    for node in startinfectednodes:
        nodestate[node]=("i",1.0)
        infecttimes[node]=0
    #dynamic graphs are bounded by their temporal duration!!
    if graphevolution in ["static","paramwise","graphwise"]:
       maxtemporaltime=globals()["maxspreadtime"]
    elif graphevolution in ["dynamic"]:
       maxtemporaltime=max(Gall.keys())+1
    for time in range(0,maxtemporaltime): #allnodes are not guaranted to be infected, since we also have spreadprob.    
        flag=False
        for node in allnodes:
            if nodestate[node][0]=="i" or nodestate[node][0]=="e":
               flag=True
               break
        if not flag:
           break 
        for node in allnodes:
            if nodestate[node][0]=="s":
               psus=nodestate[node][1] 
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
                      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:    
		      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="i":
                      timedif=time-infecttimes[spreadernode]+1
                      if timedif<=10: #for values greater than this, we assume psus have already converged to 1-spreaprob
                         psus *= nonspreadpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-psus):
                  nextnodestate[node]=("e",1.0)
                  assert node not in exposedtimes.keys()
                  exposedtimes[node]=time+1
               else:
                  nextnodestate[node]=("s",psus)
            elif nodestate[node][0]=="e":
               pexposed=nodestate[node][1] 
               timedif=time-exposedtimes[node]+1
               if timedif<=10: #for values greater than this,
                  pexposed *= exposedpertime[timedif]
               if random.random() <= (1.0-pexposed):
                  nextnodestate[node]=("i",1.0)
                  assert node not in infecttimes.keys()
                  infecttimes[node]=time+1
               else:
                  nextnodestate[node]=("e",pexposed)       
            elif nodestate[node][0]=="i":
               pinfect=nodestate[node][1] 
               timedif=time-infecttimes[node]+1
               if timedif<=10: #for values greater than this,
                  pinfect *= recoverpertime[timedif]
               if random.random() <= (1.0-pinfect):
                  nextnodestate[node]=("r",1.0)
                  assert node not in recovertimes.keys()
                  recovertimes[node]=time+1
               else:
                  nextnodestate[node]=("i",pinfect) 
            elif nodestate[node][0]=="r":    
               nextnodestate[node]=("r",1.0)
        nodestate=dict(nextnodestate)
    globals()["infectiontime"]={} #THAT IS VERY IMPORTANT!!!!!
    globals()["exposedtime"]={} #THAT IS VERY IMPORTANT!!!!!
    globals()["recovertime"]={} #THAT IS VERY IMPORTANT!!!!!
    for node in exposedtimes.keys():
        globals()["exposedtime"][node]=exposedtimes[node]    
    for node in infecttimes.keys():
        globals()["infectiontime"][node]=infecttimes[node]
    for node in recovertimes.keys():
        globals()["recovertime"][node]=recovertimes[node]          
    return


#There will be 3 extra models
#1- Infected nodes recover based on other recoveries (testing) (state changes may be loopy rather than linear like SIR or SEIR)
#2- condition 1 + recover nodes can again becomes susceptible depending on neighbours(state changes may be loopy rather than linear like SIR or SEIR)
#3- Recovered nodes has tendency to get infection again without being susceptible again
#4- The model is SIS and nodes each time becomes infected with higher and higher probability(more they stayed infected, they might also decrease chance of being infected next time(bagisiklik))
# In all of those cases, probability distributions can be arbitrary and they may change over time. We can easily model all of those
#1,4 is qute important
#asyn can also be modeled in the probability distribution.
#spreading prob and spreading time dist is different
#graph is unweighted, since weight does not mean anything. Strongness of interaction or effect can be better modeled the diffusion parameter at that time step

#1- SIS-new:
#2- SIR-new: Infected nodes recover based on other recoveries (testing) (state changes may be loopy rather than linear like SIR or SEIR)
#3- SIRS-new: condition 1 + recover nodes can again becomes susceptible depending on neighbours(state changes may be loopy rather than linear like SIR or SEIR)
#4- SIRI-new: Recovered nodes has tendency to get infection again without being susceptible again
#5- SIS-new2: The model is SIS and nodes each time becomes infected with higher and higher probability(more they stayed infected, they might also decrease chance of being infected next time(bagisiklik))


#In SIS-new, we may have problem of startnode can become susceptible without infecting any other nodes. We attack this problem by using more than one startnodes
#There is also probability r of staying at infected, 1-r and multiplication of affect of susceptible node edges we move back to susceptible state
#There will be another parameter mylambda=1-p. When lambda!=1, nodes will again pass to susceptible state even though there are no edges of suceptible nodes.
#Look at your formulation. When mylambda=1, this will reduce to only edge affect case.
#nonsusreturnpertimeedgehash -> probability of not returning back to the susceptible case.
#Number of startnodes will definitely affect the solution quality especially when mylambda=1 !!!!!
def model1_statebased_spreader(startinfectednodes,allnodes,graphevolution,nonspreadpertimeedgehash,nonsusreturnpertimeedgehash,mylambda):
    print "running model1(there is loop!!)"
    nodestate={}
    nextnodestate={}
    infecttimes={}
    sustimes={}
    for node in allnodes:
        nodestate[node]=("s",1.0)
        nextnodestate[node]=""
        sustimes[node]=[0]
    for node in startinfectednodes:
        nodestate[node]=("i",1.0)
        infecttimes[node]=[0]
        del sustimes[node]
    infectednodes=set(startinfectednodes)
    for time in range(0,globals()["maxspreadtime"]): #allnodes are not guaranted to be infected, since we also have spreadprob. 
        if len(infectednodes)==0:
           break
        for node in allnodes:
            if nodestate[node][0]=="s":
               psus=nodestate[node][1] #probability of staying at sus states 
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
                      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:    
                      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="i":
                      timedif=time-infecttimes[spreadernode][-1]+1
                      if timedif<=10: #for values greater than this, we assume psus have already converged to 1-spreaprob
                         psus *= nonspreadpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-psus):
                  assert node not in infectednodes
                  infectednodes.add(node)
                  nextnodestate[node]=("i",1.0)
                  if not infecttimes.has_key(node):
                     infecttimes[node]=[]
                  infecttimes[node].append(time+1)
               else:
                  nextnodestate[node]=("s",psus)
            elif nodestate[node][0]=="i":
               pinfect=nodestate[node][1] #probability of staying at infected states
               for spreadernode in nodestate.keys():
                   if graphevolution in ["static","paramwise","graphwise"]:
                      if not G.has_edge(spreadernode,node):
                         continue 
                   elif graphevolution in ["dynamic"]:    
                      if not Gall[time].has_edge(spreadernode,node):
                         continue
                   if nodestate[spreadernode][0]=="s":
                      timedif=time-sustimes[spreadernode][-1]+1
                      if timedif<=10: #for values greater than this,
                         pinfect *= nonsusreturnpertimeedgehash[(spreadernode,node)][timedif]
               if random.random() <= (1.0-(mylambda*pinfect)): #mylambda multiplication is here!!
                  assert node in infectednodes
                  infectednodes.remove(node)
                  nextnodestate[node]=("s",1.0)
                  if not sustimes.has_key(node):
                     sustimes[node]=[]
                  sustimes[node].append(time+1)
               else:
                  nextnodestate[node]=("i",pinfect)
        nodestate=dict(nextnodestate)

    #error checking part
    for node in infecttimes.keys():
        if not sustimes.has_key(node):
           assert len(infecttimes[node])==1 
           continue
        nodeinfecttimes=set(infecttimes[node])
        nodesustimes=set(sustimes[node])
        alltimes=set(nodeinfecttimes)
        alltimes |= set(nodesustimes)
        sortedalltimes=sorted(list(alltimes))
        lastflag=-1
        for time in sortedalltimes:
            if time in nodeinfecttimes:
               assert lastflag!=1 
               lastflag=1
            elif time in nodesustimes:
               assert lastflag!=2 
               lastflag=2
               
    for node in infecttimes.keys():
        globals()["infectiontime"][node]=list(infecttimes[node])
    for node in sustimes.keys():
        globals()["sustime"][node]=list(sustimes[node])     
    return





#NOT BEING USED RIGHT NOW        
#Minimum is 1
#In probabilistic diffusion simulation, this won't be used!!
def nodedurationassigner(dist,params,allnodes): #This will work for both static and dynamic graphs case
    durations={}
    if dist!="m":
       for node in allnodes:
           durations[node]=probassigner(dist,params)
    elif dist=="m": #multimodal
       mydists=params[1]
       sent=[]
       for distinfo in mydists:
           dist=distinfo[0]
           infoparams=tuple(distinfo[1:])
           sent.append((dist,infoparams))
       coefs=params[2]
       assert len(coefs)==len(sent)
       values=createmixturerandomvariables(sent,coefs,count=len(allnodes))
       index=0
       for node in allnodes:
           durations[node]=values[index]
           index += 1      
    return durations

#If minimum flag is 1, some of the distrubutions will be added 1 to make sure that returned value is at least 1
#generates disrete pdf from contunous ones by rounding
#NOT BEING USED RIGHT NOW SINCE WE DONT NEED IT(we are creating as a discrete way)
def generatepdf(dist,params,distlen):
    pdfvalues=[]
    if dist=="expo": #geometric
       mean=params[0]
       mylambda=1.0/mean
       for index in range(0,distlen):
           val=math.exp(index*mylambda)-math.exp(-1*(index+1)*mylambda)
           pdfvalues.append(val)
    elif dist=="powerlaw":
       a=params[0] #this will also determine coeff
       return int(round(scipy.stats.powerlaw.rvs(a,size=1)[0]))+add
    elif dist in ["zeta","zipf"]: #numpy zips generator works like zeta generator, returns discrete value(at least 1)
       a=params[0] #exponent
       return np.random.zipf(a,1)[0]
    elif dist=="pareto":
       a=params[0] #shape
       m=params[1] #location
       return int(round((np.random.pareto(a,1)+m)[0]))+add
    elif dist=="gauss":
       mean=params[0]
       stddev=params[1]
       gausstime=int(round(random.gauss(mean,stddev)))+add
       if gausstime<1:
          gausstime=1
       return gausstime
    elif dist=="weibull":
       scale=params[0]
       shape=params[1]
       return int(round(random.weibullvariate(scale,shape)))+add
    elif dist=="rayleigh":
       paramscale=params[0]
       return int(round(np.random.rayleigh(scale=paramscale,size=None)))+add
    elif dist=="lognormal":
       mu=params[0]
       sigma=params[1]
       return int(round(np.random.lognormal(mu,sigma,1)[0]))+add
    else:
        print "dist {0} is unknown!! error".format(dist)
    return pdfvalues

#If minimum flag is 1, some of the distrubutions will be added 1 to make sure that returned value is at least 1
#Then (0,1) -> interval will correspond to 1(shifting)
#powerlaw,zeta and pareto are all related!!
def probassigner(dist,params,minimumflag=1):
    add=0
    if minimumflag==1:
       add += 1 
    if dist=="expo":
       mean=params[0]
       return int(math.floor(random.expovariate(1.0/mean)))+add
    elif dist=="powerlaw":
       a=params[0] #this will also determine coeff
       return int(math.floor(scipy.stats.powerlaw.rvs(a,size=1)[0]))+add
    elif dist in ["zeta","zipf"]: #numpy zips generator works like zeta generator, returns discrete value(at least 1)
       a=params[0] #exponent
       return np.random.zipf(a,1)[0]
    elif dist=="pareto":
       a=params[0] #shape
       m=params[1] #location
       return int(math.floor((np.random.pareto(a,1)+m)[0]))+add
    elif dist=="gauss":
       mean=params[0]
       stddev=params[1]
       gausstime=int(math.floor(random.gauss(mean,stddev)))+add
       if gausstime<1:
          gausstime=1
       return gausstime
    elif dist=="weibull":
       scale=params[0]
       shape=params[1]
       return int(math.floor(random.weibullvariate(scale,shape)))+add
    elif dist=="rayleigh":
       paramscale=params[0]
       return int(math.floor(np.random.rayleigh(scale=paramscale,size=None)))+add
    elif dist=="uniform":
       val=params[0]
       return random.randrange(0,val)+add
    elif dist=="lognormal":
       mu=params[0]
       sigma=params[1]
       return int(math.floor(np.random.lognormal(mu,sigma,1)[0]))+add
    else:
        print "dist {0} is unknown!! error".format(dist)

def createmixturerandomvariables(dists,coefs,count=1):
    retlist=[]
    for varindex in range(0,count):
        cum=[]
        cum.append(coefs[0])
        for index in range(1,len(coefs)):
            tot=cum[index-1]+coefs[index]
            cum.append(tot)
        p=random.random()
        found=-1
        for index in range(0,len(cum)):
            if cum[index]>=p:
               found=index
               break
        assert found!=-1    
        dist=dists[found]
        var=probassigner(dist[0],dist[1])
        retlist.append(var)
    return retlist

def discrete_rayleigh(param,uplimit):
    probhash={}
    for index in range(1,uplimit+1):
        temp=(float(index)/(param**2))
        probhash[index]=temp*math.exp(-0.5*temp*index)
    mysum=0.0
    for index in range(1,uplimit+1):
        mysum += probhash[index]
    for index in range(1,uplimit+1):    
        probhash[index] = float(probhash[index])/mysum
    return probhash

def discrete_powerlaw(exponent,uplimit):
    probhash={}
    for index in range(1,uplimit+1):
        probhash[index]=-1.0*(exponent+1)*(index**exponent)
    mysum=0.0
    for index in range(1,uplimit+1):
        mysum += probhash[index]
    for index in range(1,uplimit+1):    
        probhash[index] = float(probhash[index])/mysum
    tempsum=0.0
    for val in probhash.values():
        tempsum += val
    #print "inside sum is {0}".format(tempsum)     
    return probhash

def discrete_exponential(mean,uplimit):
    mylambda=1.0/mean
    probhash={}
    for index in range(1,uplimit+1):
        probhash[index]=mylambda*math.exp(-1.0*mylambda*index)
    mysum=0.0
    for index in range(1,uplimit+1):
        mysum += probhash[index]
    for index in range(1,uplimit+1):    
        probhash[index] = float(probhash[index])/mysum
    tempsum=0.0
    for val in probhash.values():
        tempsum += val
    #print "inside sum is {0}".format(tempsum)     
    return probhash

def discrete_weibull(scale,shape,uplimit):
    probhash={}
    for index in range(1,uplimit+1):
        probhash[index]=(float(shape)/scale)*((float(index)/scale)**(shape-1))*(math.exp(-1.0*((float(index)/scale)**shape)))
    mysum=0.0
    for index in range(1,uplimit+1):
        mysum += probhash[index]
    for index in range(1,uplimit+1):    
        probhash[index] = float(probhash[index])/mysum
    tempsum=0.0
    for val in probhash.values():
        tempsum += val
    #print "inside sum is {0}".format(tempsum)    
    return probhash

def discrete_lognormal(mu,sigma,uplimit):
    probhash={}
    for index in range(1,uplimit+1):
        probhash[index]=(float(shape)/scale)*((float(index)/scale)**(shape-1))*(math.exp(-1.0*((float(index)/scale)**shape)))
    mysum=0.0
    for index in range(1,uplimit+1):
        mysum += probhash[index]
    for index in range(1,uplimit+1):    
        probhash[index] = float(probhash[index])/mysum
    tempsum=0.0
    for val in probhash.values():
        tempsum += val
    #print "inside sum is {0}".format(tempsum)    
    return probhash

def discrete_multimodal(mydists,coefs,uplimit):
    print mydists
    print coefs
    sent=[]
    for distinfo in mydists:
        dist=distinfo[0]
        infoparams=tuple(distinfo[1:])
        sent.append((dist,infoparams))
    assert len(coefs)==len(sent)
    tempstore=[]
    for index in range(0,len(sent)):
        tempdist,tempdistparam=sent[index]
        tempcoef=coefs[index]
        if tempdist=="rayleigh":
           tempprobhash=discrete_rayleigh(tempdistparam[0],uplimit)
        elif tempdist=="expo":
           tempprobhash=discrete_exponential(tempdistparam[0],uplimit) 
        elif tempdist=="powerlaw":
           tempprobhash=discrete_powerlaw(tempdistparam[0],uplimit)
        elif tempdist=="weibull":
           tempprobhash=discrete_weibull(tempdistparam[0],tempdistparam[1],uplimit)
        for key in tempprobhash.keys():
           tempprobhash[key]*=tempcoef
        tempstore.append(dict(tempprobhash))
    probhash={}
    for mykey in tempstore[0].keys():
        probhash[mykey]=tempstore[0][mykey]
    for elem in tempstore[1:]:
        for mykey in elem.keys():
            probhash[mykey]+=elem[mykey]
    tempsum=0.0
    for val in probhash.values():
        tempsum += val
    #print "this is multimodal"    
    #print "inside sum is {0}".format(tempsum)    
    return probhash
            
#All infected nodes must be recovered for models SIR and SEIR
#IMPORTANT: Since iteration number is bounded for dynamic graphs, some nodes may not be recovered!!!!!!
#This will be a problem whenever diffusion is cutoff before it is completed
def map2spreaddata(model):
    spreaddata={}
    if model=="si":       
       for node in globals()["infectiontime"].keys():
           infecttime=globals()["infectiontime"][node]
           spreaddata[node]={} 
           spreaddata[node]["infect"]=infecttime
    elif model=="sir":       
       for node in globals()["infectiontime"].keys():
           if not spreaddata.has_key(node):
              spreaddata[node]={}
           infecttime=globals()["infectiontime"][node]
           spreaddata[node]["infect"]=infecttime
       for node in globals()["recovertime"].keys():
           if not spreaddata.has_key(node):
              spreaddata[node]={}
           recovertime=globals()["recovertime"][node]
           spreaddata[node]["recover"]=recovertime          
    elif model=="seir":
       for node in globals()["exposedtime"].keys():
           exposedtime=globals()["exposedtime"][node]
           if not spreaddata.has_key(node):
              spreaddata[node]={}
           spreaddata[node]["exposed"]=exposedtime
       for node in globals()["infectiontime"].keys():
           infecttime=globals()["infectiontime"][node]
           if not spreaddata.has_key(node):
              spreaddata[node]={}
           spreaddata[node]["infect"]=infecttime
       for node in globals()["recovertime"].keys():
           if not spreaddata.has_key(node):
              spreaddata[node]={}
           recovertime=globals()["recovertime"][node]
           spreaddata[node]["recover"]=recovertime
    elif model=="sis": #infecttime and sustime will be list rather than single value        
       allnodes=set(globals()["infectiontime"].keys())
       allnodes|=set(globals()["sustime"].keys())
       for node in allnodes:
           spreaddata[node]={}
           if globals()["infectiontime"].has_key(node):
              infecttime=globals()["infectiontime"][node]
           else:
              infecttime=[]
           if globals()["sustime"].has_key(node):
              sustime=globals()["sustime"][node]
           else:
              sustime=[]
           spreaddata[node]["infect"]=infecttime
           spreaddata[node]["sus"]=sustime
    return spreaddata


def paramassigner(distparam):
    if distparam[0]=="m": #keep like this. will be handled during random variable assignment
       return ["m",distparam]
    else:    
       dist=distparam[0]
       if dist=="expo":
          mean=distparam[1]
          params=[mean]
       elif dist=="gauss": #? negative value handling, symettric this can be modeled better with certain webibull and rayleigh parameters
          mean=distparam[1]
          variance=distparam[2]
          stddev=math.sqrt(variance)
          params=[mean,stddev]
       elif dist=="weibull":
          scale=distparam[1]
          shape=distparam[2]
          params=[scale,shape]
       elif dist=="rayleigh":
          paramscale=distparam[1]
          params=[paramscale]
       elif dist=="lognormal":
          mu=distparam[1]
          sigma=distparam[2]
          params=[mu,sigma]
       elif dist=="uniform":
          val=distparam[1]
          params=[val]
       elif dist=="powerlaw":
          val=distparam[1]
          params=[val]      
       else:
          print "dist {0} is unknown!! error".format(dist)
    return [dist,params]


#mode can be normal, reverse cumulative, normal cumulative
#uplimitcount is 10 due to 1-this is enough and 2-precision errors!!
def partitiondist(dist,distparam,mode,spreadprob=-1.0):
    print dist
    print distparam
    if spreadprob==-1.0:
       spreadprob=1.0 
    uplimitcount=10
    if dist=="rayleigh":
       probhash=discrete_rayleigh(distparam[0],uplimitcount)
    elif dist=="expo":
       probhash=discrete_exponential(distparam[0],uplimitcount) 
    elif dist=="powerlaw":
       probhash=discrete_powerlaw(distparam[0],uplimitcount)
    elif dist=="weibull":
       probhash=discrete_weibull(distparam[0],distparam[1],uplimitcount)
    elif dist=="m":
       probhash=discrete_multimodal(distparam[1],distparam[2],uplimitcount)
    else:
       print "this dist is currently unknonwn {0}".format(dist)
       exit(1)

    #!!!THIS NORMALIZATION STEP IS MANDATORY.Even though we normalize inside probability function, it is not SAME and give error
    #Some multimodal distributions may also create very low values so righsize asser part below can give errror. BE CAREFUL for parameter of multimodal distribution   
    normsum=0
    for val in probhash.values():
        normsum += val
    for mykey in probhash.keys():
        probhash[mykey] /= float(normsum)
    maxtime=max(probhash.keys())
    for index in range(1,maxtime+1):
        probhash[index] *= spreadprob
    if mode=="normal":
       rightsize=dict(probhash) 
    elif mode=="reverse cumulative":
       rightsize={}
       mysum=0.0
       for index in range(1,maxtime+1):
           mysum += probhash[index]
           rightsize[index]=1.0-mysum
           assert rightsize[index]>=0
    elif mode=="normal cumulative":
       rightsize={}
       mysum=0.0
       for index in range(1,maxtime+1):
           mysum += probhash[index]
           rightsize[index]=mysum
    else:
       print "mode {0} is unknown".format(mode)
       exit(1)
    returnhash={}
    returnhash[1]=rightsize[1] 
    for index in range(2,uplimitcount+1):
        #print index
        #print rightsize[index]
        #print rightsize[index-1]
        ratio=float(rightsize[index])/rightsize[index-1] 
        returnhash[index]=ratio
    return returnhash

 
def nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,G,graphtype,graphevolution):
    nonspreadpertimeedgehash={}
    if tracetype=="edge": 
       s2idist,s2iparam=paramassigner(temps2iparam)
       onceforall=partitiondist(s2idist,s2iparam,"reverse cumulative",tempspreadprob)
       if graphevolution in ["static","graphwise","paramwise"]:
          for node1,node2 in G.edges():
              nonspreadpertimeedgehash[(node1,node2)]=dict(onceforall)
              if graphtype=="undirected":
                 nonspreadpertimeedgehash[(node2,node1)]=dict(onceforall)
       elif graphevolution in ["dynamic"]: #this is important
          for time in G.keys(): 
              for node1,node2 in G[time].edges():
                  nonspreadpertimeedgehash[(node1,node2)]=dict(onceforall)
                  if graphtype=="undirected":
                     nonspreadpertimeedgehash[(node2,node1)]=dict(onceforall)          
    elif tracetype=="difprobpartial": #only for static case
       for node1,node2 in G.edges():
           s2iparam,rule=temps2iparam[0]
           saveedgeresult=temps2iparam[1]
           distparamperedge=saveedgeresult[(node1,node2)]
           paramlist=[s2iparam[0]]
           paramlist.extend(distparamperedge)
           s2idist,s2iparam=paramassigner(paramlist)
           nonspreadpertime=partitiondist(s2idist,s2iparam,"reverse cumulative",tempspreadprob)
           nonspreadpertimeedgehash[(node1,node2)]=dict(nonspreadpertime)
           if graphtype=="undirected":
              nonspreadpertimeedgehash[(node2,node1)]=dict(nonspreadpertime)
    elif tracetype=="spreadprobunknown": #only for static case
       s2idist,s2iparam=paramassigner(temps2iparam)
       for node1,node2 in G.edges():
           saveedgeresult=tempspreadprob[1]
           edgespreadprob=saveedgeresult[(node1,node2)]
           nonspreadpertime=partitiondist(s2idist,s2iparam,"reverse cumulative",edgespreadprob)
           nonspreadpertimeedgehash[(node1,node2)]=dict(nonspreadpertime)
           if graphtype=="undirected":
              nonspreadpertimeedgehash[(node2,node1)]=dict(nonspreadpertime) 
    else:
       print "tracetype {0} is unknown!!".format(tracetype)
       exit(1)
    return nonspreadpertimeedgehash


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
                     for elem2 in elem:
                         if type(elem2)==tuple:
                            params.extend(elem2)
                         else:
                            params.append(elem2) 
                  else:   
                     params.append(elem)   
           tempparamlist=[]
           for param in params:
               tempparamlist.append(str(param))
           allstrs.append("_".join(tempparamlist))
       specificfolder="-".join(allstrs)
    return specificfolder

    
#Global variables
exposedtime={}
infectiontime={}
recovertime={}
sustime={} #this will especially be for looping states
maxspreadtime=-1
G=nx.Graph()
Gall={}
staticgraphevolution=["static","graphwise","paramwise"]
dynamicgraphevolution=["dynamic"]
missingcompletemode="-1"
maxtracecountperparam="-1"

if __name__ == "__main__":
    dumppath=sys.argv[1]
    firstfolder=dumppath.split("/")[0]
    realdata=firstfolder.split("_")[1]
    graphevolution=firstfolder.split("_")[2]
           
    infile=gzip.open(dumppath,"rb")
    if graphevolution=="dynamic":
       Gall=cPickle.load(infile)
    elif graphevolution in ["static","graphwise","paramwise"]:
       G=cPickle.load(infile)
    tracedir=cPickle.load(infile) #not global
    filename=cPickle.load(infile) #not global
    graphtype=cPickle.load(infile) #not global
    tracesamplenum=cPickle.load(infile) #not global
    tracestartnodenum=cPickle.load(infile) #not global
    maxspreadtime=cPickle.load(infile) 
    tracetype=cPickle.load(infile) #not global
    spreadmodel=cPickle.load(infile) #not global 
    spreadmodelparams=cPickle.load(infile) #not global
    missingcompletemode=cPickle.load(infile)
    maxtracecountperparam=cPickle.load(infile) 
    infile.close()

    allpossibleedges=set()    
    if graphevolution in ["static","graphwise","paramwise"]:
       allnodes=G.nodes()
    elif graphevolution=="dynamic":
       allnodes=set()
       for time in Gall.keys():
           G=Gall[time]
           allnodes|=set(G.nodes())
       allnodes=list(allnodes)   
    if graphtype=="undirected":
       for index1 in range(0,len(allnodes)):
           node1=allnodes[index1] 
           for index2 in range(index1+1,len(allnodes)):
               node2=allnodes[index2]
               allpossibleedges.add((node1,node2))  
    elif graphtype=="directed":   
       for node1 in allnodes:
           for node2 in allnodes:
               if node1!=node2:
                  allpossibleedges.add((node1,node2))

    if spreadmodel=="si":
       [spreadprobparams,sents2iparams]=spreadmodelparams
    elif spreadmodel=="sis":
       [spreadprobparams,sents2iparams,i2sparams]=spreadmodelparams
    elif spreadmodel=="sir":
       [spreadprobparams,sents2iparams,i2rparams]=spreadmodelparams
    elif spreadmodel=="seir":
       [spreadprobparams,sents2iparams,e2iparams,i2rparams]=spreadmodelparams

    uniquecount={}
    if tracetype=="difprobpartial":
       for distinfo,rule in sents2iparams:
           edgeparaminfo={}
           dist=distinfo[0]
           distparam=distinfo[1:]
           for node1,node2 in G.edges():
               if rule=="uniform":
                  varparam=[]
                  for tempparam in distparam:
                      lower=tempparam[0]
                      upper=tempparam[1]
                      vartemp=random.uniform(lower,upper)
                      varparam.append(vartemp)
               else:
                  print "rule {0} is unknonw!!".format(rule)
                  exit(1)
               edgeparaminfo[(node1,node2)]=varparam
           uniquecount[(distinfo,rule)]=dict(edgeparaminfo)
    elif tracetype=="spreadprobunknown":
       for spreadprobdist in spreadprobparams:
           #continous random variable not discrete(must be between 0 and 1)!! THIS IS IMPORTANT
           dist,params=paramassigner(spreadprobdist)
           edgeparaminfo={}
           for node1,node2 in G.edges():
               #we cant use probassigner method since it is for rounded discrete values. We need new method for continous values
               if dist=="uniform":
                  varvalue=random.uniform(0,1)    
               else:
                   print "currently dist {0} not implemented".format(dist)
                   print "we cant use probassigner method since it is for rounded discrete values. We need new method for continous values"
                   exit(1)
               #if dist!="m":
               #   varvalue=probassigner(dist,params)
               #else:
               #   mydists=params[1]
               #   sent=[]
               #   for distinfo in mydists:
               #       dist=distinfo[0]
               #       infoparams=tuple(distinfo[1:])
               #       sent.append((dist,infoparams))
               #   coefs=params[2]
               #   assert len(coefs)==len(sent) 
               #   varvalue=createmixturerandomvariables(sent,coefs,count=1)[0]
               edgeparaminfo[(node1,node2)]=varvalue 
           uniquecount[spreadprobdist]=dict(edgeparaminfo)

    if tracetype=="edge":
       runspreadprobparams=spreadprobparams
       runsents2iparams=sents2iparams
    elif tracetype=="difprobpartial":
       runspreadprobparams=spreadprobparams
       runsents2iparams=[]
       for dist,rule in uniquecount.keys():
           edgeelem=uniquecount[(dist,rule)]
           runsents2iparams.append(((dist,rule),edgeelem))
    elif tracetype=="spreadprobunknown":
       runspreadprobparams=[]
       for dist in uniquecount.keys():
           edgeelem=uniquecount[dist]
           runspreadprobparams.append(((dist),edgeelem))
       runsents2iparams=sents2iparams
       
    if spreadmodel=="si":
       runspreadmodelparams=[runspreadprobparams,runsents2iparams]
    elif spreadmodel=="sis":
       runspreadmodelparams=[runspreadprobparams,runsents2iparams,i2sparams]
    elif spreadmodel=="sir":
       runspreadmodelparams=[runspreadprobparams,runsents2iparams,i2rparams]
    elif spreadmodel=="seir":
       runspreadmodelparams=[runspreadprobparams,runsents2iparams,e2iparams,i2rparams]
    
    paramlist=list(itertools.product(*runspreadmodelparams)) 
    alltraces={}
    for samplenum in range(0,tracesamplenum):
        print "running trace for {0}".format(samplenum)
        if graphevolution in ["static","graphwise","paramwise"]:
           temp=G.nodes()
        elif graphevolution=="dynamic":
           firsttime=sorted(Gall.keys())[0]
           temp=Gall[firsttime].nodes()
        random.shuffle(temp)
        nodesamples=temp[0:tracestartnodenum]
        if graphevolution in ["static","graphwise","paramwise"]:
           allnodes=G.nodes()
        elif graphevolution=="dynamic":
           allnodes=set()
           for time in Gall.keys():
               allnodes=allnodes.union(Gall[time].nodes())

        #param part
        if spreadmodel=="si":
           tempspreadprob,temps2iparam=random.choice(paramlist)
           specificparam=[("spreadprob",tempspreadprob),("s2i",temps2iparam)]
           spreadfolder=returnspecificfolder(specificparam)
        elif spreadmodel=="sis":
           tempspreadprob,temps2iparam,i2sparam=random.choice(paramlist)
           specificparam=[("spreadprob",tempspreadprob),("s2i",temps2iparam),("i2s",i2sparam)]
           spreadfolder=returnspecificfolder(specificparam)
           i2sdist,i2sparam=paramassigner(i2sparam)
        elif spreadmodel=="sir":
           tempspreadprob,temps2iparam,i2rparam=random.choice(paramlist)
           specificparam=[("spreadprob",tempspreadprob),("s2i",temps2iparam),("i2r",i2rparam)]
           spreadfolder=returnspecificfolder(specificparam)
           i2rdist,i2rparam=paramassigner(i2rparam)
        elif spreadmodel=="seir":
           tempspreadprob,temps2iparam,e2iparam,i2rparam=random.choice(paramlist)
           specificparam=[("spreadprob",tempspreadprob),("s2i",temps2iparam),("e2i",e2iparam),("i2r",i2rparam)]
           spreadfolder=returnspecificfolder(specificparam)
           e2idist,e2iparam=paramassigner(e2iparam)
           i2rdist,i2rparam=paramassigner(i2rparam)
           
        lasttracefolder="{0}/{1}".format(tracedir,spreadfolder)
        if not os.path.exists(lasttracefolder):
           os.makedirs(lasttracefolder) 
        if missingcompletemode:
           filecount=0
           for fname in myutil.listfiles(lasttracefolder):
               filecount+=1
           if filecount>=maxtracecountperparam: #file count check
              continue
        tracefilename="-1"
        tracefilepath="-1"
        while True:
           fileindex=random.randint(0,10000000)
           print fileindex
           tracefilename="{0}_{1}.pkl".format(nodesamples[0],fileindex)
           tracefilepath="{0}/{1}".format(lasttracefolder,tracefilename)
           if not os.path.exists(tracefilepath):
              break    
        if spreadmodel=="si":
           if graphevolution in ["static","graphwise","paramwise"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,G,graphtype,graphevolution)
           elif graphevolution in ["dynamic"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,Gall,graphtype,graphevolution)     
           si_statebased_spreader(nodesamples,allnodes,graphevolution,nonspreadpertimeedgehash)
        elif spreadmodel=="sis":
           if graphevolution in ["static","graphwise","paramwise"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,G,graphtype,graphevolution)
           elif graphevolution in ["dynamic"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,Gall,graphtype,graphevolution)   
           notspread2ipertime=partitiondist(i2sdist,i2sparam,"reverse cumulative")
           print "running trace for {0} {1} {2}".format(temps2iparam,i2sdist,i2sparam)
           sis_statebased_spreader(nodesamples,allnodes,graphevolution,nonspreadpertimeedgehash,notspread2ipertime)
        elif spreadmodel=="sir":
           if graphevolution in ["static","graphwise","paramwise"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,G,graphtype,graphevolution)
           elif graphevolution in ["dynamic"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,Gall,graphtype,graphevolution)   
           recoverpertime=partitiondist(i2rdist,i2rparam,"reverse cumulative")
           sir_statebased_spreader(nodesamples,allnodes,graphevolution,nonspreadpertimeedgehash,recoverpertime)
        elif spreadmodel=="seir":
           if graphevolution in ["static","graphwise","paramwise"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,G,graphtype,graphevolution)
           elif graphevolution in ["dynamic"]:
              nonspreadpertimeedgehash=nonspreadpertimeedgereturner(tracetype,tempspreadprob,temps2iparam,Gall,graphtype,graphevolution)   
           exposedpertime=partitiondist(e2idist,e2iparam,"reverse cumulative")
           recoverpertime=partitiondist(i2rdist,i2rparam,"reverse cumulative")
           seir_statebased_spreader(nodesamples,allnodes,graphevolution,nonspreadpertimeedgehash,exposedpertime,recoverpertime)
           
        spreaddata=map2spreaddata(spreadmodel)
        assert tracefilepath!="-1"
        print tracefilepath
        file=gzip.open(tracefilepath,"wb")
        if tracetype=="edge":
           cPickle.dump(spreaddata,file)
        elif tracetype=="spreadprobunknown":
          cPickle.dump(alltraces,file)
          cPickle.dump(uniquecount,file)
        elif tracetype=="difprobpartial":
          cPickle.dump(uniquecount,file)
        file.close()
        continue
        
        if tracetype=="edge":
           keyspreadprob=tempspreadprob
           keys2iparam=temps2iparam
        elif tracetype=="spreadprobunknown":
           keyspreadprob=tempspreadprob[0]
           keys2iparam=temps2iparam
        elif tracetype=="difprobpartial":
           keyspreadprob=tempspreadprob
           keys2iparam=temps2iparam[0]
        #print "writing"   
        #print keyspreadprob
        #print keys2iparam
        #print uniquecount.keys()
        #print keys2iparam in uniquecount.keys()
           
        if spreadmodel=="si":
           tracekey=(tuple(sorted(nodesamples)),keyspreadprob,keys2iparam)
        elif spreadmodel=="sis":
           #convert other params to required format(tuple)
           if len(i2sparam)==1:
              i2stuple=(i2sdist,i2sparam[0])
           elif len(i2sparam)==2:
              i2stuple=(i2sdist,i2sparam[0],i2sparam[1])
           elif len(i2sparam)==3:
              i2stuple=(i2sparam[0],i2sparam[1],i2sparam[2])
           tracekey=(tuple(sorted(nodesamples)),keyspreadprob,keys2iparam,i2stuple)
        elif spreadmodel=="sir":
            #convert other params to required format(tuple) 
           if len(i2rparam)==1:
              i2rtuple=(i2rdist,i2rparam[0])
           elif len(i2rparam)==2:
              i2rtuple=(i2rdist,i2rparam[0],i2rparam[1])
           elif len(i2rparam)==3:
              i2rtuple=(i2rparam[0],i2rparam[1],i2rparam[2])
           tracekey=(tuple(sorted(nodesamples)),keyspreadprob,keys2iparam,i2rtuple)
        elif spreadmodel=="seir":
            #convert other params to required format(tuple) 
           if len(e2iparam)==1:
              e2ituple=(e2idist,e2iparam[0])
           elif len(e2iparam)==2:
              e2ituple=(e2idist,e2iparam[0],e2iparam[1])
           elif len(e2iparam)==3:
              e2ituple=(e2iparam[0],e2iparam[1],e2iparam[2])
           if len(i2rparam)==1:
              i2rtuple=(i2rdist,i2rparam[0])
           elif len(i2rparam)==2:
              i2rtuple=(i2rdist,i2rparam[0],i2rparam[1])
           elif len(i2rparam)==3:
              i2rtuple=(i2rparam[0],i2rparam[1],i2rparam[2])
           tracekey=(tuple(sorted(nodesamples)),keyspreadprob,keys2iparam,i2rtuple,e2ituple)
        if not alltraces.has_key(tracekey):
           alltraces[tracekey]=[]
        alltraces[tracekey].append(spreaddata)
               
    #output hash of all traces on the same file. The datatype will be of hash indexed by array of startnodes
    #file=gzip.open(tracefilepath,"wb")
    #if tracetype=="edge":
    #   cPickle.dump(alltraces,file)
    #elif tracetype=="spreadprobunknown":
    #   cPickle.dump(alltraces,file)
    #   cPickle.dump(uniquecount,file)
    #elif tracetype=="difprobpartial":
    #   cPickle.dump(uniquecount,file)
    #file.close()

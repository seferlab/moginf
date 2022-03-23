import networkx as nx
import numpy as np
import scipy as sp
import random
import os
import io
import sys
import math
import myutilities as myutil
import operator
import heapq
import inferenceoptimizationalgos as edgeinferenceoptalgos
import community
import gzip
import cPickle
import threading


#dynamic graphs this will be hash
#for dynamic case we assume each step might have different number of nodes!!!
def estimateedgepredictionscore(G,ourG,graphtype,graphevolution):
    print "comparison"
    #print G.number_of_edges()
    #print ourG.number_of_edges()
    #print G.edges()[0:5]
    #print ourG.edges()[0:5]
    #print set(ourG.edges()).intersection(set(G.edges()))
    if graphevolution in ["static","paramwise","graphwise"]:
       print "info pre"
       print G.number_of_nodes()
       print G.number_of_edges()
       print ourG.number_of_nodes()
       print ourG.number_of_edges()
       if graphtype=="directed":
          assert nx.DiGraph==type(G) and nx.DiGraph==type(ourG)
          tempG=nx.DiGraph(G)
          tempourG=nx.DiGraph(ourG)
          G={}
          ourG={}
          G[0]=nx.DiGraph(tempG)
          ourG[0]=nx.DiGraph(tempourG)
       elif graphtype=="undirected":
          assert nx.Graph==type(G) and nx.Graph==type(ourG)
          tempG=nx.Graph(G)
          tempourG=nx.Graph(ourG)
          G={}
          ourG={}
          G[0]=nx.Graph(tempG)
          ourG[0]=nx.Graph(tempourG)
          #sentG[0]=nx.Graph(G)
          #sentourG[0]=nx.Graph(ourG)
       #G=sentG
       #ourG=sentourG
    if graphevolution in ["static","paramwise","graphwise"]:
      print "after info"
      print G[0].number_of_nodes()
      print G[0].number_of_edges()
      print ourG[0].number_of_nodes()
      print ourG[0].number_of_edges()   
    
    #If nothing is predicted, make sure keys are assigned so there won't be any error!!
    for temptime in G.keys():
       if not ourG.has_key(temptime):
          if graphtype=="directed":
             ourG[temptime]=nx.DiGraph()
          elif graphtype=="undirected":
             ourG[temptime]=nx.Graph()
          else:
             print "this graphtype {0} is unknonw!!".format(graphtype)
             exit(1)
    assert len(set(G.keys()).difference(set(ourG.keys())))==0 and len(set(ourG.keys()).difference(G.keys()))==0
    alltimes=G.keys()
    #allnodes=set()
    #for time in G.keys():
    #    allnodes|=set(G[time].nodes())
    tp=0
    fp=0
    fn=0
    tp2=0
    tn=0  
    for time in alltimes:
       localtp=0
       localfp=0
       localfn=0
       localtp2=0
       localtn=0 
       for node1,node2 in G[time].edges():
           if ourG[time].has_edge(node1,node2):
              localtp+=1 
           else:
              localfn+=1 
       for node1,node2 in ourG[time].edges():
           if G[time].has_edge(node1,node2):
              localtp2+=1 
           else:
              localfp+=1
       assert localtp==localtp2
       #if graphtype=="directed":
       #   allcount=len(allnodes)*(len(allnodes)-1)
       #elif graphtype=="undirected":
       #   allcount=len(allnodes)*(len(allnodes)-1)/2
       if graphtype=="directed":
          allcount=G[time].number_of_nodes()*(G[time].number_of_nodes()-1)
       elif graphtype=="undirected":
          allcount=G[time].number_of_nodes()*(G[time].number_of_nodes()-1)/2
       else:
          print "unknown graph type {0}".format(graphtype)
          exit(1)
       localtn=allcount-localtp-localfn-localfp
       assert localtn>=0
       tp+=localtp
       fp+=localfp
       fn+=localfn
       tn+=localtn
       
    sen=float(tp)/(tp+fn)
    fpr=float(fp)/(fp+tn)
    recall=sen
    if (tp+fp)!=0:
       precision=float(tp)/(tp+fp)
    else:
       precision=0.0
    acc=float(tp+tn)/(tp+tn+fn+fp)
    spec=1.0-fpr #true negative rate

    #there is no true negative rate on f formula, but there is on mcc
    if precision+recall!=0.0:
       f1score=float(2.0*precision*recall)/(precision+recall)
       f2score=float(5.0*precision*recall)/((4.0*precision)+recall)
       avgnodenum=0.0
       for time in G.keys():
           avgnodenum+=G[time].number_of_nodes()
       avgnodenum/=float(len(G.keys()))    
       beta=1.0/float(avgnodenum)
       betasq=beta**2
       f1overnscore=float((1.0+betasq)*precision*recall)/((betasq*precision)+recall)
    else:
       f1score=0.0
       f2score=0.0
       f1overnscore=0.0

    if ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))!=0:
       mcc=float((tp*tn)-(fp*fn))/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))  
    else:
       mcc=0.0
       
    scores={}
    scores["sen"]=sen
    scores["fpr"]=fpr
    scores["recall"]=recall
    scores["precision"]=precision
    scores["acc"]=acc
    scores["spec"]=spec
    scores["f1"]=f1score
    scores["f2"]=f2score
    scores["f1overn"]=f1overnscore
    scores["mcc"]=mcc
    return scores

#this will be for difprobunknown
#both G and ourG are hash
#retvalues will be hash, keyed by edges and value will be probability vector
#NOT DONE YET
def estimatevectorscore(G,retvalues,graphtype):
    if graphtype=="directed":
       assert nx.DiGraph==type(G) and nx.DiGraph==type(ourG) 
    elif graphtype=="undirected":
       assert nx.Graph==type(G) and nx.Graph==type(ourG)     
    avgjsdiv=0.0
    #we might extend this similar to tf,fp,fn,tn
    
    alltimes=G.keys()
    alledges=set()
    for time in alltimes:
        alledges |= G[time].edges()
    for node1,node2 in alledges: 
        vec1=[]
        vec2=[]
        for time in alltimes:
            if not G.has_key(time):
               vec1.append(0.0)
            else:
               vec1.append(G[time][node1,node2])
            if not ourG.has_key(time):
               vec2.append(0.0)
            else:
               vec2.append(ourG[time][node1,node2])    
        val=js_divergence(vec1,vec2)
        avgjsdiv+=val
    avgjsdiv/=float(len(alledges))
    scores={}
    scores["jsdiv"]=avgjsdiv
    return scores
    
#this will be for difprobpartial and spreadprobunknown
#there are 6 types of errors
#scale dependent: mae,mse,rmse and scale independent: nmae,nmse,nrmse
def estimateerrorscore(G,ourG,graphtype):
    if graphtype=="directed":
       assert nx.DiGraph==type(G) and nx.DiGraph==type(ourG) 
    elif graphtype=="undirected":
       assert nx.Graph==type(G) and nx.Graph==type(ourG)
     
    sqe1=0.0 #squared error
    sqe2=0.0
    sqe3=0.0
    sqe4=0.0
    abse1=0.0 #abs error
    abse2=0.0
    abse3=0.0
    abse4=0.0
    nabse1=0.0 #normalized abs error
    nabse2=0.0
    nabse3=0.0
    nabse4=0.0 
    for node1,node2 in G.edges():
        if ourG.has_edge(node1,node2):
           sqe1+=(ourG[node1][node2]["weight"]-G[node1][node2]["weight"])**2
           abse1+=abs(ourG[node1][node2]["weight"]-G[node1][node2]["weight"])
           nabse1+=float(abs(ourG[node1][node2]["weight"]-G[node1][node2]["weight"]))/G[node1][node2]["weight"]
        else:
           sqe2+=(G[node1][node2]["weight"])**2
           abse2+=abs(G[node1][node2]["weight"])
           nabse2+=float(abs(G[node1][node2]["weight"]))/G[node1][node2]["weight"]
    for node1,node2 in ourG.edges():
        if G.has_edge(node1,node2):
           sqe3+=(ourG[node1][node2]["weight"]-G[node1][node2]["weight"])**2
           abse3+=abs(ourG[node1][node2]["weight"]-G[node1][node2]["weight"])
           nabse3+=float(abs(ourG[node1][node2]["weight"]-G[node1][node2]["weight"]))/G[node1][node2]["weight"]
        else:
           sqe4+=(ourG[node1][node2]["weight"])**2
           abse4+=abs(ourG[node1][node2]["weight"])
           nabse4+=float(abs(ourG[node1][node2]["weight"]))/G[node1][node2]["weight"]
    assert sqe1==sqe3
    assert abse1==abse3
    assert nabse1==nabse3
    if graphtype=="directed":
       allcount=G.number_of_nodes()*(G.number_of_nodes()-1)
    elif graphtype=="undirected":
       allcount=G.number_of_nodes()*(G.number_of_nodes()-1)/2
      
    totalsqe=sqe1+sqe2+sqe4
    totalabse=abse1+abse2+abse4
    totalnabse=nabse1+nabse2+nabse4
    
    mse=float(totalsqe)/allcount
    rmse=math.sqrt(mse)
    mabse=float(abse)/allcount
    nmabse=float(totalnabse)/allcount
    
    scores={}
    scores["mse"]=mse
    scores["rmse"]=rmse
    scores["mabse"]=mabse
    scores["nmse"]=-1 #not done
    scores["nrmse"]=-1 #not done
    scores["nmabse"]=nmabse
    return scores

def kl_divergence(dist1,dist2):
    mysum=0.0
    for index in range(0,len(dist1)):
        val=dist1[index]*math.log(float(dist1[index])/dist2[index],2)
        mysum += val
    return mysum
def js_divergence(probdist1,probdist2):
    #make both prob distributions of equal size
    len1=len(probdist1)
    len2=len(probdist2)
    maxlen=max(len1,len2)
    for index in range(0,maxlen):
        if index>=probdist1:
           probdist1.append(0.0)  
    for index in range(0,maxlen):
        if index>=probdist2:
           probdist2.append(0.0)
    half=[]
    for index in range(0,maxlen):
        val=(probdist1[index]+probdist2[index])/2.0
        half.append(val)
    return 0.5*kl_divergence(probdist1,half)+0.5*kl_divergence(probdist2,half)
def normalizedist(dist):
    mysum=0.0
    for deg in dist:
        mysum += deg
    normdist=[]
    for elem in dist:
        val=float(elem)/mysum
        normdist.append(val)
    return normdist
def estimatelouvainmodularity(G):
    part=community.best_partition(G)
    modvalue=community.modularity(part,G)
    return modvalue

#def estimatemodularity(G):
##    g=Graph.GRG(100, 0.2)
#    cl=g.community_fastgreedy()
#    print cl.membership

def assignproperties(G):
    graphproperties={}
    mygraphproperties["modularity"]=estimatemodularity(G)
    mygraphproperties["clusteringcoef"]=nx.average_clustering(G)
    mygraphproperties["degreedist"]=sorted(nx.degree(G).values(),reverse=True)  
    return graphproperties

def graphpropertiesconverge(G1,G2,Gerror):
    for stat in G1.keys():
        if stat=="modularity":
           dif=abs(G1[stat]-G2[stat])
           if dif<=Gerror[stat]:
              continue
           else:
              return False
        if stat=="clusteringcoef":
           dif=abs(G1[stat]-G2[stat])
           if dif<=Gerror[stat]:
              continue
           else:
              return False
        if stat=="degreedist":
           dif=math.sqrt(js_divergence(G1[stat],G2[stat]))
           if dif<=Gerror[stat]:
              continue
           else:
              return False
    return True     


#Discretize given scores in order to obtain the graph
#retvalues might also involve time information: (node1,node2) or (node1,node2,time)
#when algo returns directed, but we want to return undirected, if one of them is rounded, than that is ok.
#when algo return undirected but we want to return directed, we just copy edges(make sure as input to this function)
def score2graph(retvalues,method,methodparams,algo,graphtype,graphevolution):
    #we first normalize to (0,1) range if they are not
    rangevalid=True
    if graphevolution in ["static","paramwise","graphwise"]:
       for node1,node2 in retvalues.keys():
           if retvalues[(node1,node2)]>1:
              rangevalid=False
              break
    elif graphevolution in ["dynamic"]:
       for node1,node2,time in retvalues.keys():
           if retvalues[(node1,node2,time)]>1:
              rangevalid=False
              break
    if not rangevalid:
       if graphevolution in ["static","paramwise","graphwise"]:
          maxvalue=max(retvalues.values())       
          for node1,node2 in retvalues.keys():
              retvalues[(node1,node2)] /= float(maxvalue)
              if retvalues[(node1,node2)]<0: #we will round negative values to 0
                 retvalues[(node1,node2)]=0
       elif graphevolution in ["dynamic"]:
          maxvalue=max(retvalues.values())
          for node1,node2,time in retvalues.keys():
              retvalues[(node1,node2,time)] /= float(maxvalue)
              if retvalues[(node1,node2,time)]<0: #we will round negative values to 0
                 retvalues[(node1,node2,time)]=0
                 
    if graphevolution in ["static","paramwise","graphwise"]:
       if graphtype=="directed":
          retG=nx.DiGraph()
       elif graphtype=="undirected":
          retG=nx.Graph()
          
       if method=="graphproperties":
          graphpropertieserror=methodparams
          graphproperties=graphpropertieserror.keys()
          sortedkeys=[]
          retvalues2={}
          for node1,node2 in retvalues.keys():
              retvalues2[(node1,node2)]=min(retvalues[(node1,node2)],1.0-retvalues[(node1,node2)])
          for key,value in sorted(retvalues2.iteritems(),key=operator.itemgetter(1),reverse=True):
              sortedkeys.append(key)
          for node1,node2 in sortedkeys:
              if retvalues[(node1,node2)]>=0.5: 
                 retG.add_edge(node1,node2)
                 mygraphproperties=assignproperties(retG)
                 if graphpropertiesconverge(mygraphproperties,graphproperties,graphpropertieserror):
                    break
       elif method=="epsilon":
          roundingepsilon=methodparams
          for node1,node2 in retvalues.keys():
              if retvalues[(node1,node2)] >= 1.0-roundingepsilon:
                 retG.add_edge(node1,node2)
       else:
          print "rounding method {0} is unknown!!".format(method)
          exit(1)
       return retG   
    elif graphevolution in ["dynamic"]:
       alltimes=set()
       for node1,node,time in retvalues.keys():
           alltimes.add(time)
       retGall={}
       if graphtype=="directed":
          for time in alltimes:
              retG=nx.DiGraph()
              retGall[time]=retG
       elif graphtype=="undirected":
          for time in alltimes:
              retG=nx.Graph()
              retGall[time]=retG 
       if method=="graphproperties":
          graphpropertieserror=methodparams
          graphproperties=graphpropertieserror.keys()
          for curtime in alltimes:    
              sortedkeys=[]
              retvalues2={}
              for node1,node2,time in retvalues.keys():
                  if time!=curtime:
                    continue
                  retvalues2[(node1,node2)]=min(retvalues[(node1,node2,time)],1.0-retvalues[(node1,node2,time)])
              for key,value in sorted(retvalues2.iteritems(),key=operator.itemgetter(1),reverse=True):
                  sortedkeys.append(key)
              for node1,node2 in sortedkeys:
                  if retvalues[(node1,node2,curtime)]>=0.5: 
                     retGall[curtime].add_edge(node1,node2)
                     mygraphproperties=assignproperties(retGall[curtime])
                     if graphpropertiesconverge(mygraphproperties,graphproperties,graphpropertieserror):
                        break
       elif method=="epsilon":
          roundingepsilon=methodparams
          for curtime in alltimes: 
              for node1,node2,time in retvalues.keys():
                  if time!=curtime:
                     continue
                  if retvalues[(node1,node2,curtime)] >= 1.0-roundingepsilon:
                     retGall[curtime].add_edge(node1,node2)    
       else:
          print "rounding method {0} is unknown!!".format(method)
          exit(1)
       return retGall


def makedatapartial(nodestates,partialdegree,spreadmodel):
    rannode=random.choice(nodestates.keys())
    sortedalltimes=sorted(nodestates[rannode].keys())
    if spreadmodel=="si":
       for time in sortedalltimes:
           if random.random()<=partialdegree:
              continue 
           for node in nodestates.keys():
               if random.random()<=partialdegree:
                  continue
               if nodestates[node][time]["s"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["s"] -= corrupt
                  nodestates[node][time]["i"] += corrupt
               elif nodestates[node][time]["i"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["s"] += corrupt
                  nodestates[node][time]["i"] -= corrupt
    elif spreadmodel=="sir":
       for time in sortedalltimes:
           if random.random()<=partialdegree:
              continue 
           for node in nodestates.keys():
               if random.random()<=partialdegree:
                  continue
               if nodestates[node][time]["s"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["s"] -= corrupt
                  nodestates[node][time]["i"] += corrupt
               elif nodestates[node][time]["i"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  if random.random()<=0.5:
                     nodestates[node][time]["s"] += corrupt
                     nodestates[node][time]["i"] -= corrupt
                  else:
                     nodestates[node][time]["r"] += corrupt
                     nodestates[node][time]["i"] -= corrupt 
               elif nodestates[node][time]["r"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["r"] -= corrupt
                  nodestates[node][time]["i"] += corrupt                  
    elif spreadmodel=="seir":
       for time in sortedalltimes:
           if random.random()<=partialdegree:
              continue 
           for node in nodestates.keys():
               if random.random()<=partialdegree:
                  continue
               if nodestates[node][time]["s"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["s"] -= corrupt
                  nodestates[node][time]["e"] += corrupt
               elif nodestates[node][time]["e"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  if random.random()<=0.5:
                     nodestates[node][time]["s"] += corrupt
                     nodestates[node][time]["e"] -= corrupt
                  else:
                     nodestates[node][time]["i"] += corrupt
                     nodestates[node][time]["e"] -= corrupt    
               elif nodestates[node][time]["i"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  if random.random()<=0.5:
                     nodestates[node][time]["e"] += corrupt
                     nodestates[node][time]["i"] -= corrupt
                  else:
                     nodestates[node][time]["r"] += corrupt
                     nodestates[node][time]["i"] -= corrupt 
               elif nodestates[node][time]["r"]==1:
                  corrupt=random.uniform(0,1.0-partialdegree)
                  nodestates[node][time]["r"] -= corrupt
                  nodestates[node][time]["i"] += corrupt  
    return nodestates


#NOT BEING USED
#Noise will be added to both of the 3 states. After noise susceptible state noise will be kind of decreasing, infected state will be increasing and then decreasing, recovery state will be of the form decreasing.
#ts+ti/2 ->
#ts+ti+tr ->
#ti/2+tr ->       
#Noise means we can't exactly determine which time a node is infected or recovered. Therefore, we will assign infection state with probability. This will also affect recovery times.
def addnoise(nodestates):
    firstnode=nodestates.keys()[0]
    firsttimes=set(nodestates[firstnode].keys())
    for node in nodestates.keys()[1:]:
        temptimes=set(nodestates[node].keys())
        assert len(temptimes.difference(firsttimes))==0
    sortedalltimes=sorted(firsttimes)
    for timeindex in range(1,len(sortedalltimes)):
        assert sortedalltimes[timeindex]==sortedalltimes[timeindex-1]+1
    maxtime=max(sortedalltimes)
    for node in nodestates.keys():
        infecttime=-1
        for time in sortedalltimes:
            if nodestates[node][time]["i"]==1:
               infecttime=time
               break
        recovertime=-1
        for time in sortedalltimes:
            if nodestates[node][time]["r"]==1:
               recovertime=time
               break
        stime=infecttime
        itime=recovertime-infecttime
        rtime=maxtime-recovertime+1
        snoiselen=stime+((itime+1)/2)
        inoiselen=stime+itime+rtime
        rnoiselen=rtime+((itime+1)/2)
        snoisestart=0
        snoiseend=snoisestart+snoiselen-1
        inoisestart=0
        inoiseend=inoisestart+inoiselen-1
        inoisemiddle=infecttime+((itime+1)/2)
        inoisefirstlen=inoisemiddle-inoisestart
        inoisesecondlen=inoiseend+1-inoisemiddle
        rnoiseend=maxtime
        rnoisestart=rnoiseend-rnoiselen+1
        lowest=0.05
        svector=[lowest]*(maxtime+1) #0.05??
        m=(1.0-0.1)/((snoiselen)**2)
        for time in range(snoisestart,snoiseend+1):
            val=1.0-(m*(time-snoisestart+1)**2)
            if False:
             print val
             print time
             print snoisestart
             print m
             print snoiselen
            assert val>=0
            svector[time]=val
        ivector=[lowest]*(maxtime+1) #we will have piecewise function which has 2 parts
        m1=(1.0-0.1)/((inoisefirstlen)**2)
        m2=(1.0-0.1)/((inoisesecondlen)**2)
        for time in range(inoisestart,inoisemiddle):
            val=1.0-(m1*(inoisemiddle-time)**2)
            ivector[time]=val
            if False:
             print val
             print time
             print inoisestart
             print m1
             print inoisemiddle
            assert val>=0
        for time in range(inoisemiddle,inoiseend+1):
            val=1.0-(m2*(time-inoisemiddle+1)**2)
            ivector[time]=val
            if False:
             print val
             print time
             print inoiseend
             print m2
             print inoisemiddle
            assert val>=0
        rvector=[lowest]*(maxtime+1)
        m=(1.0-0.1)/((rnoiselen)**2)
        for time in range(rnoisestart,rnoiseend+1):
            val=1.0-(m*(rnoiseend-time+1)**2)
            rvector[time]=val
            if False:
             print val
             print time
             print rnoisestart
             print m
             print rnoiseend
            assert val>=0
        for time in range(0,len(svector)):
            assert svector[time]!=0.0
            assert ivector[time]!=0.0
            assert rvector[time]!=0.0
            assert svector[time]+ivector[time]!=0.0
            assert svector[time]+rvector[time]!=0.0
            assert rvector[time]+ivector[time]!=0.0
            tot=float(svector[time]+ivector[time]+rvector[time])
            assert tot >=0
            svector[time] /= tot
            ivector[time] /= tot
            rvector[time] /= tot
        for time in range(0,len(svector)):
            nodestates[node][time]["s"]=svector[time]    
        for time in range(0,len(ivector)):
            nodestates[node][time]["i"]=ivector[time]
        for time in range(0,len(rvector)):
            nodestates[node][time]["r"]=rvector[time]
        for time in range(0,len(svector)):
            nodetimesum=nodestates[node][time]["s"]+nodestates[node][time]["i"]+nodestates[node][time]["r"]
            #print "node: {0}, time {1} , val: {2}".format(node,time,valsum)
            assert round(nodetimesum,2)==1.0
    return nodestates


#NOT BEING USED
#This is uncertainity
#Do not confuse probablity of being infected at time t and probability of being at infected state at time t. They are different.
def smoothtransitions(nodestates,transitionlevel): #smooths transitions(square pulses)
    firstnode=nodestates.keys()[0]
    firsttimes=set(nodestates[firstnode].keys())
    for node in nodestates.keys()[1:]:
        temptimes=set(nodestates[node].keys())
        assert len(temptimes.difference(firsttimes))==0
    sortedalltimes=sorted(firsttimes)
    for timeindex in range(1,len(sortedalltimes)):
        assert sortedalltimes[timeindex]==sortedalltimes[timeindex-1]+1
    maxtime=max(sortedalltimes)
    for node in nodestates.keys():
        infecttime=-1
        for time in sortedalltimes:
            if nodestates[node][time]["i"]==1:
               infecttime=time
               break
        recovertime=-1
        for time in sortedalltimes:
            if nodestates[node][time]["r"]==1:
               recovertime=time
               break

        #transitionlevel=1.0 high, 0.0 low
        #our smooht transition is gaussian kind of transition
        #f=\lambdae^-(x-i)^2/q
        # We have 2 boundary conditions. One when x=t,f=a and x=infected, f=1
        # lambda=1 and q=-(t-i)^2/ln(a)    
        transition1duration=int(round(infecttime*transitionlevel))
        transition1start=infecttime-transition1duration
        transition1end=infecttime
        a=0.1
        for time in range(transition1start,transition1end):
            svalue=(-1.0*(time-infecttime)**2)/math.log(a)
            nodestates[node][time]["s"]=svalue
            nodestates[node][time]["i"]=1.0-svalue
        #    
        transition2duration=(recovertime-infecttime)/transitionlevel
        transition2start=infecttime-transition2duration
        transition2end=infecttime
        for time in range(transition2start,transition2end):
            dif=1.0
            nodestates[node][time]["s"]=1.0-dif
            nodestates[node][time]["i"]=dif
        for time in sortedalltimes:
            nodetimesum=nodestates[node][time]["s"]+nodestates[node][time]["i"]+nodestates[node][time]["r"]
            assert round(nodetimesum,2)==1.0
    return nodestates


def convert2seirstate(nodestates): #e and i will be combined into single i state
    newnodestates={}
    for node in nodestates.keys():
        newnodestates[node]={}
        for time in nodestates[node].keys():
            newnodestates[node][time]={}
            newnodestates[node][time]["s"]=nodestates[node][time]["s"]
            newnodestates[node][time]["i"]=nodestates[node][time]["e"]+nodestates[node][time]["i"]
            newnodestates[node][time]["r"]=nodestates[node][time]["r"]
    return newnodestates

#extract for given statename
def extractstatetimes(spreaddata,statename):
    assert statename in ["sus","infect","recover","exposed"]
    returntimes={}
    for node in spreaddata.keys():
        mytime=spreaddata[node][statename]
        returntimes[node]=mytime
    return returntimes
 

#seir will lose its exposed time information during conversion??
#BE CAREFUL each infected node do not have to be recovered in some cases(when we use limited horizon)
#THIS HAS BEEN MODIFIED!!
#SIS PART HAS NOT BEEN DONE YET!!!!!
def assignnodestates(spreaddata,spreadmodel):
    nodestates={}
    if spreadmodel in ["sir","seir","si"]:
       alltimes=set()
       for node in spreaddata.keys():
           for state in spreaddata[node].keys():
               alltimes.add(spreaddata[node][state])
       sortedalltimes=sorted(list(alltimes))
    elif spreadmodel in ["sis"]:
       alltimes=set()
       for node in spreaddata.keys():
           for state in spreaddata[node].keys():
               alltimes |= set(spreaddata[node][state])
       sortedalltimes=sorted(list(alltimes))
    else:
       print "this spreadmodel {0} is UNKNOWN!!".format(spreadmodel)
       exit(1)
       
    if spreadmodel=="sis": 
       globalmaxtime=sortedalltimes[-1]   
       for node in spreaddata.keys():
           nodestates[node]={}
           if spreaddata[node].has_key("sus"):
              sustimes=set(spreaddata[node]["sus"])
           else:
              sustimes=set()
           if spreaddata[node].has_key("infect"):
              infecttimes=set(spreaddata[node]["infect"])
           else:
              infecttimes=set()    
           for time in sortedalltimes:
               nodestates[node][time]={}
               nodestates[node][time]["s"]=-1
               nodestates[node][time]["i"]=-1
           if len(sustimes)==0:
               assert len(infecttimes)==1
               for time in sortedalltimes:
                   nodestates[node][time]["i"]=1
                   nodestates[node][time]["s"]=0
               continue
           if len(infecttimes)==0:
               assert len(sustimes)==1
               for time in sortedalltimes:
                   nodestates[node][time]["i"]=0
                   nodestates[node][time]["s"]=1
               continue
           for timeindex in range(1,len(sortedalltimes)):
               curtime=sortedalltimes[timeindex]
               prevtime=sortedalltimes[timeindex-1]
               if curtime in infecttimes:
                   for time in sortedalltimes:
                       if time<curtime and time>=prevtime:
                          nodestates[node][time]["i"]=0
                          nodestates[node][time]["s"]=1
               elif curtime in sustimes:
                   for time in sortedalltimes:
                      if time<curtime and time>=prevtime:
                         nodestates[node][time]["i"]=1
                         nodestates[node][time]["s"]=0
           #do not forget to assign the last state
           nodelasttime=sortedalltimes[-1]
           if nodelasttime in infecttimes:
              for time in sortedalltimes:
                  if time>=nodelasttime and time<=globalmaxtime:
                     nodestates[node][time]["i"]=1
                     nodestates[node][time]["s"]=0
           elif nodelasttime in sustimes:
              for time in sortedalltimes:
                  if time>=nodelasttime and time<=globalmaxtime:
                     nodestates[node][time]["i"]=0
                     nodestates[node][time]["s"]=1
       #nodestates assignment error checking part!!       
       for node in nodestates.keys():
           for time in nodestates[node].keys():
               print node
               print time
               print nodestates[node][time]["i"]
               print nodestates[node][time]["s"]
               assert nodestates[node][time]["i"]!=-1 and nodestates[node][time]["s"]!=-1
    elif spreadmodel=="si":
       for node in spreaddata.keys():
           nodestates[node]={}
           infecttime=spreaddata[node]["infect"]
           for time in sortedalltimes:
               nodestates[node][time]={}
               nodestates[node][time]["s"]=1
               nodestates[node][time]["i"]=0
               if time>=infecttime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["i"]=1
           for time in sortedalltimes:    
               assert (nodestates[node][time]["s"]+nodestates[node][time]["i"])==1
    elif spreadmodel=="sir":       
       for node in spreaddata.keys():
           nodestates[node]={}
           if spreaddata[node].has_key("infect"):
              infecttime=spreaddata[node]["infect"]
           else:
              infecttime=10000000000 #very large number
           if spreaddata[node].has_key("recover"):
              recovertime=spreaddata[node]["recover"]
           else:
              recovertime=10000000004 #very large number   
           for time in sortedalltimes:
               nodestates[node][time]={}
               nodestates[node][time]["s"]=1
               nodestates[node][time]["i"]=0
               nodestates[node][time]["r"]=0
               if time>=infecttime and time<recovertime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["i"]=1
               elif time>=recovertime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["r"]=1
           for time in alltimes:    
               assert (nodestates[node][time]["s"]+nodestates[node][time]["i"]+nodestates[node][time]["r"])==1      
    elif spreadmodel=="seir":
       for node in spreaddata.keys():
           nodestates[node]={}
           if spreaddata[node].has_key("exposed"):
              exposedtime=spreaddata[node]["exposed"]
           else:
              exposedtime=10000000000 #very large number
           if spreaddata[node].has_key("infect"):
              infecttime=spreaddata[node]["infect"]
           else:
              infecttime=10000000002 #very large number
           if spreaddata[node].has_key("recover"):
              recovertime=spreaddata[node]["recover"]
           else:
              recovertime=10000000004 #very large number
           for time in sortedalltimes:
               nodestates[node][time]={}
               nodestates[node][time]["s"]=1
               nodestates[node][time]["e"]=0
               nodestates[node][time]["i"]=0
               nodestates[node][time]["r"]=0
               if time>=exposedtime and time<infecttime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["e"]=1
               elif time>=infecttime and time<recovertime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["i"]=1   
               elif time>=recovertime:
                  nodestates[node][time]["s"]=0
                  nodestates[node][time]["r"]=1
           for time in alltimes:    
               assert (nodestates[node][time]["s"]+nodestates[node][time]["e"]+nodestates[node][time]["i"]+nodestates[node][time]["r"])==1
       nodestates=convert2seirstate(nodestates) #e and i will be combined into single i state
    else:
       print "nodestate converge has not been implemented for this spreadmodel {0}".format(spreadmodel)
       exit(1)
    return nodestates
   
    
def roundsolution(algo,retvalues,graphtype,temproundingmethod,graphevolution):
    roundingmethod=temproundingmethod[0]
    roundingparam=temproundingmethod[1]
    if algo in fixededgestaticcovermyalgos: #assuming directed and undirected difference would also have been explained by setcover algo
       if graphtype=="undirected": 
          ourG=nx.Graph()
          for node1,node2 in retvalues.keys():
              ourG.add_edge(node1,node2)
       elif graphtype=="directed": 
          ourG=nx.DiGraph()
          for node1,node2 in retvalues.keys():
              ourG.add_edge(node1,node2)       
    elif algo in fixededgestaticcsalgos or algo in fixededgestaticoptmyalgos:
       if roundingmethod=="epsilon":
          ourG=score2graph(retvalues,roundingmethod,roundingparam,algo,graphtype,graphevolution)
       elif roundingmethod=="graphproperties":
          ourG=score2graph(retvalues,roundingmethod,roundingparam,algo,graphtype,graphevolution)
       else:
	  print "round method {0} is unknown!".format(roundingmethod)
          exit(1)
    elif algo in fixededgestaticotheralgos:
       if graphtype=="undirected": 
          ourG=nx.Graph()
       elif graphtype=="directed": 
          ourG=nx.DiGraph()
       if algo in ["multitree","netinf"]:
          for node1,node2 in retvalues.keys():
              ourG.add_edge(node1,node2)
       elif algo in ["netrate","connie"]:
          assert roundingmethod=="epsilon"
          ourG=score2graph(retvalues,roundingmethod,roundingparam,algo,graphtype,graphevolution) 
    elif algo in fixededgedynamicotheralgos or algo in fixededgedynamicoptmyalgos: #each time point will be rounded independently
       if roundingmethod=="epsilon":
          ourG=score2graph(retvalues,roundingmethod,roundingparam,algo,graphtype,graphevolution)
       elif roundingmethod=="graphproperties": #we may not even use this for dynamic graphs!!
          ourG=score2graph(retvalues,roundingmethod,roundingparam,algo,graphtype,graphevolution)
       else:
	  print "round method {0} is unknown!".format(roundingmethod)
          exit(1)
    else:
       print "rounding strategy for algo {0} is unknown!!".format(algo)
       exit(1)
    return ourG


def graphinference(alltraces,alltracesinfo,uniquecount,timestepslist):   
    print "starting inference for type: {0}!!",format(infertype)
    print "number of traces {0}".format(len(alltraces))
    
    if spreadmodel=="si":
       #spreadparam do not need to be spreadprob(depend on infertype)
       #startnodes,spreadparam,s2iparam=specificparam
       spreadparam,s2iparam=specificparam
    elif spreadmodel=="sis":
       spreadparam,s2iparam,i2sparam=specificparam
    elif spreadmodel=="sir":
       spreadparam,s2iparam,i2rparam=specificparam
    elif spreadmodel=="seir":
       spreadparam,s2eparam,i2rparam,e2iparam=specificparam

    for tracefraction in tracefractions:
        print "fraction is {0}".format(tracefraction)
        for index in range(0,algoruncount):
            print "current algoruncount is {0}".format(index)
            mylist=range(0,len(alltraces))
            print len(mylist)
            random.shuffle(mylist)
            if graphevolution=="dynamic":   
               allposnodes=set()
               for key in Gall.keys():
                   allposnodes|=set(Gall[key].nodes())
               tracecount=int(round(tracefraction*len(allposnodes)))
            elif graphevolution in ["static","paramwise","graphwise"]:
               tracecount=int(round(tracefraction*G.number_of_nodes()))
               print "tracecount is {0}"
               print tracecount
            else:
               print "graphevolution is unknown!! big error".format(graphevolution)
               exit(1)
            assert tracecount<=len(alltraces)   
            print "i will run the algorithm with {0} many traces".format(tracecount)
            selectedindices=mylist[0:tracecount]
            allnodestates=[]
            allnodestatesinfo=[]
            for traceindex in selectedindices:
                spreaddata=alltraces[traceindex]
                nodestates=assignnodestates(spreaddata,spreadmodel) #convert spread data into nodestates hash, 
                if partialdegree!=1.0:             
                   if partialcorrelationfunction=="randomnoise":
                      nodestates=makedatapartial(nodestates,partialdegree,spreadmodel)
                allnodestates.append(nodestates)
                #allnodestatesinfo.append(alltracesinfo[traceindex],specificparam)
                allnodestatesinfo.append(specificparam)
            
            for sampleinterval in sampleintervals:
                #if infertype in ["edge"]:
                #   tempscorenames=edgepredictionscorenames
                #elif infertype in ["spreadprobunknown","difprobpartial"]:
                #   tempscorenames=errorscorenames
                #elif infertype in ["difprobunknown"]:
                #   tempscorenames=vectorscorenames
                #for scorename in tempscorenames:
                #    avgscore[sampleinterval][scorename]=0.0

                realsentresultfolder="{0}/{1}/{2}".format(sentresultfolder,sampleinterval,tracefraction)
                try:
                   os.makedirs(realsentresultfolder)
                except OSError:
                   pass
                #lock = lockfile.FileLock(realsentresultfolder)
                #lock.acquire()
                #os.makedirs(realsentresultfolder)
                #lock.release()
                
                #lock=threading.Lock()
                #lock.acquire()
                #if not os.path.exists(realsentresultfolder):
                #   os.makedirs(realsentresultfolder)
                #lock.release()

                #if any rounding method is missing run algorithm again
                usefulroundingmethods=[]
                for roundingmethod in sentroundingmethods:
                    if roundingmethod[0]=="epsilon":
                       roundingmethodstr="{0}-{1}".format(roundingmethod[0],roundingmethod[1])
                    elif roundingmethod[0]=="graphproperties":
                       roundingmethodstr="{0}-{1}".format(roundingmethod[0],"graphround")
                    elif roundingmethod=="-1":
                       roundingmethodstr="none"
                    else:
                       print "round method {0} is unknown!!".format(roundingmethod)
                       exit(1)
                    
                    tempresultfilenameprefix="result_{0}_{1}_{2}_{3}_{4}".format(sampleselection,graphtype,partialstr,roundingmethodstr,parameterstr)
                    indices=set()
                    for tempfilename in myutil.listfiles(realsentresultfolder):
                        if tempfilename.find(tempresultfilenameprefix+"_")!=-1:
                           myindex=int(tempfilename.replace(tempresultfilenameprefix+"_",""))
                           indices.add(myindex)
                    if missingcompletemode:
                       if len(indices)>=missingcompleteuptocount:
                          print "currently existing file count {0}(>=upperlimit {1}) so not running for this case".format(len(indices),missingcompleteuptocount)
                          continue
                    usefulroundingmethods.append(roundingmethod)
                if len(usefulroundingmethods)==0:
                   continue 
                    
                print "sample interval is {0}".format(sampleinterval)
                if sampleinterval==1:
                   sampleallnodestates=list(allnodestates)
                else:   
                   sampleallnodestates=[]
                   for nodestates in allnodestates:
                       samplenodestates=dict(nodestates)
                       alltimes=set()
                       for node in nodestates.keys():
                           alltimes |= set(nodestates[node].keys()) #same for all nodes so one is enough
                           break
                       sampletimes=set()
                       for time in alltimes:
			   if time%sampleinterval==0:
                              sampletimes.add(time)         
                       for node in nodestates.keys():
                           for time in nodestates[node].keys():
                               if time not in sampletimes:
                                  del samplenodestates[node][time]
                       sampleallnodestates.append(samplenodestates)

                #specificrundir="myrun_{0}_{1}_{2}_{3}_{4}_{5}_{6}".format(tracefraction,sampleinterval,index,partialstr,roundingmethodstr,sampleselection,parameterstr)
                specificrundir="myrun_{0}_{1}_{2}_{3}_{4}_{5}".format(tracefraction,sampleinterval,index,partialstr,sampleselection,parameterstr) 
                specificrunfolder="{0}/{1}".format(sentrunfolder,specificrundir)
                try:
                   os.makedirs(specificrunfolder)
                except OSError:
                   pass
                #lock = lockfile.FileLock(specificrunfolder)
                #lock.acquire()
                #os.makedirs(specificrunfolder)
                #lock.release()

                #lock=threading.Lock()
                #lock.acquire()
                #if not os.path.exists(specificrunfolder):
                #   os.makedirs(specificrunfolder)
                #lock.release()
                print len(sampleallnodestates)
                print "done"

                #timestepslist assignment
                if Gall=="":
                   timesteplist=[]
                elif type(Gall)==dict:
                   timesteplist=Gall.keys() 

                if (algo in fixededgestaticcovermyalgos) or (algo in fixededgedynamiccovermyalgos) or (algo in fixedspreadprobunknowncovermyalgos) or (algo in fixeddifprobpartialcovermyalgos) or (algo in fixeddifprobunknowncovermyalgos): #NOT WORKING!!
                      [temptraceroundingmethod,weightfunction]=algoparams
                      traceroungindmethod=temptraceroundingmethod[0]
                      traceroundingparam=temptraceroundingmethod[1]
                      #check here for some evolution and spread models which we cant infer by set cover
                      retvalues=edgeinferenceoptalgos.greedysetcoverinference(sampleallnodestates,allnodestatesinfo,graphtype,graphevolution,spreadmodel,weightfunction,traceroundingmethod,traceroundingparam)
                elif (algo in fixededgestaticoptmyalgos) or (algo in fixededgedynamicoptmyalgos) or (algo in fixedspreadprobunknownoptmyalgos) or (algo in fixeddifprobpartialoptmyalgos) or (algo in fixeddifprobunknownoptmyalgos):
                      print "running one of my algos {0}".format(algo)
                      if algo in fixededgestaticoptmyalgos:
                         [lambda1,lambda2,logbase,rightmethod,infectedprobabilityassignermethod]=algoparams
                         fusedlambda="-1"
                      elif algo in fixededgedynamicoptmyalgos:
                         [lambda1,lambda2,fusedlambda,logbase,rightmethod,infectedprobabilityassignermethod]=algoparams 
                      print lambda1
                      print lambda2
                      print "fusedlambda is {0}".format(fusedlambda)
                      retvalues=edgeinferenceoptalgos.lseinference(sampleallnodestates,allnodestatesinfo,specificrunfolder,graphevolution,spreadmodel,Gname,algo,graphtype,logbase,lambda1,lambda2,fusedlambda,rightmethod,timesteplist,infectedprobabilityassignermethod,infertype,partialdegree,mapintervals)
                elif algo in fixededgestaticcsalgos: #might not be working properly after modifications!!!! NOT WORKING IF SHOULD BE EXPANDED
                      print "running cs algo {0}".format(algo)
                      [lambda1,logbase,epsilon]=algoparams
                      retvalues=edgeinferenceoptalgos.csinference(sampleallnodestates,allnodestatesinfo,specificrunfolder,graphevolution,spreadmodel,Gname,algo,graphtype,logbase,lambda1)     
                elif (algo in fixededgestaticotheralgos) or (algo in fixededgedynamicotheralgos) or (algo in fixedspreadprobunknownotheralgos) or (algo in fixeddifprobpartialotheralgos) or (algo in fixeddifprobunknownotheralgos):
                      if partialdegree==1.0:
                         tracenoisyornot=False
                      else:
                         tracenoisyornot=True 
                      #matlabpath="/usr/local/bin/matlab"
                      matlabpath="/opt/stow/matlab-r2012a/bin/matlab"
                      print "running other algos by leskovec"
                      #traceroundingmethod=algoparams[0][0][0]
                      #traceroundingparam=algoparams[0][0][1]
                      print algoparams[0]
                      print algoparams[0][1]
                      traceroundingmethod=algoparams[0][0]
                      traceroundingparam=algoparams[0][1]
                      print traceroundingmethod
                      print traceroundingparam 
                      if spreadmodel=="si" and graphevolution in ["static","paramwise","graphwise"]:
                         newalgoparams=[]
                         if algo=="netinf": #This is currently normal neting, not our changed version!!!
                            degreeratio=algoparams[1]
                            iterationcount=int(round(degreeratio*G.number_of_nodes()))
                            print "here info"
                            print iterationcount
                            print degreeratio
                            newalgoparams.append(iterationcount)
                         elif algo=="connie":
                            newalgoparams.append(0.5)
                            newalgoparams.append(matlabpath)
                         elif algo=="multitree":
                            degreeratio=algoparams[1]
                            iterationcount=int(round(degreeratio*G.number_of_nodes()))
                            print iterationcount
                            print degreeratio
                            usedcascadenum=-1;
                            newalgoparams.append(iterationcount)
                            newalgoparams.append(usedcascadenum)
                         elif algo=="netrate":
                            newalgoparams.append(10)
                            newalgoparams.append(matlabpath)
                            newalgoparams.append(G.number_of_nodes())
                         retvalues=edgeinferenceoptalgos.otheralgoinference(sampleallnodestates,allnodestatesinfo,algo,newalgoparams,specificrunfolder,traceroundingmethod,traceroundingparam,tracenoisyornot)
                      else:
                         print "leskovec algorithms can not run on {0} {1}".format(spreadmodel,graphevolution)
                         exit(1)
                else:
                      print "No other comparison algorithm to run for model {0} {1}".format(spreadmodel,graphevolution)
                      exit(1)

                print "retvalue info::"
                print len(retvalues.keys())
                print retvalues.keys()[0:10]
                print "usefulroundingmethods"
                print usefulroundingmethods
                for roundingmethod in usefulroundingmethods:
                    usedretvalues=dict(retvalues)
                    print "finished algo execution"
                    print "running for {0}".format(roundingmethod)
                    if roundingmethod[0]=="epsilon":
                       roundingmethodstr="{0}-{1}".format(roundingmethod[0],roundingmethod[1])
                    elif roundingmethod[0]=="graphproperties":
                       roundingmethodstr="{0}-{1}".format(roundingmethod[0],"graphround")
                    elif roundingmethod=="-1":
                       roundingmethodstr="none"
                    else:
                       print "round method {0} is unknown!!".format(roundingmethod)
                       exit(1)
                    
                    if infertype in ["edge"]: #rounding is here
                       print "started rounding for algo {0}".format(algo)
                       if graphevolution in ["static","paramwise","graphwise"]:
                          ourG=roundsolution(algo,usedretvalues,graphtype,roundingmethod,graphevolution)
                          scores=estimateedgepredictionscore(G,ourG,graphtype,graphevolution)
                       elif graphevolution in ["dynamic"]:
                          ourGall=roundsolution(algo,usedretvalues,graphtype,roundingmethod,graphevolution)
                          print "dynamic graph times:"
                          print ourGall.keys()
                          print Gall.keys()
                          print Gallcompareversion.keys()
                          scores=estimateedgepredictionscore(Gallcompareversion,ourGall,graphtype,graphevolution)   
                    elif infertype in ["spreadprobunknown","difprobpartial"]:
                       if graphtype=="undirected": #what if retvalues are directed, this can only happen on other algos
                          ourG=nx.Graph()
                          for node1,node2 in usedretvalues.keys():
                              ourG.add_edge(node1,node2,weight=usedretvalues(node1,node2))
                       elif graphtype=="directed": 
                          ourG=nx.DiGraph()
                          for node1,node2 in usedretvalues.keys():
                              ourG.add_edge(node1,node2,weight=usedretvalues(node1,node2))     
                       scores=estimateerrorscore(G,ourG,graphtype,graphevolution) 
                    elif infertype in ["difprobunknown"]:
                       #retvalues will be hash, keyed by edges and value will be probability vector       
                       scores=estimatevectorscore(G,usedretvalues,graphtype,graphevolution) 

                    print "scoresinfo::"
                    print scores
                    resultfilenameprefix="result_{0}_{1}_{2}_{3}_{4}".format(sampleselection,graphtype,partialstr,roundingmethodstr,parameterstr)
                    #file outputting part
                    indices=set([1]) #rastgele ekle ki devam edebililelim
                    for tempfilename in myutil.listfiles(realsentresultfolder):
                        if tempfilename.find(resultfilenameprefix+"_")!=-1:   
                           myindex=int(tempfilename.replace(resultfilenameprefix+"_",""))
                           indices.add(myindex)
                    maxindex=max(indices)
                    resultfilepath="-1"
                    for tempindex in range(1,maxindex+2):
                        resultfilename="{0}_{1}".format(resultfilenameprefix,tempindex)
                        resultfilepath="{0}/{1}".format(realsentresultfolder,resultfilename)
                        if os.path.exists(resultfilepath):
                           continue
                        file=open(resultfilepath,"wb")
                        cPickle.dump(scores,file)
                        file.close()
                        break
                    
#NOT BEING USED RIGHT NOW       
#even if we got infertype as paramter, our code DOES NOT DEPEND ON IT
def eliminatetraces(specificparam,alltraces,spreadmodel,infertype):   
    if len(specificparam)==1 and specificparam[0][0]=="all":
       pass
    else:
       for tracename in alltraces.keys():
           if spreadmodel=="si":
              #spreadparam do not need to be spreadprob(depend on infertype)
              print tracename
              startnodes,spreadparam,s2iparam=tracename
           elif spreadmodel=="sis":
              startnodes,spreadparam,s2iparam,i2sparam=tracename
           elif spreadmodel=="sir":
              print tracename 
              startnodes,spreadparam,s2iparam,i2rparam=tracename
           elif spreadmodel=="seir":
              startnodes,spreadparam,s2eparam,i2rparam,e2iparam=tracename
           
           spreadprobflag=False
           s2iflag=False
           i2sflag=False
           i2rflag=False
           s2eflag=False
           e2iflag=False
           spreadprobseen=False
           s2iseen=False
           i2sseen=False
           i2rseen=False
           s2eseen=False
           e2iseen=False
          
           for param in specificparam:
               type,distinfo=param
               if type=="spreadprob":
                  spreadprobseen=True
                  if distinfo==spreadparam:    
                     spreadprobflag=True
               elif type=="s2e":
                  s2eseen=True
                  if distinfo==s2eparam:
                     s2eflag=True        
               elif type=="s2i":
                  s2iseen=True
                  if distinfo==s2iparam:    
                     s2iflag=True
               elif type=="i2s":
                  i2sseen=True
                  if distinfo==i2sparam:    
                     i2sflag=True      
               elif type=="i2r":
                  i2rseen=True
                  if distinfo==i2rparam:    
                     i2rflag=True
               elif type=="e2i":
                  e2iseen=True
                  if distinfo==e2iparam:    
                     e2iflag=True
               else:
                  print "this type is unknown!! {0}".format(type)
                  exit(1)
           if False:
             print "flag seen info!"
             print s2iflag
             print i2rflag
             print e2iflag
             print s2iseen
             print i2rseen
             print e2iseen
           if (spreadprobseen and not spreadprobflag) or (s2iseen and not s2iflag) or (i2rseen and not i2rflag) or (e2iseen and not e2iflag):
              del alltraces[tracename]
              
    return alltraces


#NOE BGINE USED RIGHT NOW
#minus sign can also be in parameter(as for powerlaw)
def specificfolder2spreaddist(specificfolder,spreadmodel):
    if spreadmodel=="si":
       #spreadprob_0.05-s2i_m_expo_2.0_expo_4.0_0.5_0.5
       spreadprobpart,s2ipart=specificfolder.split("-s2i_")
       spreadprob=float(spreadprobpat.replace("spreadprob_",""))
       s2iparts=s2ipart.split("_")
       return (spreadprob,s2istr)
    elif spreadmodel=="sir":
       return
    else:
      print "spreadmodel {0} is unknown!!".format(spreadmodel)
      exit(1)
 
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


edgepredictionscorenames=["sen","fpr","recall","precision","acc","spec","f1","f2","f1overn","mcc"]
errorscorenames=["mse","rmse","mabse","nmse","nrmse","nmabse"]
vectorscorenames=["jsdiv"]
G=""
Gall=""
Gallcompareversion=""
Gname=""
graphtype=""              
algo=""
algoparams=""
parameterstr=""
sentroundingmethods=""
sampleselection=""
tracefractions=""
algoruncount=""
spreadmodel=""
infertype=""
realdata=""
graphevolution=""
partialdegree=""
partialcorrelationfunction=""
partialstr=""
sentresultfolder=""
sentrunfolder=""
tracefilepaths=""
fixededgestaticoptmyalgos=""
fixededgestaticcovermyalgos=""
fixededgestaticotheralgos=""
fixededgestaticcsalgos=""
fixededgedynamicoptmyalgos=""
fixededgedynamiccovermyalgos=""
fixededgedynamicotheralgos=""
fixededgedynamiccsalgos=""
fixeddifprobpartialoptmyalgos=""
fixeddifprobpartialcovermyalgos=""
fixeddifprobpartialotheralgos=""
fixeddifprobpartialcsalgos=""
fixeddifprobunknownoptmyalgos=""
fixeddifprobunknowncovermyalgos=""
fixeddifprobunknownotheralgos=""
fixeddifprobunknowncsalgos=""
sampleintervals=""
missingcompletemode=""
missingcompleteuptocount=""
mapintervals="-1"
if __name__ == "__main__":
    dumppath=sys.argv[1]
    infile=gzip.open(dumppath,"rb")
    G=cPickle.load(infile)
    Gname=cPickle.load(infile)
    graphtype=cPickle.load(infile)

    infertype=cPickle.load(infile)
    algo=cPickle.load(infile)
    algoparams=cPickle.load(infile) #algorithms params
    parameterstr=cPickle.load(infile)
    sentroundingmethods=cPickle.load(infile)
    
    sampleselection=cPickle.load(infile)
    tracefractions=cPickle.load(infile)
    algoruncount=cPickle.load(infile)
    spreadmodel=cPickle.load(infile)
    realdata=cPickle.load(infile)
    graphevolution=cPickle.load(infile)
    
    partialdegree=cPickle.load(infile)
    partialcorrelationfunction=cPickle.load(infile)
    partialstr=cPickle.load(infile)
    sentresultfolder=cPickle.load(infile)
    sentrunfolder=cPickle.load(infile)
    senttracefolder=cPickle.load(infile)

    fixededgestaticoptmyalgos=cPickle.load(infile)
    fixededgestaticcovermyalgos=cPickle.load(infile)
    fixededgestaticotheralgos=cPickle.load(infile)
    fixededgestaticcsalgos=cPickle.load(infile)
    fixededgedynamicoptmyalgos=cPickle.load(infile)
    fixededgedynamiccovermyalgos=cPickle.load(infile)
    fixededgedynamicotheralgos=cPickle.load(infile)
    fixededgedynamiccsalgos=cPickle.load(infile)

    fixeddifprobpartialoptmyalgos=cPickle.load(infile)
    fixeddifprobpartialcovermyalgos=cPickle.load(infile)
    fixeddifprobpartialotheralgos=cPickle.load(infile)
    fixeddifprobpartialcsalgos=cPickle.load(infile)
    fixeddifprobunknownoptmyalgos=cPickle.load(infile)
    fixeddifprobunknowncovermyalgos=cPickle.load(infile)
    fixeddifprobunknownotheralgos=cPickle.load(infile)
    fixeddifprobunknowncsalgos=cPickle.load(infile)
    fixedspreadprobunknownoptmyalgos=cPickle.load(infile)
    fixedspreadprobunknowncovermyalgos=cPickle.load(infile)
    fixedspreadprobunknownotheralgos=cPickle.load(infile)
    fixedspreadprobunknowncsalgos=cPickle.load(infile)
                
    specificparam=cPickle.load(infile)
    specificfolder=cPickle.load(infile)
    sampleintervals=cPickle.load(infile)
    missingcompletemode=cPickle.load(infile)
    missingcompleteuptocount=cPickle.load(infile)
    infile.close()
  
    if graphevolution=="dynamic":
       Gall=G
       print Gall.keys()
       print "myinfo:"
       for time in Gall.keys():
           print "{0},{1},{2}".format(time,G[time].number_of_nodes(),G[time].number_of_edges())
       timestepslist=Gall.keys()
       if realdata=="synthetic":
          dynamicinfo=sentrunfolder.split("/")[1]
          nodenum,temp=dynamicinfo.replace("datan","").split("s")
          stepsize,tcount=temp.split("t")
          nodenum=int(nodenum)
          stepsize=int(stepsize)
          tcount=int(tcount)
          mapintervals={} #useful for dynamic graph inference
          for temptime in range(0,tcount):
              divratio=int(temptime)/stepsize
              mapintervals[temptime]=divratio*stepsize
       elif realdata=="real":
          print "here" 
          print sentrunfolder.split("/")
          mygraphinfo=sentrunfolder.split("/")[1].split("-")[-1]
          tcount,stepsize=mygraphinfo.replace("t","").split("s") #this stepsize has different meaning
          tcount=int(tcount)
          stepsize=int(stepsize)
          mapintervals={} #useful for dynamic graph inference 
          for temptime in range(0,tcount):
              mapintervals[temptime]=temptime%stepsize
       Gallcompareversion={}
       uniquetimes=set(mapintervals.values())
       for time in uniquetimes:
           Gallcompareversion[time]=Gall[time]
    elif graphevolution in ["static","paramwise","graphwise"]:
       timestepslist="-1"
    
    #currently only works for edge!! modify for others!!
    #currently only works for one parameter, not using multiple distribution selection mechanism!!!
    #sys.path.append("../tracepartitioner/")
    #import tracegeneratorpart     
    #specificfolder2=tracegeneratorpart.returnspecificfolder(specificparam)
    #print specificfolder2
    #exir(1)
    newtracefolder="{0}/{1}".format(senttracefolder,specificfolder)
    uniquecount="-1"
    alltraces=[]
    alltracesinfo=[] #start node info
    print "trace folder dir: {0}".format(newtracefolder)
    for tracefilename in myutil.listfiles(newtracefolder):
        #print "trace file is {0}".format(tracefilename)
        tracefilepath="{0}/{1}".format(newtracefolder,tracefilename)
        startnode,randindex=tracefilename.replace(".pkl","").split("_")
        startnode=int(startnode)
        
        #try reading file as zipped, if file is not zipped then read as normal
        mytraces="-1"
        try:
           #print "gzip trace"
           infile=gzip.open(tracefilepath,"rb")
           mytraces=cPickle.load(infile) #hash keyed by nodenumber and params.value will be list since there might be more than one spread from same node with same parameter
           if infertype in ["spreaprobunknown","difprobpartial"]:
              uniquecount=cPickle.load(infile) #this will only be used when estimating the scores 
           infile.close()
        except IOError:
           #print "normal trace" 
           infile=open(tracefilepath,"rb")
           mytraces=cPickle.load(infile) #hash keyed by nodenumber and params.value will be list since there might be more than one spread from same node with same parameter
           if infertype in ["spreaprobunknown","difprobpartial"]:
              uniquecount=cPickle.load(infile) #this will only be used when estimating the scores 
           infile.close()
        except EOFError:
           #print "normal trace" 
           infile=open(tracefilepath,"rb")
           mytraces=cPickle.load(infile) #hash keyed by nodenumber and params.value will be list since there might be more than one spread from same node with same parameter
           if infertype in ["spreaprobunknown","difprobpartial"]:
              uniquecount=cPickle.load(infile) #this will only be used when estimating the scores 
           infile.close()
           
        #startnodes=set()
        #for node in mytraces.keys():
        #    curtime=mytraces[node]["infect"]
        #    if curtime==0:
        #       startnodes.add(node)
        #print startnodes       
        #assert len(startnodes)==1
        
        #infile=open(tracefilepath,"rb")
        #mytraces=cPickle.load(infile) #hash keyed by nodenumber and params.value will be list since there might be more than one spread from same node with same parameter
        #if infertype in ["spreaprobunknown","difprobpartial"]:
        #   uniquecount=cPickle.load(infile) #this will only be used when estimating the scores 
        #infile.close()
        assert mytraces!="-1"   
        alltraces.append(mytraces)
        alltracesinfo.append(startnode)
    print "p1"    
    print len(alltraces)
    print timestepslist
    print "p2"
    graphinference(alltraces,alltracesinfo,uniquecount,timestepslist)
    

#This is package for optimization algorithms for static graph inference
import networkx as nx
import numpy as np
import scipy as sp
import random
import os
import math
import myutilities as myutil
import operator
import cPickle
import heapq
import scipy.io
import scipy.optimize
import sys
import sys
sys.path.append("../tracegenerator/")
import tracegeneratorpart

#For static case:
#There are 2 type of solvers for each case. One is our cplex implemented solver with our restricted variables constraints. Other is generic solvers which does not take our restricted set of constraints into account.
#For dynamic case:
#We will also have additional constraint for relation between the variables at various temporal steps. Therefore, none of the generic solvers can be used in our case. We extend our 8 cplex algorithms with additional temporal smoothness constraint in order to model smooth temporal change over time. We can use various smoothness constraints in this case.
#For instance, kernel based smoothing, fused lasso, or grouped lasso in order to model blockwise changes over time. Maybe we can add exponential smoothing and also our communitybased smoothing algorithm??   

#weightfunction can be
#expodecay -> 1.0/(1-s)^t
#logsum -> -log(1-s)*t
#difference -> 1.0-(1-s)^t
#sets are graph edges. We want to cover elements that are transitions between between time t-1 and time t for each node that are newly infected at time t #when there are multiple costs for an edge, take the average
#this algorithm can only be run when the states are perfectly known
#min set cover will return a parsimonious explanation of the graph. What about comparing this parsimonious approach with maximum likelihood approach????
#graphtype will directly affect the algorithm


#trace rounding
#there can be two types of rounding 1-critical point based rounding(round to determine critical points such as infection time,recovery time, etc) , 2-round each time point independently
#si and sir rounding is the most important one since other algorithms work on si.
def timewiseroundtrace(allnodestate,traceroundingmethod,traceroundingparam,spreadmodel):
    #if spreadmodel=="si":
    #   if traceroundingmethod=="randomtime":        
    return

def criticalpointroundtrace(allnodestate,traceroundingmethod,traceroundingparam,spreadmodel):
    cascadenodes=allnodestate.keys()
    randnode=cascadenodes[0]
    sortedtimes=sorted(allnodestate[randnode].keys())
    if spreadmodel=="si":
       infecttimes={}
       if traceroundingmethod=="error":
          for node in cascadenodes: #there will be only single infection time for each node
              probhash={}
              posinfected=set()
              for time in sortedtimes:
                  if allnodestate[node][time]["i"]>0.0000001:
                     posinfected.add(time)
              if len(posinfected)==0:
                 continue 
              foundone="-1"
              for infecttime in posinfected:
                  error=0.0
                  for time in sortedtimes:
                      if time<infecttime:
                         error+=allnodestate[node][time]["i"]**2
                      elif time>=infecttime:
                         error+=(1.0-allnodestate[node][time]["i"])**2
                  if error==0.0:
                     foundone=infecttime 
                     break 
                  probhash[infecttime]=1.0/error
              if foundone!="-1":
                 infecttimes[node]=foundone 
                 continue
              totprob=0.0
              for value in probhash.values():
                  totprob+=value
              for key in probhash.keys():
                  probhash[key]/=float(totprob)
              probinterval={}
              cursum=0.0
              for time in sortedtimes:
                  if probhash.has_key(time):
                     probinterval[(cursum,cursum+probhash[time])]=time
                     cursum += probhash[time]
              p=random.random()
              roundtime=-1
              for (left,right) in probinterval.keys():
                   if left<=p and p<=right:
                      roundtime=probinterval[(left,right)]
                      break
              assert roundtime!=-1
              infecttimes[node]=roundtime    
       elif traceroundingmethod=="categorical": #categorical dist rounding
          for node in cascadenodes: #there will be only single infection time for each node
               probhash={}
               avginfected=0.0
               infectedcount=0
               for time in sortedtimes:
                   if allnodestate[node][time]["i"]!=0:
                      probhash[time]=allnodestate[node][time]["i"]
                      avginfected += allnodestate[node][time]["i"]
                      infectedcount += 1
               if len(probhash.keys())==0:
                  continue        
               avginfected /= float(infectedcount)       
               if random.random()<=(1.0-avginfected): #this node will be assumed not infected at all
                  continue
               mysum=0.0
               for time in probhash.keys():
                   mysum += probhash[time]
               for time in probhash.keys():
                   probhash[time] /= float(mysum)
               probinterval={}
               cursum=0.0
               for time in sortedtimes:
                   if probhash.has_key(time):
                      probinterval[(cursum,cursum+probhash[time])]=time
                      cursum += probhash[time]
               p=random.random()
               roundtime=-1
               for (left,right) in probinterval.keys():
                   if left<=p and p<=right:
                      roundtime=probinterval[(left,right)]
                      break
               assert roundtime!=-1
               infecttimes[node]=roundtime
       elif traceroundingmethod=="random": #random dist rounding
          print "random rounding has not been implemented yet, error!!"
          exit(1)
       elif traceroundingmethod=="numeric":#numeric rounding, round to sth above threshold and best time point
          threshold=traceroundingparam  
          for node in cascadenodes: #there will be only single infection time for each node 
              besttimeval=-1
              besttime=-1
              for time in sortedtimes:
                  if allnodestate[node][time]["infected"]>=threshold and allnodestate[node][time]["infected"]>besttimeval :
                     besttime=time
                     besttimeval=allnodestate[node][time]["infected"]
              if besttime!=-1:
                 infecttimes[node]=besttime
       else:
          print "error, this traceroundingmethod is unknow!! {0}".format(traceroundingmethod)
          exit(1)
       return [infecttimes]    
    elif spreadmodel=="sir": #avginfected will be important for determining whether round sth eventually to infected state
       #we dont need rounding for other methods since other algos can only run for si 
       infecttimes={}
       recovertimes={}
       if traceroundingmethod=="categorical": #categorical dist rounding, round both states independently and then check infected<recovered, otherwise repeat
          for node in cascadenodes: #there will be only single infection time for each node
              tempinfected=-1
              temprecovered=-1
              while True: #this while is USELESS right now, i switched to do only single trial
                 #determining infection time
                 probhash={}
                 avginfected=0.0
                 infectedcount=0
                 for time in sortedtimes:
                     if allnodestate[node][time]["i"]!=0:
                        probhash[time]=allnodestate[node][time]["i"]
                        avginfected += allnodestate[node][time]["i"]
                        infectedcount += 1
                 if len(probhash.keys())==0:
                    break       
                 avginfected /= float(infectedcount)       
                 if random.random()<=(1.0-avginfected): #this node will be assumed not infected at all
                    break
                 mysum=0.0
                 for time in probhash.keys():
                     mysum += probhash[time]
                 for time in probhash.keys():
                     probhash[time] /= mysum
                 probinterval={}
                 cursum=0.0
                 for time in sortedtimes:
                     if probhash.has_key(time):
                        probinterval[(cursum,cursum+probhash[time])]=time
                        cursum+=probhash[time]
                 p=random.random()
                 roundtime=-1
                 for (left,right) in probinterval.keys():
                     if left<=p and p<=right:
                        roundtime=probinterval[(left,right)]
                        break
                 assert roundtime!=-1
                 tempinfected=roundtime

                 #determining recovery time
                 probhash={}
                 avgrecovered=0.0
                 recoveredcount=0
                 for time in sortedtimes:
                     if allnodestate[node][time]["i"]!=0:
                        probhash[time]=allnodestate[node][time]["i"]
                        avgrecovered += allnodestate[node][time]["i"]
                        recoveredcount += 1
                 if len(probhash.keys())==0:
                    break       
                 avgrecovered /= float(recoveredcount)       
                 if random.random()<=(1.0-avgrecovered): #this node will be assumed not infected at all
                    break
                 mysum=0.0
                 for time in probhash.keys():
                     mysum += probhash[time]
                 for time in probhash.keys():
                     probhash[time] /= mysum
                 probinterval={}
                 cursum=0.0
                 for time in sortedtimes:
                     if probhash.has_key(time):
                        probinterval[(cursum,cursum+probhash[time])]=time
                        cursum += probhash[time]
                 p=random.random()
                 roundtime=-1
                 for (left,right) in probinterval.keys():
                     if left<=p and p<=right:
                        roundtime=probinterval[(left,right)]
                        break
                 assert roundtime!=-1
                 temprecovered=roundtime
                 
                 if temprecovered!=-1 and tempinfected!=-1:
                    if temprecovered>=tempinfected:
                       break
                    else:
                       infecttimes[node]=tempinfected
                       recovertimes[node]=temprecovered 
                 break     
       elif traceroundingmethod=="random": #random dist rounding
          print "random rounding has not been implemented yet, error!!"
          exit(1)
       elif traceroundingmethod=="numeric":#numeric rounding, round to sth above threshold and best time point
          threshold=traceroundingparam  
          for node in cascadenodes: #there will be only single infection time for each node 
              bestinfectedtimeval=-1
              bestinfectedtime=-1
              for time in sortedtimes:
                  if allnodestate[node][time]["infected"]>=threshold and allnodestate[node][time]["infected"]>bestinfectedtimeval :
                     bestinfectedtime=time
                     bestinfectedtimeval=allnodestate[node][time]["infected"]
              bestrecoveredtimeval=-1
              bestrecoveredtime=-1
              for time in sortedtimes:
                  if allnodestate[node][time]["recovered"]>=threshold and allnodestate[node][time]["recovered"]>bestrecoveredtimeval :
                     bestrecoveredtime=time
                     bestrecoveredtimeval=allnodestate[node][time]["recovered"]       
              if bestinfectedtime!=-1 and bestrecoveredtime!=-1 and bestinfectedtime<bestrecoveredtime:
                 infecttimes[node]=bestinfectedtime
                 recovertimes[node]=bestrecoveredtime
       else:
          print "error, this traceroundingmethod is unknow!! {0}".format(traceroundingmethod)
          exit(1)
       return [infecttimes,recovertimes]   
    else:
       print "rounding mechanism for {0} has not yet implemented!!".format(spreadmodel)
       exit(1)

       
#A node that has infected value above threshold will be assumed to be infected and its infection time will be time when it has maximum infected value.
#threshold can be either number, "categorical", "random"
#nodemap is mapping nodes to 0-n-1 interval for netrate algorithm       
def createtracefile(allnodestates,tracefilepath,traceroundingmethod,traceroundingparam,seperator,tracenoisyornot,maptype="nomap"):
    print "trace file creation"
    print len(allnodestates)
    print traceroundingmethod
    print traceroundingparam
    node2newnode={}
    if maptype=="nodemap":
       tempnodes=set()
       for allnodestate in allnodestates:
           tempnodes |= set(allnodestate.keys())
       sortednodes=sorted(list(tempnodes))
       for index in range(0,len(sortednodes)):
           node2newnode[sortednodes[index]]=index
    elif maptype=="nomap":
       tempnodes=set()
       for allnodestate in allnodestates:
           tempnodes |= set(allnodestate.keys())
       tempnodes=list(tempnodes)    
       for index in range(0,len(tempnodes)):
           node2newnode[tempnodes[index]]=tempnodes[index]
    else:
       print "maptype {0} is unknown:".format(maptype)
       exit(1)
       
    allnodes=set()
    for allnodestate in allnodestates:
        allnodes |= set(allnodestate.keys())
    print "inside trace generation"
    print len(allnodes)
    file=open(tracefilepath,"w")
    for node in allnodes:
        #file.write("{0},{1}\n".format(node,node))
        file.write("{0},{1}\n".format(node2newnode[node],node2newnode[node]))
    file.write("\n")

    for allnodestate in allnodestates:
        if not tracenoisyornot: #that means trace is not noisy
           infecttimes={}
           for node in allnodestate.keys():
               for time in sorted(allnodestate[node].keys()):
                   if allnodestate[node][time]["i"]==1:
                      infecttimes[node]=time
                      break
        else:        
           infecttimes=criticalpointroundtrace(allnodestate,traceroundingmethod,traceroundingparam,"si")[0]
        time2node={}
        for tempnode in infecttimes.keys():
            mytime=infecttimes[tempnode]
            if not time2node.has_key(mytime):
               time2node[mytime]=set()
            time2node[mytime].add(tempnode)   
        usedtimes=sorted(time2node.keys())
        outlist=[]
        for usedtime in usedtimes:
            for usednode in time2node[usedtime]:
                partstr="{0}{1}{2}".format(node2newnode[usednode],",",usedtime)
                #partstr="{0}{1}{2}".format(usednode,",",usedtime)
                outlist.append(partstr)
        tempstr=seperator.join(outlist)
        file.write("{0}\n".format(tempstr))
    file.close()
    return node2newnode #this will only be useful for netrate kind of algorithms


#NOT BEING USED RIGHT NOW. THIS OUR CHANGED CUTOFF VERSION
#codedirname will be given as input
#this will return output assuming graph is directed
#truegraphfilepath will be used as param but won't be useful
#netinf assumes graph is directed
def si_static_netinf_cutoff(allnodestates,outrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,cutoffratio):
    print "running netinf"
    spreaddist=s2idist[0]
    spreadparam=s2idist[1:]
    if spreaddist=="expo":
       disttype=0
       alpha=1.0/spreadparam[0]
    elif spreaddist=="powerlaw":
       disttype=1
       alpha=-1.0*spreadparam[0]
    elif spreaddist=="rayleigh": #not sure how parameter is converted!!
       disttype=2
       alpha=spreadparam[0]
    else:
       print "netinf can not be run with dist type {0}".format(spreaddist)
       exit(1)
    print "outrunfolder"
    print outrunfolder
    outputprefix="-".join(outrunfolder.split("/"))
    outputprefix +="-{0}-{1}".format(spreaddist,spreadparam[0])
    outputprefix="{0}-{1}".format("netinf",outputprefix)
    tracefilepath="{0}/{1}.netinf".format(codedirname,outputprefix)
    print "prefix info!"
    print outputprefix
    allnodes=set()
    for nodestate in allnodestates:
        allnodes |= set(nodestate.keys())
    print "my node num {0}".format(len(allnodes))
    s=2
    seperator=";"
    createtracefile(allnodestates,tracefilepath,traceroundingmethod,traceroundingparam,seperator)
    codefile="netinf"
    codepath="{0}/{1}".format(codedirname,codefile)
    maxiterationnum=len(allnodes)*(len(allnodes)-1)/2
    code="{0} -i:{1} -o:{2} -e:{3} -a:{4} -m:{5} -s:{6} -f:{7} -b:{8}".format(codepath,tracefilepath,outputprefix,maxiterationnum,alpha,disttype,s,cutoffratio,100)
    os.system(code)
    edge2value={}
    outputfilename="{0}.txt".format(outputprefix)
    file=open(outputfilename,"r") #read solution file
    #count=0
    #for line in file:
    #    if count!=0:
    #        line=line.rstrip()
    #       node1,node2,volume,marginalgain,temp1,temp2=line.split("/")
    #       node1=int(node1)
    #       node2=int(node2)
    #       marginalgain=float(marginalgain)
    #       edge2value[(node1,node2)]=marginalgain
    #    else:
    #       count+=1
    #file.close()
    for line in file:
        line=line.rstrip()
        if line=="":
           continue 
        node1,node2=line.split(",")
        node1=int(node1)
        node2=int(node2)
        if node1==node2:
           continue
        edge2value[(node1,node2)]=1.0
    file.close()
    
    #extensions=["tab","plt","eps"]
    #extensions2=["plt","eps"]
    #deletefilename="{0}.txt".format(outputprefix)
    #deletefilenames=[deletefilename]
    #for extension in extensions:
    #    deletefilename="{0}-objective.{1}".format(outputprefix,extension)
    #    deletefilenames.append(deletefilename)
    #for extension in extensions2:    
    #    deletefilename="{0}-obj.{1}".format(outputprefix,extension)
    #    deletefilenames.append(deletefilename)
    #print "deletefilenames"
    #print deletefilenames
    #for deletefilename in deletefilenames:
    #    assert os.path.exists(deletefilename)
    #    code="rm -rf {0}".format(deletefilename)
    #    os.system(code)
    code="rm -rf {0}".format(tracefilepath)
    os.system(code)
    print "num of edges {0}".format(len(edge2value.keys()))
    return edge2value



#codedirname will be given as input
#this will return output assuming graph is directed
#truegraphfilepath will be used as param but won't be useful
#netinf assumes graph is directed
def si_static_netinf(allnodestates,outrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,iterationcount,tracenoisyornot,probtype):
    print "inside netinf!!"
    print iterationcount
    print "running netinf"
    spreaddist=s2idist[0]
    spreadparam=s2idist[1:]
    if probtype=="discrete":
       if spreaddist=="expo":
          disttype=0
          alpha=1.0/spreadparam[0]
       elif spreaddist=="powerlaw":
          disttype=1
          alpha=-1.0*spreadparam[0]
       elif spreaddist=="rayleigh":
          disttype=2
          alpha=1.0/(spreadparam[0]**2)
       elif spreaddist=="weibull":
          disttype=3
          alpha=spreadparam[0] 
       else:
          print "netinf can not be run with dist type {0}".format(spreaddist)
          exit(1)
    elif probtype=="interval":
       if spreaddist=="expo":
          disttype=4
          alpha=1.0/spreadparam[0]
       elif spreaddist=="powerlaw":
          disttype=5
          alpha=-1.0*spreadparam[0]
       elif spreaddist=="rayleigh":
          disttype=6
          alpha=1.0/(spreadparam[0]**2)
       elif spreaddist=="weibull":
          disttype=7
          alpha=spreadparam[0] 
       else:
          print "netinf can not be run with dist type {0}".format(spreaddist)
          exit(1)
    else:
         print "netinf can not be run with dist type {0}".format(spreaddist)
         exit(1)     
    print "outrunfolder"
    print outrunfolder
    outputprefix="-".join(outrunfolder.split("/"))
    outputprefix +="-{0}-{1}-{2}".format(spreaddist,spreadparam[0],probtype)
    outputprefix="{0}-{1}".format("netinf",outputprefix)
    tracefilepath="{0}/{1}.netinf".format(codedirname,outputprefix)
    print "prefix info!"
    print outputprefix
    allnodes=set()
    for nodestate in allnodestates:
        allnodes |= set(nodestate.keys())
    print "my node num {0}".format(len(allnodes))
    s=1
    seperator=";"
    createtracefile(allnodestates,tracefilepath,traceroundingmethod,traceroundingparam,seperator,tracenoisyornot)
    codefile="netinf"
    codepath="{0}/{1}".format(codedirname,codefile)
    code="{0} -i:{1} -o:{2} -e:{3} -a:{4} -m:{5} -s:{6}".format(codepath,tracefilepath,outputprefix,iterationcount,alpha,disttype,s)
    print "code is:"
    print code
    #exit(1)
    os.system(code)
    edge2value={}
    outputfilename="{0}.txt".format(outputprefix)
    #outputfilename="{0}-edge.info".format(outputprefix)
    file=open(outputfilename,"r") #read solution file
    #count=0
    #for line in file:
    #    if count!=0:
    #       line=line.rstrip()
    #       node1,node2,volume,marginalgain,temp1,temp2=line.split("/")
    #       node1=int(node1)
    #       node2=int(node2)
    #       marginalgain=float(marginalgain)
    #       edge2value[(node1,node2)]=marginalgain
    #    else:
    #       count+=1
    #file.close()
    for line in file:
        line=line.rstrip()
        if line=="":
           continue 
        node1,node2=line.split(",")
        node1=int(node1)
        node2=int(node2)
        if node1==node2:
           continue
        edge2value[(node1,node2)]=1.0
    file.close()
    
    extensions=["-edge.info",".txt"]
    for extension in extensions:
        deletefilename="{0}{1}".format(outputprefix,extension)
        code="rm -rf {0}".format(deletefilename)
        os.system(code)
        
    #extensions=["tab","plt","eps"]
    #extensions2=["plt","eps"]
    #deletefilename="{0}.txt".format(outputprefix)
    #deletefilenames=[deletefilename]
    #for extension in extensions:
    #    deletefilename="{0}-objective.{1}".format(outputprefix,extension)
    #    deletefilenames.append(deletefilename)
    #for extension in extensions2:    
    #    deletefilename="{0}-obj.{1}".format(outputprefix,extension)
    #    deletefilenames.append(deletefilename)
    #print "deletefilenames"
    #print deletefilenames
    #for deletefilename in deletefilenames:
    #    assert os.path.exists(deletefilename)
    #    code="rm -rf {0}".format(deletefilename)
    #    os.system(code)
    code="rm -rf {0}".format(tracefilepath)
    os.system(code)
    print "num of edges {0}".format(len(edge2value.keys()))
    return edge2value
       
def matrix2matlabstr(mat):
    matstr=""
    for index1 in range(0,np.shape(mat)[0]):
        for index2 in range(0,np.shape(mat)[1]):
            matstr+=" {0} ".format(mat[index1,index2])
        matstr+=";"
    matstr=matstr[0:-1]
    matstr ="[{0}]".format(matstr)
    return matstr

#due to snopt, maximum diffusion number and nodenumber is 300
#not recent one, needs to be improved
def si_static_connie(allnodestates,outrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,lambda1,matlabpath,tracenoisyornot,probtype):
    spreaddist=s2idist[0]
    spreadparam=s2idist[1:]
    if spreaddist=="expo": #only with lambda 1
       disttype=2
       assert spreadparam[0]==1.0
    elif spreaddist=="powerlaw": #only alpha -2
       disttype=1
       assert spreadparam[0]==-2.0
    elif spreaddist=="uniform": 
       disttype=3
       assert spreadparam==1
    elif spreaddist=="weibull": #lambda = 9.479 and k = 2.3494
       disttype=4
       assert spreadparam[0]==9.749 and spreadparam[1]==2.3494
    else:
       print "connie can not be run with dist type {0}".format(spreaddist)
       exit(1)

    difnum=len(allnodestates)
    allnodes=set()
    for allnodestate in allnodestates:
        tempnodes=allnodestate.keys()
        allnodes=allnodes.union(tempnodes)
    print len(allnodes)    
    node2index={}
    nodeindex=0
    for node in allnodes:
        node2index[node]=nodeindex
        nodeindex += 1
    index2node={}
    for node in node2index.keys():
        index2node[node2index[node]]=node
    difmatrix=np.zeros((difnum,len(allnodes)),dtype=np.int)
    for difindex in range(0,difnum):
        #if difindex==2:
        #   break
        allnodestate=allnodestates[difindex]
        for allnodestate in allnodestates:
            infecttimes=criticalpointroundtrace(allnodestate,traceroundingmethod,traceroundingparam,"si")[0]
            #myi=0
            for node in allnodes:
                #myi+=1
                #if myi==200:
                #   break 
                nodeindex=node2index[node]
                if infecttimes.has_key(node):
                   difmatrix[difindex,nodeindex]=infecttimes[node]
                else:
                   difmatrix[difindex,nodeindex]=-1
    difmatrixstr=matrix2matlabstr(difmatrix)
    olddirname=os.getcwd()
    newdirname=codedirname
    suboptimal_tol=0.1
    outputprefix="-".join(outrunfolder.split("/"))
    outputprefix +="-{0}-{1}-{2}".format(spreaddist,spreadparam[0],probtype)
    outputprefix="{0}-{1}".format("connie",outputprefix)
    print "prefix info!"
    print outputprefix
    print difmatrixstr[0:100]
    print olddirname
    print newdirname
    outfilename="{0}.mat".format(outputprefix)
    outfilepath="{0}/{1}/{2}".format(olddirname,outrunfolder,outfilename)
    print outfilepath
    codepath="connie"
    #code="opt/stow/matlab-r2012a/bin/matlab -r \" cd('./algos/connie'); connie(0.2,2,[1,2,3;2,3,4;1,2,2],0.2,'../myoutfilename.txt') ; cd('../../'); quit; \"" 
    code="{0} -r \" cd('{1}'); {2}({3},{4},{5},{6},'{7}') ; cd('{8}'); quit; \" ".format(matlabpath,newdirname,codepath,lambda1,disttype,difmatrixstr,suboptimal_tol,outfilepath,olddirname)
    os.system(code)
    
    data=scipy.io.loadmat(outfilepath)
    difmat=data['A_mle']
    edge2value={}
    for index1 in range(0,np.shape(difmat)[0]):
        node1=index2node[index1]
        for index2 in range(0,np.shape(difmat)[1]):
            if index1==index2:
               continue
            node2=index2node[index2]
            edge2value[(node1,node2)]=difmat[index1,index2]
    print "key count {0}".format(len(edge2value.keys()))        
    code="rm -rf {0}".format(outfilepath)
    os.system(code)
    return edge2value        

#netrate assumes nodes are sorted starting from index 0 to n-1
def si_static_netrate(allnodestates,outrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,windowsize,matlabpath,allnodenum,tracenoisyornot,probtype):
    spreaddist=s2idist[0]
    spreadparam=s2idist[1:]
    if spreaddist=="expo":
       disttype='exp'
       alpha=spreadparam[0]
    elif spreaddist=="powerlaw":
       disttype='pl'
       alpha=-1.0-(1.0*spreadparam[0])
    elif spreaddist=="rayleigh":
       disttype='rayleigh'
       alpha=0.5/(spreadparam[0]**2)
    else:
       print "Netrate can not be run with dist type {0}".format(spreaddist)
       exit(1)
           
    olddirname=os.getcwd()
    newdirname=codedirname
    outputprefix="-".join(outrunfolder.split("/"))
    outputprefix +="-{0}-{1}-{2}".format(spreaddist,spreadparam[0],probtype)
    outputprefix="{0}-{1}".format("netrate",outputprefix)
    print "prefix info!"
    print "outrunfolder is: {0}".format(outrunfolder)
    print outputprefix
    print olddirname
    print newdirname
    outdirname2="{0}/{1}".format(olddirname,outrunfolder)
    if not os.path.exists(outdirname2):
       os.makedirs(outdirname2) 
    outfilename="{0}.mat".format(outputprefix)
    outfilepath="{0}/{1}".format(outdirname2,outfilename)
    print "outfilepath is {0}".format(outfilepath)
    codepath="netrate"
    tracefilename="{0}.netratetrace".format(outputprefix)
    tracefilepath="{0}/{1}".format(codedirname,tracefilename)
    seperator=","
    node2newnode=createtracefile(allnodestates,tracefilepath,traceroundingmethod,traceroundingparam,seperator,tracenoisyornot,maptype="nodemap")
    print tracefilepath
    
    cvxfolder="cvx"
    cvxreturnfolder=".."
    code="{0} -r \" cd('{1}'); cd('{2}'); cvx_setup; cd('{3}'); {4}('{5}','{6}',{7},'{8}',{9},'{10}'); cd('{11}'); quit; \" ".format(matlabpath,newdirname,cvxfolder,cvxreturnfolder,codepath,"[]",tracefilename,windowsize,disttype,allnodenum,outfilepath,olddirname)
    print code
    os.system(code)
    
    data=scipy.io.loadmat(outfilepath)
    print np.shape(data)
    difmat=data['A_hat']
    print "info about matrix"
    print np.shape(difmat)
    print difmat[0,2]
    print difmat[0,4]
    print difmat[2,5]
    print difmat[1,8]
    edge2value={}
    for node1 in range(0,np.shape(difmat)[0]):
        for node2 in range(0,np.shape(difmat)[1]):
            if node1==node2:
               continue
            edge2value[(node1,node2)]=difmat[node1,node2]
    print "key count {0}".format(len(edge2value.keys()))        
    code="rm -rf {0}".format(outfilepath)
    os.system(code)
    code="rm -rf {0}".format(tracefilepath)
    os.system(code)
    return edge2value

    
def si_static_multitree(allnodestates,outrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,iterationnum,usedcascadenum,tracenoisyornot,probtype):
    spreaddist=s2idist[0]
    spreadparam=s2idist[1:]
    delta=1.0 #not important, just random
    if probtype=="discrete":
       if spreaddist=="expo":
          disttype=0
          alpha=1.0/spreadparam[0]
       elif spreaddist=="powerlaw":
          disttype=1
          alpha=-1.0*spreadparam[0]
       elif spreaddist=="rayleigh":
          disttype=2
          alpha=1.0/(spreadparam[0]**2)
       elif spreaddist=="weibull":
          disttype=3
          alpha=spreadparam[0] 
       else:
          print "multitree can not be run with dist type {0}".format(spreaddist)
          exit(1)
    elif probtype=="interval":
       if spreaddist=="expo":
          disttype=4
          alpha=1.0/spreadparam[0]
       elif spreaddist=="powerlaw":
          disttype=5
          alpha=-1.0*spreadparam[0]
       elif spreaddist=="rayleigh":
          disttype=6
          alpha=1.0/(spreadparam[0]**2)
       elif spreaddist=="weibull":
          disttype=7
          alpha=spreadparam[0]
       else:
          print "multitree can not be run with dist type {0}".format(spreaddist)
          exit(1)
    else:
         print "multitree can not be run with dist type {0}".format(spreaddist)
         exit(1)
    outputprefix="-".join(outrunfolder.split("/"))
    outputprefix += "-{0}-{1}-{2}".format(spreaddist,spreadparam[0],probtype)
    outputprefix="{0}-{1}".format("multitree",outputprefix)
    tracefilepath="{0}/{1}.netratetrace".format(codedirname,outputprefix)
    allnodes=set()
    for nodestate in allnodestates:
        allnodes |= set(nodestate.keys())
    s=1
    seperator=","
    outputfilename="{0}-edge.info".format(outputprefix)
    createtracefile(allnodestates,tracefilepath,traceroundingmethod,traceroundingparam,seperator,tracenoisyornot)
    codefile="network-inference-multitree"
    codepath="{0}/{1}".format(codedirname,codefile)
    code="{0} -i:{1} -o:{2} -e:{3} -a:{4} -d:{5} -nc:{6} -m:{7} -s:{8}".format(codepath,tracefilepath,outputprefix,iterationnum,alpha,delta,usedcascadenum,disttype,s)
    #./network-inference-multitree -i example-cascades.txt -o myout -e 5 -a 1 -nc -1 -m 1 -s 1
    os.system(code)
    edge2value={}
    file=open(outputfilename,"r") #read solution file
    count=0
    for line in file:
        if count!=0:
           line=line.rstrip()
           node1,node2,vol,marginalgain,temp3,temp4=line.split("/")
           node1=int(node1)
           node2=int(node2)
           marginalgain=float(marginalgain)
           edge2value[(node1,node2)]=marginalgain
        else:
           count+=1 
    file.close()
    code="rm -rf {0}".format(outputfilename)
    os.system(code)
    code="rm -rf {0}".format(tracefilepath)
    os.system(code)
    return edge2value

    
#this will only make sense if we will be using same distribution parameters for all of them
#those algos can only be run for SI model
def otheralgoinference(allnodestates,allnodestatesinfo,algo,algoparams,specificrunfolder,traceroundingmethod,traceroundingparam,tracenoisyornot,probtype):
    print "inside other algoss"
    print allnodestatesinfo[0]
    [spreadparam,s2idist]=allnodestatesinfo[0] #spreaddist and spreadparam will be obtained here
    for traceindex in range(1,len(allnodestates)):
        nodestates=allnodestates[traceindex]
        [spreadparam1,s2idist1]=allnodestatesinfo[traceindex]
        assert spreadparam==spreadparam1
        assert s2idist==s2idist1
    print spreadparam
    print s2idist
    spreadparam=spreadparam[1]
    s2idist=s2idist[1]
    print spreadparam
    print s2idist
    # NETINF CUTOFFF VERSION S NOT BEIN USED RIGHT NOW
    #if algo=="netinf":
    #   cutoffratio=algoparams[0]
    #   codedirname="./algos/netinf-unix"
    #   retvalues=si_static_netinf(allnodestates,specificrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,cutoffratio)
    if algo=="netinf":
       iterationcount=algoparams[0]
       codedirname="./algos/netinf-unix"
       retvalues=si_static_netinf(allnodestates,specificrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,iterationcount,tracenoisyornot,probtype)   
    elif algo=="netrate":
       codedirname="./algos/netrate"
       windowsize=algoparams[0]
       matlabpath=algoparams[1]
       allnodenum=algoparams[2] #allnodenum must be given as a parameter since netrate assumes nodes are sorted from 0 to n-1!!
       print "prenetrate info"
       print windowsize
       print matlabpath
       print allnodenum
       retvalues=si_static_netrate(allnodestates,specificrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,windowsize,matlabpath,allnodenum,tracenoisyornot,probtype)
    elif algo=="connie": #this is implemented but won't be using due to snopt limitations
       lambda1=algoparams[0]
       matlabpath=algoparams[1]
       codedirname="./algos/connie"
       retvalues=si_static_connie(allnodestates,specificrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,lambda1,matlabpath,tracenoisyornot,probtype)
    elif algo=="multitree":
       codedirname="./algos/multitree"
       iterationnum=algoparams[0]
       usedcascadenum=algoparams[1]
       retvalues=si_static_multitree(allnodestates,specificrunfolder,s2idist,traceroundingmethod,traceroundingparam,codedirname,iterationnum,usedcascadenum,tracenoisyornot,probtype)
    return retvalues



#return given probability distribution to our format, probdist must sum up to 1. If not, make it happen
#probdist should start from 1 not 0
#not being used right now!!!
def prob2spreaddist(probdist):
    totalsum=0.0
    for pro in probdist:
        totalsum += pro
    normprobdist=[]
    for pro in probdist:
        tempval=pro/float(totalsum)
        normprobdist.append(tempval)
    probdist=list(normprobdist)    
    sumsofar=0.0
    betas=[]
    rights=[]
    rights.append(1.0)
    for prob in probdist:    
        sumsofar += prob
        rights.append(1.0-sumsofar)
    for index in range(1,len(rights)):
        ratio=rights[index]/rights[index-1]
        curbeta=1.0-ratio
        betas.append(curbeta)
    #test
    assert len(betas)==len(probdist)    
    betasofar=1.0
    for betaindex in range(0,len(betas)):
        beta=betas[betaindex]
        betasofar *= (1.0-beta)
        cumsum=0.0
        for index in range(0,betaindex+1):
            cumsum += probdist[index]
        right=1.0-cumsum    
        assert betasofar==right
    return betas


#node becomes infected when they first become greater than firstthres and node becomes recovered when node drops below secondthres after becoming greater than firstthres
#for certain cases, firsthres will be 1 and secondthres will be 0
#NOT BEING USED RIGHT NOW!!
def discreteinfecttimeassigner(nodestates,firstthres,secondthres):
    infecttimes={}
    randnode=nodestates.keys()[0]
    alltimes=nodestates[randnode].keys()
    for node in nodestates.keys():
        infecttime=-1
        for time in sorted(alltimes):
            if nodestates[node][time]["i"] >= firstthres:
               infecttime=time
               break
        if infecttime!=-1: #infected
           infecttimes[node]=infecttime
    recovertimes={}
    for node in nodestates.keys():
        if not infecttimes.has_key(node):
           continue 
        recovertime=-1
        for time in sorted(alltimes):
            if time<=infecttimes[node]:
               continue 
            if nodestates[node][time]["i"] <= secondthres:
               recovertime=time
               break
        recovertimes[node]=recovertime       
    difset=set(infecttimes.keys()).difference(set(recovertimes.keys()))
    assert len(difset)==0
    return [infecttimes,recovertimes]



#Current version assumes the knowledge has not been corrupted by noise. if corrupted, maybe we should take infection duration distribution into account on weights???
#this can also handle fractional states via threshold based rounding    
def sir_static_greedysetcover(allnodestates,allnodestatesinfo,graphtype,weightfunction,spreadtimemode,lagdepth,infinitylag,traceroundingmethod,traceroundingparam):
    if spreadtimemode=="asyn": #this has not yet been implemeneted for asyn case
       return
    elem2set=[]
    set2elem={}
    cost={} #estimate average cost for each element(set)
    count={}
    for traceindex in range(0,len(allnodestates)):
        nodestates=allnodestates[traceindex]
        [spreadprob,spreadist,spreadparam,infecteddist,intectedparam]=allnodestatesinfo[traceindex]
        infecttimes,recovertimes=criticalpointroundtrace(nodestates,traceroundingmethod,traceroundingparam,"sir")
        allnodes=nodestates.keys()
        randnode=nodestates.keys()[0]
        alltimes=nodestates[randnode].keys()
        sortedtimes=sorted(alltimes)
        if graphtype=="undirected":
           sortednodes=sorted(allnodes)
           edge2index={}
           for index1 in range(0,len(sortednodes)):
               node1=sortednodes[index1]
               for index2 in range(index1+1,len(sortednodes)):
                   node2=sortednodes[index2]
                   edge2index[(node1,node2)]=(node1,node2)
                   edge2index[(node2,node1)]=(node1,node2)
        elif graphtype=="directed": #why we need edge2index here, just for indexing??
           sortednodes=sorted(allnodes)
           edge2index={}
           for node1 in sortednodes:
               for node2 in sortednodes:
                   if node1!=node2:
                      edge2index[(node1,node2)]=(node1,node2)
                   
        allinfect={}
        for time in sortedtimes:
            myset=set()
            for node in allnodes:
                if not infecttimes.has_key(node):
                   continue
                if time>=infecttimes[node] and time<recovertimes[node]:
                   myset.add(node) 
            allinfect[time]=myset
        newinfect={}    
        newinfect[sortedtimes[0]]=allinfect[sortedtimes[0]]
        for index in range(1,len(sortedtimes)):
            before=sortedtimes[index-1]
            current=sortedtimes[index]
            newinfect[current]=set(allinfect[current].difference(allinfect[before]))
           
        alledges=[]
        allcons=[]
        for index in range(1,len(sortedtimes)):
            before=sortedtimes[index-1]
            current=sortedtimes[index]
            for node1 in newinfect[current]:
                cons=[]
                for node2 in allinfect[before]:
                    assert node1!=node2
                    cons.append(edge2index[(node1,node2)])
                if len(cons)==0: #important
                   continue 
                allcons.append(cons)
                alledges.extend(cons)
        assert len(set(alledges))==len(alledges)
        varnum=len(alledges)
        consnum=len(allcons)
        
        for cons in allcons:
            elem2set.append(cons)
        for node1,node2 in alledges:
            if not cost.has_key((node1,node2)):
               cost[(node1,node2)]=0.0
            if not count.has_key((node1,node2)):
               count[(node1,node2)]=0
            diftime=infecttimes[node2]-infecttimes[node1]
            #diftime=infecttimes[node2]-infecttimes[node1]-1
            if diftime<0:
               diftime *= -1
            assert diftime>0   
            if weightfunction=="expodecay":    
               mycost=1.0/((1.0-spreadprob)**diftime)
            elif weightfunction=="logsum":
               mycost=-1.0*math.log(1.0-spreadprob,2)*diftime 
            elif weightfunction=="difference":
               mycost=1.0-((1.0-spreadprob)**diftime)
            assert mycost>=0   
            cost[(node1,node2)] += mycost
            count[(node1,node2)] += 1
          
    for node1,node2 in cost.keys():
        cost[(node1,node2)] /= float(count[(node1,node2)])
        
    #convert elem2set list into set2elem hashs
    for index in range(0,len(elem2set)):
        elem=elem2set[index]
        for node1,node2 in elem:
            if not set2elem.has_key((node1,node2)):
               set2elem[(node1,node2)]=set()
            set2elem[(node1,node2)].add(index)
            
    #run the algorithm
    uncovered=set(range(0,len(elem2set))) #universal elements
    solution=[]
    while len(uncovered)!=0:
       bestelem=-1
       bestscore=sys.maxint 
       assert len(cost.keys())>0
       for node1,node2 in cost.keys():
           tempcovered=set2elem[(node1,node2)] # we don't explicity check whether length of tempcovered is 0 or whether it has been included in the solution. This is done below by deleting the key from cost hash
           setscore=float(cost[(node1,node2)])/len(tempcovered)
           if setscore <= bestscore:
              bestscore=setscore
              bestelem=(node1,node2)
       solution.append(bestelem)
       del cost[bestelem] 
       deletenext=set2elem[bestelem]
       uncovered -= deletenext
       for node1,node2 in set2elem.keys():
           set2elem[(node1,node2)] -= deletenext
       for node1,node2 in set2elem.keys():    
           if len(set2elem[(node1,node2)])==0 and cost.has_key((node1,node2)): #set2elem wont be deleted even value list has 0 element
              del cost[(node1,node2)] 
    
    retvalues={}
    for node1,node2 in solution:
        retvalues[(node1,node2)]=1.0
    if graphtype=="undirected":    
       for node1,node2 in solution:
           assert (node2,node1) not in solution
    return retvalues




#constraint generator for interdependent traces only for si
def edge_interdependent_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    return 
    

def edge_othermodels_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    s2iprobdistcoef=sentstatechangeparams[0]
    if spreadmodel=="model1":
       i2rprobdistcoef=sentstatechangeparams[0]
       senti2sparam=sentstatechangeparams[1]
    elif spreadmodel=="model2": #?
       i2rprobdistcoef=sentstatechangeparams[0]
       senti2sparam=sentstatechangeparams[1]
    elif spreadmodel=="model3":
       print "algo for spreadmodel {0} not implemented".format(spreadmodel)
       exit(1)
    elif spreadmodel=="model4":
       print "algo for spreadmodel {0} not implemented".format(spreadmodel)
       exit(1)
    elif spreadmodel=="model5": #?
       i2rprobdistcoef=sentstatechangeparams[0]
       senti2sparam=sentstatechangeparams[1]    
    else:
       print "error, unknown model {0}".format(spreadmode)
       exit(1)

    localsubstr=""
    localconsnum=0
    #only consraint will come from s-i interaction part
    if allnodestate[node][pretime]["s"]==0: #everything involving edges will be 0
       return ["",0]
    #if timedif is not 1, we should somehow fill the gap!
    #right now assume i of pretime for all states in between, but this can easily be modified
    tempstr=""
    for temptime in range(curtime,pretime,-1):
        #this part will determine the affect of temptime
        activenode2expo={}
        for neighnode in allnodestate.keys():
            if neighnode==node:
               continue
            if allnodestate[neighnode][pretime]["i"]==0: #exact 0 olmasi lazim!!! dikkat
               continue
            activenode2expo[neighnode]=allnodestate[neighnode][pretime]["i"]
        for neighnode in activenode2expo.keys():    
            basesum=0.0
            added=0
            for infectiontime in estimatedstatechangetime["i"][neighnode].keys():
                if estimatedstatechangetime["i"][neighnode][infectiontime]==0:
                   continue
                if not probdistcoef.has_key(temptime-infectiontime) or probdistcoef[temptime-infectiontime]==0:
                   continue 
                basesum+=(estimatedstatechangetime["i"][neighnode][infectiontime]*probdistcoef[temptime-infectiontime])
                added +=1
            if added==0: #neighnode has no diffusion affect on node
               continue 
            neighnodecoef=math.log(basesum,logbase)
            neighnodecoef*=activenode2expo[neighnode]
            if graphevolution in ["dynamic"]: #we are assuming all times are already in timestepslist!
               varname=edge2var[(neighnode,node,pretime)]
            elif graphevolution in ["static","graphwise","paramwise"]:
               varname=edge2var[(neighnode,node)]
            if neighnodecoef>0:
               tempstr+=" + {0} {1} ".format(neighnodecoef,varname)
            elif neighnodecoef<0:
               tempstr+=" {0} {1} ".format(neighnodecoef,varname)
            else:
               pass
    if tempstr=="":
       return ["",0]
    localconsnum+=1
    plusvarname="pluse{0}{1}".format(con,curconsnum+localconsnum)
    minusvarname="minuse{0}{1}".format(con,curconsnum+localconsnum)
            
    if rightmethod=="rightnoise":
       print "rightnoise log normal noise addition has not been implemented!!"
       print "exiting"
       exit(1)
    elif rightmethod=="epsilon":
       epsilon=rightmethodparam
  
    #for sis model, there will only be single right side(Both left and right will be equal)        
    if spreadmodel in ["sir","seir","si"]:
       firstright=float(allnodestate[node][curtime]["s"])/allnodestate[node][pretime]["s"]
    elif spreadmodel in ["sis"]: #right not works only for exponential i2s distribution(not s2i)!!
       assert senti2sparam[0]=="expo"
       mylambda=float(senti2sparam[1])
       mu=1.0-math.exp(-1.0*mylambda)
       flow1=((1.0-mu)**(curtime-pretime-1))*mu 
       firstright=float(allnodestate[node][curtime]["s"]-flow1)/allnodestate[node][pretime]["s"]
    if spreadmodel in ["sir","seir"]:
       dif=allnodestate[node][curtime]["r"]-allnodestate[node][pretime]["r"]
       recovflow=allnodestate[node][pretime]["i"]-dif
       secondright=allnodestate[node][pretime]["s"]-allnodestate[node][curtime]["i"]+recovflow
    elif spreadmodel in ["si"]:
       secondright=1.0-(float(allnodestate[node][curtime]["i"])/allnodestate[node][pretime]["s"])
    elif spreadmodel in ["sis"]: #same as leftright due to s=1-i equality
       secondright=firstright 
    #first add both equations right sides and then take the log   
    rightsum=(firstright+secondright)/2.0
    if rightsum<=0: #make sure log is defined
       rightsum=epsilon
    right=math.log(rightsum,logbase)
    tempstr += " + {0} - {1} = {2}\n".format(plusvarname,minusvarname,right)
    localsubstr += tempstr
            
    return [localsubstr,localconsnum]   


# lag will also be implemented here??
def edge_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum,sumormul):
    probdistcoef=sentstatechangeparams[0]
    if spreadmodel=="sis":
       senti2sparam=sentstatechangeparams[1]
    localsubstr=""
    localconsnum=0
    #only consraint will come from s-i interaction part
    if allnodestate[node][pretime]["s"]==0: #everything involving edges will be 0
       return ["",0]
    #if timedif is not 1, we should somehow fill the gap!
    #right now assume i of pretime for all states in between, but this can easily be modified
    tempstr=""
    for temptime in range(curtime,pretime,-1):
        #this part will determine the affect of temptime
        activenode2expo={}
        for neighnode in allnodestate.keys():
            if neighnode==node:
               continue
            if allnodestate[neighnode][pretime]["i"]==0: #exact 0 olmasi lazim!!! dikkat
               continue
            activenode2expo[neighnode]=allnodestate[neighnode][pretime]["i"]
        for neighnode in activenode2expo.keys():
            #neighnodecoef assignment part
            if sumormul=="sum":
               basesum=0.0
               added=0
               for infectiontime in estimatedstatechangetime["i"][neighnode].keys():
                   if estimatedstatechangetime["i"][neighnode][infectiontime]==0:
                      continue
                   if not probdistcoef.has_key(temptime-infectiontime) or probdistcoef[temptime-infectiontime]==0:
                      continue 
                   basesum+=(estimatedstatechangetime["i"][neighnode][infectiontime]*probdistcoef[temptime-infectiontime])
                   added +=1
               if added==0: #neighnode has no diffusion affect on node
                  continue 
               neighnodecoef=math.log(basesum,logbase)
               neighnodecoef*=activenode2expo[neighnode]
            elif sumormul=="mul":
               neighnodecoef=0.0
               added=0
               for infectiontime in estimatedstatechangetime["i"][neighnode].keys():
                   if estimatedstatechangetime["i"][neighnode][infectiontime]==0:
                      continue
                   if not probdistcoef.has_key(temptime-infectiontime) or probdistcoef[temptime-infectiontime]==0:
                      continue
                   mypart=math.log(probdistcoef[temptime-infectiontime],logbase)
                   mypart*=estimatedstatechangetime["i"][neighnode][infectiontime]
                   neighnodecoef+=mypart
                   added +=1
               if added==0: #neighnode has no diffusion affect on node
                  continue 
               neighnodecoef*=activenode2expo[neighnode]
            else:
               print "this type of sumormul is unknown!!: {0}".format(sumormul)
               exit(1)
                
            if graphevolution in ["dynamic"]: #we are assuming all times are already in timestepslist!
               varname=edge2var[(neighnode,node,pretime)]
            elif graphevolution in ["static","graphwise","paramwise"]:
               varname=edge2var[(neighnode,node)]
            if neighnodecoef>0:
               tempstr+=" + {0} {1} ".format(neighnodecoef,varname)
            elif neighnodecoef<0:
               tempstr+=" {0} {1} ".format(neighnodecoef,varname)
            else:
               pass
    if tempstr=="":
       return ["",0]
    localconsnum+=1
    plusvarname="pluse{0}{1}".format(con,curconsnum+localconsnum)
    minusvarname="minuse{0}{1}".format(con,curconsnum+localconsnum)
            
    if rightmethod=="rightnoise":
       print "rightnoise log normal noise addition has not been implemented!!"
       print "exiting"
       exit(1)
    elif rightmethod=="epsilon":
       epsilon=rightmethodparam
  
    #for sis model, there will only be single right side(Both left and right will be equal)        
    if spreadmodel in ["sir","seir","si"]:
       firstright=float(allnodestate[node][curtime]["s"])/allnodestate[node][pretime]["s"]
    elif spreadmodel in ["sis"]: #right not works only for exponential i2s distribution(not s2i)!!
       assert senti2sparam[0]=="expo"
       mylambda=float(senti2sparam[1])
       mu=1.0-math.exp(-1.0*mylambda)
       flow1=((1.0-mu)**(curtime-pretime-1))*mu 
       firstright=float(allnodestate[node][curtime]["s"]-flow1)/allnodestate[node][pretime]["s"]
    if spreadmodel in ["sir","seir"]:
       dif=allnodestate[node][curtime]["r"]-allnodestate[node][pretime]["r"]
       recovflow=allnodestate[node][pretime]["i"]-dif
       secondright=allnodestate[node][pretime]["s"]-allnodestate[node][curtime]["i"]+recovflow
    elif spreadmodel in ["si"]:
       secondright=1.0-(float(allnodestate[node][curtime]["i"])/allnodestate[node][pretime]["s"])
    elif spreadmodel in ["sis"]: #same as leftright due to s=1-i equality
       secondright=firstright 
    #first add both equations right sides and then take the log   
    rightsum=(firstright+secondright)/2.0
    if rightsum<=0: #make sure log is defined
       rightsum=epsilon
    right=math.log(rightsum,logbase)
    tempstr += " + {0} - {1} = {2}\n".format(plusvarname,minusvarname,right)
    localsubstr += tempstr
            
    return [localsubstr,localconsnum]


def lag_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    return

def spreadprobunknown_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    return

def difprobunknown_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    return

#impl. both si and sir
def macropartial_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    return

#This can only be solved explicitly when state information is perfect
#we need only distribution type parameter right now
#This can also work for any generic distribution but right now we assume dist is one of those given types
#piecewise models given distribution as piecewise those 4 distributions
def difprobpartial_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,curconsnum):
    if spreadmodel in ["si","sir","seir","sis"]:
       pass
    else:
       print "difprobpartial has not been impleneted for spreadmodel: {}".format(spreadmodel)
       exit(1)
    disttype=sentstatechangeparams[0]
    assert disttype in ["expo","powerlaw","rayleigh","weibull","piecewise"]
   
    if disttype=="weibull":
       shapeparam=sentstatechangeparams[1]
    elif disttype=="piecewise":
       distinfos=sentstatechangeparams[1]
       
    timedistinfo={}
    if disttype in ["expo","powerlaw","rayleigh","weibull"]:
       for temptime in range(0,200):
           if disttype in ["weibull"]:
              timedistinfo[temptime]=disttype,shapeparam
           else:
              timedistinfo[temptime]=disttype   
    elif disttype in ["piecewise"]:
       for distinfo in distinfos:
           if len(distinfo)==2:
              (distname,distinterval)=distinfo
              assert distname in ["weibull"] 
           elif len(distinfo)==3:
              (distname,distparam,distinterval)=distinfo
           else:
              print "wrong length of distinfom{0}".format(len(distinfo))
              exit(1)
           left,right=distinterval
           for temptime in range(left,right): #right is not included
               if distname in ["weibull"]:
                  timedistinfo[temptime]=distname,distparam
               else:
                  timedistinfo[temptime]=distname 
           
    localsubstr=""
    localconsnum=0
    #only consraint will come from s-i interaction part
    if allnodestate[node][pretime]["s"]==0: #everything involving edges will be 0
       return ["",0] 
    #if timedif is not 1, we should somehow fill the gap!
    #right now assume i of pretime for all states in between, but this can easily be modified
    tempstrlist=[]
    for temptime in range(curtime,pretime,-1):
        #this part will determine the affect of temptime
        activenode2expo={}
        for neighnode in allnodestate.keys():
            if neighnode==node:
               continue
            if allnodestate[neighnode][pretime]["i"]==0: #exact 0 olmasi lazim!!! dikkat
               continue
            activenode2expo[neighnode]=allnodestate[neighnode][pretime]["i"]
        for neighnode in activenode2expo.keys():    
            basesum=0.0
            added=0
            for infectiontime in estimatedstatechangetime["i"][neighnode].keys():
                if estimatedstatechangetime["i"][neighnode][infectiontime]==0:
                   continue
                if probdistcoef[temptime-infectiontime]==0 or not probdistcoef.has_key(temptime-infectiontime):
                   continue 
                basesum+=(estimatedstatechangetime["i"][neighnode][infectiontime]*probdistcoef[temptime-infectiontime])
                added +=1
            assert added in [0,1]
            if added==0: #neighnode has no diffusion affect on node
               continue
            neighnodecoef=math.log(basesum,logbase)
            neighnodecoef*=activenode2expo[neighnode]
            if graphevolution in ["dynamic"]: #we are assuming all times are already in timestepslist!
               varname=edge2var[(neighnode,node,pretime)]
            elif graphevolution in ["static","grahwise","paramwise"]:
               varname=edge2var[(neighnode,node)]   
            mystr = " {0} {1} ".format(neighnodecoef,varname)
            tempstrlist.append(mystr)
    if len(tempstrlist)==0:
       return ["",0]
    tempstr="+".join(tempstrlist)
    localconsnum+=1
    plusvarname="pluse{0}{1}".format(con,consnum)
    minusvarname="minuse{0}{1}".format(con,consnum)
            
    if rightmethod=="rightnoise":
       print "rightnoise log normal noise addition has not been implemented!!"
       print "exiting"
       exit(1)
    elif rightmethod=="epsilon":
       epsilon=rightmethodparam
    
    #for sis model, there will only be single right side(Both left and right will be equal)        
    if spreadmodel in ["sir","seir","si"]:
       firstright=float(allnodestate[node][curtime]["s"])/allnodestate[node][pretime]["s"]
    elif spreadmodel in ["sis"]: #right not works only for exponential i2s distribution(not s2i)!!
       assert senti2sparam[0]=="expo"
       mylambda=float(senti2sparam[1])
       mu=1.0-math.exp(-1.0*mylambda)
       flow1=((1.0-mu)**(curtime-pretime-1))*mu 
       firstright=float(allnodestate[node][curtime]["s"]-flow1)/allnodestate[node][pretime]["s"]
    if spreadmodel in ["sir","seir"]:
       dif=allnodestate[node][curtime]["r"]-allnodestate[node][pretime]["r"]
       recovflow=allnodestate[node][pretime]["i"]-dif
       secondright=allnodestate[node][pretime]["s"]-allnodestate[node][curtime]["i"]+recovflow
    elif spreadmodel in ["si"]:
       secondright=1.0-(float(allnodestate[node][curtime]["i"])/allnodestate[node][pretime]["s"])
    elif spreadmodel in ["sis"]: #same as leftright due to s=1-i equality
       secondright=firstright 
    #first add both equations right sides and then take the log   
    rightsum=(firstright+secondright)/2.0
    if rightsum<=0: #make sure log is defined
       rightsum=epsilon
    right=math.log(rightsum,logbase)
    tempstr += " + {0} - {1} = {2}\n".format(plusvarname,minusvarname,right)
    localsubstr += tempstr
            
    print "consnum {0}".format(localconsnum)
    return [localsubstr,localconsnum]


#INFINITYTIME prob is remaininng prob=1.0-summmation
#IMPORTANT: For nodes that are never infected we assume they are infected at infinitytime and assign it to 1!!
#We introduce INFINITYTIME for those
def recoverlseprobassigner(allnodestate,sortedtimes,recoverprobdistcoef):
    INFINITYTIME=1000000
    timeduration=len(sortedtimes)
    retprobs={}
    for node in allnodestate.keys():
        A=np.zeros((timeduration-1,timeduration-1),dtype=np.float)
        b=np.zeros((timeduration-1),dtype=np.float)
        for timeindex in range(1,timeduration):
            curtime=sortedtimes[timeindex]
            pretime=sortedtimes[timeindex-1]
            rightval=allnodestate[node][curtime]["r"]-allnodestate[node][pretime]["r"]
            b[timeindex-1]=rightval
            for timeindex2 in range(1,timeindex+1):
                if recoverprobdistcoef.has_key(timeindex-timeindex2+1):
                   A[timeindex-1,timeindex2-1]=recoverprobdistcoef[timeindex-timeindex2+1]
        optx=scipy.optimize.nnls(A,b)[0] #nonnegative one, better!
        assert len(optx)==timeduration-1
        retprobs[node]={}
        for timeindex in range(0,timeduration-1):
            mytime=sortedtimes[timeindex]
            retprobs[node][mytime]=optx[timeindex]
        optsum=0.0
        for elem in optx:
            optsum+=elem
        retprobs[node][INFINITYTIME]=max(0.0,1.0-optsum)    
        weightedsum=0.0
        for temptime in retprobs[node].keys():
            weightedsum+=retprobs[node][temptime]
        for temptime in retprobs[node].keys():
            retprobs[node][temptime]/=float(weightedsum)    
    return retprobs

def errorprobassigner(allnodestate,sortedtimes):
    INFINITYTIME=1000000
    retprobs={}
    for node in allnodestate.keys():
        retprobs[node]={}
        probhash={}
        posinfected=set()
        for time in sortedtimes:
            if allnodestate[node][time]["i"]>0.0000001:
               posinfected.add(time)
        if len(posinfected)==0:
           for time in sortedtimes:
               retprobs[node][time]=0.0 
           retprobs[node][INFINITYTIME]=1.0
           continue
        for time in sortedtimes:
            if time not in posinfected:
               retprobs[node][time]=0.0  
        foundone="-1"        
        for infecttime in posinfected:
            error=0.0
            for time in sortedtimes:
                if time<infecttime:
                   error+=(allnodestate[node][time]["i"]**2)
                elif time>=infecttime:
                   error+=(1.0-allnodestate[node][time]["i"])**2
            if error==0.0:
               foundone=infecttime
               break
            probhash[infecttime]=1.0/error  
        if foundone!="-1":
           for time in posinfected:
               if time!=foundone:
                  retprobs[node][time]=0.0
               else:    
                  retprobs[node][time]=1.0
           continue
        #error for not rounded case
        notrounderror=0.0
        for time in sortedtimes:
            notrounderror+=allnodestate[node][time]["i"]**2
        probhash[INFINITYTIME]=1.0/notrounderror    
        totprob=0.0
        for value in probhash.values():
            totprob+=value
        for key in probhash.keys():
            probhash[key]/=float(totprob)
        for time in probhash.keys():
            retprobs[node][time]=probhash[time]
    return retprobs

#IMPORTANT: For nodes that are never infected we assume they are infected at infinitytime and assign it to 1!!
#We introduce INFINITYTIME for those
#This method is biased towards assigning infection probability to state. It can not assign score for case of not being infected at any time!!.
def flowdifferenceprobassigner(allnodestate,sortedtimes):
    INFINITYTIME=1000000
    timeduration=len(sortedtimes)
    retprobs={}
    for node in allnodestate.keys():
        retprobs[node]={}
        mytime=sortedtimes[0]
        retprobs[node][mytime]=allnodestate[node][mytime]["i"]
        for timeindex in range(1,timeduration):
            curtime=sortedtimes[timeindex]
            pretime=sortedtimes[timeindex-1]
            rightval=allnodestate[node][pretime]["s"]-allnodestate[node][curtime]["s"]
            if rightval<=0:
               rightval=0.0
            retprobs[node][curtime]=rightval
        weightedsum=0.0
        for temptime in retprobs[node].keys():
            weightedsum+=retprobs[node][temptime]
        if round(weightedsum,6)==0:
           retprobs[node][INFINITYTIME]=1.0
           weightedsum+=retprobs[node][INFINITYTIME]
        for temptime in retprobs[node].keys():
            retprobs[node][temptime]/=float(weightedsum)
    return retprobs


def constgenerate(allnodestates,allnodestatesinfo,graphtype,logbase,graphevolution,spreadmodel,rightmethod,rightmethodparam,timesteplist,infectedprobabilityassignermethod,infertype,partialdegree,sumormul,temporalsparsetype,mapintervals,probtype): 
    con="?"
    
    allnodes=set()
    for allnodestate in allnodestates:
        allnodes |= set(allnodestate.keys())
    allnodes=sorted(list(allnodes))
    edge2var={}
    if graphevolution in ["static","paramwise","graphwise"]:
       if graphtype=="undirected":
          for index1 in range(0,len(allnodes)):
              node1=allnodes[index1]
              for index2 in range(index1+1,len(allnodes)):
                  node2=allnodes[index2]
                  edge2var[(node1,node2)]="x{0}{1}{2}{3}".format(con,node1,con,node2)
                  edge2var[(node2,node1)]="x{0}{1}{2}{3}".format(con,node1,con,node2)
       elif graphtype=="directed":
          for node1 in allnodes:
              for node2 in allnodes:
                  if node1!=node2:
                     edge2var[(node1,node2)]="x{0}{1}{2}{3}".format(con,node1,con,node2)
    elif graphevolution in ["dynamic"]:
       if graphtype=="undirected":
          for index1 in range(0,len(allnodes)):
              node1=allnodes[index1]
              for index2 in range(index1+1,len(allnodes)):
                  node2=allnodes[index2]
                  for time in timesteplist:
                      mappedtime=mapintervals[time]
                      edge2var[(node1,node2,time)]="x{0}{1}{2}{3}{4}{5}".format(con,node1,con,node2,con,mappedtime)
                      edge2var[(node2,node1,time)]="x{0}{1}{2}{3}{4}{5}".format(con,node1,con,node2,con,mappedtime)
       elif graphtype=="directed":
          for node1 in allnodes:
              for node2 in allnodes:
                  if node1!=node2:
                     for time in timesteplist:
                         mappedtime=mapintervals[time]
                         edge2var[(node1,node2,time)]="x{0}{1}{2}{3}{4}{5}".format(con,node1,con,node2,con,mappedtime)
           
    #there are 2 types of constraints(s(t) and i(t)), but we will handle them in the same constraint
    substr="Subject To\n"
    consnum=0
    for traceindex in range(0,len(allnodestates)):
        print "running trace {0}".format(traceindex)
        allnodestate=allnodestates[traceindex]
        if spreadmodel=="si":
           #spreadparam do not need to be spreadprob(depend on infertype)
           spreadparam,s2iparam=allnodestatesinfo[traceindex]
           spreadparam=spreadparam[1]
           s2iparam=s2iparam[1]
        elif spreadmodel=="sis":
           spreadparam,s2iparam,i2sparam=allnodestatesinfo[traceindex]
           spreadparam=spreadparam[1]
           s2iparam=s2iparam[1]
           i2sparam=i2sparam[1]
        elif spreadmodel=="sir":
           spreadparam,s2iparam,i2rparam=allnodestatesinfo[traceindex]
           spreadparam=spreadparam[1]
           s2iparam=s2iparam[1]
           i2rparam=i2rparam[1]
        elif spreadmodel=="seir":
           spreadparam,s2eparam,i2rparam,e2iparam=allnodestatesinfo[traceindex]
           spreadparam=spreadparam[1]
           s2eparam=s2eparam[1]
           i2rparam=i2rparam[1]
           e2iparam=e2iparam[1]

        print "info about distribution"
        print spreadparam
        print s2iparam
        if infertype in ["edge"]:
           if spreadmodel in ["seir"]:
              print "writing!!" 
              print s2eparam
              probdistcoef=tracegeneratorpart.partitiondist(s2eparam[0],s2eparam[1:],probtype,"reverse cumulative",spreadparam)
           elif spreadmodel in ["si","sir","sis"]:    
              probdistcoef=tracegeneratorpart.partitiondist(s2iparam[0],s2iparam[1:],probtype,"reverse cumulative",spreadparam)
           else:
              print "probdistcoef assignment has not been defined for model {0}".format(spreadmodel)
              exit(1)
           if spreadmodel in ["sir","seir"]:
              recoverprobdistcoef=tracegeneratorpart.partitiondist(i2rparam[0],i2rparam[1:],probtype,"normal")
        elif infertype=="spreadprobunknown": #??
           probdistcoef=tracegeneratorpart.partitiondist(s2iparam[0],s2iparam[1:],probtype,"reverse cumulative",spreadparam)
        elif infertype=="difprobpartial": #??
           probdistcoef=tracegeneratorpart.partitiondist(s2iparam[0],s2iparam[1:],probtype,"reverse cumulative",spreadparam)
        elif infertype=="difprobunknown": # ??
           print "cons generation for {0} not done!".format(infertype)
           exit(1)   
        elif infertype=="macropartial": # ??
           print "cons generation for {0} not done!".format(infertype)
           exit(1)   
        elif infertype=="lag": # ??    
           print "cons generation for {0} not done!".format(infertype)
           exit(1)    
        else:
           print "infertype {0} is unknonw".format(infertype)
           exit(1)
         
        randnode=allnodestate.keys()[0]
        sortedtimes=sorted(allnodestate[randnode].keys())
        estimatedstatechangetime={}
        estimatedstatechangetime["i"]={}
        print "partial degree is {0}".format(partialdegree)
        if partialdegree!=1.0:
           #Partial case can also be run for types other than edges especially with mul model(not sum)!!!!
           if infectedprobabilityassignermethod=="flowdifference":
              print "entering flowdifference" 
              assert spreadmodel in ["si"]
              estimatedstatechangetime["i"]=dict(flowdifferenceprobassigner(allnodestate,sortedtimes))
           elif infectedprobabilityassignermethod=="recoverlse":
              print "entering recoverlse" 
              assert spreadmodel in ["sir"]
              estimatedstatechangetime["i"]=dict(recoverlseprobassigner(allnodestate,sortedtimes,recoverprobdistcoef))
           elif infectedprobabilityassignermethod=="recoverlse":
              print "entering error infectedprobability assignment" 
              assert spreadmodel in ["si"]
              estimatedstatechangetime["i"]=dict(errorprobassigner(allnodestate,sortedtimes))
           else:
              print "This infectedprobabilityassignermethod is unknown!!".format(infectedprobabilityassignermethod)
              exit(1)
        else:
           for node in allnodestate.keys():
               estimatedstatechangetime["i"][node]={}
               sortedtimes=sorted(allnodestate[node].keys())
               for time in sortedtimes:
                   if allnodestate[node][time]["i"]==1:
                      estimatedstatechangetime["i"][node][time]=1.0
                      break
        #infected state checking probability
        for tempnode in estimatedstatechangetime["i"].keys():
            probsum=0.0
            myadded=0
            for temptime in estimatedstatechangetime["i"][tempnode].keys():
                probsum += estimatedstatechangetime["i"][tempnode][temptime]
                myadded+=1
            if myadded!=0:    
               assert round(probsum,3)==1.0
           
        for node in allnodestate.keys():    
            for timeindex in range(1,len(sortedtimes)):
                curtime=sortedtimes[timeindex]
                pretime=sortedtimes[timeindex-1]
                if infertype=="edge":
                   if spreadmodel in ["sir","seir","si","sis"]:
                      sentstatechangeparams=[]
                      sentstatechangeparams.append(probdistcoef)
                      if spreadmodel=="sis":
                         sentstatechangeparams.append(i2sparam) 
                      localsubstr,localconsnum=edge_sifamily_constgenerate(allnodestate,node,pretime,curtime,estimatedstatechangetime,graphevolution,logbase,rightmethod,rightmethodparam,edge2var,spreadmodel,sentstatechangeparams,con,consnum,sumormul) #for dynamic graphs, mapinterval info will be retrieved on edge2var
                   elif spreadmodel in ["model1","model2","model5"]:
                      localsubstr,localconsnum=edge_othermodels_constgenerate()
                   elif spreadmodel in ["model3","model4"]:   
                      print "Not impl yet!! {0}".format(infertype)
                      exit(1)
                   elif spreadmodel in ["singlemulti-si","multimulti-si"]:   
                      print "Not impl yet!! {0}".format(infertype)
                      exit(1)
                   else:
                      print "spreadmodel is unknown {0}".format(spreadmodel)
                      exit(1)
                elif infertype=="spreadprobunknown":
                   print "Not impl yet!! {0}".format(infertype)
                   exit(1) 
                elif infertype=="difprobpartial":
                   print "Not impl yet!! {0}".format(infertype)
                   exit(1) 
                elif infertype=="difprobunknown":     
                   print "Not impl yet!! {0}".format(infertype)
                   exit(1) 
                elif infertype=="macropartial":     
                   print "Not impl yet!! {0}".format(infertype)
                   exit(1) 
                elif infertype=="lag":     
                   print "Not impl yet!! {0}".format(infertype)
                   exit(1) 
                else:
                   print "No inference algo for {0}".format(infertype)
                   exit(1)
                substr+=localsubstr
                consnum+=localconsnum

    #temporal sparsity objective function constraints(absolute value constraints)
    temporalconstr=""
    temporalconsnum=0
    if graphevolution in ["static","graphwise","paramwise"]:
       assert temporalsparsetype=="-1"
    elif graphevolution in ["dynamic"]:
       if temporalsparsetype in ["grouplasso","sparsegrouplasso"]: #this will especially be important for dynamic graph inference
          print "group sparseness has not been implemented yet!!"
          exit(1)
       elif temporalsparsetype=="fused":
          sortedtimes=sorted(timesteplist)
          alledges=set()
          if graphtype=="undirected":
             for index1 in range(0,len(allnodes)):
                 node1=allnodes[index1]
                 for index2 in range(index1+1,len(allnodes)):
                     node2=allnodes[index2]
                     alledges.add((node1,node2))
          elif graphtype=="directed":
             for node1 in allnodes:
                 for node2 in allnodes:
                     if node1!=node2:
                        alledges.add((node1,node2))
          for node1,node2 in alledges:
              for timeindex in range(1,len(sortedtimes)):
                  pretime=sortedtimes[timeindex-1]
                  curtime=sortedtimes[timeindex]
                  varname1=edge2var[(node1,node2,curtime)]
                  varname2=edge2var[(node1,node2,pretime)]
                  if varname1==varname2: #no temporal sparsity over same edges over time!!
                     continue 
                  varname3="tempabs{0}".format(edge2var[(node1,node2,curtime)]) #temporal absolute variable for dynamic case
                  cons1="{0} - {1} + {2} >=0 ".format(varname3,varname1,varname2)
                  temporalconstr+="{0}\n".format(cons1)
                  cons2="{0} - {1} + {2} >= 0".format(varname3,varname2,varname1)
                  temporalconstr+="{0}\n".format(cons2)
                  temporalconsnum+=2
       else:
          print "this temporalsparse type is unknonw!! ERROR {0}".format(temporalsparsetype)
          exit(1)
    
    if temporalconstr!="": #dynamic graph case
       consnum+=temporalconsnum 
       substr+="{0}".format(temporalconstr)
    return [substr,consnum,edge2var]       


#THIS IS NOT BEING USED RIGHT NOW???
#Converts trace information to lp usable format 
def preparevariables(allnodestates):
    for nodestates in allnodestates:
       #estimate all possible edges and estimate their possible duration, then assign score based on this possible duration 
       randnode=nodestates.keys()[0]
       alltimes=nodestates[randnode].keys()
       infecttimes={}
       for node in nodestates.keys():
           infecttime=-1
           for time in sorted(alltimes):
               if nodestates[node][time]["i"]==1:
                  infecttime=time
                  break
           infecttimes[node]=infecttime
       recovertimes={}
       for node in nodestates.keys():
           recovertime=-1
           for time in sorted(alltimes):
               if nodestates[node][time]["r"]==1:
                  recovertime=time
                  break
           recovertimes[node]=recovertime    
       allnodes=nodestates.keys()       
       #print "Number of nodes affected by spread is {0}".format(len(allnodes))
       sortedtimes=sorted(alltimes)
       if graphtype=="undirected":
          sortednodes=sorted(allnodes)
          posedges=[]
          for index1 in range(0,len(sortednodes)):
              node1=sortednodes[index1]
              for index2 in range(index1+1,len(sortednodes)):
                  node2=sortednodes[index2]
                  posedges.append((node1,node2))
       
       allinfect={}
       for time in sortedtimes:
           myset=set()
           for node in allnodes:
               if time>=infecttimes[node] and time<recovertimes[node]:
                  myset.add(node) 
           allinfect[time]=myset
       newinfect={}    
       newinfect[sortedtimes[0]]=allinfect[sortedtimes[0]]
       for index in range(1,len(sortedtimes)):
           before=sortedtimes[index-1]
           current=sortedtimes[index]
           newinfect[current]=set(allinfect[current].difference(allinfect[before]))
           
       alledges=[]
       allcons=[]
       for index in range(1,len(sortedtimes)):
           before=sortedtimes[index-1]
           current=sortedtimes[index]
           for node1 in newinfect[current]:
               cons=[]
               for node2 in allinfect[before]:
                   assert node1!=node2
                   if graphtype=="undirected":
                      if (node1,node2) in posedges:
                         cons.append((node1,node2))
                      elif (node2,node1) in posedges:
                         cons.append((node2,node1))
                      else:
                         print "there is no such possibility ERROR"
                         exit(1)
                   elif graphtype=="directed":
                      cons.append((node2,node1))
               allcons.append(cons)
               alledges.extend(cons)
       assert len(set(alledges))==len(alledges)
       varnum=len(alledges)
       consnum=len(allcons)
       
       for cons in allcons:
           elem2set.append(cons)
       for node1,node2 in alledges:
           if not cost.has_key((node1,node2)):
              cost[(node1,node2)]=0.0
           if not count.has_key((node1,node2)):
              count[(node1,node2)]=0
           diftime=infecttimes[node2]-infecttimes[node1]
           if diftime<0:
              diftime *= -1
           assert diftime>0   
           #diftime=infecttimes[node2]-infecttimes[node1]-1
           if weightfunction=="expodecay":    
              mycost=1.0/((1.0-spreadprob)**diftime)
           elif weightfunction=="logsum":
              mycost=-1.0*math.log(1.0-spreadprob,2)*diftime 
           elif weightfunction=="difference":
              mycost=1.0-((1.0-spreadprob)**diftime)
           assert mycost>=0   
           cost[(node1,node2)] += mycost
           count[(node1,node2)] += 1
          
    for node1,node2 in cost.keys():
        cost[(node1,node2)] /= float(count[(node1,node2)])
        
    #convert elem2set list into set2elem hash
    for index in range(0,len(elem2set)):
        elem=elem2set[index]
        for node1,node2 in elem:
            if not set2elem.has_key((node1,node2)):
               set2elem[(node1,node2)]=set()
            set2elem[(node1,node2)].add(index)

 
#All of those algorithms will have same set of constraints but their objective function will be different(Their variable boundaries might also be different)
#epsilon will only be used when the state value we will use is either 1 or 0 but not sth fractional
#Name is lseinference but we can also model absolute error.            
def mainoptinference(allnodestates,allnodestatesinfo,specificrunfolder,graphevolution,spreadmodel,Gname,algo,graphtype,logbase,lambda1,lambda2,fusedlambda,temprightmethod,timesteplist,infectedprobabilityassignermethod,infertype,partialdegree,mapintervals,probtype):
    rightmethod=temprightmethod[0]
    rightmethodparam=temprightmethod[1]
    if rightmethod=="epsilon":
       rightmethodstr="{0}-{1}".format(rightmethod,rightmethodparam)
    else:
       print "this rightmethod is unknown ".format(rightmethod)
       exit(1)

    #there is no such check at macropartial   
    if infertype in ["difprobpartial","spreadprobunknown","difprobunknown","lag"]:
       #perfect state information check
       for nodestate in allnodestates:
           for node in nodestate.keys():
               for time in nodestate[node].keys():
                   for state in nodestate[node][time].keys():
                       if nodestate[node][time][state] not in [0,1]:
                          print "ERROR:Trace info is not perfect but noisy!!" 
                          exit(1)
                          
    if graphevolution in ["static","paramwise","graphwise"]:                
       errortype,sparsetype,sumormul=algo
       temporalsparsetype="-1"
       assert mapintervals=="-1"
    elif graphevolution in ["dynamic"]:                      
       errortype,sparsetype,sumormul,temporalsparsetype=algo   
    print "errortype is {0}".format(errortype)
    print "sparsetype is {0}".format(sparsetype)
    print "sumormul is {0}".format(sumormul)
    con="?"
    consstr,consnum,edge2name=constgenerate(allnodestates,allnodestatesinfo,graphtype,logbase,graphevolution,spreadmodel,rightmethod,rightmethodparam,timesteplist,infectedprobabilityassignermethod,infertype,partialdegree,sumormul,temporalsparsetype,mapintervals,probtype) #timessteplist will be used for dynamic graph inference
    print "constaints are written!!"
    print mapintervals
    print edge2name.keys()[0:10]
    print len(edge2name.keys())
    print consnum
    
    #objective function will have 2 parts. One part minimize the error and other part will force sparsity
    objstr="Minimize\n obj: "
    if errortype=="abse":
       for index in range(0,consnum):
           coef=1.0
           varname="pluse{0}{1}".format(con,index+1)
           objstr += " + {0} {1}".format(coef,varname)
           varname="minuse{0}{1}".format(con,index+1)
           objstr += " + {0} {1}".format(coef,varname)
    elif errortype=="lse":
       objstr += " [ " 
       for index in range(0,consnum):
           coef=2.0
           varname="pluse{0}{1}".format(con,index+1)
           objstr += " + {0} {1} * {2}".format(coef,varname,varname)
           varname="minuse{0}{1}".format(con,index+1)
           objstr += " + {0} {1} * {2}".format(coef,varname,varname)
       objstr += " ] / 2"  
    else:
        print "errortype {0} is unknown!!".format(errortype)
        exit(1)

    if sparsetype in ["l1","l1l2"]:
       for edgeinfokey in edge2name.keys(): #edgeinfokey is either (node1,node2) or (node1,node2,time)
           coef=lambda1
           varname=edge2name[edgeinfokey]
           objstr += " + {0} {1}".format(coef,varname)
    if sparsetype in ["l2","l1l2"]:
       objstr += " [ "  
       for edgeinfokey in edge2name.keys():
           coef=2.0*lambda2
           varname=edge2name[edgeinfokey]
           objstr += " + {0} {1} * {2}".format(coef,varname,varname)
       objstr += " ] / 2"

    #temporal sparsity objective function
    #WARNING!!: Since tempabs is guaranteed to be greater than 0(due to constraints put), we don't add any other boundry constraint below for that   
    if graphevolution in ["dynamic"]:
       if temporalsparsetype in ["grouplasso","sparsegrouplasso"]: #dynamic graph inference case
          print "group sparseness has not been implemented yet!!"
          exit(1)
       elif temporalsparsetype=="fused":
          sortedtimes=sorted(timesteplist)
          alledges=set()
          for edgeinfokey in edge2name.keys():
              node1,node2,time=edgeinfokey
              alledges.add((node1,node2))
          uniquetimes=set(mapintervals.values()) #all real times on dynamic graph
          for node1,node2 in alledges:
              for timeindex in range(1,len(sortedtimes)):
                  curtime=sortedtimes[timeindex]
                  if curtime not in uniquetimes:
                     continue 
                  curvarname="tempabs{0}".format(edge2name[(node1,node2,curtime)]) #temporal absolute variable for dynamic case
                  objstr += " + {0} {1} ".format(fusedlambda,curvarname) #trial purpose 0.01 
       else:
          print "this temporalsparse type is unknonw!! for objective function ERROR {0}".format(temporalsparsetype)
          exit(1)
          
    #define boundary string
    boundstr="Bounds\n"
    for edgeinfokey in edge2name.keys():
        varname=edge2name[edgeinfokey]
        boundstr+="0 <= {0} <= 1\n".format(varname)         
    for index in range(0,consnum):
        varname="pluse{0}{1}".format(con,index+1)
        boundstr+="0 <= {0} \n".format(varname)
        varname="minuse{0}{1}".format(con,index+1)
        boundstr+="0 <= {0} \n".format(varname)

    #write cplex script file
    outlpfile="{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.lp".format(Gname,graphtype,errortype,sparsetype,logbase,lambda1,lambda2,rightmethodstr)
    outlppath="{0}/{1}".format(specificrunfolder,outlpfile)
    file=open(outlppath,"w")
    file.write(objstr+"\n")
    file.write(consstr)
    file.write(boundstr+"\n")
    file.write("End\n")
    file.close()
    cplexoutfile="{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.lpout".format(Gname,graphtype,errortype,sparsetype,logbase,lambda1,lambda2,rightmethodstr)
    cplexscriptfile="{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.script".format(Gname,graphtype,errortype,sparsetype,logbase,lambda1,lambda2,rightmethodstr)
    cplexoutpath="{0}/{1}".format(specificrunfolder,cplexoutfile)
    cplexscriptpath="{0}/{1}".format(specificrunfolder,cplexscriptfile)
    processornum=4
    file=open(cplexscriptpath,"w")
    file.write("read {0}\n".format(outlppath))
    file.write("set threads {0}\n".format(processornum))
    file.write("optimize\n") 
    file.write("display solution objective\n")
    file.write("display solution variables -\n")  
    file.close()
            
    #run cplex script file 
    code="cplex < {0} > {1}".format(cplexscriptpath,cplexoutpath)
    os.system(code)

    #make temp text file
    temptextfilename="mytext.txt"
    temptextfilepath="{0}/{1}".format(specificrunfolder,temptextfilename)
    file=open(temptextfilepath,"w")
    file.write("done\n")
    file.close()
    
    #read cplex output and return values
    retvalues=readcplexoutput(cplexoutpath,specific=["x"])
    retvalues2={}
    if graphevolution in ["dynamic"]:
       for varname in retvalues.keys():
           temp,node1,node2,time=varname.replace("x","").split(con)
           node1=int(node1)
           node2=int(node2)
           time=int(time)
           retvalues2[(node1,node2,time)]=retvalues[varname]
    elif graphevolution in ["static","graphwise","paramwise"]:
       for varname in retvalues.keys():
           temp,node1,node2=varname.replace("x","").split(con)
           node1=int(node1)
           node2=int(node2)
           retvalues2[(node1,node2)]=retvalues[varname]
    print "done lse inference part"
    #del codes, don't delete cplex output
    delcode="rm -rf {0}".format(cplexscriptpath)
    os.system(delcode)
    delcode="rm -rf {0}".format(outlppath)
    os.system(delcode)
    #delcode="rm -rf {0}".format(cplexoutpath)
    #os.system(delcode)
    return retvalues2



def readmatlaboutput(outfilepath,varsymbol):
    myflag=0
    file=open(outfilepath,"r")
    values=[]
    for line in file:
        line=line.rstrip()
        if myflag==0 and line.find("{0} =".format(varsymbol))!=-1:
           myflag=1
           continue
        elif myflag==1:
           myflag=2
           continue
        elif myflag==2:
           print line 
           replaced=line.replace(" ","")
           if replaced!="":
              val=float(line.replace(" ",""))
              values.append(val)
           else:
              break 
    file.close()    
    return values   

def gpsr_runner(datafilepath,lambda1,matlabpath,outfilepath,codepath="GPSR_Basic"):
    code="{0} -r \"{1}({2},{3}) ; quit ;\" > {4}".format(matlabpath,codepath,datafilepath,lambda1,outfilepath)
    os.system(code)
    values=readmatlaboutput(outfilepath,"x")
    return values
def nuirls_runner(datafilepath,matlabpath,outfilepath,codepath="NUIRLS"): #not working right now
    p=1
    nmfiter=10
    irlsiter=10
    showflag=0
    code="{0} -r \"{1}({2},{3},{4},{5},{6},{7}) ; quit ; \" > {8}".format(matlabpath,codepath,y,A,p,nmfiter,irlsiter,showflag,outfilepath)
    os.system(code)
    values=readmatlaboutput(outfilepath,"X")
    return values
def l1ls_runner(datafilepath,lambda1,matlabpath,outfilepath,codepath="l1_ls"):
    code="{0} -r \"{1}({2},{3}) ; quit ;\" > {4}".format(matlabpath,codepath,matlabpath,lambda1,outfilepath)
    os.system(code)
    values=readmatlaboutput(outfilepath,"x")
    return values
def nonl1ls_runner(datafilepath,lambda1,matlabpath,outfilepath,codepath="l1_ls_nonneg"):
    code="{0} -r \"{1}({2},{3}) ; quit ;\" > {4}".format(matlabpath,codepath,matlabpath,lambda1,outfilepath)
    os.system(code)
    values=readmatlaboutput(outfilepath,"x")
    return values
def elasticnet_runner(datafilepath,lambda1,matlabpath,outfilepath,codepath="glmnet"):
    code="{0} -r \"".format(matlabpath)
    algo="{0}".format(codepath)
    code+=" fit1={0}({1}) ; ".format(algo,datafilepath)
    algo="{0}Coef".format(codepath)
    code+=" {0}(fit1,{1}) ; ".format(algo,lambda1) #% extract coefficients at a single value of lambda
    code+=" quit ;\" > {0}".format(outfilepath)
    os.system(code)

#def arrayOfArrayToMatlabString(array):
#    return '[' + "\n ".join(" ".join("%6g" % val for val in line) for line in array) + ']'

#compressed sensing is here to come up with solutions with least samples. When the system is undetermined, however this might now work since we don't know whether matrix satisfy isometry property ??
#asyn done????
#NOT BEING USED RIGHT NOW, NOT SYCHRNOZIED WITH LAST IMPLEMENTATION
def sir_static_csalgo(allnodestates,allnodestatesinfo,specificrunfolder,Gname,algo,graphtype,spreadtimemode,lagdepth,infinitylag,epsilon,logbase,lambda1=1.0):    
    A,y,pos2edge=sir_static_cons2matrices(allnodestates,allnodestatesinfo,graphtype,spreadtimemode,lagdepth,infinitylag,epsilon,logbase) #directed undirected handled here
    print np.shape(A)
    print "constraints generated!! {0} {1}".format(np.shape(A)[0],np.shape(A)[1])
    assert np.shape(A)[1]==len(pos2edge.keys())
  
    #Astr=""
    #for index1 in range(0,np.shape(A)[0]):
    #    for index2 in range(0,np.shape(A)[1]):
    #        Astr+=" {0} ".format(A[index1,index2])
    #    Astr+=";"
    #Astr=Astr[0:-1]
    #Astr ="[{0}]".format(allstr)
    #ystr=""
    #for elem in y:
    #    ystr+="{0};".format(elem)
    #ystr=ystr[0:-1]
    #ystr="[{0}]".format(ystr)
    
    #alllist=[]
    #for index1 in range(0,np.shape(A)[0]):
    #    rowstr=""
    #    for index2 in range(0,np.shape(A)[1]):
    #        templist.append(str(A[index1,index2]))
    #        rowstr+=" {0} ".format(A[index1,index2]))
    #    rowstr=" ".join(templist)
    #    rowstr+=";"
    #    alllist.append(rowstr)
    #allstr=";".join(alllist)    
    #Astr ="[{0}]".format(allstr)
    #ylist=[]
    #for elem in y:
    #    ylist.append(str(elem))
    #ystr=";".join(ylist)
    #ystr="[{0}]".format(ystr)
    
    matlabpath="/opt/stow/matlab-r2012a/bin/matlab"
    #matlabpath="/Applications/MATLAB_R2012a.app/bin/matlab"
    #matlabpath="/usr/local/bin/matlab"

    #"myMatrix = " + arrayOfArrayToMatlabString(array)
    
    if not os.path.exists(specificrunfolder):
       os.makedirs(specificrunfolder)
    datafile="{0}_{1}_{2}_{3}_{4}_{5}.mat".format(Gname,algo,graphtype,epsilon,logbase,lambda1)
    datafilepath="{0}/{1}".format(specificrunfolder,datafile)
    data={}
    data['A']=A
    data['y']=y
    scipy.io.savemat(datafilepath,data)
    outfile="{0}_{1}_{2}_{3}_{4}_{5}.csout".format(Gname,algo,graphtype,epsilon,logbase,lambda1)
    outfilepath="{0}/{1}".format(specificrunfolder,outfile)
    if algo=="elasticnet":
       values=elasticnet_runner(datafilepath,lambda1,matlabpath,outfilepath) 
    elif algo=="gpsr":
       values=gpsr_runner(datafilepath,lambda1,matlabpath,outfilepath) 
    elif algo=="l1ls":
       values=l1ls_runner(datafilepath,lambda1,matlabpath,outfilepath) 
    elif algo=="nonl1ls":
       values=nonl1ls_runner(datafilepath,lambda1,matlabpath,outfilepath) 
    elif algo=="nuirls": #BE CAREFUL: This won't work when the system is undetermined. there is no lambda either, so screw this
       values=nuirls_runner(datafilepath,matlabpath,outfilepath) 
    else:
        print "error cs algo {0}".format(csalgo)
        exit(1)

    code="rm -rf {0}".format(datafilepath)
    os.system(code)
    
    retvalues={}
    for pos in range(0,len(values)):
        retvalues[pos2edge[pos]]=values[pos]
    return retvalues


#Converts constraints obtained from allnodestates information to matrix form which will be used by CS algorithms
#NOT USED RIGHT NOW, NOT SYCNHRONIZED WITH LAST IMPLEMENTATION!!
def sir_static_cons2matrices(allnodestates,allnodestatesinfo,graphtype,spreadtimemode,lagdepth,infinitylag,epsilon,logbase): 
    allnodes=set()
    for allnodestate in allnodestates:
        tempnodes=allnodestate.keys()
        allnodes=allnodes.union(tempnodes)
    sortedallnodes=sorted(allnodes)    
    edge2pos={} #maps each edge to its position on matrix
    if graphtype=="undirected":
       pos=0 
       for index1 in range(0,len(sortedallnodes)):
           node1=sortedallnodes[index1]
           for index2 in range(index1+1,len(sortedallnodes)):
               node2=sortedallnodes[index2]
               for lag in range(1,lagdepth+1):
                   edge2pos[(node1,node2,lag)]=pos
                   edge2pos[(node2,node1,lag)]=pos
                   pos += 1
    elif graphtype=="directed":
       pos=0 
       for node1 in sortedallnodes:
           for node2 in sortedallnodes:
               if node1!=node2:
                  for lag in range(1,lagdepth+1):
                      edge2pos[(node1,node2,lag)]=pos
                      pos += 1
    edgenum=pos
    
    print "starting!"
    A=[]
    y=[]
    seenposmap={}
    #there are 2 types of constraints(s(t) and i(t)), Handle them in a single constraint since they are both the same constraints with different values
    for traceindex in range(0,len(allnodestates)):
        print "matrix filling for trace {0}".format(traceindex)
        allnodestate=allnodestates[traceindex]
        spreadprob,recoverdist,recoverparam=allnodestatesinfo[traceindex]
        randnode=allnodestate.keys()[0]
        sortedtimes=sorted(allnodestate[randnode].keys())
        for timeindex in range(1,len(sortedtimes)):
            time=sortedtimes[timeindex]
            pretime=sortedtimes[timeindex-1]
            for node in allnodestate.keys():
                posslags=[] #all possible lags
                for lag in range(1,lagdepth+1):
                    if allnodestate[node].has_key(time-lag):
                       if allnodestate[node][time-lag]["s"]==0: #exact 0 olmasi lazim
                          continue
                       posslags.append(lag)
                if len(posslags)==0:
                   continue
                
                row=[0.0]*edgenum
                scoefs=1.0
                sthdoneflag=False
                for lag in posslags:    
                    varnodes=[]
                    varcoefs=[]
                    for neighnode in allnodestate.keys():
                        if neighnode==node:
                           continue
                        if allnodestate[neighnode][time-lag]["i"]==0: #exact 0 olmasi lazim!!! dikkat
                           continue 
                        varnodes.append(neighnode)
                        mycoef=allnodestate[neighnode][time-lag]["i"]
                        varcoefs.append(mycoef)
                    if len(varcoefs)==0:
                       continue
                    sthdoneflag=True
                    scoefs*=allnodestate[node][time-lag]["s"]  
                    for neighindex in range(0,len(varnodes)):
                        neighnode=varnodes[neighindex]
                        neighcoef=varcoefs[neighindex]
                        pos=edge2pos[(neighnode,node,lag)]
                        if graphtype=="undirected":
                           node1,node2=sorted([neighnode,node])
                           if not seenposmap.has_key((node1,node2)):
                              seenposmap[(node1,node2)]=set()
                           seenposmap[(node1,node2)].add(lag) 
                        elif graphtype=="directed":
                           if not seenposmap.has_key((neighnode,node)):
                              seenposmap[(neighnode,node)]=set()
                           seenposmap[(neighnode,node)].add(lag) 
                        row[pos]=neighcoef
                if not sthdoneflag:
                   continue

                #rightval must be between 0 and 1. If greater than 1, we will get negative values after taking the log. closest will be 0 so making right 0 makes sense.
                #but do not make 0 quickly add them up and then make 0 if necessary
                firstright=allnodestate[node][time]["s"]/scoefs
                if firstright<=0:
                   firstright=epsilon
                firstright=math.log(firstright,1.0-spreadprob)
                #will secondright depend on probability of state allnodestate[node][pretime]["s"]???
                dif=allnodestate[node][time]["r"]-allnodestate[node][pretime]["r"]
                recovflow=allnodestate[node][pretime]["i"]-dif
                secondright=allnodestate[node][pretime]["s"]-allnodestate[node][time]["i"]+recovflow
                #here what if secondright<0, log won't be defined so, make it epsilon(closest possible)
                #make sure seconddright is between above epsilon
                #being 1 causes negative value after log and this will be handled after adding both values
                if secondright<=0:
                   secondright=epsilon
                secondright=math.log(secondright,1.0-spreadprob)
                #for both first and secondright, <0 means it is less possible
                #assert firstright>=0
                #assert secondright>=0
                right=(firstright+secondright)/2.0
                #if right<0: #why make right 0, keep like this!! 
                #   right=0
                A.append(row)
                y.append(right)   

    #column deletion way 1
    delindices1=set(range(0,edgenum))
    for node1,node2 in seenposmap.keys():
        for lag in seenposmap[(node1,node2)]:
            pos=edge2pos[(node1,node2,lag)]
            delindices1.remove(pos)
    
    #column deletion way 2                
    print "column deletion part"
    delindices=set()
    consnum=len(A)
    varnum=len(A[0])
    for varindex in range(0,varnum):
        flag=True
        for consindex in range(0,consnum):
            if A[consindex][varindex]!=0.0:
               flag=False
               break
        if flag:
           delindices.add(varindex)

    assert len(delindices1.difference(delindices))==0

    tempmat=np.array(A)
    print np.shape(tempmat)
    sortededges=[(-1,-1,-1)]*edgenum
    for node1,node2,lag in edge2pos.keys():
        pos=edge2pos[(node1,node2,lag)]
        sortededges[pos]=(node1,node2,lag)
    for offset,delindex in enumerate(sorted(delindices)):
        delindex -= offset
        del sortededges[delindex]
        #for myindex in range(0,len(A)):
        #    del A[myindex][delindex]
    tempmat=sp.delete(tempmat,list(delindices),1)
    print np.shape(tempmat)
    print "number of remaining edges left {0}".format(len(seenposmap.keys()))
    myzeros=np.zeros((np.shape(tempmat)[0], len(seenposmap.keys())), dtype=tempmat.dtype)
    tempmat=np.concatenate((tempmat,myzeros), axis=1)
    del myzeros
    A=tempmat.tolist()
    print len(delindices)
    print len(A),len(A[0])
    
    tempedge2pos={}
    for index in range(0,len(sortededges)):
        node1,node2,lag=sortededges[index]
        if graphtype=="undirected":
           tempedge2pos[(node1,node2,lag)]=index
           tempedge2pos[(node2,node1,lag)]=index
        elif graphtype=="directed":
           tempedge2pos[(node1,node2,lag)]=index
    edge2pos=dict(tempedge2pos)
    
    #add infinitylag for seen edges
    if graphtype=="undirected":
       pos=len(edge2pos.keys())/2
       pos2=edgenum-len(delindices)
       print pos
       print pos2
       assert pos==pos2
    elif graphtype=="directed":
       pos=len(edge2pos.keys())
       pos2=edgenum-len(delindices)
       print pos
       print pos2
       assert pos==pos2
    for node1,node2 in seenposmap.keys():
        if graphtype=="undirected":
           edge2pos[(node1,node2,infinitylag)]=pos
           edge2pos[(node2,node1,infinitylag)]=pos
        elif graphtype=="directed":
           edge2pos[(node1,node2,infinitylag)]=pos
        pos += 1

    #This constraint is summation constraint for asyn case and also for syn case
    print "cons writing part!!"
    print "remaining var count {0}".format(len(seenposmap.keys()))
    edgenum=pos
    print edgenum
    print "len before adding cons!!"
    print len(A),len(A[0])
    assert edgenum==len(A[0])
    for node1,node2 in seenposmap.keys():
        row=[0.0]*edgenum
        for lag in seenposmap[(node1,node2)]:
            pos=edge2pos[(node1,node2,lag)] 
            row[pos]=1.0
        pos=edge2pos[(node1,node2,infinitylag)] 
        row[pos]=1.0    
        A.append(row)
        y.append(1.0)
         
    Amat=np.array(A)
    print np.shape(Amat)
    pos2edge={}
    for node1,node2,lag in edge2pos.keys():
        pos=edge2pos[(node1,node2,lag)]
        pos2edge[pos]=(node1,node2,lag)
    return [Amat,y,pos2edge]

 
def readcplexoutput(outfile,specific=[]):
    retvalues={}
    varflag=False
    file=open(outfile,"r")
    for line in file:
        line=line.rstrip()
        if not varflag and line.find("CPLEX> Variable Name")!=-1 and line.find("Solution Value")!=-1:
           varflag=True
           continue
        if varflag:
           for varname in specific: 
               if line.startswith(varname):
                  key,value=line.split()
                  value=float(value)
                  retvalues[key]=value
                  break
    file.close()
    return retvalues



#We will be using graoh6 or sparse6 for storing undirected graphs when node names are numbers(they don't have labels). When there are labels on the nodes and directed graphs, we will use gml file format and then compress it by gzip since networkx can also read compressed graph files with extension .gz
import networkx as nx
import os
import math
import numpy as np
import itertools
import random
import sys
from optparse import OptionParser
from collections import deque
import myutilities as myutil


def pbsrunner(code2run,pbsfolder,pbsfilename):
    pbsfilepath="{0}/{1}.pbs".format(pbsfolder,pbsfilename)
    errorpath="{0}/{1}.err".format(pbsfolder,pbsfilename)
    outpath="{0}/{1}.out".format(pbsfolder,pbsfilename)
    print pbsfilepath
    print errorpath
    print outpath
    file=open(pbsfilepath,"w")
    file.write("#PBS -q pool1\n")
    file.write("#PBS -r n\n")
    file.write("#PBS -V\n")
    file.write("#PBS -o {0}\n".format(outpath))
    file.write("#PBS -e {0}\n".format(errorpath))
    file.write("cd $PBS_O_WORKDIR\n")
    file.write(code2run+"\n")    
    file.close()
    code="qsub {0}".format(pbsfilepath)
    os.system(code)


#kronecker and directed_ff will be tough to use here!!
def algoparamassigner(n,m):
    algosparameter={}
    
    start=0.1
    end=1.0
    inc=0.2
    p4list=myutil.myrange(start,end,inc)
    p3list=myutil.myrange(start,end,inc)
    p2list=myutil.myrange(start,end,inc)
    p1list=myutil.myrange(start,end,inc)
    algosparameter["kronecker"]={}
    algosparameter["kronecker"]["p4"]=p4list
    algosparameter["kronecker"]["p3"]=p3list
    algosparameter["kronecker"]["p2"]=p2list
    algosparameter["kronecker"]["p1"]=p1list
    algosparameter["kronecker"]["n"]=[n]
    
    dmcstart=0.2
    dmcinc=0.2
    dmcend=1.0
    qcon1=myutil.myrange(dmcstart,dmcend,dmcinc)    
    dmcstart=0.01
    dmcinc=0.01
    dmcend=1.00
    qdel1=myutil.myrange(dmcstart,dmcend,dmcinc)
    #dmcstart=0.80
    #dmcinc=0.02
    #dmcend=1.0
    #qdel2=myutil.myrange(dmcstart,dmcend,dmcinc)
    #qdel1.extend(qdel2)
    algosparameter["dmc"]={}
    algosparameter["dmc"]["qcon"]=qcon1
    algosparameter["dmc"]["qdel"]=qdel1
    algosparameter["dmc"]["n"]=[n]

    smwstart=0.1
    smwend=1.0
    smwinc=0.1  
    algosparameter["smw"]={}
    algosparameter["smw"]["qrewire"]=myutil.myrange(smwstart,smwend,smwinc)
    algosparameter["smw"]["n"]=[n]
    algosparameter["smw"]["m"]=[m]

    algosparameter["rds"]={}
    algosparameter["rds"]["n"]=[n]
    algosparameter["rds"]["m"]=[m]
     
    algosparameter["lpa"]={}
    algosparameter["lpa"]["n"]=[n]
    algosparameter["lpa"]["m"]=[m] 

    algosparameter["undirected_ff"]={}
    algosparameter["undirected_ff"]["n"]=[n]
    ffstart=0.05
    ffend=0.36
    ffinc=0.05
    temp1=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.35
    ffend=0.61
    ffinc=0.01
    temp2=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.61
    ffend=1.0
    ffinc=0.02
    temp3=myutil.myrange(ffstart,ffend,ffinc)
    temp1.extend(temp2)
    temp1.extend(temp3)
    algosparameter["undirected_ff"]["p"]=temp1

    algosparameter["directed_ff"]={}
    algosparameter["directed_ff"]["n"]=[n] 
    ffstart=0.05
    ffend=0.36
    ffinc=0.05
    temp1=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.35
    ffend=0.62
    ffinc=0.02
    temp2=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.62
    ffend=1.0
    ffinc=0.05
    temp3=myutil.myrange(ffstart,ffend,ffinc)
    temp1.extend(temp2)
    temp1.extend(temp3)
    algosparameter["directed_ff"]["pforward"]=temp1
    ffstart=0.05
    ffend=0.36
    ffinc=0.05
    temp1=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.35
    ffend=0.62
    ffinc=0.02
    temp2=myutil.myrange(ffstart,ffend,ffinc)
    ffstart=0.62
    ffend=1.0
    ffinc=0.05
    temp3=myutil.myrange(ffstart,ffend,ffinc)
    temp1.extend(temp2)
    temp1.extend(temp3)
    algosparameter["directed_ff"]["pbackward"]=temp1

    geostart=0.02
    geoend=1.0
    geoinc=0.02
    algosparameter["geo"]={}
    algosparameter["geo"]["n"]=[n]
    algosparameter["geo"]["radius"]=myutil.myrange(geostart,geoend,geoinc)
    
    dmcstart=0.2
    dmcinc=0.2
    dmcend=1.0
    qnew1=myutil.myrange(dmcstart,dmcend,dmcinc)
    dmcstart=0.005
    dmcinc=0.005
    dmcend=1.00
    qdel1=myutil.myrange(dmcstart,dmcend,dmcinc)
    #dmcstart=0.90
    #dmcinc=0.02
    #dmcend=1.0
    #qdel2=myutil.myrange(dmcstart,dmcend,dmcinc)
    #qdel1.extend(qdel2)
    algosparameter["dmr"]={}
    algosparameter["dmr"]["qnew"]=qnew1
    algosparameter["dmr"]["qdel"]=qdel1
    algosparameter["dmr"]["n"]=[n]
    
    #algosparameter["agv"]={}
    #algosparameter["agv"]["rho"]=np.arange(0,1,0.1)
    #algosparameter["agv"]["alpha"]=np.arange(0,1,0.1)
    
    #algosparameter["bpa"]={}
    #algosparameter["bpa"]["n"]=[n]
    #algosparameter["bpa"]["m"]=[m]

    plcstart=0.02
    plcinc=0.02
    plcend=1.0
    algosparameter["plc"]={}
    algosparameter["plc"]["p"]=myutil.myrange(plcstart,plcend,plcinc)
    algosparameter["plc"]["n"]=[n]
    algosparameter["plc"]["m"]=[m]
    
    return algosparameter

#algo parameter assigner for paramwise case!!
#!!For RTG model, parameter will be edgenum rather than nodenum
def paramwisealgoparamassigner(nodenum):
    algosparameters={}
     
    paramdict1={}
    paramdict1["edgenum"]=nodenum
    paramdict1["beta"]=0.5
    paramdict1["k"]=4
    paramdict1["q"]=0.2
    paramdict2={}
    paramdict2["edgenum"]=nodenum
    paramdict2["beta"]=0.5
    paramdict2["k"]=4
    paramdict2["q"]=0.3
    paramdict3={}
    paramdict3["edgenum"]=nodenum
    paramdict3["beta"]=0.2
    paramdict3["k"]=5
    paramdict3["q"]=0.2
    paramdict4={}
    paramdict4["edgenum"]=nodenum
    paramdict4["beta"]=0.4
    paramdict4["k"]=5
    paramdict4["q"]=0.3
    algosparameters["rtg"]=[paramdict1,paramdict2,paramdict3,paramdict4]
    
    paramdict1={}
    paramdict1["n"]=nodenum
    paramdict1["qcon"]=0.5
    paramdict1["qmod"]=0.6
    paramdict2={}
    paramdict2["n"]=nodenum
    paramdict2["qcon"]=0.5
    paramdict2["qmod"]=0.8
    paramdict3={}
    paramdict3["n"]=nodenum
    paramdict3["qcon"]=0.3
    paramdict3["qmod"]=0.6
    paramdict4={}
    paramdict4["n"]=nodenum
    paramdict4["qcon"]=0.4
    paramdict4["qmod"]=0.8
    algosparameters["dmc"]=[paramdict1,paramdict2,paramdict3,paramdict4]

    paramdict1={}
    k1=nodenum/10
    paramdict1["m"]=k1*(nodenum-k1)
    paramdict1["n"]=nodenum
    paramdict2={}
    k2=nodenum/50
    paramdict2["m"]=k2*(nodenum-k2)
    paramdict2["n"]=nodenum
    algosparameters["lpa"]=[paramdict1,paramdict2]

    paramdict1={}
    p1=0.05
    paramdict1["m"]=p1*nodenum*(nodenum-1)/2
    paramdict1["n"]=nodenum
    paramdict2={}
    p2=0.1
    paramdict2["m"]=p2*nodenum*(nodenum-1)/2
    paramdict2["n"]=nodenum
    algosparameters["rds"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["m"]=nodenum*10
    paramdict1["n"]=nodenum
    paramdict1["qrewire"]=0.2
    paramdict2={}
    paramdict2["m"]=nodenum*20
    paramdict2["n"]=nodenum
    paramdict2["qrewire"]=0.3
    algosparameters["smw"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["n"]=nodenum
    paramdict1["p"]=0.2
    paramdict2={}
    paramdict2["n"]=nodenum
    paramdict2["p"]=0.4
    algosparameters["undirected_ff"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["n"]=nodenum
    paramdict1["pbackward"]=0.1
    paramdict1["pforward"]=0.3
    paramdict2={}
    paramdict2["n"]=nodenum
    paramdict2["pbackward"]=0.2
    paramdict2["pforward"]=0.4
    algosparameters["directed_ff"]=[paramdict1,paramdict2]

    paramdict1={}
    paramdict1["n"]=nodenum
    paramdict1["p1"]=0.4
    paramdict1["p2"]=0.8
    paramdict1["p3"]=0.8
    paramdict1["p4"]=0.4
    paramdict2={}
    paramdict2["n"]=nodenum
    paramdict2["p1"]=0.5
    paramdict2["p2"]=0.7
    paramdict2["p3"]=0.7
    paramdict2["p4"]=0.4
    algosparameters["kronecker"]=[paramdict1,paramdict2]
    
    return algosparameters

#only nodenumbers will be fixed.
def data_generate_paramwise(n,samplecount,underfoldername):
    algosparameter=paramwisealgoparamassigner(n)
    for algo in algos:
        if algo in directedalgos:
           algotype="directed"
        else:
           algotype="undirected"
           
        algodirname="{0}/{1}".format(underfoldername,algo)
        if not os.path.exists(algodirname):
           os.makedirs(algodirname)

        params=algosparameter[algo]
        paramnames=sorted(algosparameter[algo][0].keys())

        for param in params:
            strparamlist=[]
            for index in range(0,len(paramnames)):
                paramname=paramnames[index]
                elem="{0}{1}".format(paramname,param[paramname])
                strparamlist.append(elem)
            strsent="_".join(strparamlist)

            strparam=[str(param[paramname]) for paramname in paramnames]
            paramstr=" ".join(strparam)

            code='python graphgeneratorpart.py {0} {1} {2} {3} {4} {5} {6} {7}'.format(strsent,algodirname,samplecount,algo,algotype,-1,n,paramstr)
            pbsfilename="{0}-{1}-{2}-{3}-{4}-{5}-{6}".format(strsent,samplecount,algo,algotype,-1,n,paramstr)
            pbsfilename=pbsfilename.replace(" ","")
            pbsrunner(code,pbsfolder,pbsfilename)
            #code="runCmd -c \""+code+"\" --nowait"
            #os.system(code)
               

#where node and edge numbers are fixed!Same density graphs
def data_generate(n,m,samplecount,underfoldername):
    algosparameter=algoparamassigner(n,m)
    for algo in algos:
        if algo in directedalgos:
           algotype="directed"
           edgecount=2*int(m)
        else:
           algotype="undirected"
           edgecount=int(m)
           
        algodirname="{0}/{1}".format(underfoldername,algo)
        if not os.path.exists(algodirname):
           os.makedirs(algodirname)
        
        paramnames=sorted(algosparameter[algo].keys())
        templist=[]
        for key in paramnames:
            templist.append(algosparameter[algo][key])
        params=list(itertools.product(*templist))

        if algo not in ["undirected_ff"]: #algorithms other than ffs
           print algo 
           for param in params:
               print param
               strparamlist=[]
               for index in range(0,len(paramnames)):
                   elem="{0}{1}".format(paramnames[index],param[index])
                   strparamlist.append(elem)
               strsent="_".join(strparamlist)
               
               strparam=[str(elem) for elem in param]
               paramstr=" ".join(strparam)

               code='python graphgeneratorpart.py {0} {1} {2} {3} {4} {5} {6} {7}'.format(strsent,algodirname,samplecount,algo,algotype,edgecount,n,paramstr)
               pbsfilename="{0}-{1}-{2}-{3}-{4}-{5}-{6}".format(strsent,samplecount,algo,algotype,edgecount,n,paramstr)
               pbsfilename=pbsfilename.replace(" ","")
               pbsrunner(code,pbsfolder,pbsfilename)
               #code="runCmd -c \""+code+"\" --nowait"
               #os.system(code)
               
        else: #We will execute serially and once a parameter satisfy the graph. After satistication, if next 5 samples do not create a graph with satisfying properties, execution stops!
           myqu=deque() 
           myqu.append(True)
           myqu.append(True)
           myqu.append(True)
           #myqu.append(True)
           #myqu.append(True)            
           
           fastffflag=False 
           for param in params:
               strparamlist=[]
               for index in range(0,len(paramnames)):
                   elem="{0}{1}".format(paramnames[index],param[index])
                   strparamlist.append(elem)
               strsent="_".join(strparamlist)
               
               strparam=[str(elem) for elem in param]
               paramstr=" ".join(strparam)
               
               code='python graphgeneratorpart.py {0} {1} {2} {3} {4} {5} {6} {7}'.format(strsent,algodirname,samplecount,algo,algotype,edgecount,n,paramstr)
               retval=int(os.system(code))
               myqu.popleft()
               if retval==256:
                  fastffflag=True
                  myqu.append(True)
               elif retval==65280:
                  myqu.append(False)
               else:
                  print "ERROR:Return value is {0} but this is not valid!!".format(retval)
                  exit(1)
                     
               if fastffflag:
                  flag=False
                  for elem in list(myqu):
                      if elem==True:
                         flag=True
                         break
                      elif elem==False:
                         pass 
                      else:
                         print "ERROR:Elem value is {0} but it should have been boolean!!".format(elem)
                         exit(1)
                  
                  if not flag:
                     break
                  

def nodeparamassigner(mode):
    paramhash={}
    if mode=="degree":
       #paramhash[126]=[6,10,26]
       paramhash[250]=[6,10,26]
       #paramhash[500]=[10,26,50]
       #paramhash[1000]=[10,50,100]
       #paramhash[2000]=[20,100,250]

       #paramhash[500]=[50]
       #paramhash[1000]=[100]
       #paramhash[2000]=[100,250]
    elif mode=="density":
       paramhash[126]=[0.05,0.1,0.2]
       paramhash[250]=[0.05,0.1,0.2]
       paramhash[500]=[0.05,0.1,0.2]
       paramhash[1000]=[0.05,0.1,0.2]
       paramhash[2000]=[0.05,0.1,0.2]
    return paramhash


pbsfolder="pbsfolder"
samplecount=20
#mode="density" #mode becoming density or degree does not change much
mode="degree"
createmode="graphwise" #This will create graphs from all models that have the same node and edgenumbers without considering specific parameters for that model
#createmode="paramwise" #This will create graphs from each model for given model parameters. Graphs do not need to have same node and edge numbers
algos=["dmc","smw","rds","lpa","undirected_ff","directed_ff","kronecker"] #geo,dmr,plc,opa,bpa is not being used right now
algos=["dmc","lpa","smw","rds"]
directedalgos=["directed_ff","kronecker"]
allnodenums=[126,500,1000,2000] #those are node numbers for paramwise case

#In order to create paramwise graphs, we don't need to create graphs from stratch. We can use graphwise graphs created.
if __name__ == "__main__":
    if not os.path.exists(pbsfolder):
       os.makedirs(pbsfolder)
    
    graphfolder="{0}_{1}_{2}".format("synthetic",createmode,"graphs")
    
    if createmode=="graphwise":
       paramhash=nodeparamassigner(mode)
       params=[]
       for key in paramhash.keys():
           for value in paramhash[key]:
               params.append((key,value)) 
       if mode=="degree":
          for nodenum,degree in params:
              underfoldername="{0}/datan{1}deg{2}".format(graphfolder,nodenum,degree)
              edgenum=(nodenum*degree)/2 
              data_generate(nodenum,edgenum,samplecount,underfoldername)
       elif mode=="density":
          for nodenum,density in params:
              underfoldername="{0}/datan{1}den{2}".format(graphfolder,nodenum,density)
              edgenum=int(density*(nodenum*(nodenum-1))/2) #this edgenum is for undirected case
              data_generate(nodenum,edgenum,samplecount,underfoldername)
    elif createmode=="paramwise": #paramwise will fix the number of nodes for all algorithms
       for nodenum in allnodenums:
           underfoldername="{0}/datan{1}".format(graphfolder,nodenum)  
           data_generate_paramwise(nodenum,samplecount,underfoldername) 
           

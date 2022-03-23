#When there are more than one parameters that satisfy the given node and edge numbers, random sampling is done.
import networkx as nx
import os
import math
import numpy as np
import itertools
import random
import sys
from optparse import OptionParser
import myutilities as myutil
import operator

mode="degree"
#mode="density"
extensions=[".gml",".sparse6",".graph6"]
diflimit={}
diflimit["dmc"]=0.25
diflimit["undirected_ff"]=0.2
diflimit["directed_ff"]=0.5
diflimit["kronecker"]=1.0
newgraphdir="newgraphs"
directedalgos=["directed_ff","kronecker"]
algos2select=["dmc","dmr","undirected_ff","directed_ff","kronecker"]
algosnot2select=["lpa","rds","smw"]

if __name__ == "__main__":
    assert len(sys.argv)==3
    samplecount=int(sys.argv[1]) #number of samples
    maingraphdir=sys.argv[2] #graph folder directory
    print "samplecount is {0}".format(samplecount)
    
    for data in myutil.listdirectories(maingraphdir):
        print data
        if data.find("yedek")!=-1:
           continue
        if mode=="degree":
           tempsplit=data.replace("datan","").split("deg")
           n1=int(tempsplit[0])
           deg=int(tempsplit[1])
           originaledge=n1*deg/2 #directed or undirected??
        elif mode=="density":
           tempsplit=data.replace("datan","").split("den")
           n1=int(tempsplit[0])
           den=float(tempsplit[1])
           originaledge=((n1*(n1-1))/2)*den
           
        datapath="{0}/{1}".format(maingraphdir,data)
        for algo in myutil.listdirectories(datapath):
            if algo in directedalgos:
               edgecount=2*originaledge
            else:
               edgecount=originaledge
            algopath="{0}/{1}".format(datapath,algo)
            valuehash={}
            newalgopath="{0}/{1}/{2}".format(newgraphdir,data,algo)
            if not os.path.exists(newalgopath):
               os.makedirs(newalgopath)
            for param in myutil.listdirectories(algopath):
                parampath="{0}/{1}".format(algopath,param)
                files=[]
                for myfile in myutil.listfiles(parampath):
                    for extension in extensions: 
                        if myfile.find(extension)!=-1:
                           files.append(myfile)
                           break       
                assert len(files)==samplecount

                if algo in  algosnot2select:
                   pass
                elif algo in algos2select:
                   splitted=param.split("_")
                   if algo=="dmc":
                      n2=int(splitted[0].replace("n",""))
                      qcon=float(splitted[1].replace("qcon",""))
                      qdel=float(splitted[2].replace("qdel",""))
                   elif algo=="dmr":
                      n2=int(splitted[0].replace("n",""))
                      qdel=float(splitted[1].replace("qdel",""))
                      qnew=float(splitted[2].replace("qnew",""))
                   elif algo=="undirected_ff":
                      n2=int(splitted[0].replace("n",""))
                      p=float(splitted[1].replace("p",""))
                   elif algo=="directed_ff":
                      n2=int(splitted[0].replace("n",""))
                      pbackward=float(splitted[1].replace("pbackward",""))
                      pforward=float(splitted[1].replace("pforward",""))
                   elif algo=="kronecker":
                      n2=int(splitted[0].replace("n",""))
                      p1=float(splitted[1].replace("p1",""))
                      p2=float(splitted[2].replace("p2",""))
                      p3=float(splitted[3].replace("p3",""))
                      p4=float(splitted[4].replace("p4",""))
               
                   assert n1==n2
                   originalnode=n1
                   total=0.0
                   for graphfile in files:
                       graphpath="{0}/{1}".format(parampath,graphfile)
                       if graphfile.find(".gml")!=-1:
                          G=nx.read_gml(graphpath)
                       elif graphfile.find(".sparse6")!=-1:
                          G=nx.read_sparse6(graphpath)
                       elif graphfile.find(".graph6")!=-1:
                          G=nx.read_graph6(graphpath)
                       else:
                          print "{0} extension is unknown!!".format(graphfile)
                          exit(1)
                       total += abs(originalnode-G.number_of_nodes()) + abs(edgecount-G.number_of_edges())
                   valuehash[param]=total
                else:
                   print "This algo {0} is not known!!".format(algo)
                   exit(1)
                
            if algo in algos2select:
               print data 
               print algo 
               minvalue=10000000000  
               minparam=-1
               for param in valuehash.keys():
                   if valuehash[param] < minvalue:
                      minvalue=valuehash[param]
                      minparam=param
               assert minparam!=-1
    
               difhash={}
               splitted=minparam.split("_")   
               if algo=="dmc":
                  minqcon=float(splitted[1].replace("qcon",""))
                  minqdel=float(splitted[2].replace("qdel",""))
                  for param in valuehash.keys():
                      if param!=minparam:
                         splitted=param.split("_")
                         paramqcon=float(splitted[1].replace("qcon",""))
                         paramqdel=float(splitted[2].replace("qdel","")) 
                         dif=abs(paramqdel-minqdel)
                         difhash[param]=dif
               elif algo=="dmr":
                  minqdel=float(splitted[1].replace("qdel",""))
                  minqnew=float(splitted[2].replace("qnew",""))
                  for param in valuehash.keys():
                      if param!=minparam:
                         splitted=param.split("_")
                         paramqdel=float(splitted[1].replace("qdel",""))
                         paramqnew=float(splitted[2].replace("qnew","")) 
                         dif=abs(paramqdel-minqdel)
                         difhash[param]=dif
               elif algo=="undirected_ff":
                  minp=float(splitted[1].replace("p",""))
                  for param in valuehash.keys():
                      if param!=minparam:
                         splitted=param.split("_")
                         paramp=float(splitted[1].replace("p","")) 
                         dif=abs(paramp-minp)
                         difhash[param]=dif
               elif algo=="kronecker":
                  minp1=float(splitted[1].replace("p1",""))
                  minp2=float(splitted[2].replace("p2",""))
                  minp3=float(splitted[3].replace("p3",""))
                  minp4=float(splitted[4].replace("p4",""))
                  for param in valuehash.keys():
                      if param!=minparam:
                         splitted=param.split("_")
                         paramp1=float(splitted[1].replace("p1",""))
                         paramp2=float(splitted[2].replace("p2",""))
                         paramp3=float(splitted[3].replace("p3",""))
                         paramp4=float(splitted[4].replace("p4",""))
                         dif=abs(paramp1-minp1)+abs(paramp2-minp2)+abs(paramp3-minp3)+abs(paramp4-minp4)
                         difhash[param]=dif
               elif algo=="directed_ff":
                  minpbackward=float(splitted[1].replace("pbackward",""))
                  minpforward=float(splitted[2].replace("pforward",""))
                  for param in valuehash.keys():
                      if param!=minparam:
                         splitted=param.split("_")
                         parampbackward=float(splitted[1].replace("pbackward",""))
                         parampforward=float(splitted[2].replace("pforward",""))
                         dif=abs(parampforward-minpforward)+abs(parampbackward-minpbackward)
                         difhash[param]=dif

               keepparams=[minparam]
               sorteditems = sorted(difhash.iteritems(), key=operator.itemgetter(1),reverse=True)
               for param,value in sorteditems:
                   if value >= diflimit[algo]:
                      keepparams.append(param)
                   else:
                      break
               for param in myutil.listdirectories(algopath):
                   if param not in keepparams:
                      parampath="{0}/{1}".format(algopath,param)
                      os.system("rm -rf {0}".format(parampath))
            
            #Among all valid parameters left, we randomly select samplecount and renumerate their names.
            allfilepaths=[]
            for param in myutil.listdirectories(algopath):
                parampath="{0}/{1}".format(algopath,param)
                for myfile in myutil.listfiles(parampath):
                    for extension in extensions: 
                        if myfile.find(extension)!=-1:
                           filepath="{0}/{1}".format(parampath,myfile)
                           allfilepaths.append(filepath)
                           break
            random.shuffle(allfilepaths)    
            ourfilepaths=allfilepaths[0:samplecount]
            for index in range(0,len(ourfilepaths)):
                ourfilepath=ourfilepaths[index]
                extension=ourfilepath.split(".")[-1]
                newfilepath="{0}/{1}.{2}".format(newalgopath,index+1,extension)
                code="cp -r {0} {1}".format(ourfilepath,newfilepath)
                os.system(code)
             

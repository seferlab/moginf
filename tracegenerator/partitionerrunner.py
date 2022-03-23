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
import gzip
import cPickle


def returnallfilesrecursive(rootdir):
    import os
    import sys
    filelist=[]
    folderCount=0
    for root, subFolders, files in os.walk(rootdir):
        folderCount += len(subFolders)
        for filename in files:
            f = os.path.join(root,filename)
            filelist.append(f)
    return filelist


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


tracefolderprefix="../tracegenerator"
newtracefolderprefix="newone"
tracefolder="traces_real_static_si_edge"
#tracefolder="traces_real_dynamic_si_edge"
tracefolder="traces_synthetic_graphwise_si_edge"
spreadmodel="si"
if __name__ == "__main__":
    newtracefolder="{0}_{1}".format(newtracefolderprefix,tracefolder)
    tracefolderpath="{0}/{1}".format(tracefolderprefix,tracefolder)
    filenames=returnallfilesrecursive(tracefolderpath)
    selectedones=set()
    for filename in filenames:
        if filename.endswith(".pkl"):
           selectedones.add(filename)
    print len(selectedones)
    
    for filepath in filenames:
        parts=filepath.split("/")
        predirname="/".join(parts[2:-1])
        tracefilename=parts[-1]
        outfilename=tracefilename.replace(".pkl","").split("_")[0]
        tracenum=int(tracefilename.replace(".pkl","").split("_")[2])
        newdirname="{0}/{1}/{2}".format(newtracefolderprefix,predirname,outfilename)
        if not os.path.exists(newdirname):
           os.makedirs(newdirname)
        
        infile=gzip.open(filepath,"rb")
        mytraces=cPickle.load(infile) #hash keyed by nodenumber and params.value will be list since there might be more than one spread from same node with same parameter
        #print mytraces.keys()[0:10]
        #if infertype in ["spreaprobunknown","difprobpartial"]:
        #    uniquecount=cPickle.load(infile) #this will only be used when estimating the scores 
        infile.close()

        #count=0
        #mykeys=set()
        #for elem in mytraces.keys():
        #    count+=len(mytraces[elem])
        #    mykey=elem[1:]
        #    mykeys.add(mykey)
        for keyname in mytraces.keys():
            startnode=int(keyname[0][0])
            spreadprob=float(keyname[1])
            sentparam=[("spreadprob",spreadprob),("s2i",keyname[2])]
            print sentparam
            paramfolder=returnspecificfolder(sentparam)
            paramfolder2="{0}/{1}".format(newdirname,paramfolder)
            #print len(mytraces[keyname])
            #print type(mytraces[paramkey][0])
            #print mytraces[keyname][0].keys()
            if not os.path.exists(paramfolder2):
               os.makedirs(paramfolder2)
            for mytrace in mytraces[keyname]:
                while True:
                   fileindex=random.randint(0,10000000)
                   newfilename="{0}_{1}.pkl".format(startnode,fileindex)
                   newfilepath="{0}/{1}".format(paramfolder2,newfilename)
                   if not os.path.exists(newfilepath):
                      outfile=open(newfilepath,"wb")
                      cPickle.dump(mytrace,outfile)
                      outfile.close()
                      break 

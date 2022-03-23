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

extensions=[".gml",".sparse6",".graph6"]
newgraphdir="newgraphs"
directedalgos=["directed_ff","kronecker"]

if __name__ == "__main__":
    assert len(sys.argv)==3
    samplecount=int(sys.argv[1]) #number of samples
    maingraphdir=sys.argv[2] #graph folder directory
    print "samplecount is {0}".format(samplecount)
    
    for data in myutil.listdirectories(maingraphdir):
        print data
        if data.find("yedek")!=-1:
           continue
        datapath="{0}/{1}".format(maingraphdir,data)
        for algo in myutil.listdirectories(datapath):
            algopath="{0}/{1}".format(datapath,algo)
            newalgopath="{0}/{1}/{2}".format(newgraphdir,data,algo)
            if not os.path.exists(newalgopath):
               os.makedirs(newalgopath)
            allfilepaths=[]
            for param in myutil.listdirectories(algopath):
                parampath="{0}/{1}".format(algopath,param)
                allfiles=[]
                for myfile in myutil.listfiles(parampath):
                    for extension in extensions: 
                        if myfile.find(extension)!=-1:
                           filepath="{0}/{1}".format(parampath,myfile)
                           allfiles.append(filepath)
                           break
                assert len(allfiles)==samplecount
                allfilepaths.extend(allfiles)
            random.shuffle(allfilepaths)    
            ourfilepaths=allfilepaths[0:samplecount]
            for index in range(0,len(ourfilepaths)):
                ourfilepath=ourfilepaths[index]
                extension=ourfilepath.split(".")[-1]
                newfilepath="{0}/{1}.{2}".format(newalgopath,index+1,extension)
                code="cp -r {0} {1}".format(ourfilepath,newfilepath)
                os.system(code)             

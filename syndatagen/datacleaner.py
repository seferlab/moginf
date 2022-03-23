import networkx as nx
import os
import math
import numpy as np
import itertools
import random
import sys
import myutilities as myutil


def learndatacontent(maingraphdir):  #print number of files of each algo and parameter
    datapaths=[]
    for datadir in myutil.listdirectories(maingraphdir):
        if datadir.startswith("datan"):
           datapath="{0}/{1}".format(maingraphdir,datadir)
           datapaths.append(datapath) 
       
    for datapath in datapaths:
        print "Data {0}".format(datapath)
        algodirs=myutil.listdirectories(datapath)
        for algo in algodirs:
            algopath="{0}/{1}".format(datapath,algo)
            paramdirs=myutil.listdirectories(algopath)
            print "{0} : {1}".format(algo,len(paramdirs))


def datacleaner(samplecount,maingraphdir):
    datapaths=[]
    for dataname in myutil.listdirectories(maingraphdir):
        if dataname.startswith("datan"):
           datapath="{0}/{1}".format(maingraphdir,dataname) 
           datapaths.append(datapath) 
    
    for datapath in datapaths:
        print datapath 
        algodirs=myutil.listdirectories(datapath)
        for algo in algodirs:
            print algo
            algopath="{0}/{1}".format(datapath,algo)
            params=myutil.listdirectories(algopath)
            for param in params:
                parampath="{0}/{1}".format(algopath,param)
                files=[]
                for myfile in myutil.listfiles(parampath):
                    for extension in extensions: 
                        if myfile.find(extension)!=-1:
                           files.append(myfile)
                           break
                        
                if len(files)!=samplecount:
                   os.system("rm -rf {0}".format(parampath))             


extensions=[".gml",".sparse6",".graph6"]
if __name__ == "__main__":
    print "Looking for file extensions: {0}".format(extensions)
    assert len(sys.argv)==3
    samplecount=int(sys.argv[1]) #number of samples
    maingraphdir=sys.argv[2] #graph folder directory
    print "samplecount is {0}".format(samplecount)
    
    datacleaner(samplecount,maingraphdir)
    learndatacontent(maingraphdir)    
    

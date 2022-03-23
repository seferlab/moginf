import networkx as nx
import math
import myutilities as myutil
import pickle

for filename in myutil.listfiles("."):
    if filename.find(".pkl")!=-1:
       print filename 
       filepath="{0}".format(filename) 
       infile=open(filepath,"rb")
       alltraces=pickle.load(infile)
       infile.close()
       
       print len(alltraces.keys())
       print alltraces.keys()[3:4]
       print alltraces.values()[3:4]
       val=alltraces.values()[3][0]
       print val
       for node in val.keys():
           infect=val[node]["infect"]
           recover=val[node]["recover"]
           if infect!=recover:
              print "node is {0}".format(node) 
              exit(1)

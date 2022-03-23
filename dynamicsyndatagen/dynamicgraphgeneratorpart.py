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

#Fused Lasso based Smoothing evolution over time
#Edge existence and edge nonexistence interval lenghts will both be obtained from interval distribution that is geometric with mean tcount/5. Interval lengths are related to over length of evolution(tcount).
#Starting graph will be obtained from erdos-renyi with prob 0.1. The model of starting graph won't affect the performance that much.
#there will be 4 parameters determining how to progress graph over time
def temporal_fusedlasso(intervaldist,intervaldistparam,startgraphalgo,startgraphalgoparam):
    Gall=[]
    if intervaldist=="uniform":
       changeprob=intervaldistparam
       if startgraphalgo=="rds":
          G=nx.erdos_renyi_graph(startnodenum,startgraphalgoparam)
       else:
          print "this startgraph algo is unknonw!!".format(startgraphalgo)
          exit(1)
       allnodes=G.nodes()   
       intervalcount=int(math.ceil(float(tcount)/stepsize))
       curG=nx.Graph(G)
       for interval in range(0,intervalcount):
           if interval==0:
              for timeindex in range(0,stepsize):
                  time=(interval*stepsize)+timeindex
                  Gall.append(nx.Graph(curG))
              continue    
           for timeindex in range(0,stepsize):
               #time=(interval*stepsize)+timeindex
               if timeindex==0:
                  for index1 in range(0,len(allnodes)):
                      node1=allnodes[index1]
                      for index2 in range(index1+1,len(allnodes)):
                          node2=allnodes[index2]
                          if random.random()<=(1.0-changeprob):
                             continue 
                          if curG.has_edge(node1,node2):
                             curG.remove_edge(node1,node2)
                          elif not curG.has_edge(node1,node2):
                             curG.add_edge(node1,node2)
               Gall.append(nx.Graph(curG))
    elif intervaldist=="expo": 
       meaninterval=tcount/float(intervaldistparam)
       Gall=[]
       G=nx.erdos_renyi_graph(startnodenum,startgraphalgoparam)
       allnodes=G.nodes()
       for time in range(0,tcount):
           Gall.append(nx.Graph())
       edgeintervals={} #edge intervals will be assigned to all possible edges
       for index1 in range(0,len(allnodes)):
           node1=allnodes[index1]
           for index2 in range(index1+1,len(allnodes)):
                   node2=allnodes[index2]
                   durations=[]
                   if G.has_edge(node1,node2):
                      nextmode=1
                   else:
                      nextmode=0
                   start=0   
                   while start<tcount:
                      while True:
                         interval=int(round(random.expovariate(1.0/meaninterval)))
                         if interval!=0:
                            break
                      end=start+interval
                      if end>tcount:
                         end=tcount
                      if nextmode==1:
                         durations.append((start,end))
                         nextmode=0
                      elif nextmode==0:
                         nextmode=1
                      start=end
                   edgeintervals[(index1,index2)]=durations
           
       for node1,node in edgeintervals.keys():
           for start,end in edgeintervals[(node1,node2)]:
               for time in range(start,end):
                   Gall[time].add_edge(node1,node2)
    else:
        print "intervaldist {0} is unknown".format(intervaldist)
        exit(1)
        
    return Gall
 

#Snap implementation. Currently, we run seperately for each nodenumbers(not iteratively)!!
def temporal_directed_ff(pbackward,pforward):
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        realtime=startnodenum+interval
        intervalhash[realtime]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[realtime]+=1
  
    for nodenum in sorted(intervalhash.keys()):
        gcount=intervalhash[nodenum]    
        myrand=random.random()
        tempfile="tempout{0}".format(myrand)
        code="./forestfire -o:{0} -n:{1} -f:{2} -b:{3}".format(tempfile,nodenum,pforward,pbackward)
        os.system(code)
        G=nx.DiGraph()
        file=open(tempfile,"r")
        for line in file:
            line=line.rstrip()
            if line.startswith("#"):
               continue
            node1,node2=line.split("\t")
            node1=int(node1)
            node2=int(node2)
            G.add_edge(node1,node2)
        code="rm -rf {0}".format(tempfile)
        os.system(code)
        
        assert G.number_of_nodes()==nodenum
        if nx.number_weakly_connected_components(G)!=1:#if disconnected make it connected
           components=nx.weakly_connected_components(G)
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)
        for mycount in range(0,gcount):
            Gall.append(nx.DiGraph(G))      

    assert len(Gall)==tcount
    return Gall
   
    
#The SNAP implementation
#It always creates directed graph
#NOT WORKING RIGHT NOW, SINCE WE ARE INCREMENTING ONE BY ONE THIS CAN NOT DO IT
def temporal_kronecker(p1,p2,p3,p4):
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        realtime=startnodenum*(2**interval)
        intervalhash[realtime]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[realtime]+=1
    for nodenum in sorted(intervalhash.keys()):
        gcount=intervalhash[nodenum]
        myrand=random.random()
        itcount=int(math.ceil(math.log(nodenum,2)))
        tempfile="tempout{0}".format(myrand) 
        code="./krongen -o:{0} -m:\"{1} {2}; {3} {4}\" -i:{5}".format(tempfile,p1,p2,p3,p4,itcount)
        os.system(code)
        G=nx.DiGraph()
        file=open(tempfile,"r")
        for line in file:
            line=line.rstrip()
            if line.startswith("#"):
               continue
            node1,node2=line.split("\t")
            node1=int(node1)
            node2=int(node2)
            G.add_edge(node1,node2)
        code="rm -rf {0}".format(tempfile)
        os.system(code)

        #if disconnected make it connected
        #but graph must have the required number of nodes.
        if (G.number_of_nodes()>=int(nodenum)) and nx.number_weakly_connected_components(G)!=1:
           components=nx.weakly_connected_components(G)
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)           
        for mycount in range(0,gcount):
            Gall.append(nx.DiGraph(G))
            
    assert len(Gall)==tcount
    return Gall

#Continously evolves and returns the evolved graphs at various steps
def temporal_dmc(qcon,qmod):
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        intervalhash[interval]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[interval]+=1
    
    #Create the seed graph (node indices starts from 0)
    G=nx.Graph()
    G.add_edge(0,1)

    maxtime=startnodenum+intervalcount+1
    currentnode=G.number_of_nodes()
    for count in range(currentnode,maxtime):
        G.add_node(count)
         
        nodenum=random.randrange(0,count)
        if G[nodenum]:
           neighlist=G[nodenum].keys()
        else:
           neighlist=[]
            
        for elem in neighlist:
            G.add_edge(count,elem)
         
        for elem in neighlist:
            coin=random.random() 
            if coin<0.5:
               num=random.random()
               if num <= qmod:
                  G.remove_edge(nodenum,elem)                   
            else:
               num=random.random()
               if num <= qmod:
                  G.remove_edge(count,elem)   
             
        num=random.random()
        if num <= qcon:
             G.add_edge(count,nodenum)

        if (count-startnodenum+1) not in intervalhash.keys():
           continue
        
        for addinex in range(0,intervalhash[count-startnodenum+1]):
            Gall.append(nx.Graph(G))
               
    assert len(Gall)==tcount
    #for index in range(1,len(Gall)):
    #    G1=Gall[index-1]
    #    G2=Gall[index]
    #    interedges=set(G1.edges()).difference(set(G2.edges()))
    #    print "{0} -> {1}, {2}, {3}".format(index,len(interedges),G1.number_of_nodes(),G1.number_of_edges())
    #print "last: -> {0}, {1}".format(G2.number_of_nodes(),G2.number_of_edges())   
    #exit(1)
    
    #For all graphs; if disconnected make it connected after evolution is completely done!!
    for G in Gall:
        if nx.number_connected_components(G)!=1:
           components=nx.connected_components(G)  
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)
           
    return Gall 


def temporal_dmr(floatn,qdel,qnew): #NOT implemented!!!
    n=int(floatn)
    G=nx.cycle_graph(5)
    currentnode=G.number_of_nodes()
    for count in range(currentnode,n):
        G.add_node(count)
        alpha=float(qnew)/(count+1)
        n0=random.randrange(0,count)
        if G[n0]:
             neighlist=G[n0].keys()
        else:
             neighlist=[]
        for elem in neighlist:
             G.add_edge(count,elem)
             if random.random() <= qdel:
                G.remove_edge(count,elem)
                
        remainall=G.nodes()
        if G[count]:
             neighlist=G[count].keys()
        else:
             neighlist=[]
        neighlist.append(count)     
        remainall=list(set(remainall) - set(neighlist))     
        
        for remain in remainall:
            if random.random() <= alpha:
               G.add_edge(remain,count)

    #if disconnected make it connected
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)
           
    return G

def plc(floatm,floatn,p): #NOT implemented
    n=int(floatn)
    m=int(floatm)
    G=nx.powerlaw_cluster_graph(n,m,p)
    return G

    
def temporal_rds(p): #our own iterative implementation
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        intervalhash[interval]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[interval]+=1
      
    G=nx.Graph()
    G.add_node(0)
    maxtime=startnodenum+intervalcount+1
    currentnode=G.number_of_nodes()
    for count in range(currentnode,maxtime):
        G.add_node(count)
        for node in G.nodes():
            if node!=count:
               if random.random()<p:
                  G.add_edge(count,node) 

        if (count-startnodenum+1) not in intervalhash.keys():
           continue
           
        for addinex in range(0,intervalhash[count-startnodenum+1]):
            Gall.append(nx.Graph(G))

    assert len(Gall)==tcount
     
    #For all graphs; if disconnected make it connected after evolution is completely done!!
    for G in Gall:
        if nx.number_connected_components(G)!=1:
           components=nx.connected_components(G)  
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)
    return Gall  

        
def temporal_lpa(k): #our iterative implementation of barabasi albert graph
    k=int(k)
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        intervalhash[interval]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[interval]+=1
      
    #Create the seed graph (node indices starts from 0)
    G=nx.Graph()
    for node in range(0,k):
        G.add_node(node)
        
    maxtime=startnodenum+intervalcount+1    
    currentnode=G.number_of_nodes()
    for count in range(currentnode,maxtime):
        if count==k:
           for index in range(0,k):
               randnode=random.choice(G.nodes())
               G.add_edge(count,randnode) 
           continue
        
        deghash={}
        total=0.0
        cumlist=[]
        denom=2.0*G.number_of_edges()
        for node in G.nodes():
            deghash[node]=float(G.degree(node))/denom
            cumlist.append((node,total,total+deghash[node])) 
            total += deghash[node]
        for index in range(0,k):
            randnum=random.random()
            node2=-1
            for node,start,end in cumlist:
                if start<=randnum and randnum<=end:  
                   node2=node
                   break
            assert node2!=-1    
            G.add_edge(count,node2)    

        if (count-startnodenum+1) not in intervalhash.keys():
           continue
        
        for addinex in range(0,intervalhash[count-startnodenum+1]):
            Gall.append(nx.Graph(G))
               
    assert len(Gall)==tcount
    
    #For all graphs; if disconnected make it connected after evolution is completely done!!
    for G in Gall:
        if nx.number_connected_components(G)!=1:
           components=nx.connected_components(G)  
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)
           
    return Gall


def opa(n): #NOT IMPLEMENTED
    #original preferential attachment
    """ Creates a random network of size 'n' that follows the standard
        preferential attachment model of evolution. No parameters, no daddy.
    """

    # Starting configuration, dumbbell.
    G=nx.Graph()
    G.add_node(0)
    G.add_node(1)
    G.add_edge(0,1)

    zeros = 0
    for v in xrange(2,n):
        G.add_node(v)

        # Add each to existing node u with probability deg(u)/2m
        num_edges = G.num_edges # value will otherwise change in the loop.
        for u in G:
            if u == v: continue
            if random.random() < (len(G.neighbors(u)) / (2*num_edges)):
                G.add_edge(v,u)

        # If no edges were added, just add one random edge (component).
        if len(G.neighbors(v)) == 0:
            zeros += 1
            rand_node = v
            while (rand_node == v): rand_node = G.random_node()
            G.add_edge(v,rand_node)

    assert G.num_nodes == n
    return G


def bpa(n,k): #NOT implemented
    """ Creates a random network of size 'n' that follows the Barabasi preferential
        attachment model of evolution. In each step, a node enters and links to
        'k' existing nodes.
    """
    assert n >= k
    a = 0 # smooth parameter.

    G = Graph("bpa-%i" %(n))
    Order = {} # Mapping from node to its arrival time.

    # Starting configuration: Erdos-Renyi random graph with k nodes.
    #histogram = G.make_random(k)
    histogram = G.make_clique(k) # k+1 nodes added. so start with k+2.
    assert len(histogram) > 0
    #assert G.num_nodes == k
    assert G.num_nodes == k+1

    #for v in xrange(k+1,n+1):
    for v in xrange(k+2,n+1):
        assert v not in G
        G.add_node(v)
        for ai in xrange(a): histogram.append(v)

        # histogram contains node v 'degree(v)'-times, essentially creating a histogram.
        added = set([v])
        while len(added) < k+1:
            u = random.sample(histogram,1)[0]
            if u in added: continue
            G.add_edge(v,u)
            histogram.append(v)
            histogram.append(u)
            added.add(u)

    assert G.num_nodes == n

    return (G,Order)


def geometric(p):
    return int(round(math.log(random.random()) / math.log(1-p)))
def temporal_undirected_ff(p):
    """ Creates a random network of size 'n' tht follows the basic forest
        fire model of evolution, with burning probability 'p'.
        UNDIRECTED VERSION.
    """
    Gall=[]
    intervalhash={} 
    intervalcount=int(math.ceil(float(tcount)/stepsize))
    for interval in range(0,intervalcount):
        intervalhash[interval]=0
        for timeindex in range(0,stepsize):
            time=(interval*stepsize)+timeindex
            if time<tcount:
               intervalhash[interval]+=1
      
    #Create the seed graph (node indices starts from 0)
    G=nx.Graph()
    G.add_edge(0,1)    
    maxtime=startnodenum+intervalcount+1    
    currentnode=G.number_of_nodes()
    for v in range(currentnode,maxtime):#count
        u = random.choice(G.nodes()) #pick node first, so random node != new node.
        G.add_node(v)
        G.add_edge(v,u) # link to anchor.
        
        visited = set([u,v]) # visited nodes, to prevent cycling.
        queue = deque([u])
        
        # start the forest fire.
        while len(queue) > 0:
            u = queue.popleft()

            # Flip a geometric coin with mean 1-p and link to that many neighbors.
            x = geometric(1-p)
              
            neighbors = G.neighbors(u)
            to_visit = random.sample(neighbors,x) if len(neighbors) > x else neighbors

            for candidate in to_visit:
                if candidate in visited:
                    continue
                visited.add(candidate)
                queue.append(candidate)
                if v!= candidate:
                   G.add_edge(v,candidate)

        if (v-startnodenum+1) not in intervalhash.keys():
           continue
        
        for addinex in range(0,intervalhash[v-startnodenum+1]):
            Gall.append(nx.Graph(G))
               
    assert len(Gall)==tcount
    
    #For all graphs; if disconnected make it connected after evolution is completely done!!
    for G in Gall:
        if nx.number_connected_components(G)!=1:
           components=nx.connected_components(G)  
           for index in range(1,len(components)):
               node1=random.choice(components[index-1])
               node2=random.choice(components[index])
               G.add_edge(node1,node2)
           
    return Gall


#global parameters startnodenum,stepsize,tcount
startnodenum=-1
stepsize=-1
tcount=-1

def isfloat(mystr):
    try:
       float(mystr)
    except ValueError:
       return False
    return True

#WARNING: After stepsize modification, Algorithms other than dmc and fusedlasso does not work!!
#IMPORTANT::

if __name__ == "__main__":
    argcount=len(sys.argv)
    #print "initial info::"
    #for temp in sys.argv:
    #    print temp
    strsent=sys.argv[1]
    algodirname=sys.argv[2]
    samplecount=int(sys.argv[3])
    algo=sys.argv[4]
    algotype=sys.argv[5]
    startnodenum=int(sys.argv[6])
    stepsize=int(sys.argv[7]) 
    tcount=int(sys.argv[8])
    param=[]
    for index in range(9,argcount):
        if isfloat(sys.argv[index]):
           param.append(float(sys.argv[index]))
        else:    
           param.append(sys.argv[index])
     
    paramdirname="{0}/{1}".format(algodirname,strsent)
    if not os.path.exists(paramdirname):
       os.makedirs(paramdirname)
    
    for i in range(0,samplecount):
        sampledirname="{0}/{1}".format(paramdirname,i+1)
        if not os.path.exists(sampledirname):
           os.makedirs(sampledirname)
        algofuncname="temporal_{0}".format(algo)
        Gall=globals()[algofuncname](*param)
        for time in range(0,len(Gall)):
            G=Gall[time]
            print "dynamic algo {0}; param {1}; node: {2}; edge: {3}".format(algo,param,G.number_of_nodes(),G.number_of_edges())
            if algotype=="directed":
               graphfilename="{0}.gml".format(time)
               graphfilepath="{0}/{1}".format(sampledirname,graphfilename)
               nx.write_gml(G,graphfilepath)
               code="gzip -9 -v -c {0} > {0}.gz".format(graphfilepath)
               os.system(code)
               code="rm -rf {0}".format(graphfilepath)
               os.system(code)
            elif algotype=="undirected": 
               #outputs the graph either in graph6 or sparse6 format(which size is minimum, it just outputs!!)
               graphfilename1="{0}.graph6".format(time)
               graphfilepath1="{0}/{1}".format(sampledirname,graphfilename1)
               amtogpath="./amtog"
               tempstr1=myutil.convertgraph6(G,amtogpath)
               graphfilename2="{0}.sparse6".format(time)
               graphfilepath2="{0}/{1}".format(sampledirname,graphfilename2)
               tempstr2=myutil.convertsparse6(G,amtogpath)
               if len(tempstr1)>len(tempstr2):
                  file=open(graphfilepath2,"w")
                  file.write("{0}\n".format(tempstr2))
                  file.close()
               else:    
                  file=open(graphfilepath1,"w")
                  file.write("{0}\n".format(tempstr1))
                  file.close()
                  

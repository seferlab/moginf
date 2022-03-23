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
import string

#Graph returned won't be number of edges since input edgenum is number of multiedges
def rtg(beta,edgenum,k,q):
    matlabpath="/opt/stow/matlab-r2012a/bin/matlab"
    #matlabpath="/Applications/MATLAB_R2012a.app/bin/matlab"
    outfiledir="{0}/".format(''.join(random.choice(string.ascii_uppercase) for x in range(20)))
    codepath="RTG"
    if not os.path.exists(outfiledir):
       os.makedirs(outfiledir) 
    isbipartite=0
    isselfloop=0 
    timeticks=10
    code="{0} -r \"{1}({2},{3},{4},{5},{6},\'{7}\',{8},{9}) ; quit ;\" ".format(matlabpath,codepath,edgenum,k,beta,q,isbipartite,outfiledir,isselfloop,timeticks)
    os.system(code)
    filename="edgefile_{0}_{1}_{2}_{3}.txt".format(edgenum,k,beta,q)
    filepath="{0}{1}".format(outfiledir,filename)
    G=nx.Graph()
    file=open(filepath,"r")
    for line in file:
        line=line.rstrip()
        time,node1,node2=line.split()
        G.add_edge(node1,node2)
    file.close()	
    
    code="rm -rf {0}".format(filepath)
    os.system(code)
    code="rm -rf {0}".format(outfiledir)
    os.system(code)
    
    #if disconnected make it connected
    if nx.number_weakly_connected_components(G)!=1:
       components=nx.weakly_connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)           
    return G
    
    
#Snap implementation
def directed_ff(floatn,pbackward,pforward):
    myrand=random.random()
    tempfile="tempout{0}".format(myrand)
    code="./forestfire -o:{0} -n:{1} -f:{2} -b:{3}".format(tempfile,floatn,pforward,pbackward)
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
    assert G.number_of_nodes()==int(floatn)
    if nx.number_weakly_connected_components(G)!=1:
       components=nx.weakly_connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)           
    return G

#The SNAP implementation
#the graph is directed
#It is hard to use kronecker when we want to evolve graph that has required number of edges and nodes since it is tough to control the parameters.
def kronecker(floatn,p1,p2,p3,p4):
    myrand=random.random()
    itcount=int(math.ceil(math.log(floatn,2)))
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
    if (G.number_of_nodes()>=int(floatn)) and nx.number_weakly_connected_components(G)!=1:
       components=nx.weakly_connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)           
    return G

def dmc(floatn,qcon,qmod):
    #Create the seed graph (node indices starts from 0)
    n=int(floatn)
    G=nx.Graph()
    G.add_edge(0,1)

    currentnode=G.number_of_nodes()
    for count in range(currentnode,n):
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

    #if disconnected make it connected
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)  
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)
           
    return G  


def dmr(floatn,qdel,qnew):
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

#powerlawcluster graph
def plc(floatm,floatn,p):
    n=int(floatn)
    m=int(floatm)
    G=nx.powerlaw_cluster_graph(n,m,p)
    return G

def rds(floatm,floatn):
    n=int(floatn)
    m=int(floatm)
    p=float(2*m)/((n*n)-n)
    G=nx.erdos_renyi_graph(n,p)
    #if disconnected make it connected
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)
           
    return G  

def lpa(floatm,floatn):
    n=int(floatn)
    m=int(floatm)
    #delta=(n+1)**2-(4*m)
    delta=(n)**2-(4*m)
    if delta>=0:
       delta=math.sqrt(delta) 
       #k1=int(math.ceil((n+1-delta)/2.0))
       #k2=int(math.ceil((n+1+delta)/2.0))
       #k3=int(math.floor((n+1-delta)/2.0))
       #k4=int(math.floor((n+1+delta)/2.0))
       k1=int(math.ceil((n-delta)/2.0))
       k2=int(math.ceil((n+delta)/2.0))
       k3=int(math.floor((n-delta)/2.0))
       k4=int(math.floor((n+delta)/2.0))
       G1=nx.barabasi_albert_graph(n,k1)
       G2=nx.barabasi_albert_graph(n,k2)
       G3=nx.barabasi_albert_graph(n,k3)
       G4=nx.barabasi_albert_graph(n,k4)
       
       #if disconnected make it connected
       if nx.number_connected_components(G1)!=1:
          components=nx.connected_components(G1)
          for index in range(1,len(components)):
              node1=random.choice(components[index-1])
              node2=random.choice(components[index])
              G1.add_edge(node1,node2)
       if nx.number_connected_components(G2)!=1:
          components=nx.connected_components(G2)
          for index in range(1,len(components)):
              node1=random.choice(components[index-1])
              node2=random.choice(components[index])
              G2.add_edge(node1,node2)
       if nx.number_connected_components(G3)!=1:
          components=nx.connected_components(G3)
          for index in range(1,len(components)):
              node1=random.choice(components[index-1])
              node2=random.choice(components[index])
              G3.add_edge(node1,node2)
       if nx.number_connected_components(G4)!=1:
          components=nx.connected_components(G4)
          for index in range(1,len(components)):
              node1=random.choice(components[index-1])
              node2=random.choice(components[index])
              G4.add_edge(node1,node2)
              
       x1=abs(G1.number_of_edges()-m)
       x2=abs(G2.number_of_edges()-m)
       x3=abs(G3.number_of_edges()-m)
       x4=abs(G4.number_of_edges()-m)
       xlist=[x1,x2,x3,x4]
       mapG={}
       mapG[x1]=G1
       mapG[x2]=G2
       mapG[x3]=G3
       mapG[x4]=G4
    
       return mapG[min(xlist)]
    else:
       G=nx.Graph()
       return G
    
def smw(floatm,floatn,qrewire):
    n=int(floatn)
    m=int(floatm)
    k=int(math.ceil(float(2*m)/n))
    G=nx.watts_strogatz_graph(n,k,qrewire)
    
    #if disconnected make it connected
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)

    return G

def geo(floatn,radius):
    n=int(floatn)
    G=nx.random_geometric_graph(n,radius)

    #if disconnected make it connected
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)
    return G

def opa(n): #original preferential attachment
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

def bpa(n,k):
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
def undirected_ff(floatn,p):
    """ Creates a random network of size 'n' tht follows the basic forest
        fire model of evolution, with burning probability 'p'.
        UNDIRECTED VERSION.
    """
    n=int(floatn)
    G=nx.Graph()
    G.add_node(0)
    G.add_node(1)
    G.add_edge(0,1) # starting configuration, dumbell.
  
    for v in range(2,n):
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

    assert G.number_of_nodes()==n
    if nx.number_connected_components(G)!=1:
       components=nx.connected_components(G)  
       for index in range(1,len(components)):
           node1=random.choice(components[index-1])
           node2=random.choice(components[index])
           G.add_edge(node1,node2)
    return G


degreedeviate=2.0
krondegreedeviate=8.0
generaliterthres=40
ffiterthres=100

if __name__ == "__main__":
    argcount=len(sys.argv)
    strsent=sys.argv[1]
    algodirname=sys.argv[2]
    samplecount=int(sys.argv[3])
    algo=sys.argv[4]
    algotype=sys.argv[5]
    m=int(sys.argv[6])
    n=int(sys.argv[7])
    param=[]
    for index in range(8,argcount):
        param.append(float(sys.argv[index]))
    
    if m==-1:
       createmode="paramwise"
    else:
       createmode="graphwise"
    
    paramdirname="{0}/{1}".format(algodirname,strsent)
    if not os.path.exists(paramdirname):
       os.makedirs(paramdirname)
    
    if createmode=="graphwise":
       if algotype=="directed": #Directed graphs will have 2*m edges
          m*=2
       
       for i in range(0,samplecount):
           itnumber=0
           myflag=False
           while True:
              itnumber += 1
              G = globals()[algo](*param)
              print "algo {0}; param {1}; node: {2}; edge: {3}".format(algo,param,G.number_of_nodes(),G.number_of_edges())
              if (algo=="kronecker" and abs(G.number_of_edges()-m) <= ((krondegreedeviate/2.0)*n) and nx.number_connected_components(G)==1) or (algo!="kronecker" and abs(G.number_of_edges()-m) <= ((degreedeviate/2.0)*n) and nx.number_connected_components(G)==1 and G.number_of_nodes()==n):    
                 if algotype=="directed":
                    graphfilename="{0}.gml".format(i+1)
                    graphfilepath="{0}/{1}".format(paramdirname,graphfilename)
                    nx.write_gml(G,graphfilepath)
                    code="gzip -9 -v -c {0} > {0}.gz".format(graphfilepath)
                    os.system(code)
                    code="rm -rf {0}".format(graphfilepath)
                    os.system(code)
                 elif algotype=="undirected": 
                    #outputs the graph either in graph6 or sparse6 format(which size is minimum, it just outputs!!)
                    graphfilename1="{0}.graph6".format(i+1)
                    graphfilepath1="{0}/{1}".format(paramdirname,graphfilename1)
                    amtogpath="./amtog"
                    tempstr1=myutil.convertgraph6(G,amtogpath)
                    graphfilename2="{0}.sparse6".format(i+1)
                    graphfilepath2="{0}/{1}".format(paramdirname,graphfilename2)
                    tempstr2=myutil.convertsparse6(G,amtogpath)
                    if len(tempstr1)>len(tempstr2):
                       file=open(graphfilepath2,"w")
                       file.write("{0}\n".format(tempstr2))
                       file.close()
                    else:    
                       file=open(graphfilepath1,"w")
                       file.write("{0}\n".format(tempstr1))
                       file.close()
                    break
                   
              if (algo!="ff" and itnumber > generaliterthres) or (algo=="ff" and itnumber > ffiterthres):
                 myflag=True
                 break
        
           if myflag:
              #clean previously generated files because this parameter combination is not good enough to create samplecount number of samples of that graph generation process
              for graphfilename in myutil.listfiles(paramdirname):
                  if graphfilename.endswith(".graph6") or graphfilename.endswith(".sparse6") or graphfilename.endswith(".gml") or graphfilename.endswith(".gml.gz"):
                     graphfilepath="{0}/{1}".format(paramdirname,graphfilename) 
                     os.system("rm -rf {0}".format(graphfilepath)) 
              break

       if myflag:
          exit(-1)
       else:
          exit(1)
 
    elif createmode=="paramwise":
       if algotype=="directed": #Directed graphs will have 2*m edges
          m*=2
       
       for i in range(0,samplecount):
           G=globals()[algo](*param)
           print "paramwise algo {0}; param {1}; node: {2}; edge: {3}".format(algo,param,G.number_of_nodes(),G.number_of_edges())
           if algotype=="directed":
              graphfilename="{0}.gml".format(i+1)
              graphfilepath="{0}/{1}".format(paramdirname,graphfilename)
              nx.write_gml(G,graphfilepath)
              code="gzip -9 -v -c {0} > {0}.gz".format(graphfilepath)
              os.system(code)
              code="rm -rf {0}".format(graphfilepath)
              os.system(code)
           elif algotype=="undirected": 
              #outputs the graph either in graph6 or sparse6 format(which size is minimum, it just outputs!!)
              graphfilename1="{0}.graph6".format(i+1)
              graphfilepath1="{0}/{1}".format(paramdirname,graphfilename1)
              amtogpath="./amtog"
              tempstr1=myutil.convertgraph6(G,amtogpath)
              graphfilename2="{0}.sparse6".format(i+1)
              graphfilepath2="{0}/{1}".format(paramdirname,graphfilename2)
              tempstr2=myutil.convertsparse6(G,amtogpath)
              if len(tempstr1)>len(tempstr2):
                 file=open(graphfilepath2,"w")
                 file.write("{0}\n".format(tempstr2))
                 file.close()
              else:    
                 file=open(graphfilepath1,"w")
                 file.write("{0}\n".format(tempstr1))
                 file.close()
                   
       exit(1)
 

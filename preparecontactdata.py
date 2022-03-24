import os
import random
import networkx as nx

def makeGraph(fpath):
    """ makes graph out of file
    """
    G = nx.Graph()
    with open(fpath, "r") as infile:
        for line in infile:
            line = line.rstrip()
            splitted = line.split("\t")
            node1, node2 = int(splitted[0]), int(splitted[1])
            G.add_edge(node1, node2)

    return G


def generateHumanContactNetwork():
    """
    """
    # Make contact graph out of raw data
    fpath = os.path.join("rawcontactdata","staticWeightedEdgeList_at=1350_min=540_max=2159.txt")
    # Store sparsified version of graph
    Graw = makeGraph(fpath)
    G = nx.Graph()
    for node1,node2 in Graw.edges():
        if random.random() < 0.35:
            G.add_edge(node1,node2)
    print(G.number_of_nodes())
    print(G.number_of_edges())
    
    outpath = "human_contact_network.edgelist"
    nx.write_edgelist(G, outpath, data=False)

    
def generateLFRBenchmark(count):
    """ generates count LFR synethetic networks
    """
    node_count = 1000
    tau1=2
    tau2=1.5
    mus = [0.1, 0.2]
    min_degree_lower = 20
    min_degree_upper = 30
    max_degree_lower = 50
    max_degree_upper = 55
    
    for index in range(1,count+1):
        min_degree = random.randint(min_degree_lower, min_degree_upper)
        max_degree = random.randint(max_degree_lower, max_degree_upper)
        mu = random.choice(mus)
        G = nx.LFR_benchmark_graph(node_count, tau1, tau2, mu, min_degree=min_degree, max_degree=max_degree)
        
        outpath = "lfr_benchmark/lfr{0}.edgelist".format(index)
        nx.write_edgelist(G, outpath, data=False)

    
#generateHumanContactNetwork()
#generateLFRBenchmark(20)


